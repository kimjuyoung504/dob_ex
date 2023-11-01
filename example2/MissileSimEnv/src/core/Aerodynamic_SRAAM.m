classdef Aerodynamic_SRAAM < handle
    properties
        aeroTable
        refLength = 0.1524  % % Reference Length for Moment Derivatives [m]
        refArea =  0.01824  % Reference Area for Aero Coefficients [m^2]
        Xcg_ref = 1.536     % Launch CG aft of Vehicle Nose   [m]    
    end
    properties (Hidden)
        MachLim
        AlphaLim
        noiseStd = 0
        biasRate = 0
        WindShearFcn
    end
    methods
        function obj = Aerodynamic_SRAAM(coeffdatafilename)
            splitedName = split(coeffdatafilename, '.');
            ext = splitedName{end};
            if strcmp(ext, 'mat')
                load(coeffdatafilename, 'Mach','Alpha','ca0','caa','cad','caoff','cyp',...
                'cydr','cn0','cnp','cndq','cllap','cllp','clldp','clm0','clmp',...
                'clmq','clmdq','clnp','clnr','clndr');
            elseif strcmp(ext, 'm')
                run(coeffdatafilename);
            end
            
            samplePoints = {Mach(:), Alpha(:)};
            obj.aeroTable.ca0 = DataTable(ca0(:), samplePoints{1});
            obj.aeroTable.caa = DataTable(caa(:), samplePoints{1});
            obj.aeroTable.cad = DataTable(cad(:), samplePoints{1});
            obj.aeroTable.caoff = DataTable(caoff(:), samplePoints{1});
            obj.aeroTable.clmq = DataTable(clmq(:), samplePoints{1});
            obj.aeroTable.clnr = DataTable(clnr(:), samplePoints{1});
            obj.aeroTable.cyp = DataTable(cyp, samplePoints{:});
            obj.aeroTable.cydr = DataTable(cydr, samplePoints{:});
            obj.aeroTable.cn0 = DataTable(cn0, samplePoints{:});
            obj.aeroTable.cnp = DataTable(cnp, samplePoints{:});
            obj.aeroTable.cndq = DataTable(cndq, samplePoints{:});
            obj.aeroTable.cllap = DataTable(cllap, samplePoints{:});
            obj.aeroTable.cllp = DataTable(cllp, samplePoints{:});
            obj.aeroTable.clldp = DataTable(clldp, samplePoints{:});
            obj.aeroTable.clm0 = DataTable(clm0, samplePoints{:});
            obj.aeroTable.clmdq = DataTable(clmdq, samplePoints{:});
            obj.aeroTable.clnp = DataTable(clnp, samplePoints{:});
            obj.aeroTable.clndr = DataTable(clndr, samplePoints{:});
            obj.aeroTable.clmp = DataTable(clmp, samplePoints{:});        
            
            obj.MachLim = [min(Mach); max(Mach)];
            obj.AlphaLim = [min(Alpha); max(Alpha)];
        end
        
        function obj = setParams(obj, biasRate, noiseStd)
            obj.biasRate = biasRate;
            obj.noiseStd = noiseStd;
        end
        function [F_aero,M_aero] = aeroForceAndMoment(obj, alpha_t, phi_t, pqr, delta, Vm, Xcg, mprop, mach, rho)
            r2d = 180/pi;
            biasRate_ = obj.biasRate;
            noiseStd_ = obj.noiseStd;
            delta = rad2deg(delta); % unit [rad] -> [deg]
            del_p = delta(1);
            del_q = delta(2);
            del_r = delta(3);
            
            p = pqr(1); % [rad/s]
            q = pqr(2); % [rad/s]
            r = pqr(3); % [rad/s]
            
            Q = 1/2 * rho * Vm^2;
            alpha_t = rad2deg(alpha_t); % [deg]
            phi_t = rad2deg(phi_t); % [deg]
            mach = min(max(mach, obj.MachLim(1)), obj.MachLim(2)); 
            alpha_t = min(max(alpha_t, obj.AlphaLim(1)), obj.AlphaLim(2)); 
            
            ca0 = obj.aeroTable.ca0.interpolate(mach);
            caa = obj.aeroTable.caa.interpolate(mach);
            cad = obj.aeroTable.cad.interpolate(mach);
            caoff = obj.aeroTable.caoff.interpolate(mach);
            clmq = obj.aeroTable.clmq.interpolate(mach);
            clnr = obj.aeroTable.clnr.interpolate(mach); 
            
            cyp = obj.aeroTable.cyp.interpolate(mach, alpha_t); 
            cydr = obj.aeroTable.cydr.interpolate(mach, alpha_t);
            cn0 = obj.aeroTable.cn0.interpolate(mach, alpha_t);
            cnp = obj.aeroTable.cnp.interpolate(mach, alpha_t);
            cndq = obj.aeroTable.cndq.interpolate(mach, alpha_t);
            cllap = obj.aeroTable.cllap.interpolate(mach, alpha_t);
            cllp = obj.aeroTable.cllp.interpolate(mach, alpha_t);
            clldp = obj.aeroTable.clldp.interpolate(mach, alpha_t);
            clm0 = obj.aeroTable.clm0.interpolate(mach, alpha_t); 
            clmdq = obj.aeroTable.clmdq.interpolate(mach, alpha_t);
            clnp = obj.aeroTable.clnp.interpolate(mach, alpha_t); 
            clndr = obj.aeroTable.clndr.interpolate(mach, alpha_t);
            clmp = obj.aeroTable.clmp.interpolate(mach, alpha_t); 
            
            %.. Angular Velocity & Control Deflection Transformation Body to Aeroballistic coordiate system
            
            qa              =       q     * cosd ( phi_t ) - r     * sind ( phi_t ) ;                                     % Pitch Angular Velocity in Aeroballistic Axes                                  (rad/s)
            ra              =       q     * sind ( phi_t ) + r     * cosd ( phi_t ) ;                                     % Yaw Angular Velocity in Aeroballistic Axes                                    (rad/s)
            del_qa          =       del_q * cosd ( phi_t ) - del_r * sind ( phi_t ) ;                                     % Pitch Control Deflection in Aeroballistic Axes                                (deg)
            del_ra          =       del_q * sind ( phi_t ) + del_r * cosd ( phi_t ) ;                                     % Yaw Control Deflection in Aeroballistic Axes                                  (deg)
     
            %.. Compute Force Coefficient
            
            deff            =       0.5 * ( abs( del_qa ) + abs( del_ra ) ) ;                                           % Effective Fin Deflection                                                      (Deg)
            ca              =       ca0 + caa * ( alpha_t ) + cad * deff^2 + (1 - mprop ) * caoff ;                     % Total Axial Force Coefficient                                                 (ND)
            cya             =       cyp * sind( 4 * phi_t ) + cydr * del_ra ;                                            % Total Side Force Coefficient in Aeroballistic Axes                            (ND)
            cna             =       cn0 + cnp * ( sind( 2 * phi_t ) )^2 + cndq * del_qa ;                                % Total Normal Force Coefficient in Aeroballistic Axes                          (ND)
            
            cza = -cna;
            
            %.. Compute Moment Coefficient
            
            cll             =       cllap * ( alpha_t )^2 * sind( 4 * phi_t ) + ...                                      % Total Rolling Moment Coefficient                                              (ND)
                                    cllp * ( p * r2d ) * obj.refLength / ( 2 * Vm ) + clldp * del_p ;
            clma            =       clm0 + clmp * ( sind( 2 * phi_t ) )^2 + ...                                          % Total Pitching Moment Coefficient in AeroBallistic Axes                       (ND)
                                    clmq * ( qa * r2d ) * obj.refLength / ( 2 * Vm ) + clmdq * del_qa - cna * ( obj.Xcg_ref - Xcg ) / obj.refLength ;
            clna            =       clnp * sind( 4 * phi_t ) + clnr * ( ra * r2d ) * obj.refLength / ( 2 * Vm ) + ...             % Total Yawing Moment Coefficient in AeroBallistic Axes                (ND)
                                    clndr * del_ra - cya * ( obj.Xcg_ref - Xcg ) / obj.refLength ;
            
            %.. Force & Moment Coefficient Transformation Aeroballistic to Body Coordinate system
            
            cy              =       cya  * cosd ( phi_t ) + cza  * sind ( phi_t ) ;                                       % Total Side Force Coefficient Body Axes                                        (ND)
            cz              =      -cya  * sind ( phi_t ) + cza  * cosd ( phi_t ) ;                                       % Total Normal Force Coefficient Body Axes                                      (ND)
            clm             =       clma * cosd ( phi_t ) + clna * sind ( phi_t ) ;                                       % Total Pitching Moment Coefficient Body Axes                                   (ND)
            cln             =      -clma * sind ( phi_t ) + clna * cosd ( phi_t ) ;                                       % Total Yawing Moment Coefficient Body Axes                                     (ND)
            
            %.. Aerodynamic bias and uncertainties
            ca              =       ca * (1 + biasRate_ + noiseStd_ * randn);
            cy              =       cy * (1 + biasRate_ + noiseStd_ * randn);
            cz              =       cz * (1 + biasRate_ + noiseStd_ * randn);
            cll             =       cll * (1 + biasRate_ + noiseStd_ * randn);
            clm             =       clm * (1 + biasRate_ + noiseStd_ * randn);
            cln             =       cln * (1 + biasRate_ + noiseStd_ * randn);

            %.. Output to Other Modules
            
            F_aero     =      Q * obj.refArea *[-ca; cy; cz];
            M_aero     =      Q * obj.refArea * obj.refLength * [cll; clm; cln];

        end

        function out = get_aerocoeff(obj, alpha, beta, mach, delta, t, airframe)
            r2d = 180/pi;
            veh = obj.vehicle;
            del_q = delta(2); % [rad]
            del_r = delta(3); % [rad]
            
            Xcg = airframe.getCG_data(t);
            mprop = veh.propulsion.isCombust(t);
            alpha_t = acos(cos(alpha)*cos(beta));
            phi_t = atan2d(tan(beta),sin(alpha));
            alpha_t = rad2deg(alpha_t); % [deg]
            mach = min(max(mach, obj.MachLim(1)), obj.MachLim(2)); 
            alpha_t = min(max(alpha_t, obj.AlphaLim(1)), obj.AlphaLim(2)); 
            Xcg_r = obj.Xcg_ref;
            d = obj.refLength;
            
            Ca_0 = obj.aeroTable.ca0.interpolate(mach);
            Ca_alp = obj.aeroTable.caa.interpolate(mach);
            Ca_deff = obj.aeroTable.cad.interpolate(mach);
            Ca_off = obj.aeroTable.caoff.interpolate(mach);
            CMa_q = obj.aeroTable.clmq.interpolate(mach) * r2d;
            CNa_r = obj.aeroTable.clnr.interpolate(mach) * r2d; 
            
            Cya_phi = obj.aeroTable.cyp.interpolate(mach, alpha_t); 
            Cya_dr = obj.aeroTable.cydr.interpolate(mach, alpha_t)*r2d;
            
            Cna_0 = obj.aeroTable.cn0.interpolate(mach, alpha_t);
            Cna_phi = obj.aeroTable.cnp.interpolate(mach, alpha_t);
            Cna_dq = obj.aeroTable.cndq.interpolate(mach, alpha_t)*r2d;

            CLa_phi = obj.aeroTable.cllap.interpolate(mach, alpha_t);
            CLa_p = obj.aeroTable.cllp.interpolate(mach, alpha_t) * r2d;
            CLa_dp = obj.aeroTable.clldp.interpolate(mach, alpha_t) * r2d;

            CMa_0 = obj.aeroTable.clm0.interpolate(mach, alpha_t); 
            CMa_phi = obj.aeroTable.clmp.interpolate(mach, alpha_t); 
            CMa_dq = obj.aeroTable.clmdq.interpolate(mach, alpha_t)*r2d;

            CNa_phi = obj.aeroTable.clnp.interpolate(mach, alpha_t); 
            CNa_dr = obj.aeroTable.clndr.interpolate(mach, alpha_t)*r2d;
            
            %.. Compute Coefficient
            del_qa          =       del_q * cosd ( phi_t ) - del_r * sind ( phi_t ) ;                                     % Pitch Control Deflection in Aeroballistic Axes                                (deg)
            del_ra          =       del_q * sind ( phi_t ) + del_r * cosd ( phi_t ) ;                                     % Yaw Control Deflection in Aeroballistic Axes                                  (deg)
            deff            =       0.5 * ( abs( del_qa * r2d ) + abs( del_ra * r2d) ) ;                                           % Effective Fin Deflection                                                      (Deg)
            Cna_0           =       Cna_0 + Cna_phi * ( sind( 2 * phi_t ) )^2;
            CMa_0 = CMa_0 + CMa_phi * ( sind( 2 * phi_t ) )^2 - Cna_0 * ( Xcg_r - Xcg ) / d;
            CNa_0 = CNa_phi * sind( 4 * phi_t ) - Cya_phi * sind( 4 * phi_t ) * ( Xcg_r - Xcg ) / d;
            
            Cx = -(Ca_0 + Ca_alp * alpha_t + Ca_deff * deff^2 + (1 - mprop ) * Ca_off);
            Cy = sind(phi_t) * (-Cna_0) + cosd(phi_t) * Cya_phi*sind(4*phi_t) ...
                 + (Cya_dr - Cna_dq)*cosd(phi_t)*sind(phi_t)*del_q + (Cna_dq*sind(phi_t)^2 + Cya_dr*cosd(phi_t)^2)*del_r;
            Cn = cosd(phi_t)*Cna_0  - sind(phi_t) * Cya_phi*sind(4*phi_t)...
                 - (Cna_dq*cosd(phi_t)^2 + Cya_dr*sind(phi_t)^2)*del_q - (Cya_dr - Cna_dq)*cosd(phi_t)*sind(phi_t)*del_r; 
            Cz = -Cn;
            CL_0 = CLa_phi * ( alpha_t )^2 * sind( 4 * phi_t );
            CL_p = CLa_p;
            CL_dp = CLa_dp;
            CL_dq = 0*r2d;
            CL_dr = 0*r2d;
            CM_0 = cosd(phi_t)*CMa_0 + sind(phi_t)*CNa_0;
            CM_q = CMa_q*cosd(phi_t)^2 + CNa_r*sind(phi_t)^2;
            CM_dq = ( CMa_dq - (Xcg_r -  Xcg)/d*Cna_dq )* cosd(phi_t)^2 + (CNa_dr-(Xcg_r -  Xcg)/d*Cya_dr)*sind(phi_t)^2;
            CN_0 = -sind(phi_t) * CMa_0 + cosd(phi_t)*CNa_0;
            CN_r = CMa_q*sind(phi_t)^2 + CNa_r*cosd(phi_t)^2; 
            CN_dr = ( CMa_dq - (Xcg_r - Xcg)/d*Cna_dq )*sind(phi_t)^2 + ( CNa_dr - (Xcg_r - Xcg)/d*Cya_dr)*cosd(phi_t)^2;
            out = [Cx, Cy, Cz, CL_0, CL_p, CM_0, CM_q, CN_0, CN_r, CL_dp, CL_dq, CL_dr, CM_dq, CN_dr];

        end

         function out = get_aerocoeff2(obj, t, alpha, beta, mach, delta, Xcg, d)
            r2d = 180/pi;
            del_q = delta(2); % [rad]
            del_r = delta(3); % [rad]
                        
            mprop = 1;
            alpha_t = acos(cos(alpha)*cos(beta));
            phi_t = atan2d(tan(beta),sin(alpha));
            alpha_t = rad2deg(alpha_t); % [deg]
            mach = min(max(mach, obj.MachLim(1)), obj.MachLim(2)); 
            alpha_t = min(max(alpha_t, obj.AlphaLim(1)), obj.AlphaLim(2)); 
            Xcg_r = obj.Xcg_ref;
            
            Ca_0 = obj.aeroTable.ca0.interpolate(mach);
            Ca_alp = obj.aeroTable.caa.interpolate(mach);
            Ca_deff = obj.aeroTable.cad.interpolate(mach);
            Ca_off = obj.aeroTable.caoff.interpolate(mach);
            CMa_q = obj.aeroTable.clmq.interpolate(mach) * r2d;
            CNa_r = obj.aeroTable.clnr.interpolate(mach) * r2d; 
            
            Cya_phi = obj.aeroTable.cyp.interpolate(mach, alpha_t); 
            Cya_dr = obj.aeroTable.cydr.interpolate(mach, alpha_t)*r2d;
            
            Cna_0 = obj.aeroTable.cn0.interpolate(mach, alpha_t);
            Cna_phi = obj.aeroTable.cnp.interpolate(mach, alpha_t);
            Cna_dq = obj.aeroTable.cndq.interpolate(mach, alpha_t)*r2d;

            CLa_phi = obj.aeroTable.cllap.interpolate(mach, alpha_t);
            CLa_p = obj.aeroTable.cllp.interpolate(mach, alpha_t) * r2d;
            CLa_dp = obj.aeroTable.clldp.interpolate(mach, alpha_t) * r2d;

            CMa_0 = obj.aeroTable.clm0.interpolate(mach, alpha_t); 
            CMa_phi = obj.aeroTable.clmp.interpolate(mach, alpha_t); 
            CMa_dq = obj.aeroTable.clmdq.interpolate(mach, alpha_t)*r2d;

            CNa_phi = obj.aeroTable.clnp.interpolate(mach, alpha_t); 
            CNa_dr = obj.aeroTable.clndr.interpolate(mach, alpha_t)*r2d;
            
            %.. Compute Coefficient
            del_qa          =       del_q * cosd ( phi_t ) - del_r * sind ( phi_t ) ;                                     % Pitch Control Deflection in Aeroballistic Axes                                (deg)
            del_ra          =       del_q * sind ( phi_t ) + del_r * cosd ( phi_t ) ;                                     % Yaw Control Deflection in Aeroballistic Axes                                  (deg)
            deff            =       0.5 * ( abs( del_qa * r2d ) + abs( del_ra * r2d) ) ;                                           % Effective Fin Deflection                                                      (Deg)
            Cna_0           =       Cna_0 + Cna_phi * ( sind( 2 * phi_t ) )^2;
            CMa_0 = CMa_0 + CMa_phi * ( sind( 2 * phi_t ) )^2 - Cna_0 * ( Xcg_r - Xcg ) / d;
            CNa_0 = CNa_phi * sind( 4 * phi_t ) - Cya_phi * sind( 4 * phi_t ) * ( Xcg_r - Xcg ) / d;
            
            Cx = -(Ca_0 + Ca_alp * alpha_t + Ca_deff * deff^2 + (1 - mprop ) * Ca_off);
            Cy = sind(phi_t) * (-Cna_0) + cosd(phi_t) * Cya_phi*sind(4*phi_t) ...
                 + (Cya_dr - Cna_dq)*cosd(phi_t)*sind(phi_t)*del_q + (Cna_dq*sind(phi_t)^2 + Cya_dr*cosd(phi_t)^2)*del_r;
            Cn = cosd(phi_t)*Cna_0  - sind(phi_t) * Cya_phi*sind(4*phi_t)...
                 - (Cna_dq*cosd(phi_t)^2 + Cya_dr*sind(phi_t)^2)*del_q - (Cya_dr - Cna_dq)*cosd(phi_t)*sind(phi_t)*del_r; 
            Cz = -Cn;
            CL_0 = CLa_phi * ( alpha_t )^2 * sind( 4 * phi_t );
            CL_p = CLa_p;
            CL_dp = CLa_dp;
            CL_dq = 0*r2d;
            CL_dr = 0*r2d;
            CM_0 = cosd(phi_t)*CMa_0 + sind(phi_t)*CNa_0;
            CM_q = CMa_q*cosd(phi_t)^2 + CNa_r*sind(phi_t)^2;
            CM_dq = ( CMa_dq - (Xcg_r -  Xcg)/d*Cna_dq )* cosd(phi_t)^2 + (CNa_dr-(Xcg_r -  Xcg)/d*Cya_dr)*sind(phi_t)^2;
            CN_0 = -sind(phi_t) * CMa_0 + cosd(phi_t)*CNa_0;
            CN_r = CMa_q*sind(phi_t)^2 + CNa_r*cosd(phi_t)^2; 
            CN_dr = ( CMa_dq - (Xcg_r - Xcg)/d*Cna_dq )*sind(phi_t)^2 + ( CNa_dr - (Xcg_r - Xcg)/d*Cya_dr)*cosd(phi_t)^2;
            out = [Cx, Cy, Cz, CL_0, CL_p, CM_0, CM_q, CN_0, CN_r, CL_dp, CL_dq, CL_dr, CM_dq, CN_dr];
         end

        
        function [alpha_t, phi_t] = TotalAoAandPhi(~, airVelocity)
            [ua, va, wa] = disperse(airVelocity);
            Va = vecnorm(airVelocity,2);
            alpha_t = acos(ua/Va); % Eq. 3.22
            phi_t = atan2(va, wa); % Eq. 3.23
        end

        function airVel = getAirVelocity(obj, h)
            if isempty(obj.WindShearFcn)                
                airVel = obj.vehicle.velBody;
                return
            end
            vel = obj.vehicle.velLocal(:);
            [windSpd, windDir] = obj.WindShearFcn(h);
            W_x = windSpd * sin(windDir);
            W_y = -windSpd * cos(windDir);

            airVel = vel + [W_x; W_y; 0];
            airVel = obj.vehicle.TLocal2Body * airVel;            
        end
        
        function airVel = anglesToAirvel(obj, alpha, beta, Vm)
            airVel = [Vm*cos(alpha)*cos(beta);
                      Vm*sin(beta);
                      Vm*sin(alpha)*cos(beta)];
        end
        function setWindShear(obj, seaVel, tropVel, seaDir, tropDir)             
            obj.WindShearFcn = @(h) windShear(h, seaVel, tropVel, seaDir, tropDir);
        end
    end
end