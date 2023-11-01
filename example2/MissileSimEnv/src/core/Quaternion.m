classdef Quaternion < DynSystems
    properties
    end
    methods 
        function obj = Quaternion(name, q0, logOn)
            obj = obj@DynSystems(name, q0, logOn);
        end

        function dsdt = dynEqns(obj, ~, s, pqr)
            [p, q, r] = disperse(pqr);
            [q0, q1, q2, q3] = disperse(s);
            
            % Compute Quaternion Derivatives
            ck  = 400;                                    % Quaternion Orthogonalizing Factor
            erq = 1.0 - ( q0^2 + q1^2 + q2^2 + q3^2 ) ;   % Error of Quaternion Nonorthonomarlity 

            % "erq*ck{q}" is added to maintain the unit norm (See pp. 126, Zipfel "Modeling and Simulation of Aerospace Vehicle Dynamics")
            dq0dt = 0.5 * ( - p * q1 - q * q2 - r * q3 ) + ck * erq * q0 ;
            dq1dt = 0.5 * ( p * q0 + r * q2 - q * q3 )+ ck * erq * q1 ;
            dq2dt = 0.5 * ( q * q0 - r * q1 + p * q3 )+ ck * erq * q2 ;
            dq3dt = 0.5 * ( r * q0 + q * q1 - p * q2 )+ ck * erq * q3 ;
            dsdt = [dq0dt; dq1dt; dq2dt; dq3dt];
            if obj.logData               
            end
        end
    end

    methods (Static)
        function [Phi, Theta, Psi] = quat2euler(quat)
            [q0, q1, q2, q3] = disperse(quat);
            Phi = atan2(2*(q2.*q3 + q0.*q1), q0.^2 - q1.^2 -q2.^2 +q3.^2);
            Theta = asin(-2*(q1.*q3 - q0.*q2));
            Psi = atan2(2*(q1.*q2 + q0.*q3), q0.^2 + q1.^2 - q2.^2 - q3.^2);
        end

        function TBL = quat2dcm(quat)
            [q0, q1, q2, q3] = disperse(quat);
            %.. Compute DCM Matrix Using Quaternion: Quaternion -> Direction Cosine Matrix (Local coord. sys. to Body coord. sys.)            
            TBL(1,1)        =       q0^2 + q1^2 - q2^2 - q3^2 ; 
            TBL(1,2)        =       2 * ( q1 * q2 + q0 * q3 ) ;
            TBL(1,3)        =       2 * ( q1 * q3 - q0 * q2 ) ;
            
            TBL(2,1)        =       2 * ( q1 * q2 - q0 * q3 ) ;
            TBL(2,2)        =       q0^2 - q1^2 + q2^2 - q3^2 ;
            TBL(2,3)        =       2 * ( q2 * q3 + q0 * q1 ) ;
            
            TBL(3,1)        =       2 * ( q1 * q3 + q0 * q2 ) ;
            TBL(3,2)        =       2 * ( q2 * q3 - q0 * q1 ) ;
            TBL(3,3)        =       q0^2 - q1^2 - q2^2 + q3^2 ;
        end
    end
end