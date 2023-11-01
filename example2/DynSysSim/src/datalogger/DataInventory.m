classdef DataInventory < handle
    properties
        fieldName
        dataName        
        data
        dataNum
        dataLen        
        lastAppenedIdx = 0
        indepVar
    end
    
    properties (Hidden)
        interval
        startVal
    end
    methods
        function [obj, str] = DataInventory(fieldNames, dataNames, indepVarSpan)           
            if nargin < 1
                return
            end
            if ischar(fieldNames)
                    fieldNames = {fieldNames};
            end
            numFieldNames = numel(fieldNames);
            if numFieldNames == 1
                obj.fieldName = fieldNames;
                obj.dataLen = numel(indepVarSpan);
                obj.startVal = indepVarSpan(1);
                obj.interval = diff(indepVarSpan(1:2));
                obj.dataNum = length(dataNames);
                obj.dataName = dataNames;
                obj.data = nan(length(indepVarSpan), obj.dataNum);
                obj.indepVar = indepVarSpan;
                str.(fieldNames{1}) = obj;
            else                
                objs = cell(1, numFieldNames);
                for i = 1:numFieldNames
                    objs{i} = DataInventory(fieldNames(i), dataNames{i}, indepVarSpan);
                end
                arguments = [fieldNames(:).';objs];
                arguments = arguments(:).';
                str = struct(arguments{:});
            end
        end
             
        function append(obj, data)
            startIdx = obj.lastAppenedIdx + 1;           
            obj.data(startIdx, :) = data;
            obj.lastAppenedIdx = startIdx;
            if obj.lastAppenedIdx > obj.dataLen
                obj.data(obj.lastAppenedIdx+1:2*obj.lastAppenedIdx, :) = nan;
                obj.dataLen = 2*obj.lastAppenedIdx;
            end
        end
        
        function p = subplots(obj, fignum, Title, idx)
            switch nargin
                case 1
                    figure();
                    idx = 1:1:obj.dataNum;
                case 2
                    figure(fignum);
                    idx = 1:1:obj.dataNum;
                case 3
                    fig = figure(fignum);
                    fig.Name = Title;
                    idx = 1:1:obj.dataNum;
                case 4
                    fig = figure(fignum);
                    fig.Name = Title;
            end
            obj.postProcessData;
            p = cell(length(idx) , 1); % preallocation
            for i = 1:1:length(idx)
                k = idx(i);
                subplot(length(idx), 1, i); hold on;
                if i == 1 && nargin > 2
                    title(Title);
                end
                if isempty(obj.dataName{k})
                    obj.dataName{k} = '';
                end
                p{i} = plot(obj.indepVar, obj.data(:, k),'DisplayName',obj.dataName{k});
                ylabel(obj.dataName{k});
                box on; grid on;
            end      
             xlabel('Time [sec]');
        end
        
        function p = plot(obj, fignum, Title, idx)
            switch nargin
                case 1
                    figure();
                    idx = 1:1:obj.dataNum;
                    Title = '';
                case 2
                    figure(fignum);
                    Title = '';
                    idx = 1:1:obj.dataNum;
                case 3
                    fig = figure(fignum);
                    fig.Name = Title;
                    idx = 1:1:obj.dataNum;
                case 4
                    fig = figure(fignum);
                    fig.Name = Title;
            end
            obj.postProcessData;
            p = cell(length(idx) , 1); % preallocation
            hold on;
            title(Title);
            xlabel('Time [sec]');box on; grid on;
            for i = 1:length(idx)                
                p{i} = plot(obj.indepVar, obj.data(:, idx(i)),'DisplayName', obj.dataName{idx(i)});
            end           
        end
        
        function plot2(obj, fignum, i, j)
            if nargin < 2
                figure('Color', 'white');
            else
                fig = figure(fignum);
                fig.Color = 'white';
            end
            if nargin < 3
                i = 1;
                j = 2;
            end
            obj.postProcessData;
            plot(obj.data(:,i), obj.data(:,j))
            xlabel(obj.dataName{i});
            ylabel(obj.dataName{j});
            grid on; box on; 
        end
        
        function plot3(obj, fignum, i, j, k)
            if nargin < 2
                figure('Color', 'white');
            else
                fig = figure(fignum);
                fig.Color = 'white';
            end
            if nargin < 3
                i = 1;
                j = 2;
                k = 3;
            end
            obj.postProcessData;
            plot3(obj.data(:,i), obj.data(:,j), obj.data(:,k))
            xlabel(obj.dataName{i});
            ylabel(obj.dataName{j});
            zlabel(obj.dataName{k})
            grid on; box on; 
        end
        
        function erase(obj, idx) 
            if isnan(idx)
                obj.data(obj.lastAppenedIdx+1:end, :) = [];
                return
            end
            obj.indepVar(idx) = [];
            obj.data(idx, :) = []; 
            obj.dataLen = size(obj.data, 1);
        end           
        
        function str = saveobj(obj)
            str.fieldName = obj.fieldName;
            str.dataName = obj.dataName;
            str.indepVar = obj.indepVar;
            str.data = obj.data;
            str.dataNum = obj.dataNum;
        end
        
        function postProcessData(obj)
            len = obj.lastAppenedIdx;
            obj.data(len+1:end, :) = [];
        end
    end
    
    methods (Static)
        function test()
            timeSpan = 0:0.1:10;
            % singular data case
            fieldNames = {'aa'};
            dataNames = {'a1','a2'};
            aa = DataInventory(fieldNames, dataNames, timeSpan);
            
            for i = 1:length(timeSpan)
                adata = [sin(timeSpan(i)); cos(timeSpan(i))];               
                aa.save(adata, i);
            end
            aa.multiPlot(1, 'a'); hold on;
            
            %% multiple data case
            fieldNames = {'aa','bb','cc'};
            dataNames = {{'a1','a2'},{'b1','b2'},{'c1','c2'}};
            [~, log] = DataInventory(fieldNames, dataNames, timeSpan);
            
            for i = 1:length(timeSpan)
                adata = [sin(timeSpan(i)); cos(timeSpan(i))];
                bdata = 2*adata;
                cdata = 3*adata;
                log.aa.save(adata, i);
                log.bb.save(bdata, i);
                log.cc.save(cdata, i);
            end
           
            log.aa.multiPlot(2, 'a'); hold on;
            log.bb.multiPlot(2, 'b');
            log.cc.multiPlot(2, 'c')
            log.aa.plot2(3);
        end
        
        function str = obj2str(datastr)
            fieldNames = fieldnames(datastr);
            for i = 1:length(fieldNames)
                object = datastr.(fieldNames{i});
                dataStr = object.saveobj;
                str.(fieldNames{i}) =  dataStr;
            end            
        end
        function strs = correctName(str)
            fieldNames = fieldnames(str);
            nfield = length(fieldNames);
            for i = 1:nfield
                fieldName = fieldNames{i}; 
                strs.(fieldName).fieldName = fieldName;
                strs.(fieldName).dataName = str.(fieldName).dataName;
                strs.(fieldName).indepVar = str.(fieldName).indepVar;
                strs.(fieldName).data = str.(fieldName).data;
                strs.(fieldName).dataNum = str.(fieldName).dataLen;
            end
        end
        function objs = str2obj(str)
            fieldNames = fieldnames(str);
            nfield = length(fieldNames);
            for i  = 1:nfield
                fieldName = fieldNames{i};
                objs.(fieldName) = DataInventory;
                objs.(fieldName).dataName = str.(fieldName).dataName;
                objs.(fieldName).indepVar = str.(fieldName).indepVar;
                objs.(fieldName).data = str.(fieldName).data;
                objs.(fieldName).dataNum = str.(fieldName).dataNum;
            end
            
        end   
             
        function str = addData(obj1, obj2, operator, fieldName, dataName)
            assert(obj1.dataNum == obj2.dataNum, 'the length of entered two data are not equal');
            if nargin == 3
                fieldName = 'addedData';
                dataName = cell(obj1.dataNum, 1);
            end
            if nargin == 4
                dataName = cell(obj1.dataNum, 1);
            end
            
            
            [~, str] = DataInventory(fieldName, dataName, obj1.indepVar);
            switch operator
                case '+'
                    str.(fieldName).append(obj1.data + obj2.data);
                case '-'
                    str.(fieldName).append(obj1.data - obj2.data);
            end
            
        end
        
    end
end