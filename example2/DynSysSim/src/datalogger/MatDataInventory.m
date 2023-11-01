classdef MatDataInventory < handle
    properties
        data
        dataVarNum
        matData
        fieldNames
        indepVar
        indepVarName
    end
    methods
        function obj = MatDataInventory(fieldNames, colNum, rowNum, indepVar, indepVarName)
            if nargin <= 4
                indepVarName = 'Time [sec]';
            end
            varNum = numel(fieldNames);
            nanMat = nan(colNum, rowNum);
            fieldValPairs = cell(1, 2*varNum);
            for i = 1:varNum
                fieldValPairs(2*i-1) = fieldNames(i);
                fieldValPairs(2*i) = {nanMat};
            end
            obj.data = struct(fieldValPairs{:});
            obj.matData = nan(colNum, rowNum, varNum);
            obj.dataVarNum = varNum;
            obj.fieldNames = fieldNames;
            obj.indepVar = indepVar;
            obj.indepVarName = indepVarName;
        end
        
        function save(obj, fieldName, data, idxArray)
            obj.data.(fieldName)(idxArray) = data;                    
        end
        
        function mat2str(obj)
            for i = 1:obj.dataVarNum
                obj.data.(obj.fieldNames{i})= obj.matData(:, :, i);
            end
        end
        
        function plot(obj, fieldNames, fignum)
           if nargin <= 2
               figure()
           else
               figure(fignum);
           end
           plot(obj.indepVar, obj.data.(fieldNames))
           grid on;
           xlabel(obj.indepVarName);
           ylabel(fieldNames);
        end
    end
end