classdef Data_Tools
    methods(Static)
        %% > Tools
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [Dir_1,Dir_2] = Set_Directories(ij)
            % >> Set...
            %  > ...destination directory.
            Dir_1 = convertCharsToStrings(['../[DataFiles]/Test_',ij]);
            if ~exist(Dir_1,'dir')
                mkdir(Dir_1)
            end
            %  > ...main directory.
            Dir_2 = '../../src';
        end
        % >> 1.2. ---------------------------------------------------------
        function [] = SaveData_1(Dir_1,Dir_2,ik,np,h,Variable)
            %  > Change directory.
            cd(Dir_1);
            %  > Save data and print to terminal.
            save(['T','_',sprintf('%.0f',np),'_',sprintf('%.0f',ik)],'np','h','Variable');
            fprintf('Loop #%d complete.\n',ik);
            %  > Change directory.
            cd(Dir_2);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [X] = LoadData(Dir_1,Dir_2)
            %  > Change directory.
            cd(Dir_1);
            %  > Load data...
            Data = dir('*.mat');
            for i = 1:length(Data)
                X_Mat(i) = load(Data(i).name);
                X  {1,i} = X_Mat(i).np;
                X  {2,i} = X_Mat(i).h;
                X  {3,i} = X_Mat(i).Variable;
            end
            %  > Change directory.
            cd(Dir_2);
        end
    end
end