classdef Tools
    methods (Static)
        %% > Set directories.
        function [] = Set_Directories()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('C_1D'));
            addpath(genpath('D_1D'));
            addpath(genpath('[Tools - Post-processing]'));
            addpath(genpath('../[Tools - Data]'));
            addpath(genpath('../[Tools - Numerical]'));            
        end
    end
end


