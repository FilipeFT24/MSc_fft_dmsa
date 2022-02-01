classdef Tools_1D
    methods (Static)
        %% > Set directories.
        function [] = Set_Directories()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('C_1D'));
            addpath(genpath('D_1D'));
            addpath(genpath('../[Tools]/[Tools - Data]'));
            addpath(genpath('../[Tools]/[Tools - Numerical]'));
            addpath(genpath('../[Tools]/[Tools - Post-processing]'));           
        end
    end
end


