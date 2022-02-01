classdef Tools
    methods (Static)
        %% > Set directories.
        function [] = Set_Directories()
            addpath(genpath('A'));
            addpath(genpath('B'));
            addpath(genpath('C'));
            addpath(genpath('D'));
            addpath(genpath('../[Tools - Data]'));
            addpath(genpath('../[Tools - Numerical]'));
            addpath(genpath('[Tools - Post-processing]'));
        end
        %% > Measure elapsed time.
        function [T] = Time_AB(inp,msh)
            %  > Class A.
            Class_A = @() A.WrapUp_A(); 
            TA      = timeit(Class_A);
            %  > Class B.
            Class_B = @() B.WrapUp_B(inp,msh); 
            TB      = timeit(Class_B);    
            %  > Elpased time.
            T = TA+TB;
        end
    end
end


