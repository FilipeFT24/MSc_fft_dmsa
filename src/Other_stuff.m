classdef Other_stuff
    methods (Static)
        %% > Set directories.
        function [] = Add_ABC()
            addpath(genpath('A'));
            addpath(genpath('B'));
            addpath(genpath('C'));
        end
        %% > Measure elapsed time.
        function [TA,TB] = Time_AB(inp,msh)
            %  > Class A.
            Class_A = @() A.WrapUp_A(); 
            TA      = timeit(Class_A);
            %  > Class B.
            Class_B = @() B.WrapUp_B(inp,msh); 
            TB      = timeit(Class_B);          
        end
    end
end


