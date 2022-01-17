classdef Tools
    methods (Static)
        %% > Set directories.
        function [] = Set_Directory(str)
            addpath(genpath(str));
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


