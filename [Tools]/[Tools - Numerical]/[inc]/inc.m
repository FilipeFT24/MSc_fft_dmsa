classdef inc < handle
    properties
        v
        i
        k
    end
    methods
        function [obk] = inc(i,k)
            obk.v = i;
            obk.k = k; 
        end
        function [v] = get.i(obk)
            obk.v = obk.v+obk.k;
            v     = obk.v;
        end        
    end
end