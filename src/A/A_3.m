classdef A_3
    methods (Static)
        %% > Wrap-up A_3.
        % >> --------------------------------------------------------------
        % >> 1. Set 'struct' structure.
        % >> 2. Determine grid properties.
        % >> 3. Set up stencil.
        % >> --------------------------------------------------------------
        function [msh] = WrapUp_A_3(inp)
            % >> Local variables.
            nt = inp.fr.nt;
            p  = inp.fr.np;
            nl = 1./2.*(p+1);  
            et = inp.fr.et;
            
            % >> 1.
            [struct,msh] = A_2.WrapUp_A_2(inp);
            % >> 2.
            msh = A_3_1.WrapUp_A_3_1(struct,msh);
            % >> 3.
            msh = A_3_2.WrapUp_A_3_2(msh,nt,nl,p,et);
        end
    end
end