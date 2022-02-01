classdef A_3_2
    methods (Static)
        %% > Wrap-up A_3_2.
        % >> --------------------------------------------------------------  
        % >> 1. Stencil setup.
        % >> 2. Stencil extension(s).
        % >> --------------------------------------------------------------
        function [msh] = WrapUp_A_3_2(msh,nt,nl,np_x,np_y,et)
            %  > Auxiliary arrays.
            i         = 1:size(msh.bnd.c,2);
            j         = 1:size(msh.bnd.f,2);
            bnd_cc(i) = [msh.bnd.c{2,i}];
            bnd_ff(j) = [msh.bnd.f{2,j}];
            bnd_fc(j) = [msh.bnd.f{3,j}];
               
            % >> 1.
            msh = A_3_2_1.Set_Stencil(msh,bnd_cc,bnd_ff,bnd_fc,nt,nl);           
            % >> 2.
            msh = A_3_2_2.Extend_Stencil(msh,bnd_cc,bnd_ff,bnd_fc,nt,np_x,np_y,et);
        end
    end
end