classdef A_3_1_1D
    methods (Static)
        %% > Wrap-up A_3_1 (1D).
        function [msh]   = WrapUp_A_3_1_1D(msh)
            i            = 1:msh.c.NC;
            %  > Xc.
            msh.c.Xc (i) = 1./2.*(msh.f.Xv(i)+msh.f.Xv(i+1));
            %  > Cell volume.
            msh.c.Vol(i) = msh.f.Xv(i+1)-msh.f.Xv(i);
            %  > Reference length.
            msh.d.H_ref  = (msh.f.Xv(msh.f.NF)-msh.f.Xv(1))./msh.c.NC; 
        end
    end
end