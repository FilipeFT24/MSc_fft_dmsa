classdef A_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        function [inp,msh] = Set_A(h)
            %  > Set 'inp' structure.
            inp  = A_1_1D.Set_inp(h);
            %  > Auxiliary variables.
            u    = inp.msh.Uniform;
            XLim = inp.msh.XLim;
            h    = inp.msh.h;
            NC   = round(1./h.*(XLim(2)-XLim(1)));
            
            switch u
                case true
                    Xv = A_2_1D.msh_1(XLim,NC);
                case false
                    Xv = A_2_1D.msh_2(XLim,NC,inp.msh.A,inp.msh.c);
                otherwise
                    return;
            end
            msh = A_2_1D.Set_Grid(XLim,Xv);
        end
    end
end