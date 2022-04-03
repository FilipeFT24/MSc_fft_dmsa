classdef A_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [msh] = Set_A1(h)
            %  > Auxiliary variables.
            inp  = A_1_1D.Set_inp_1(h);
            u    = inp.Uniform;
            XLim = inp.XLim;
            h    = inp.h;
            NC   = round(1./h.*(XLim(2)-XLim(1)));
            
            switch u
                case true
                    Xv = A_2_1D.msh_1(XLim,NC);
                case false
                    Xv = A_2_1D.msh_2(XLim,NC,inp.A,inp.c);
                otherwise
                    return;
            end
            msh     = A_2_1D.Set_Grid(XLim,Xv);
            msh.d.h = h;
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set 'inp' structure.
        function [inp] = Set_A2(v)
            inp = A_1_1D.Set_inp_2(v);
        end 
    end
end