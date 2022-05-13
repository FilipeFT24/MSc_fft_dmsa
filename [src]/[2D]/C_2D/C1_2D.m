classdef C1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [V] = Run_p(run)
            %  > Load/run...
            if ~run(1)
                load('[.mat Files]/V3.mat');
            else
                V = C1_2D.Check_Decay;
            end
            Fig_2D_3.Plot(V.msh,V.obj);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [V] = Check_Decay()
            %  > "inp".
            h_lim = [5.0E-2,2.5E-2];
            n     = 3;
            h     = exp(1).^(linspace(log(h_lim(1)),log(h_lim(2)),n));
            inp   = A1_2D.Set_inp_2(1,[0.5,0.5,100]);   %  > f_type/xc/yc/i.
            
            %  > Set up "P-standard" run.
            if ~inp.p_adapt.allow
                for i = 1:numel(h)
                    %  > Execute...
                    msh(i) = A2_2D.Set_msh(h(i));       %  > h.
                    obj(i) = B3_2D.Run_p  (inp,msh(i)); %  > inp/msh.
                    %  > Print to terminal...
                    fprintf("Cycle #%d\n",i);
                end
            end
            V.msh = msh;
            V.obj = obj;
            Tools_1.Save_mat('3','C_2D/[.mat Files]/',V);
        end
    end
end