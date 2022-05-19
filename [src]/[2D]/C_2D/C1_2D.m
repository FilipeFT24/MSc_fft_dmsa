classdef C1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [V] = Run_p(run,T)
            if ~run(1)
                load(strjoin(["[.mat Files]/V",T,".mat"],''));
            else
                V = C1_2D.Check_Decay(T);
            end
            Fig_2D_3.Plot(V.msh,V.obj);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [V] = Check_Decay(T)
            %  > "inp".
            h_lim    = [5.0E-2,2.0E-2];
            n        = 5;
            h        = exp(1).^(linspace(log(h_lim(1)),log(h_lim(2)),n));
            inp      = A1_2D.Set_inp([1,1],[100,0.5,0.5]);     %  > c/f.
            inp.plot = [0,0];
            
            %  > Set up "P-standard" run.
            if ~inp.p_adapt.allow
                for i = 1:numel(h)
                    %  > Execute...
                    msh  (i)   = A2_2D.Set_msh(h(i));       %  > h.
                    obj  (i)   = B3_2D.Run_p  (inp,msh(i)); %  > inp/msh.
                    %  > Assign...
                    V.msh(i).d = msh(i).d;
                    V.obj(i).e = obj(i).e;
                    for j = 1:numel(obj(i).m)
                        V.obj(i).m{j}.nnz = obj(i).m{j}.nnz;
                    end
                    %  > Print to terminal...
                    fprintf("Cycle #%d\n",i);
                end
            end
            Tools_1.Save_mat(T,"C_2D/[.mat Files]/",V);
        end
    end
end