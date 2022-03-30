classdef C_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Check_p(run,load_V)
            if run(1)
                C_1D.T2('1','C_1D/[.mat Files]/',load_V(1));
            end 
        end
               
        %% > 2. ----------------------------------------------------------- 
        % >> 2.1. ---------------------------------------------------------
        %  > Save .mat file.
        function [] = Save_mat(td,wd,V)
            save(join([wd,'V',td,'.mat']),'V');
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Set...
        %  > nc  : Number of cycles.
        %  > h(i): Reference length(s).
        function [msh] = Set_msh(nc,h)
            lin_h = linspace(log(h(1)),log(h(2)),nc);
            for i = 1:length(lin_h)
                msh(i) = A_1D.Set_A2(exp(1).^(lin_h(i)));
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Compute error slope.
        function [s] = Slope(h,e)
            [m,n]  = size(e);
            i      = 1:n;
            j      = 1:m-1;
            s(j,i) = log(e(j+1,i)./e(j,i))./log(h(j+1)./h(j)); 
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [] = T2(td,wd,load_V)
            if ~load_V
                %  > Set 'inp' and 'msh' structures.
                V.inp = A_1D.Set_A1;
                V.msh = C_1D.Set_msh(20,[0.5E-1,2.5E-3]);
                %  > Set up 'p-standard' run.
                if ~V.inp.pa.adapt
                    for i = 1:size(V.msh,2)
                        %  > Execute...
                        [V.obj(i),V.pde(i)] = B_1D.Initialize    (V.inp,V.msh(i));
                        [V.obj(i),V.msh(i)] = B_2_2_1D.p_standard(V.inp,V.obj(i),V.msh(i),V.pde(i));
                        %  > Print to terminal...
                        fprintf("Cycle #%d\n",i);
                    end
                end
            else
                load(['C_1D/[.mat Files]/V',td,'.mat']);
            end
            %  > Save structures...
            if ~load_V
                C_1D.Save_mat(td,wd,V);
            end
            %  > Plot...
            if V.inp.pl.all
                Fig_V2_1_1D.Plot(V);
            end
        end        
    end
end