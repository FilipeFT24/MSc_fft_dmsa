classdef C_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Check_p(load_V)
            % >> Set working directories.
            td = '1';
            wd = 'C_1D/[.mat Files]/';
            % >> Load(?).
            if ~load_V
                %  > Set 'inp' and 'msh' structures.
                V.inp = A_1D.Set_A1;
                V.msh = C_1D.Set_msh(25,[1E-1,2.5E-3]);
                %  > Set up 'p-standard' and 'p-adaptative' runs.
                switch V.inp.pa.adapt
                    case false
                        %  > 'p-standard' run.
                        for i = 1:size(V.msh,2)
                            %  > Execute...
                            [V.obj(i),V.pde(i)] = B_1D.Initialize    (V.inp,V.msh(i));
                            [V.obj(i),V.msh(i)] = B_2_2_1D.p_standard(V.inp,V.obj(i),V.msh(i),V.pde(i));
                            %  > Print to terminal...
                            fprintf("Cycle #%d\n",i);
                        end
                    otherwise
                        return;
                end
            else
                load(['C_1D/[.mat Files]/V',td,'.mat']);
            end
            %  > Save structures...
            if ~load_V
                C_1D.Save_struct(td,wd,V);
            end
            %  > Plot...
            if V.inp.pl.all
               Fig_V2_1_1D.Plot(V);
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > nc  : Number of cycles.
        %  > h(i): Reference length(s).
        function [msh] = Set_msh(nc,h)
            lin_h = linspace(log(h(1)),log(h(2)),nc);
            for i = 1:length(lin_h)
                msh(i) = A_1D.Set_A2(exp(1).^(lin_h(i)));
            end
        end
        % >> 1.3. ---------------------------------------------------------
        function [] = Save_struct(td,wd,V)
            save(join([wd,'V',td,'.mat']),'V');
        end
    end
end