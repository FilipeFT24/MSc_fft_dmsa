classdef C_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Check_P1(run,load_V)
            % >> Run(?).
            if run
                % >> Set working directories.
                td = '1';
                wd = 'C_1D/[.mat Files]/';
                % >> Load(?).
                if ~load_V
                    %  > Set 'inp' and 'msh' structures.
                    V.inp = A_1D.Set_A1;
                    V.msh = Tools_1D.Set_msh(25,[0.5E-1,2.5E-3]);
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
                    Tools_1D.Save_mat(td,wd,V);
                end
                %  > Plot...
                if V.inp.pl.all
                    Fig_V2_1_1D.Plot(V);
                end
            end
        end
    end
end