classdef D
    methods(Static)
        %% > Wrap-up D.
        function [] = WrapUp_D(save,load,ij)
            % >> Select...
            %  > Save (?).
            if save
                %  > Set directories...
                [Dir_1,Dir_2] = Data_Tools.Set_Directories(ij);
                %  > Save...
                h_max = 0.10;
                h_min = 0.01;
                N     = 20;
                h     = linspace(h_max,h_min,N);
                for ik = 1:N
                    %  > Compute...
                    [inp{ik},msh{ik}] = A.WrapUp_A(h(ik));
                    [pde{ik}]         = B.WrapUp_B(inp{ik},msh{ik});
                    %  > Save data...
                    np  {ik}          = inp{ik}.fr.np;
                    h_rf{ik}          = msh{ik}.d.h_ref;
                    Var {ik}          = pde{ik}.E.EN;
                    Data_Tools.SaveData_1(Dir_1,Dir_2,ik,np{ik},h_rf{ik},Var{ik});
                end
            end
            %  > Load (?).
            if load
                %  > Set directories...
                [Dir_1,Dir_2] = Data_Tools.Set_Directories(ij);
                %  > Load...
                Xi = Data_Tools.LoadData(Dir_1,Dir_2);
                Xo = D.Process_Data(Xi);
                %  > Remove duplicates.
                for i = 1:size(Xo.NR.H,2)
                    [~,~,rem_j{i}] = RunLength(Xo.NR.H{i});
                    [Xo.NR.H{i},Xo.NR.E{i}] = ...
                        deal(Xo.NR.H{i}(rem_j{i}),Xo.NR.E{i}(rem_j{i}));
                end
                %  > Convergence rate.
                for i = 1:size(Xo.NR.H,2)
                    if size(Xo.NR.E{i},2) ~= 1
                        [Xo.CR.H{i},Xo.CR.R{i}] = D.Compute_CR(Xo.NR.H{i},Xo.NR.E{i});
                    end
                end
                
                % >> Select...
                %  > Test #1.
                Test_1.WrapUp_T1(true,false,[1,2],Xo);
            end
        end
        
        %% > Auxiliary functions.
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [Xo] = Process_Data(Xi)
            %  > np.
            for i = 1:size(Xi,2)
                np(i) = Xi{1,i};
            end
            j    = unique(np); 
            Xo.n = j+1;
            
            %  > H,X.
            for i = 1:length(j)
                k{i} = find(np == j(i));
                for l = 1:size(k{i},2)
                    H{i}(l) = Xi{2,k{i}(l)};
                    E{i}{l} = Xi{3,k{i}(l)};
                end
                %  > Sort...
                [~,m{i}] = sort(H{i},'descend');
                for l = 1:size(k{i},2)
                    Xo.NR.H{i}(l) = H{i}(m{i}(l));
                    Xo.NR.E{i}{l} = [E{i}{m{i}(l)}{:}];
                end
            end
        end           
        % >> 1.2. ---------------------------------------------------------
        function [H_mean,CR] = Compute_CR(H,X)
            for i = 1:size(X,2)-1
                for j = 1:3
                    dE  (i,j) = log(X{i+1}(j)./X{i}(j));
                    CR  (i,j) = dE(i,j)./log(H(i+1)./H(i));
                    H_mean(i) = 1./2.*(H(i)+H(i+1));
                end
            end
        end
    end
end