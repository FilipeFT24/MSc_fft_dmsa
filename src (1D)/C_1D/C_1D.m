classdef C_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Check_p()
            %  > Run/load...
            run    = [1,0];
            load_V = [1,1];
            %  > Boundary fix: E vs NNZ.
            if run(1)
                plot = [0,1];
                if plot(1)
                    C_1D.T1_1('1','C_1D/[.mat Files]/T1/',load_V(1));
                end
                if plot(2)
                    C_1D.T1_2('2','C_1D/[.mat Files]/T1/',load_V(1));
                end
            end 
            %  > Error estimator: E vs NNZ.
            if run(2)
                C_1D.T2('3','C_1D/[.mat Files]/T2/',load_V(2));
            end 
        end
               
        %% > 2. ----------------------------------------------------------- 
        % >> 2.1. ---------------------------------------------------------
        %  > Save .mat file.
        function [] = Save_mat(td,wd,V)
            save(join([wd,'V',td,'.mat']),'V');
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [msh_v] = inp_variables()
            %  > Number of cycles.
            n         = 5;
            %  > (h)_[min,max].
            x_min     = 0.5E-1;
            x_max     = 2.5E-3;
            lim       = [x_min,x_max];
            %  > Set...
            msh_v.n   = n;
            msh_v.lim = lim;
        end
        %  > 2.2.2. -------------------------------------------------------
        function [V] = SetUp_p_Decay(add,p,v)
            %  > Set 'inp' and 'msh' structures.
            msh_v        = C_1D.inp_variables;
            V.msh        = C_1D.Set_msh(msh_v.n,msh_v.lim);
            V.inp        = A_1D.Set_A2 (v);
            V.inp.ps.add = add;
            V.inp.ps.p   = p;

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
        end
        % >> 2.3. ---------------------------------------------------------
        function [msh] = Set_msh(nc,h)
            lin_h = linspace(log(h(1)),log(h(2)),nc);
            for i = 1:length(lin_h)
                msh(i) = A_1D.Set_A1(exp(1).^(lin_h(i)));
            end
        end
        % >> 2.4. ---------------------------------------------------------
        %  > Compute error slope.
        function [s] = Slope(h,e)
            [m,n]  = size(e);
            i      = 1:n;
            j      = 1:m-1;
            s(j,i) = log(e(j+1,i)./e(j,i))./log(h(j+1)./h(j)); 
        end
        % >> 2.5. ---------------------------------------------------------
        %  > Set regression model (fit).
        function [RM] = Set_RM(x,y)
            % >> From...
            %    https://www.mathworks.com/matlabcentral/answers/712943-help-finding-r-2-value-for-curves
            %  > y(fit) = a*x^b.
            y_fit   = @(k,x) k(1).*x.^k(2);
            %  > Auxiliary variables.
            a       = 1;
            b       = mean(C_1D.Slope(x,y));
            %  > SSECF.
            SSECF   = @(k) sum(((y_fit(k,x)-y)./x.^b).^2);
            %  > y(predicted).
            options = optimset('Display','off');
            K       = fminsearch(SSECF,[a;b],options);
            y_pred  = y_fit(K,x);
            %  > SSE.
            SSE     = sum(((y-y_pred)./x.^K(2)).^2);
            %  > SST.
            SST     = sum(((y-mean(y))./x.^K(2)).^2);
            %  > R^2.
            Rsq     = 1-SSE./SST;
            % >> Set structure...
            RM.K    = K;
            RM.r_sq = Rsq;
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        function [] = T1_1(td,wd,load_V)
            if ~load_V
                %  > Set 'V' structures.
                V(1) = C_1D.SetUp_p_Decay(0,[1,1],[0.5,0.1]);
                V(2) = C_1D.SetUp_p_Decay(1,[1,1],[0.5,0.1]);
                V(3) = C_1D.SetUp_p_Decay(0,[1,1],[0.5, 25]);
                V(4) = C_1D.SetUp_p_Decay(1,[1,1],[0.5, 25]);
                %  > Save...
                C_1D.Save_mat(td,wd,V);
            else
                load(['C_1D/[.mat Files]/T1/V',td,'.mat']);
            end
            %  > Save structures...
            if ~load_V
                C_1D.Save_mat(td,wd,V);
            end
            %  > Plot...
            Fig_V2_1_1D.Plot(V);
        end
        %  > 3.1.2. -------------------------------------------------------
        function [] = T1_2(td,wd,load_V)
            if ~load_V
                %  > Set 'V' structures.
                V(1)  = C_1D.SetUp_p_Decay(1,[1,1],[0.5,  1]);
                V(2)  = C_1D.SetUp_p_Decay(1,[3,3],[0.5,  1]);
                V(3)  = C_1D.SetUp_p_Decay(1,[5,5],[0.5,  1]);
                V(4)  = C_1D.SetUp_p_Decay(1,[7,7],[0.5,  1]);
                V(5)  = C_1D.SetUp_p_Decay(1,[9,9],[0.5,  1]);
                V(6)  = C_1D.SetUp_p_Decay(1,[1,1],[0.5,250]);
                V(7)  = C_1D.SetUp_p_Decay(1,[3,3],[0.5,250]);
                V(8)  = C_1D.SetUp_p_Decay(1,[5,5],[0.5,250]);
                V(9)  = C_1D.SetUp_p_Decay(1,[7,7],[0.5,250]);
                V(10) = C_1D.SetUp_p_Decay(1,[9,9],[0.5,250]);
                %  > Save structures...
                C_1D.Save_mat(td,wd,V);
            else
                load(['C_1D/[.mat Files]/T1/V',td,'.mat']);
            end
            %  > Plot...
            Fig_V2_2_1D.Plot(V);
        end
        % >> 3.2. ---------------------------------------------------------
        function [] = T2(td,wd,load_V)
            if ~load_V
                %  > Set 'V' structure.
                V = C_1D.SetUp_p_Decay(1,[1,1],[0.5,1]);
                %  > Save...
                C_1D.Save_mat(td,wd,V);
            else
                load(['C_1D/[.mat Files]/T2/V',td,'.mat']);
            end
            Fig_V2_3_1D.Plot(V);
        end
    end
end