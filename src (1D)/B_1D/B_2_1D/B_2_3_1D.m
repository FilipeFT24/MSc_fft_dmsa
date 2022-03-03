classdef B_2_3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = SetUp_ee(obj,msh,pde,s)
            %  > Test flags.
            flag_1 = 1;
            
            %  > #1.
            if flag_1
                %  > Set schemes order/type (v/g) / number of truncated terms.
                stl.t1 = ["UDS","CDS"];
                stl.t2 = [1,1];
                stl.nt = [2,2];
                %  > Compute truncated terms' magnitude (w/ analytic derivatives).
                B_2_3_1D.EE_1(obj,msh,pde,s,stl);
            end
            %  > #2.
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. 
        %  > Initialize 'stl' structure.
        function [stl] = SetUp_stl(msh,stl)
            %  > stl.p, stl.s, stl.t.
            stl_pst = A_2_1D.Initialize_stl(msh,stl.t1,stl.t2);
            stl.p   = stl_pst.p;
            stl.s   = stl_pst.s;
            stl.t   = stl_pst.t;
            %  > stl.o
            for i   = 1:length(stl.t1)
                stl.o(i) = A_2_1D.Compute_p(stl.t2(i),stl.t1(i));
            end
        end

        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Compute analytic derivatives up to order n.
        function [df] = Compute_dfa(msh,pde,stl)
            %  > Auxiliary variables.
            of = stl.o;
            nt = stl.nt;
            nf = length(of);
            Xv = msh.f.Xv;
            
            syms x;
            for i = 1:nf
                for j = 0:nt(i)-1
                    n      = of(i)+j;
                    dfn{i} = diff(pde.f.f{1},x,n);
                    dfn{i} = matlabFunction(dfn{i});
                    if nargin(dfn{i}) ~= 0
                        df{i}(:,j+1) = dfn{i}(Xv)./factorial(n);
                    else
                        df{i}(:,j+1) = zeros(length(Xv),1);
                    end
                end
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > Compute magnitude of truncated terms (w/ analytic derivatives).
        function [] = EE_1(obj,msh,pde,s,stl)
            %  > Initialize. 
            [stl]     = B_2_3_1D.SetUp_stl    (msh,stl);
            %  > 'p-standard' run.
            [s,stl,~] = B_2_1_1D.Update_s     (obj,msh,pde,s,stl);
            [x.f.a,~] = B_2_1_1D.Update_pde_x (msh,s,pde.a);
            [pde.e.t] = B_2_1_1D.Update_pde_et(obj,msh,s,pde.a,x.f.a); 
            %  > Compute truncated terms.
            [df]      = B_2_3_1D.Compute_dfa  (msh,pde,stl);
            [pde.tm]  = A_2_1D.Compute_tm     (obj,msh,s,stl,df);
            %  > Plot.
            Fig_2_1D.WrapUp_Fig_2_1_1D(msh,pde,stl.o-1);
        end

        %% > 4. -----------------------------------------------------------
        % >> 4.1. ---------------------------------------------------------
        %  > Compute derivatives w/ higher-order solution up to order n.
        function [dfn] = Compute_dfn_2(pde_ho,nt,p,s_ho)
            [m,n] = size(s_ho.Inv);
            for i = 1:m
                ki = p.lo(i)+1;
                kf = p.lo(i)+nt(i);
                for j = 1:n
                    for k = ki:kf
                        l           = k+1-ki;
                        df          = zeros(1,size(s_ho.Inv{i,j},1));
                        df    (1,k) = 1;
                        dfn{i}(j,l) = df*s_ho.Inv{i,j}*pde_ho.x.v.f{i,j};
                    end
                end
            end
        end
        
        
        
        
        
%         %  > Initialize 'stl' structure.
%         function [stl] = SetUp_stl(msh,flag)
%             switch flag
%                 case false
%                     %  > Schemes order/type (v/g) /number of truncated terms.
%                     t1     = ["CDS","CDS"];
%                     t2     = [2,2];
%                     %  > stl.
%                     stl    = A_2_1D.Initialize_stl(msh,t1,t2);
%                     stl.nt = [2,2];
%                     for i  = 1:length(t1)
%                         stl.o(i) = A_2_1D.Compute_p(t2(i),t1(i));
%                     end
%                 case true
%                     %  > Schemes order/type (v/g).
%                     t1{1} = ["UDS","CDS"]; % > LO.
%                     t1{2} = ["CDS","CDS"]; % > HO.
%                     t2{1} = [1,1];         % > LO.
%                     t2{2} = [1,2];         % > HO.
%                     %  > Number of truncated terms.
%                     ns    = size(t1,2);
%                     p_lo  = 1;
%                     p_ho  = ns;
%                     for i = 1:ns
%                         p.lo(i) = A_2_1D.Compute_p(t2{p_lo}(i),t1{p_lo}(i));
%                         p.ho(i) = A_2_1D.Compute_p(t2{p_ho}(i),t1{p_ho}(i));
%                     end
%                     nt    = p.ho-p.lo;
%                     %  > stl.
%                     for i = 1:length(t1{p_lo})
%                         [stl.lo.p{i}(:,1),stl.lo.s{i}(:,1),stl.lo.t{i}(:,1)] = ...
%                             A_2_1D.Initialize_stl(msh,t1{p_lo}(i),t2{p_lo}(i));
%                     end
%                     for i = 1:length(t1{p_ho})
%                         [stl.ho.p{i}(:,1),stl.ho.s{i}(:,1),stl.ho.t{i}(:,1)] = ...
%                             A_2_1D.Initialize_stl(msh,t1{p_ho}(i),t2{p_ho}(i));
%                     end
%                 otherwise
%                     return;
%             end
%         end
        %% > 2. -----------------------------------------------------------
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        %  > Truncated terms' magnitude (w/ higher-order solution).
        function [] = EE_2(obj,msh,pde,s)
            %  > Initialize.
            [stl,nt,p]           = B_2_3_1D.Initialize_stl(msh,true);
            %  > Compute tuncated term's magnitude.
            [pde_lo,s_lo,stl.lo] = B_2_1_1D.WrapUp_B_2_1_1D(obj,msh,pde,s,stl.lo);
            [pde_ho,s_ho,stl.ho] = B_2_1_1D.WrapUp_B_2_1_1D(obj,msh,pde,s,stl.ho);
            dfn_1                = B_2_3_1D.Compute_dfn_1  (msh,pde_lo,p.lo,nt);
            ttm_1                = B_2_3_1D.Compute_ttm    (obj,msh,s_lo,stl.lo.s,p.lo,nt,dfn_1);
            dfn_2                = B_2_3_1D.Compute_dfn_2  (pde_ho,nt,p,s_ho);
            ttm_2                = B_2_3_1D.Compute_ttm    (obj,msh,s_lo,stl.lo.s,p.lo,nt,dfn_2);
            %  > Plot.
            Fig_2_1D.WrapUp_Fig_2_2_1D(msh,pde_lo,p.lo-1,dfn_1,dfn_2,ttm_1,ttm_2);    
        end       
    end
end