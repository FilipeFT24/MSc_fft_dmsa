classdef B_2_3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        %  > Initialize 'stl' structure.
        function [stl,nt,p] = Initialize_stl(msh,flag)
            switch flag
                case false
                    %  > Schemes order/type (v/g).
                    t1 = ["CDS","CDS"];
                    t2 = [2,2];
                    j  = 1:length(t1);
                    %  > Number of truncated terms.
                    nt = [2,2];
                    %  > stl.
                    for i = j
                        [stl.p{i}(:,1),stl.s{i}(:,1),stl.t{i}(:,1)] = ...
                            A_2_1D.Initialize_stl(msh,t1(i),t2(i));
                        p(i) = A_2_1D.Compute_p(t2(i),t1(i));
                    end
                case true
                    %  > Schemes order/type (v/g).
                    t1{1} = ["UDS","CDS"]; % > LO.
                    t1{2} = ["CDS","CDS"]; % > HO.
                    t2{1} = [1,1];         % > LO.
                    t2{2} = [1,2];         % > HO.
                    %  > Number of truncated terms.
                    ns    = size(t1,2);
                    p_lo  = 1;
                    p_ho  = ns;
                    for i = 1:ns
                        p.lo(i) = A_2_1D.Compute_p(t2{p_lo}(i),t1{p_lo}(i));
                        p.ho(i) = A_2_1D.Compute_p(t2{p_ho}(i),t1{p_ho}(i));
                    end
                    nt    = p.ho-p.lo;
                    %  > stl.
                    for i = 1:length(t1{p_lo})
                        [stl.lo.p{i}(:,1),stl.lo.s{i}(:,1),stl.lo.t{i}(:,1)] = ...
                            A_2_1D.Initialize_stl(msh,t1{p_lo}(i),t2{p_lo}(i));
                    end
                    for i = 1:length(t1{p_ho})
                        [stl.ho.p{i}(:,1),stl.ho.s{i}(:,1),stl.ho.t{i}(:,1)] = ...
                            A_2_1D.Initialize_stl(msh,t1{p_ho}(i),t2{p_ho}(i));
                    end
                otherwise
                    return;
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        %  > Compute analytic derivatives up to order n.
        function [dfn] = Compute_dfn_1(msh,pde,p,nt)
            syms x;
            for i = 1:length(p)
                for j = 0:nt(i)-1
                    n     = p(i)+j;
                    df{i} = diff(pde.f.f{1},x,n);
                    df{i} = matlabFunction(df{i});
                    if nargin(df{i}) ~= 0
                        dfn{i}(:,j+1) = df{i}(msh.f.Xv)./factorial(n);
                    else
                        dfn{i}(:,j+1) = zeros(length(msh.f.Xv),1);
                    end
                end
            end
        end
        %  > 2.1.2. -------------------------------------------------------
        %  > Truncated terms' magnitude (w/ analytic derivatives).
        function [] = EE_1(obj,msh,pde,s)
            %  > Initialize.
            [stl,nt,p]  = B_2_3_1D.Initialize_stl (msh,false);
            %  > Compute tuncated term's magnitude.
            [pde,s,stl] = B_2_1_1D.WrapUp_B_2_1_1D(obj,msh,pde,s,stl);
            dfn         = B_2_3_1D.Compute_dfn_1  (msh,pde,p,nt);
            ttm         = B_2_3_1D.Compute_ttm    (obj,msh,s,stl.s,p,nt,dfn);
            %  > Plot.
            Fig_2_1D.WrapUp_Fig_2_1_1D(msh,pde,p-1,ttm);
        end
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
        %  > 2.2.2. -------------------------------------------------------
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
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [] = EE_3(obj,msh,pde,s)
            
            ll = 1;
            
        end
    end
end