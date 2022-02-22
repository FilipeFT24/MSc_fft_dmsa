classdef B_2_3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [msh,pde] = Check_TE(msh,pde,s,A,B,v,g,bnd)
            %  > Schemes order/type.
            t1(1)   = "CDS";
            t1(2)   = "CDS";
            t2(1)   = 1;
            t2(2)   = 1;
            %  > Number of truncation error terms.
            nt      = 5;
            %  > stl.
            [stl,p] = B_2_3_1D.Initialize_stl(msh.f.NF,t1,t2);
            
            %  > Evaluate analytic derivatives at cell faces.
            for i = 1:length(p)
                for j = 1:nt
                    dfn{i}(:,j) = B_2_3_1D.Compute_dfn(p(i)-1+j,pde.f.f{1},msh.f.Xv);
                end
            end
            %  > Solve PDE.
            [~,s,x,ea] = B_2_1_1D.Update_stl(msh,stl,s,A,B,pde.a,bnd,v,g);
            [e,~]      = B_2_1_1D.Update_pde(msh,pde.a,s,x,ea,v,g);
            
            %  > Compute truncation error (finite number of terms: nt).
            TE = B_2_3_1D.Compute_TE(msh,stl.s,s,p,nt,dfn);
            
            hold on;
            plot(msh.f.Xv,abs(e.t.f(:,1)),'or');
            plot(msh.f.Xv,abs(TE{1})); set(gca,'YScale','log');
            
        end
        % >> 1.2. ---------------------------------------------------------
        function [stl,n] = Initialize_stl(NF,t1,t2)
            for i = 1:length(t1)
                stl.p{i}(:,1) = repelem(t2(i),NF);
                stl.s{i}(:,1) = 1:NF;
                stl.t{i}(:,1) = repelem(t1(i),NF);
                n    (i)      = B_2_2_1D.Compute_p(t2(i),t1(i));
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set symbolic functions/handles and compute analytic values of derivative of order n.
        function [dfn_v] = Compute_dfn(p,func,Xv)
            syms x;
            dfn = diff(func,x,p);
            dfn = matlabFunction(dfn);
            if nargin(dfn) ~= 0
                dfn_v(:,1) = dfn(Xv)./factorial(p);
            else
                dfn_v(:,1) = zeros(length(Xv),1);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Compute truncation error terms' magnitude (truncated terms of Df).
        function [TE] = Compute_TE(msh,stl_s,s,p,nt,dfn)
            %  > Auxiliary variables.
            Xv   = msh.f.Xv;
            trsh = 10e-12;
            
            % >> Loop through selected faces...
            for i = 1:size(stl_s,2)
                if isempty(stl_s{i})
                    continue;
                else
                    for j = 1:size(stl_s{i},1)
                        %  > Auxiliary variables.
                        k   = stl_s{i}(j);
                        xt  = s.xt{i,k};
                        len = length(xt);
                        l   = 1:nt;
                        m   = 1:len;
                        fx  = Xv(k);
                        
                        %  > Df_T.
                        Df_T (l,m) = transpose((xt(m)-fx)'.^(p(i)-1+l));
                        %  > Truncation error (TE).
                        TE{i}(k,l) = transpose(Df_T*s.xf{i,k}').*dfn{i}(k,l);
                    end
                end
                %  > Set elements below a given treshold to 0.
                TE{i}(abs(TE{i})<trsh) = 0;
            end
        end
    end
end