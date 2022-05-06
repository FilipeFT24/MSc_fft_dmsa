classdef B1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "s" (stencil cell/face indices, coordinates, etc.).
        function [s] = Initialize_1(nc,ns,Nf)
            for i = 1:ns
                s{i}.i       = cell(Nf,nc);
                s{i}.logical = cell(Nf,nc);
                s{i}.xt      = cell(Nf,nc);
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil coordinates, etc.
        function [s] = Update_1(msh,s,u)
            s        = B1_2D.Update_sc(msh,s,u);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update stencil coordinates, etc.
        function [s] = Update_sc(msh,s,u)
            for i = 1:size(u.s,2)
                if isempty(u.s{i})
                    continue;
                else
                    for j = u.s{i}'
                        p = u.p(j,i);
                        
                        %  > Face "j" vertex indices.
                        fv = msh.f.iv(j,:);
                        %  > To what cells do these vertices belong to?
                        cv = unique(cat(1,msh.v.ic{fv}));
                        
                        %  > Are these boundary cells?
                        bc = msh.c.logical(cv);
                        %  > Add boundary faces...
                        fc = unique(reshape(msh.c.f.if(cv(bc),:),1,[]));
                        bf = fc(~msh.f.logical(fc))';
                        %  > Check whether the selected boundary faces contain any of "fv"'s vertices.
                        if ~isempty(bf)
                            %  > Initialize array "l".
                            l = false(numel(bf),numel(fv));
                            %  > For each of the vertices in "fv"...
                            for k = 1:numel(fv)
                                l(:,k) = any(bsxfun(@eq,msh.f.iv(bf,:),fv(k)),2);
                            end
                            bf = bf(any(l,2));
                            %  > Remove face that shares 1 vertex w/ face "j" if it belongs to the same cell.
                            if any(bf == j)
                                fv_c             = unique (cat(1,msh.v.if{fv}));
                                b                = ismembc(bf,fv_c);
                                b   (bf == j)    = false;
                                %  > Check whether face "j" belongs to the same cell...
                                c_b              = cat(1,msh.f.ic{bf});
                                c_j              = msh.f.ic{j};
                                b   (c_b ~= c_j) = false;
                                %  > Remove...
                                if nnz(b) == 1
                                    bf(b) = [];
                                end
                            end
                        end
                        
                        %  > Compute coordinates...
                        xc = B1_2D.xt_c(msh.c.c.xy.c,cv);
                        if ~isempty(bf)
                            %  > xt.
                            xf      = B1_2D.xt_f(msh.f.xy.c,bf);
                            xt      = cat(1,xc,xf);
                            %  > logical: 0-bnd.
                            %             1-blk.
                            logical = cat(1,true(numel(cv),1),false(numel(bf),1));
                        else
                            %  > xt.
                            xt      = xc;
                            %  > logical: 1-blk.
                            logical = true(numel(cv),1);
                        end
                        %  > Update field "s".
                        if ~isempty(bf)
                            s.i  {j,i} = cat(1,cv,bf);
                        else
                            s.i  {j,i} = cv;
                        end
                        s.logical{j,i} = logical;
                        s.xt     {j,i} = xt;
                    end
                end
            end
        end
        %  > 2.1.1. -------------------------------------------------------
        %  > Compute stencil (adimensional) parameters and verify the need for stencil extension(s).
        function [f] = Extend(h,p,sc,xy)
            %  > Auxiliary variables.
            sz_s = size(xy);
            
            %  > Stencil limits/adimensional parameters(n) in the x/y-direction(s).
            for i = 1:sz_s(2)
                [L(1,i),L(2,i)] = MinMaxElem(xy(:,i),'finite');
            end
            Lt = diff(L,1);
            n  = ceil(Lt./Tools_1.mean(h(sc,:)));
            %  > Extend(?).
            f  = ~(n >= p);
        end
        %  > 2.1.2. -------------------------------------------------------
        %  > Compute stencil coordinates.
        %  > 2.1.2.1. -----------------------------------------------------
        %  > c.
        function [xt] = xt_c(xc,c)
            xt = xc(c,:);
        end
        %  > 2.1.2.2. -----------------------------------------------------
        %  > f.
        function [xt] = xt_f(xf,f)
            xt = xf(f,:);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize field "x" (stencil coefficients/nodal solution,etc.).
        function [x] = Initialize_24(f,u,nc,ns,Nc,Nf)
            for i = 1:ns
                %  > "sx".
                x{i}.Tf   = cell(Nf,nc);
                x{i}.Tf_V = cell(Nf,nc);
                x{i}.Pf   = cell(Nf,nc);
                %  > "x".
                x{i}.nv.a.f = zeros(Nf,1);
                for j = u{i}.f
                    x{i}.cf.(j) = cell(Nf,nc);
                    x{i}.vf.(j) = cell(Nf,nc);
                    x{i}.xf.(j) = cell(1 ,nc);
                    for k = 1:nc
                        switch k
                            case 1
                                x{i}.xf.(j){k} = zeros(Nf,1); %  > \phi_f.
                            case 2
                                x{i}.xf.(j){k} = zeros(Nf,2); %  > \nabla\phi_f(x,y).
                        end
                    end
                    if j == "a"
                        x{i}.nv.(j)   = f.av;
                    else
                        x{i}.nv.(j).c = zeros(Nc,1);
                    end
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        %  > Update stencil coefficients, etc.
        function [x] = Update_2(inp,msh,f,s,u,x)
            x        = B1_2D.Update_sx(inp,msh,f,s,u,x);
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Update nodal/face values, etc.
        function [x] = Update_4(f,s,u,x)
            x        = Tools_2.Update_xv(f,s,u,x);
            x        = Tools_2.Update_cf(u,x);
            x        = Tools_2.Update_xf(u,x);
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Update stencil coefficients.
        function [x] = Update_sx(inp,msh,f,s,u,x)
            %  > Auxiliary variables.
            k = ["v","g"];
            
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        %  > Polynomial regression coefficients/exponents.
                        [c,e]       = Tools_2.Terms(inp.p.p(i),0);
                        %  > Df.
                        xf_c        = msh.f.xy.c(j,:);
                        xt          = s.xt{j,i};
                        Df          = B1_2D.Assemble_Df(c.v,e.v,xf_c,xt);
                        %  > Tf.
                        xf_v        = msh.f.xy.v{j};
                        df          = B1_2D.Assemble_df(c.(k(i)),e.(k(i)),f.qd,xf_v);
                        Pf          = B1_2D.Assemble_Pf(Df,inp.wf,xf_c,xt);
                        Tf          = df*Pf;
                        %  > Update field x".
                        x.Tf  {j,i} = Tf;
                        x.Tf_V{j,i} = inp.c(i,:)'.*Tf;
                        x.Pf  {j,i} = Pf;
                    end
                end
            end
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Check boundary type and assemble matrix Df.
        %    Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
        function [Df] = Assemble_Df(c,e,xf_c,xt)
            d_ft = xt-xf_c;
            Df   = c.*d_ft(:,1).^e{1}(1,:).*d_ft(:,2).^e{1}(2,:);
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Assemble matrix Tf.
        %    Tf = [1,(x-xf),(y-yf),...]*Pf*\phi = df*Pf*\phi.
        function [df] = Assemble_df(c,e,qd,xf_v)
            %  > Auxiliary variables.
            p = qd.pw.Points;
            w = qd.pw.Weights;
            
            %  > xu(u,x).
            for i = 1:size(xf_v,2)
                j         = 1:numel(p);
                xp_v(j,i) = qd.xu(xf_v(:,i),p(j));
            end
            d_fp = xf_v-xp_v;
            %  > df.
            for i = 1:size(c,1)
                j            = 1:numel(p);
                df_i{i}(j,:) = w(j).*c(i,:).*d_fp(j,1).^e{i}(1,:).*d_fp(j,2).^e{i}(2,:);
                df     (i,:) = sum(df_i{i},1);
            end
        end
        %  > 2.2.3. -------------------------------------------------------
        %  > Assemble matrix Pf.
        function [Pf] = Assemble_Pf(Df,wf,xf_c,xt)
            %  > Compute distances to face centroid and check for nil d's.
            d   = sqrt(sum((xt-xf_c).^2,2));
            d_n = d == 0;
            if any(d_n)
                d (d_n) = min(d(~d_n),[],'all')./2;
            end
            %  > Pf = inv(Dwf_T*Df)*Dwf_T, where Dwf = W*Df.
            Dwf = diag(wf(d,max(d)))*Df;
            Pf  = inv (transpose(Dwf)*Df)*transpose(Dwf);
        end
    end
end