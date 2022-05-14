classdef B1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "s" (stencil cell/face indices, coordinates, etc.).
        function [s] = Initialize_1(nc,ns,Nf)
            for i = 1:ns
                s{i}.i       = cell(Nf,nc(2));
                s{i}.logical = cell(Nf,nc(2));
                s{i}.sc      = cell(Nf,nc(2));
                s{i}.sf      = cell(Nf,nc(2));
                s{i}.xt      = cell(Nf,nc(2));
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil coordinates, etc.
        function [s] = Update_1(inp,msh,s,u)
            s        = B1_2D.Update_sc(inp,msh,s,u);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update stencil coordinates, etc.
        function [s] = Update_sc(inp,msh,s,u)
            %  > NOTE: treat convective/diffusive terms in a unified manner...
            i = 1;
            for j = 1:numel(u.s{i})
                for k = u.s{i}{j}', p = u.p{i}(k,j);
                    %  > Initialize...
                    sc = cell(1,ceil(p./2));
                    sf = cell(1,ceil(p./2));
                    %  > Loop through levels...
                    n = 1;
                    while n <= ceil(p./2)
                        if n == 1
                            %  > Cells to which face "k"'s vertices belong to.
                            fv    = msh.f.iv (k,:);
                            sc{n} = RunLength(sort(cat(1,msh.v.ic{fv})));
                            %  > Are these boundary cells?
                            bc    = msh.c.logical(sc{n});
                            %  > Add boundary faces...
                            fc    = RunLength(sort(reshape(msh.c.f.if(sc{n}(bc),:),[],1)));
                            sf{n} = fc(~msh.f.logical(fc));
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            %  > Check whether the selected boundary faces contain any of "fv"'s vertices.
                            if ~isempty(sf{n})
                                %  > Initialize array "l".
                                m = false(numel(sf{n}),numel(fv));
                                %  > For each of the vertices in "fv"...
                                for l = 1:numel(fv)
                                    m(:,l) = any(bsxfun(@eq,msh.f.iv(sf{n},:),fv(l)),2);
                                end
                                sf{n} = sf{n}(any(m,2));
                                %  > Remove face that shares 1 vertex w/ face "j" if it belongs to the same cell.
                                if any(sf{n} == k)
                                    fv_c            = unique (cat(1,msh.v.if{fv}));
                                    b               = ismembc(sf{n},fv_c);
                                    b  (sf{n} == k) = false;
                                    %  > Check whether face "j" belongs to the same cell...
                                    c_b             = cat(1,msh.f.ic{sf{n}});
                                    c_j             = msh.f.ic{k};
                                    b  (c_b ~= c_j) = false;
                                    %  > Remove...
                                    if nnz(b) == 1
                                        sf{n}(b) = [];
                                    end
                                end
                            end
                            
                            
                            
                            
                            
                        else
                            if ~inp.p.nb_t
                                %  > Face neighbours.
                                [sc{n},sf{n}] = B1_2D.fn(msh,sc{n-1},sf{n-1}); %  > [0].
                            else
                                %  > Vertex neighbours.
                                [sc{n},sf{n}] = B1_2D.vn(msh,sc{n-1},sf{n-1}); %  > [1].
                            end
                        end
                        n = n+1;
                    end
                    
                    %  > Compute coordinates...
                    xc = msh.c.c.xy.c(cat(1,sc{:}),:);
                    if ~isempty(cat(1,sf{:}))
                        %  > xt.
                        xf      = msh.f.xy.c(cat(1,sf{:}),:);
                        xt      = cat(1,xc,xf);
                        %  > logical: 0-bnd.
                        %             1-blk.
                        logical = cat(1,true(size(xc,1),1),false(size(xf,1),1));
                    else
                        %  > xt.
                        xt      = xc;
                        %  > logical: 1-blk.
                        logical = true(size(xc,1),1);
                    end
                    %  > Update field "s".
                    for l = 1:numel(u.s)
                        if ~isempty(cat(1,sf{:}))
                            s.i  {k,l}{j} = cat(1,cat(1,sc{:}),cat(1,sf{:}))';
                        else
                            s.i  {k,l}{j} = cat(1,sc{:});
                        end
                        s.logical{k,l}{j} = logical;
                        s.sc     {k,l}{j} = sc;
                        s.sf     {k,l}{j} = sf;
                        s.xt     {k,l}{j} = xt;
                    end
                end
            end
        end
        %  > 2.1.1. -------------------------------------------------------
        %  > Compute stencil (adimensional) parameters and verify the need for stencil extension(s).
        function [f] = Extend(h,p,sc,xy)
            %  > Stencil limits/adimensional parameters(n) in the x/y-direction(s).
            for i = 1:size(xy,2)
                [L(1,i),L(2,i)] = MinMaxElem(xy(:,i),'finite');
            end
            Lt = diff(L,1);
            n  = ceil(Lt./Tools_1.mean(h(sc,:)));
            %  > Extend(?).
            f  = ~(n >= p);
        end
        %  > 2.1.2. -------------------------------------------------------
        %  > 2.1.2.1. -----------------------------------------------------
        %  > Face neighbours [0].
        function [sc,sf] = fn(msh,c,f)
            %  > Stencil cell(s).
            sc   = Tools_1.setdiff(RunLength(sort(cat(1,msh.c.c.nb.v{c}))),c);
            %  > Stencil face(s).
            fn_f = RunLength(sort(reshape(msh.c.f.if(c,:),[],1)));
            %  > Check faces to be added...
            vb_f = fn_f(~msh.f.logical(fn_f));
            if ~isempty(f)
                sf = Tools_1.setdiff(vb_f,f);
            else
                sf = vb_f;
            end
        end
        %  > 2.1.2.2. -----------------------------------------------------
        %  > Vertex neighbours [1].
        function [sc,sf] = vn(msh,c,f)
            %  > Stencil cell(s).
            sc   = Tools_1.setdiff(RunLength(sort(cat(1,msh.c.c.nb.v{c}))),c);
            %  > Stencil face(s).
            v    = RunLength(sort(reshape(msh.struct.ConnectivityList(c,:),[],1)));
            vn_f = RunLength(sort(cat(1,msh.v.if{v})));
            %  > Check faces to be added...
            vb_f = vn_f(~msh.f.logical(vn_f));
            if ~isempty(f)
                sf = Tools_1.setdiff(vb_f,f);
            else
                sf = vb_f;
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize field "x" (stencil coefficients/nodal solution,etc.).
        function [x] = Initialize_24(f,u,nc,ns,Nc,Nf)
            for i = 1:ns
                %  > "sx".
                x{i}.gf   = cell (Nf,nc(2));
                x{i}.goV  = zeros(Nf,nc(2));
                x{i}.Pf   = cell (Nf,nc(2));
                x{i}.Tf_V = cell (Nf,nc(2));
                %  > "x".
                x{i}.nv.a.f = zeros(Nf,1);
                for j = u{i}.f
                    x{i}.cf.(j) = cell(Nf,nc(2));
                    x{i}.vf.(j) = cell(Nf,nc(2));
                    x{i}.xf.(j) = cell(1 ,nc(2));
                    for k = 1:nc(2)
                        x{i}.xf.(j){k} = zeros(Nf,nc(1)); %  > \phi_f*V(x,y) and \nabla\phi_f*V(x,y).
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
            % >> gf.
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    for k = u.s{i}{j}', p = u.p{i}(k,j);
                        x.gf{k,i}{j} = B1_2D.Gauss_f(inp.c{i,j},f.qd{j},msh.f.xy.v{k});
                    end
                end
            end
            %  > goV.
            %  > NOTE:   Treat convective/diffusive terms in a unified manner...
            %  > Remark: If V=0, goV=Inf (do not try to use a Robin BC if no convection exists!)
            for i = 1:size(x.goV,1)
                for j = 1:size(x.goV,2)
                    x.goV(i,j) = x.gf{i,2}{j}.Vc./x.gf{i,1}{j}.Vc;
                end
            end
            % >> Pf.
            %  > NOTE: Treat convective/diffusive terms in a unified manner...
            i = 1;
            for j = 1:numel(u.s{i})
                for k = u.s{i}{j}', p = u.p{i}(k,j);
                    %  > Pf.
                    Pf = B1_2D.Assemble_Pf(inp,msh,x.goV(k,:),p,...
                        f.bd,s.logical{k,i}{j},cat(1,s.sf{k,i}{j}{:}),msh.f.xy.c(k,:),s.xt{k,i}{j});
                    %  > Update field "x.Pf".
                    for l = 1:numel(u.s)
                        x.Pf{k,l}{j} = Pf;
                    end
                end
            end
            % >> Tf_V.
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    for k = u.s{i}{j}', p = u.p{i}(k,j);
                        switch i
                            case 1, t = Tools_2.Terms_1(p);
                            case 2, t = Tools_2.Terms_2(p,j);
                        end
                        x.Tf_V{k,i}{j} = B1_2D.Assemble_Tf_V(x.gf{k,i}{j},t,x.Pf{k,i}{j});
                    end
                end
            end
        end
        %  > 2.3.1. -------------------------------------------------------
        %  > Check boundary type and assemble matrices Df and Pf.
        function [Pf] = Assemble_Pf(inp,msh,goV,p,f_bd,lc,sf,xy_fc,xy_t)
            % >> Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
            %  > d_ft.
            d_ft       = xy_t-xy_fc;
            %  > Df.
            T1         = Tools_2.Terms_1(p);
            Df         = zeros(size(xy_t,1),size(T1.c,2));
            Df  (lc,:) = T1.c.*d_ft(lc,1).^T1.e(1,:).*d_ft(lc,2).^T1.e(2,:);
            %  > Check whether the stencil contains boundary faces...
            if ~isempty(sf)
                %  > Identify boundary...
                [bd_i,~] = find(bsxfun(@eq,shiftdim(sf,1),f_bd.i)); bd_loc = f_bd.t(bd_i); j = find(~lc);
                %  > Fill remaining rows of Df...
                for k = 1:numel(bd_i)
                    switch inp.b.t(bd_loc(k))
                        case "Dirichlet"
                            %  > Df(j(k)).
                            Df(j(k),:) = T1.c.*d_ft(j(k),1).^T1.e(1,:).*d_ft(j(k),2).^T1.e(2,:);
                        case "Neumann"
                            %  > Sf.
                            c  = msh.f.ic  {sf(k)};
                            Sf = msh.c.f.Sf{c}(msh.c.f.if(c,:) == sf(k),:);
                            %  > Df(j(k)).
                            T2 = Tools_2.Terms_3(p);
                            for l = 1:size(T2.c,1)
                                df(l,:) = T2.c(l,:).*d_ft(j(k),1).^T2.e{l}(1,:).*d_ft(j(k),2).^T2.e{l}(2,:);
                            end
                            Df(j(k),:) = Sf*df;
                        case "Robin"
                            %  > Sf.
                            c  = msh.f.ic  {sf(k)};
                            Sf = msh.c.f.Sf{c}(msh.c.f.if(c,:) == sf(k),:);
                            %  > Df(j(k)).
                            T2 = Tools_2.Terms_3(p);
                            for l = 1:numel(Sf)
                                df_1(l,:) = T1.c     .*d_ft(j(k),1).^T1.e   (1,:).*d_ft(j(k),2).^T1.e   (2,:);
                                df_2(l,:) = T2.c(l,:).*d_ft(j(k),1).^T2.e{l}(1,:).*d_ft(j(k),2).^T2.e{l}(2,:);
                            end
                            Df(j(k),:) = Sf*(df_1+goV'.*df_2);
                    end
                end
            end
            
            % >> Pf = inv(Df'*W*Df)*Df'*W.
            if ~inp.p.wls
                DTWf = Df';
            else
                %  > Compute distances to face centroid and check for nil d's.
                d   = sqrt(sum((xy_t-xy_fc).^2,2));
                d_n = d == 0;
                if any(d_n)
                    d (d_n) = min(d(~d_n))./2;
                end
                DTWf = Df'*diag(inp.p.wf(d,max(d)));
            end
            Pf = inv(DTWf*Df)*DTWf;
        end
        %  > 2.3.2. -------------------------------------------------------
        %  > 2.3.2.1. -----------------------------------------------------
        %  > Handle (1D) face quadrature.
        function [g] = Gauss_f(c,qd,xy_fv)
            %  > (xg,yg).
            xy_c = Tools_1.mean(xy_fv,1);
            xy_g = qd.xu(qd.Points,xy_fv);
            g.d  = xy_c-xy_g;
            %  > V.
            g.Vc = c(xy_c);
            g.Vg = qd.Weights.*c(xy_g)./2;
        end
        %  > 2.3.2.2. -----------------------------------------------------
        %  > Assemble matrix Tf_V.
        function [Tf_V] = Assemble_Tf_V(g,t,Pf)
            %  > d_fV.
            df   = t.c.*g.d(:,1).^t.e(1,:).*g.d(:,2).^t.e(2,:);
            df_V = df.*g.Vg;
            %  > T_fV.
            Tf_V = sum(df_V,1)*Pf;
        end
    end
end