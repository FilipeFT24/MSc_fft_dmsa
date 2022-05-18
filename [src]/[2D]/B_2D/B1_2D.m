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
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    for k = u.s{i}{j}', p = u.p{i}(k,j);
                        %  > Initialize...
                        sc = cell(1,ceil(p./2));
                        sf = cell(1,ceil(p./2));
                        
                        %  > Loop through levels...
                        n = 1;
                        while n <= ceil(p./2)
                            if n == 1
                                %  > Face "k"'s vertices.
                                kv    = msh.f.iv (k,:);
                                %  > Cells to which "kv" belong to.
                                sc{n} = RunLength(sort(cat(1,msh.v.ic{kv})));
                                %  > Boundary faces to which "kv" belong to.
                                fv    = RunLength(sort(cat(1,msh.v.if{kv})));
                                sf{n} = fv(~msh.f.logical(fv));
                                %  > Remove faces that belong to the same cell (case "k" is a boundary face).
                                if ~msh.f.logical(k)
                                    ck    = msh.f.ic{k};
                                    sf{n} = sf{n}(cat(1,msh.f.ic{sf{n}}) ~= ck | sf{n} == k);
                                end
                            else
                                if ~inp.p.nb_t
                                    [sc{n},sf{n}] = B1_2D.fn(msh,sc{n-1},sc_t,sf_t); %  > [0]: Face   neighbours.
                                else
                                    [sc{n},sf{n}] = B1_2D.vn(msh,sc{n-1},sc_t,sf_t); %  > [1]: Vertex neighbours.
                                end
                            end
                            sc_t = cat(1,sc{:});
                            sf_t = cat(1,sf{:});
                            %  > Extend...
                            if n ~= 1
                                if ~isempty(sf_t)
                                    %  > ...until no extension is necessary.
                                    while 1
                                        e = B1_2D.Extend_1(msh,2*n-1,sc_t,sf_t);
                                        if ~all(e.f)
                                            break;
                                        end
                                        %  > Loop through x/y-directions...
                                        nd = numel(e.f);
                                        for l = 1:nd
                                            if e.f(l)
                                                %  > Direction.
                                                dir     = Tools_1.setdiff(1:nd,l);
                                                %  > Update stencil elements...
                                                sn      = B1_2D.Extend_2(inp,msh,dir,e.L,sc{n},sc_t,sf_t);
                                                sc  {n} = [sc{n};sn.c];
                                                sc_t    = cat(1,sc{:});
                                                if ~isempty(sn.f)
                                                    sf  {n} = [sf{n};sn.f];
                                                    sf_t    = cat(1,sf{:});
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            n = n+1;
                        end
                        %  > Compute coordinates...
                        xt = B1_2D.s_xt(msh,sc_t,sf_t);
                        %  > Assign logical indexing: 0-bnd.
                        %                             1-blk.
                        nc = numel (sc_t);
                        nf = numel (sf_t);
                        if ~isempty(sf_t)
                            logical = [true(nc,1);false(nf,1)];
                        else
                            logical = true(nc,1);
                        end
                        %  > Update field "s".
                        if ~isempty(sf_t)
                            s.i  {k,i}{j} = [sc_t;sf_t];
                        else
                            s.i  {k,i}{j} = sc_t;
                        end
                        s.logical{k,i}{j} = logical;
                        s.sc     {k,i}{j} = sc;
                        s.sf     {k,i}{j} = sf;
                        s.xt     {k,i}{j} = xt;
                    end
                end
            end
        end
        %  > 1.3.1. -------------------------------------------------------
        function [xt] = s_xt(msh,sc,sf)
            xc = msh.c.c.xy.c(sc,:);
            if ~isempty(sf)
                xt = [xc;msh.f.xy.c(sf,:)];
            else
                xt = xc;
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > 1.3.2.1. -----------------------------------------------------
        %  > Face neighbours [0].
        function [sc,sf] = fn(msh,c,ct,ft)
            %  > Stencil cell(s).
            sc   = Tools_1.setdiff(RunLength(sort(cat(1,msh.c.c.nb.f{c}))),ct);
            %  > Stencil face(s).
            fn_f = RunLength(sort(reshape(msh.c.f.if(c,:),[],1)));
            %  > Check faces to be added...
            vb_f = fn_f(~msh.f.logical(fn_f));
            if ~isempty(ft)
                sf = Tools_1.setdiff(vb_f,ft);
            else
                sf = vb_f;
            end
        end
        %  > 1.3.2.2. -----------------------------------------------------
        %  > Vertex neighbours [1].
        function [sc,sf] = vn(msh,c,ct,ft)
            %  > Stencil cell(s).
            sc   = Tools_1.setdiff(RunLength(sort(cat(1,msh.c.c.nb.v{c}))),ct);
            %  > Stencil face(s).
            v    = RunLength(sort(reshape(msh.struct.ConnectivityList(c,:),[],1)));
            vn_f = RunLength(sort(cat(1,msh.v.if{v})));
            %  > Check faces to be added...
            vb_f = vn_f(~msh.f.logical(vn_f));
            if ~isempty(ft)
                sf = Tools_1.setdiff(vb_f,ft);
            else
                sf = vb_f;
            end
        end
        %  > 1.3.3. -------------------------------------------------------
        %  > 1.3.3.1. -----------------------------------------------------
        %  > Compute stencil (adimensional) parameters and verify the need for stencil extension(s).
        function [e] = Extend_1(msh,p,sc,sf)
            %  > Limits.
            xt = B1_2D.s_xt(msh,sc,sf);
            for i = 1:size(xt,2)
                [e.L(1,i),e.L(2,i)] = MinMaxElem(xt(:,i),'finite');
            end
            %  > Adimensional parameters in the x/y-direction(s).
            Lt  = e.L(2,:)-e.L(1,:);
            n   = ceil(Lt./Tools_1.mean(msh.c.h.xy(sc,:),1));
            %  > Extend(?).
            e.f = ~(n >= p);
        end
        %  > 1.3.3.2. -----------------------------------------------------
        %  > Perform extension in the appropriate direction.
        function [sn] = Extend_2(inp,msh,dir,L,c,ct,ft)
            % >> Compute "new layer"...
            if ~inp.p.nb_t
                %  > Face neighbours.
                [sn.c,sn.f] = B1_2D.fn(msh,c,ct,ft); %  > [0].
            else
                %  > Vertex neighbours.
                [sn.c,sn.f] = B1_2D.vn(msh,c,ct,ft); %  > [1].
            end
            % >> Select direction...
            %  > c.
            sn.c      = Tools_1.setdiff(sn.c,ct);
            xc        = msh.c.c.xy.c(sn.c,:);
            logical_c = xc(:,dir) >= L(1,dir) & xc(:,dir) <= L(2,dir);
            sn.c      = sn.c(logical_c);
            %  > f.
            if ~isempty(sn.f)
                sn.f      = Tools_1.setdiff(sn.f,ft);
                xf        = msh.f.xy.c(sn.f,:);
                logical_f = xf(:,dir) >= L(1,dir) & xf(:,dir) <= L(2,dir);
                sn.f      = sn.f(logical_f);
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize field "x" (stencil coefficients/nodal solution,etc.).
        function [x] = Initialize_24(f,u,nc,ns,Nc,Nf)
            for i = 1:ns
                %  > "sx".
                x{i}.Df   = cell (Nf,nc(2));
                x{i}.gf   = cell (Nf,nc(2));
                x{i}.goV  = zeros(Nf,nc(2));
                x{i}.Pf   = cell (Nf,nc(2));
                x{i}.Tf_V = cell (Nf,nc(2));
                %  > "x".
                x{i}.nv.a.f = zeros(Nf,1);
                for j = ["a","x"]
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
        %  > Update stencil coefficients, etc.
        function [x] = Update_2(inp,msh,f,s,u,x)
            x        = B1_2D.Update_sx(inp,msh,f,s,u,x);
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
            %  > NOTE: If V=0, goV=Inf (do not try to use a Robin BC if no convection exists!)
            for i = 1:size(x.goV,1)
                for j = 1:size(x.goV,2)
                    x.goV(i,j) = x.gf{i,2}{j}.Vc./x.gf{i,1}{j}.Vc;
                end
            end
            % >> Pf and Tf_V.
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    for k = u.s{i}{j}', p = u.p{i}(k,j);
                        %  > Df and DWf'.
                        [Df,DTWf] = B1_2D.Assemble_Df(inp,msh,x.goV(k,:),p,...
                            f.bd,s.logical{k,i}{j},cat(1,s.sf{k,i}{j}{:}),msh.f.xy.c(k,:),s.xt{k,i}{j});
                        %  > Pf.
                        Pf = (DTWf*Df)\DTWf;
                        %  > Tf.
                        switch i
                            case 1, t = Tools_2.Terms_1(p);
                            case 2, t = Tools_2.Terms_2(p,j);
                        end
                        Tf = B1_2D.Assemble_Tf_V(x.gf{k,i}{j},t,Pf);
                        %  > Update field "x".
                        x.Df  {k,i}{j} = Df;
                        x.Pf  {k,i}{j} = Pf;
                        x.Tf_V{k,i}{j} = Tf;
                    end
                end
            end
        end
        %  > 2.3.1. -------------------------------------------------------
        %  > Check boundary type and assemble matrices Df and Pf.
        function [Df,DTWf] = Assemble_Df(inp,msh,goV,p,f_bd,lc,sf,xy_fc,xy_t)
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
                d(1,:) = sqrt(sum((xy_t-xy_fc).^2,2)); d_n = d == 0;
                if any(d_n)
                    d (d_n) = min(d(~d_n))./2;
                end
                DTWf = bsxfun(@times,Df',inp.p.wf(d));
            end
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