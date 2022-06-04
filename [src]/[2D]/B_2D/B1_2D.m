classdef B1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "s" (stencil cell/face indices, coordinates, etc.).
        function [s] = Initialize_1(inp,msh,nc,ns,Nf)
            for i = 1:ns
                s{i}.i       = cell(Nf,nc(2));
                s{i}.logical = cell(Nf,nc(2));
                if ~inp.m.WLS.t && msh.flag
                    s{i}.nc  = cell(Nf,nc(2));
                end
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
                    if ~isempty(u.s{i}{j})
                        for k = u.s{i}{j}', p = u.p{i}{j}(k,:);
                            %  > Initialize...
                            q  = max (p);
                            sc = cell(1,ceil(q./2));
                            sf = cell(1,ceil(q./2));
                            
                            %  > Loop through levels...
                            n = 1;
                            while n <= ceil(q./2)
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
                                    if ~inp.m.nb
                                        [sc{n},sf{n}] = B1_2D.fn(msh,sc{n-1},sc_t,sf_t); %  > [0]: Face   neighbours.
                                    else
                                        [sc{n},sf{n}] = B1_2D.vn(msh,sc{n-1},sc_t,sf_t); %  > [1]: Vertex neighbours.
                                    end
                                end
                                %  > If anisotropic, exclude cell(s)/face(s) (if necessary)...
                                
                                
                                
                                
                                
                                sc_t = cat(1,sc{:});
                                sf_t = cat(1,sf{:});
                                
                                %  > Extend...
                                if n ~= 1
                                    if ~isempty(sf_t)
                                        %  > ...until no extension is necessary.
                                        break_extension = false;
                                        while 1
                                            e = B1_2D.Extend_1(msh,2*n-1,sc_t,sf_t);
                                            if all(~e.f)
                                                break;
                                            end
                                            %  > Loop through x/y-directions...
                                            nd = numel(e.f);
                                            for l = 1:nd
                                                if e.f(l)
                                                    %  > Direction.
                                                    d  = src_Tools.setdiff(1:nd,l);
                                                    %  > Update stencil elements...
                                                    sn = B1_2D.Extend_2(inp,msh,d,e.L,sc{n},sc_t,sf_t);
                                                    if ~isempty(sn.c)
                                                        sc{n} = [sc{n};sn.c]; sc_t = cat(1,sc{:});
                                                    end
                                                    if ~isempty(sn.f)
                                                        sf{n} = [sf{n};sn.f]; sf_t = cat(1,sf{:});
                                                    end
                                                    %  > Update stencil length (if we need to extend in both directions).
                                                    if all(e.f)
                                                        e.L = B1_2D.Compute_L(msh,sc{n},sf{n});
                                                    end
                                                    %  > Unable to extend stencil...
                                                    if isempty(sn.c)
                                                        break_extension = true;
                                                    end
                                                end
                                            end
                                            if break_extension
                                                break;
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
                            if ~inp.m.WLS.t && msh.flag
                                s_adp         = B1_2D.Compute_adp(msh,sc_t,sf_t);
                                s.nc {k,i}{j} = s_adp.n;
                            end
                            s.sc     {k,i}{j} = sc;
                            s.sf     {k,i}{j} = sf;
                            s.xt     {k,i}{j} = xt;
                        end
                    end
                end
            end
        end
        %  > 1.3.1. -------------------------------------------------------
        %  > 1.3.1.1. -----------------------------------------------------
        %  > Face neighbours.
        function [sc,sf] = fn(msh,c,ct,ft)
            %  > Stencil cell(s).
            sc   = src_Tools.setdiff(RunLength(sort(cat(1,msh.c.c.nb.f{c}))),ct);
            %  > Stencil face(s).
            fn_f = RunLength(sort(reshape(msh.c.f.if(c,:),[],1)));
            %  > Check faces to be added...
            vb_f = fn_f(~msh.f.logical(fn_f));
            if ~isempty(ft)
                sf = src_Tools.setdiff(vb_f,ft);
            else
                sf = vb_f;
            end
        end
        %  > 1.3.1.2. -----------------------------------------------------
        %  > Vertex neighbours.
        function [sc,sf] = vn(msh,c,ct,ft)
            %  > Stencil cell(s).
            sc   = src_Tools.setdiff(RunLength(sort(cat(1,msh.c.c.nb.v{c}))),ct);
            %  > Stencil face(s).
            v    = RunLength(sort(reshape(msh.struct.ConnectivityList(c,:),[],1)));
            vn_f = RunLength(sort(cat(1,msh.v.if{v})));
            %  > Check faces to be added...
            vb_f = vn_f(~msh.f.logical(vn_f));
            if ~isempty(ft)
                sf = src_Tools.setdiff(vb_f,ft);
            else
                sf = vb_f;
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > 1.3.2.1. -----------------------------------------------------
        %  > Compute coordinates.
        function [xt] = s_xt(msh,sc,sf)
            xc = msh.c.c.xy.c(sc,:);
            if ~isempty(sf)
                xt = [xc;msh.f.xy.c(sf,:)];
            else
                xt = xc;
            end
        end
        %  > 1.3.2.2. -----------------------------------------------------
        %  > Compute stencil limits.
        function [L] = Compute_L(msh,sc,sf)
            xt = B1_2D.s_xt(msh,sc,sf);
            for i = 1:size(xt,2)
                [L(1,i),L(2,i)] = MinMaxElem(xt(:,i),'finite');
            end
        end
        %  > 1.3.2.3. -----------------------------------------------------
        %  > Exclude...
        function [sc,sf] = Exclude(msh,d,c,f,ct,ft)
            %  > Use "round" to avoid round-off errors when evaluating expressions ">=" or "<="...
            n         = 10;
            
            %  > L.
            L         = round(B1_2D.Compute_L(msh,ct,ft),n);
            %  > Select direction...
            %  > c.
            sc        = src_Tools.setdiff(c,ct);
            xc        = round(msh.c.c.xy.c(sc,:),n);
            logical_c = xc(:,d) >= L(1,d) & xc(:,d) <= L(2,d);
            sc        = sc(logical_c);
            %  > f.
            if ~isempty(f)
                sf        = src_Tools.setdiff(f,ft);
                xf        = round(msh.f.xy.c(sf,:),n);
                logical_f = xf(:,d) >= L(1,d) & xf(:,d) <= L(2,d);
                sf        = sf(logical_f);
            else
                sf        = [];
            end
        end
        %  > 1.3.3. -------------------------------------------------------
        %  > 1.3.3.1. -----------------------------------------------------
        %  > Compute stencil length/aimensional parameters in the x/y-direction(s).
        function [s] = Compute_adp(msh,sc,sf)
            s.L = B1_2D.Compute_L(msh,sc,sf);
            s.n = ceil((s.L(2,:)-s.L(1,:))./src_Tools.mean(msh.c.h.xy(sc,:),1));
        end
        %  > 1.3.3.2. -----------------------------------------------------
        %  > Return extension flag.
        function [e] = Extend_1(msh,p,sc,sf)
            %  > Adimensional parameters in the x/y-direction(s).
            s   = B1_2D.Compute_adp(msh,sc,sf);
            %  > Extend(?).
            e.f = ~(s.n >= p);
            e.L = s.L;
        end
        %  > 1.3.3.3. -----------------------------------------------------
        %  > Perform extension in the appropriate direction.
        function [sn] = Extend_2(inp,msh,d,L,c,ct,ft)
            %  > Use "round" to avoid round-off errors when evaluating expressions ">=" or "<="...
            n = 10;
            L = round(L,n);
            
            % >> Compute "new layer"...
            if ~inp.m.nb
                %  > Face neighbours.
                [sn.c,sn.f] = B1_2D.fn(msh,c,ct,ft); %  > [0]: Face   neighbours.
            else
                %  > Vertex neighbours.
                [sn.c,sn.f] = B1_2D.vn(msh,c,ct,ft); %  > [1]: Vertex neighbours.
            end
            % >> Select direction...
            %  > c.
            sn.c      = src_Tools.setdiff(sn.c,ct);
            xc        = round(msh.c.c.xy.c(sn.c,:),n);
            logical_c = xc(:,d) >= L(1,d) & xc(:,d) <= L(2,d);
            sn.c      = sn.c(logical_c);
            %  > f.
            if ~isempty(sn.f)
                sn.f      = src_Tools.setdiff(sn.f,ft);
                xf        = round(msh.f.xy.c(sn.f,:),n);
                logical_f = xf(:,d) >= L(1,d) & xf(:,d) <= L(2,d);
                sn.f      = sn.f(logical_f);
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
                x{i}.Tf_V = cell (Nf,nc(2));
                %  > "x".
                x{i}.nv.a.f = zeros(Nf,1);
                for j = ["a","x"]
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
                    if ~isempty(u.s{i}{j})
                        for k = u.s{i}{j}', p = max(u.p{i}{j}(k,:)); p = 7;
                            x.gf{k,i}{j} = B1_2D.Gauss_f(inp.c{i,j},p,f.qd.xu,msh.f.xy.v{k});
                        end
                    end
                end
            end
            % >> Tf_V.
            for i = 1:numel(u.s)
                for j = 1:numel(u.s{i})
                    if ~isempty(u.s{i}{j})
                        for k = u.s{i}{j}', p = u.p{i}{j}(k,:);
                            %  >> d_cf.
                            xc = msh.f.xy.c(k,:);
                            xf = s.xt      {k,i}{j};
                            %  > Remove (boundary) face "k" from "xt".
%                             if ~msh.f.logical(k)
%                                 l       = find(~s.logical{k,i}{j});
%                                 m       = l(cat(1,s.sf{k,i}{j}{:}) == k);
%                                 xf(m,:) = [];
%                             end
                            d_fc = xf-xc;
                            % >> Df.
                            C1       = B1_2D.pr_1       (inp,p);
                            [Df,Dwf] = B1_2D.Assemble_Df(inp,C1,d_fc);
                            %  > Add gauss points' constraints...
                            if ~msh.f.logical(k)
                                d_gc    = x.gf{k,i}{j}.xy-xc;
                                b       = x.gf{k,i}{j}.v;
                                [Cf,~]  = B1_2D.Assemble_Df(inp,C1,d_gc);
                                %Df      = cat(1,Cf,Df);
                                %Dwf     = cat(1,Cf,Dwf);
                                [Gf,Hf] = B1_2D.Assemble_Pf(b,Cf,Df,Dwf);
                            else
                                Pf     = (Dwf'*Df)\Dwf';
                            end
                            
                            Pf     = (Dwf'*Df)\Dwf';%-Gf;
                            %  > Tf_V.
                            switch i
                                case 1, C2 = C1;
                                case 2, C2 = B1_2D.pr_2(inp,p,j);
                            end
                            Tf_V = B1_2D.Assemble_Tf_V(d_gc,x.gf{k,i}{j}.v,C2,Pf);
                            %  > Update field "x".
                            x.Tf_V{k,i}{j} = Tf_V;
                        end
                    end
                end
            end
        end
        %  > 2.3.1. -------------------------------------------------------
        %  > 2.3.1.1. -----------------------------------------------------
        %  > Assemble Df (cell dependent matrix).
        function [Df,Dwf] = Assemble_Df(inp,C1,d_cf)
            %  > Df.
            Df = C1.c.*d_cf(:,1).^C1.e(1,:).*d_cf(:,2).^C1.e(2,:);
            %  > Df*W.
            if ~inp.m.WLS.t
                Dwf = Df;
            else
                w   = inp.m.WLS.Wf(sqrt(sum(d_cf.^2,2)));
                Dwf = bsxfun(@times,Df,w);
            end
            Dwf = Df;
        end
        %  > 2.3.1.2. -----------------------------------------------------
        %  > Assemble Gf and Hf.
        function [Gf,Hf] = Assemble_Pf(b,Cf,Df,Dwf)
            %  > Ef and Ff.
            A1 = Dwf'*Df;
            Ef = A1\Cf';
            Ff = A1\Dwf';
            %  > Gf and Hf.
            A2 = Cf*Ef;
            A3 = Cf'/A2;
            Gf = A1\(Dwf'-A3*(Cf*Ff));
            Hf = A1\(A3*b);
        end
        %  > 2.3.2. -------------------------------------------------------
        %  > 2.3.2.1. -----------------------------------------------------
        %  > Handle (1D) face quadrature.
        function [g] = Gauss_f(c,p,xu,xy_fv)
            %  > "Q" structure.
            Q    = A3_2D.Q_1D_2(ceil(p./2));
            %  > (xg,yg).
            xy_g = xu(Q.Points,xy_fv);
            %  > "g" structure.
            g.xy = xy_g;
            g.v  = c(xy_g);
            g.w  = Q.Weights./2;
        end
        %  > 2.3.2.2. -----------------------------------------------------
        %  > Assemble matrix Tf_V.
        function [Tf_V] = Assemble_Tf_V(dg,vg,t,Pf)
            %  > df_V.
            df   = t.c.*dg(:,1).^t.e(1,:).*dg(:,2).^t.e(2,:);
            df_V = df.*vg;
            %  > Tf_V.
            Tf_V = sum(df_V*Pf,1);
        end
        %  > 2.3.3. -------------------------------------------------------
        %  > Compute polynomial regression coefficients/exponents.
        %  > 2.3.3.1. -----------------------------------------------------
        %  > \phi_f.
        function [t] = pr_1(inp,p)
            %  > x^n*y^n.
            [n{2},n{1}] = meshgrid(0:max(p));
            c           = ones(size(n{1}));
            %  > logical.
            if inp.m.WLS.t
                logical = n{1} <= p(1) & n{2} <= p(2) & sum(cat(3,n{:}),3) <= max(p);
            end
            %  > Assign to structure "t".
            t.c(1,:) = c   (logical); %  > c.
            t.e(1,:) = n{1}(logical); %  > x.
            t.e(2,:) = n{2}(logical); %  > y.
        end
        %  > 2.3.3.2. -----------------------------------------------------
        %  > \nabla\phi_f_x or \nabla\phi_f_y.
        function [t] = pr_2(inp,p,k)
            %  > (dx)^n*y^n or x^n*(dy)^n.
            i = 0:max(p)-1; i = [0,i];
            j = 0:max(p);
            switch k
                case 1
                    %  > (dx)^n*y^n.
                    [~,n{1}] = meshgrid(i);
                    [n{2},~] = meshgrid(j);
                    c        = n{2}';
                case 2
                    %  > x^n*(dy)^n.
                    [~,n{1}] = meshgrid(j);
                    [n{2},~] = meshgrid(i);
                    c        = n{1}';
                otherwise
                    return;
            end
            %  > logical.
            if inp.m.WLS.t
                v = c; v(v ~= 0) = 1;
                switch k
                    case 1, logical = n{1} <= p(1)-1 & n{2} <= p(2)   & v.*sum(cat(3,n{:}),3) <= max(p)-1;
                    case 2, logical = n{1} <= p(1)   & n{2} <= p(2)-1 & v.*sum(cat(3,n{:}),3) <= max(p)-1;
                end
            end
            %  > Assign to structure "t".
            t.c(1,:) = c   (logical); % > c.
            t.e(1,:) = n{1}(logical); % > x.
            t.e(2,:) = n{2}(logical); % > y.
        end
        %  > 2.3.3.3. -----------------------------------------------------
        %  > \nabla\phi_f_x and \nabla\phi_f_y.
        function [t] = pr_3(inp,p)
            t{1} = B1_2D.pr_2(inp,p,1); %  > x.
            t{2} = B1_2D.pr_2(inp,p,2); %  > y.
        end
    end
end