classdef B1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "s" (stencil cell/face indices, coordinates, coefficents, etc.).
        function [s] = Initialize_s(inp,f,nc,ns,Nc,Nf,XYc)
            A = [0,2];
            for i = 1:ns
                %  > Fiels: "i", "logical", "sc", "sf" and "tfV".
                s{i}.i       = cell(Nf,nc(2));
                s{i}.logical = cell(Nf,nc(2));
                s{i}.tfV     = cell(Nf,nc(2));
                %  > For plotting...
                if inp.plot{1}(1)
                    s{i}.sc  = cell(Nf,nc(2));
                    s{i}.sf  = cell(Nf,nc(2));
                end
                %  > Field: "u".
                for j = 1:numel(inp.p.p)
                    for k = 1:size(inp.p.p{j},1)
                        s{i}.u.p{j}{k}      = repelem(inp.p.p{j}(k,:)+A(i),Nf,1);
                        s{i}.u.s{j}{k}(:,1) = 1:Nf;
                    end
                end
                %  > Field: "x".
                for j = ["a","x"]
                    s{i}.x.vf. (j) = cell(Nf,nc(2));
                    s{i}.x.xfV.(j) = cell(1 ,nc(2));
                    for k = 1:nc(1)
                        s{i}.x.xfV.(j){k} = zeros(Nf,nc(2));
                    end
                end
                s{i}.x.nv.a = f.fh.f.f(XYc);
                s{i}.x.nv.x = zeros(Nc,1);
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil coordinates, etc.
        function [s] = Update_ss(inp,msh,f,s)
            for i = 1:numel(s.u.s)
                for j = 1:numel(s.u.s{i})
                    if ~isempty(s.u.s{i}{j})
                        for k = s.u.s{i}{j}', p = s.u.p{i}{j}(k,:);
                            %  > Initialize...
                            n   = 1;
                            q   = ceil(max(p)./2);
                            sc  = cell(1,q);
                            sf  = cell(1,q);
                            tfV = cell(1,2); %  > {1}-Cf.
                                             %  > {2}-kf.
                            
                            %  > Loop through levels...
                            while n <= ceil(q)
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
                                %  > (...)
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
                                                    d  = func.setdiff(1:nd,l);
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
                            %  > Logical indexing: 0-bnd.
                            %                      1-blk.
                            if ~msh.f.logical(k) && inp.m.cls
                                sf_t = sf_t(sf_t ~= k);
                            end
                            nc = numel(sc_t);
                            nf = numel(sf_t);
                            if nf ~= 0
                                logical = [true(nc,1);false(nf,1)];
                            else
                                logical = true(nc,1);
                            end
                            %  > Compute coordinates.
                            xf = msh.f.xy.c(k,:);
                            xt = B1_2D.s_xt(msh,sc_t,sf_t);
                            
                            % >> Vdf.
                            %  > Polynomial regression coefficients/exponents.
                            c_Df  = B1_2D.C1_1(p);
                            switch i
                                case 1, c_df = c_Df;            %  > \phi.
                                case 2, c_df = B1_2D.C2_1(p,j); %  > \nabla\phi(x or y).
                                otherwise
                                    return;
                            end
                            %  > Integration point(s)' location.
                            n_ip  = q;
                            Q_ip  = A3_2D.Q_1D_2(n_ip);
                            x_ip  = f.qd.xu(Q_ip.x,msh.f.xy.v{k});
                            %  > Vdf.
                            dx_ip = x_ip-xf;
                            Vdf   = inp.c{i,j}(x_ip).*c_df.c.*dx_ip(:,1).^c_df.e(1,:).*dx_ip(:,2).^c_df.e(2,:);
                            
                            % >> tfV.
                            if ~msh.f.logical(k) && inp.m.cls
                                %  > Compute structure "D".
                                D     = B1_2D.Assemble_cls_D(inp,c_Df,xf,xt);
                                %  > Check boundary type of face "k".
                                bd_c  = B1_2D.Check_bd_t(inp,msh,f,k);
                                %  > Apply constraints...
                                cls_m = B1_2D.Assemble_cls_m(inp,msh,bd_c,f.fh.f,p,msh.f.xy.c(k,:),x_ip);
                                cls_t = func.cls_t(cls_m.b,cls_m.C,D.Dw,D.DwTD);
                                %  > tfV.
                                for i_tfV = 1:numel(cls_t)
                                    tfV{i_tfV} = Q_ip.w'./2*Vdf*cls_t{i_tfV};
                                end
                            else
                                %  > Compute structure "D".
                                D      = B1_2D.Assemble_uls_D(inp,msh,c_Df,f,logical,p,sf_t,xf,xt);
                                %  > tfV.
                                tfV{1} = Q_ip.w'./2*Vdf*func.backlash(D.DwTD,D.Dw'); %  > ...from Cf.
                                tfV{2} = [];                                         %  > ...from kf.
                            end
                            %  > Update field "s".
                            if nf ~= 0
                                s.i  {k,i}{j} = [sc_t;sf_t];
                            else
                                s.i  {k,i}{j} = sc_t;
                            end
                            s.logical{k,i}{j} = logical;
                            s.tfV    {k,i}{j} = tfV;
                            %  > For plotting...
                            if inp.plot{1}(1)
                                if ~msh.f.logical(k) && inp.m.cls
                                    sf{1} = sf{1}(sf{1} ~= k);
                                end
                                s.sf {k,i}{j} = sf;
                                s.sc {k,i}{j} = sc;
                            end
                        end
                    end
                end
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > 1.2.1.1. -----------------------------------------------------
        %  > Face neighbours.
        function [sc,sf] = fn(msh,c,ct,ft)
            %  > Stencil cell(s).
            sc   = func.setdiff(RunLength(sort(cat(1,msh.c.c.nb.f{c}))),ct);
            %  > Stencil face(s).
            fn_f = RunLength(sort(reshape(msh.c.f.if(c,:),[],1)));
            %  > Check faces to be added...
            vb_f = fn_f(~msh.f.logical(fn_f));
            if ~isempty(ft)
                sf = func.setdiff(vb_f,ft);
            else
                sf = vb_f;
            end
        end
        %  > 1.2.1.2. -----------------------------------------------------
        %  > Vertex neighbours.
        function [sc,sf] = vn(msh,c,ct,ft)
            %  > Stencil cell(s).
            sc   = func.setdiff(RunLength(sort(cat(1,msh.c.c.nb.v{c}))),ct);
            %  > Stencil face(s).
            v    = RunLength(sort(reshape(msh.struct.ConnectivityList(c,:),[],1)));
            vn_f = RunLength(sort(cat(1,msh.v.if{v})));
            %  > Check faces to be added...
            vb_f = vn_f(~msh.f.logical(vn_f));
            if ~isempty(ft)
                sf = func.setdiff(vb_f,ft);
            else
                sf = vb_f;
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > 1.2.2.1. -----------------------------------------------------
        %  > Compute coordinates.
        function [xt] = s_xt(msh,sc,sf)
            xc = msh.c.c.xy.c(sc,:);
            if ~isempty(sf)
                xt = [xc;msh.f.xy.c(sf,:)];
            else
                xt = xc;
            end
        end
        %  > 1.2.2.2. -----------------------------------------------------
        %  > Compute stencil limits.
        function [L] = Compute_L(msh,sc,sf)
            xt = B1_2D.s_xt(msh,sc,sf);
            for i = 1:size(xt,2)
                [L(1,i),L(2,i)] = MinMaxElem(xt(:,i),'finite');
            end
        end
        %  > 1.2.3. -------------------------------------------------------
        %  > 1.2.3.1. -----------------------------------------------------
        %  > Compute stencil length/aimensional parameters in the x/y-direction(s).
        function [s] = Compute_adp(msh,sc,sf)
            s.L = B1_2D.Compute_L(msh,sc,sf);
            s.n = ceil((s.L(2,:)-s.L(1,:))./func.mean(msh.c.h.xy(sc,:),1));
        end
        %  > 1.2.3.2. -----------------------------------------------------
        %  > Return extension flag.
        function [e] = Extend_1(msh,p,sc,sf)
            %  > Adimensional parameters in the x/y-direction(s).
            s   = B1_2D.Compute_adp(msh,sc,sf);
            %  > Extend(?).
            e.f = ~(s.n >= p);
            e.L = s.L;
        end
        %  > 1.2.3.3. -----------------------------------------------------
        %  > Perform extension in the appropriate direction.
        function [sn] = Extend_2(inp,msh,d,L,c,ct,ft)
            %  > Compute "new layer"...
            if ~inp.m.nb
                %  > Face neighbours.
                [sn.c,sn.f] = B1_2D.fn(msh,c,ct,ft); %  > [0]: Face   neighbours.
            else
                %  > Vertex neighbours.
                [sn.c,sn.f] = B1_2D.vn(msh,c,ct,ft); %  > [1]: Vertex neighbours.
            end
            %  > Select direction...
            n = 10;
            L = round(L,n);
            %  > c.
            sn.c      = func.setdiff(sn.c,ct);
            xc        = round(msh.c.c.xy.c(sn.c,:),n);
            logical_c = xc(:,d) >= L(1,d) & xc(:,d) <= L(2,d);
            sn.c      = sn.c(logical_c);
            %  > f.
            if ~isempty(sn.f)
                sn.f      = func.setdiff(sn.f,ft);
                xf        = round(msh.f.xy.c(sn.f,:),n);
                logical_f = xf(:,d) >= L(1,d) & xf(:,d) <= L(2,d);
                sn.f      = sn.f(logical_f);
            end
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        %  > Check boundary type.
        function [bd] = Check_bd_t(inp,msh,f,k)
            for i = 1:numel(k)
                bd_k(i,1) = inp.b.t(f.bd.t(f.bd.i == k(i)));
                c   (i)   = msh.f.ic  {k(i)};
                Sf_k(i,:) = msh.c.f.Sf{c(i)}(msh.c.f.if(c(i),:) == k(i),:);
            end
            %  > Assign to structure "bd".
            bd.t  = bd_k;
            bd.Sf = Sf_k;
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > Check boundary type and assemble matrices b and C accordingly.
        function [cls_m] = Assemble_cls_m(inp,msh,bd,fh,p,xf,x_ip)
            %  > dx.
            dx = x_ip-xf;
            %  > Check...
            switch bd.t
                case "Dirichlet"
                    %  > b.
                    cls_b = fh.f(x_ip);
                    %  > C.
                    c     = B1_2D.C1_1(p);
                    cls_C = c.c.*dx(:,1).^c.e(1,:).*dx(:,2).^c.e(2,:);
                case {"Neumann","Robin"}
                    switch bd.t
                        case "Neumann"
                            c = B1_2D.C3_1(p);
                            for j = 1:size(x_ip,2)
                                %  > v.
                                v  (:,j) = fh.d{j}(x_ip);
                                %  > Sf*C.
                                C_Sf {j} = bd.Sf(j).*c{j}.c.*dx(:,1).^c{j}.e(1,:).*dx(:,2).^c{j}.e(2,:);
                            end
                        case "Robin"
                            c = B1_2D.C1_1(p);
                            d = B1_2D.C3_1(p);
                            for j = 1:size(x_ip,2)
                                %  > v.
                                gov (:,j) = inp.c{2,j}(x_ip)./inp.c{1,j}(x_ip);
                                v   (:,j) = fh.f(x_ip)+gov(j).*fh.d{j}(x_ip);
                                %  > Sf*C(\phi).
                                C_Sfv {j} = bd.Sf(j).*c.c.*dx(:,1).^c.e(1,:).*dx(:,2).^c.e(2,:);
                                %  > Sf*C(\nabla\phi).
                                C_Sfg {j} = bd.Sf(j).*d{j}.c.*dx(:,1).^d{j}.e(1,:).*dx(:,2).^d{j}.e(2,:);
                                %  > Sf*C.
                                C_Sf  {j} = C_Sfv{j}+gov(j)*C_Sfg{j};
                            end
                    end
                    cls_b = v*bd.Sf';
                    cls_C = sum(cat(3,C_Sf{:}),3);
            end
            %  > Assign to structure "cls_m".
            cls_m.b = cls_b;
            cls_m.C = cls_C;
        end
        %  > 1.3.3. -------------------------------------------------------
        %  > 1.3.3.1. -----------------------------------------------------
        function [M] = Assemble_cls_D(inp,c,xf,xt)
            %  > dx.
            dx     = xt-xf;
            %  > dd.
            dd     = sqrt(sum(dx.^2,2)); dd(dd == 0) = min(dd(dd ~= 0));
            %  > Assemble matrices D and Dw.
            D      = c.c.*dx(:,1).^c.e(1,:).*dx(:,2).^c.e(2,:);
            Dw     = D.*inp.m.wf(dd);
            %  > Assign to structure "M".
            M.Dw   = Dw;
            M.DwTD = Dw'*D;
        end
        %  > 1.3.3.2. -----------------------------------------------------
        function [M] = Assemble_uls_D(inp,msh,c,f,logical,p,sf_t,xf,xt)
            %  > dx.
            dx = xt-xf;
            %  > dd.
            dd = sqrt(sum(dx.^2,2)); dd(dd == 0) = min(dd(dd ~= 0));
            %  > Assemble matrices D and Dw.
            n  = find(~logical);
            if ~isempty(sf_t)
                %  > Check boundary type of face(s) "sf_t".
                bd = B1_2D.Check_bd_t(inp,msh,f,sf_t);
                if any(bd.t == ["Neumann","Robin"], 'all')
                    d  = B1_2D.C3_1(p);
                end
                %  > #1.
                a      = [find(logical);n(bd.t == "Dirichlet")];
                b      = func.setdiff(1:numel(logical),a);
                D(a,:) = c.c.*dx(a,1).^c.e(1,:).*dx(a,2).^c.e(2,:);
                %  > #2.
                for i  = b
                    switch bd.t(n == i)
                        case "Neumann"
                            for j = 1:numel(d)
                                Dd (j,:) = d{j}.c.*dx(i,1).^d{j}.e(1,:).*dx(i,2).^d{j}.e(2,:);
                            end
                        case "Robin"
                            for j = 1:numel(d)
                                gov(j)   = inp.c{2,j}(xt(i,:))./inp.c{1,j}(xt(i,:));
                                Dd (j,:) = c.c.*dx(i,1).^c.e(1,:).*dx(i,2).^c.e(2,:)+...
                                    gov(j).*d{j}.c.*dx(i,1).^d{j}.e(1,:).*dx(i,2).^d{j}.e(2,:);
                            end
                    end
                    D(i,:) = bd.Sf(n == i,:)*Dd;
                end
            else
                D = c.c.*dx(:,1).^c.e(1,:).*dx(:,2).^c.e(2,:);
            end
            Dw     = D.*inp.m.wf(dd);
            %  > Assign to structure "M".
            M.Dw   = Dw;
            M.DwTD = Dw'*D;
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Compute polynomial regression coefficients/exponents.
        %  > 1.4.1. -------------------------------------------------------
        %  > WLS.
        %  > 1.4.1.1. -----------------------------------------------------
        %  > \phi_f.
        function [t] = C1_1(p)
            %  > x^n*y^n.
            [n{2},n{1}] = meshgrid(0:max(p));
            c           = ones(size(n{1}));
            %  > logical.
            logical  = n{1} <= p(1) & n{2} <= p(2) & sum(cat(3,n{:}),3) <= max(p);
            %  > Assign to structure "t".
            t.c(1,:) = c   (logical); %  > c.
            t.e(1,:) = n{1}(logical); %  > x.
            t.e(2,:) = n{2}(logical); %  > y.
        end
        %  > 1.4.1.2. -----------------------------------------------------
        %  > \nabla\phi_f(x or y).
        function [t] = C2_1(p,k)
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
            v = c; v(v ~= 0) = 1;
            switch k
                case 1, logical = n{1} <= p(1)-1 & n{2} <= p(2)   & v.*sum(cat(3,n{:}),3) <= max(p)-1;
                case 2, logical = n{1} <= p(1)   & n{2} <= p(2)-1 & v.*sum(cat(3,n{:}),3) <= max(p)-1;
            end
            %  > Assign to structure "t".
            t.c(1,:) = c   (logical); % > c.
            t.e(1,:) = n{1}(logical); % > x.
            t.e(2,:) = n{2}(logical); % > y.
        end
        %  > 1.4.1.3. -----------------------------------------------------
        %  > \nabla\phi_f(x and y).
        function [t] = C3_1(p)
            t{1} = B1_2D.C2_1(p,1); %  > x.
            t{2} = B1_2D.C2_1(p,2); %  > y.
        end
        %  > 1.4.2. -------------------------------------------------------
        %  > Direct.
    end
end