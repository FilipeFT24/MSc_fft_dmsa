classdef B1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "s" (stencil cell/face indices, coordinates, coefficents, etc.).
        function [s] = Initialize_s(inp,f,ns,Nf,XYc)
            %  > Auxiliary variables.
            A = [0,2];
            for i = 1:ns
                %  > Fiels: "i", "logical", "sc", "sf" and "tfV".
                s{i}.D        = cell(Nf,1);
                s{i}.i        = cell(Nf,1);
                s{i}.logical  = cell(Nf,1);
                s{i}.sc       = cell(Nf,1);
                s{i}.sf       = cell(Nf,1);
                %  > Field: "u".
                s{i}.u.p      = repelem(inp.p.p+A(i),Nf,1);
                s{i}.u.s      = (1:Nf)';
                %  > Field: "x".
                s{i}.x.s.tfV  = cell(Nf,2); %  > Convection/diffusion.
                s{i}.x.s.wVdf = cell(Nf,2); %  > Convection/diffusion.
                for j = ["a","x"]
                    s{i}.x.x.cf .(j)    = cell (Nf,1);
                    switch j
                        case "a", s{i}.x.x.nv.(j) = f.fh.f.f(XYc);
                        case "x", s{i}.x.x.nv.(j) = zeros   (size(XYc,1),1);
                    end
                    s{i}.x.x.vf. (j)    = cell (Nf,1);
                    s{i}.x.x.xfV.(j){1} = zeros(Nf,2); %  > \phi.
                    s{i}.x.x.xfV.(j){2} = zeros(Nf,2); %  > \nabla\phi.
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil cell/face indices, coordinates, coefficents, etc.
        function [s] = Update_ss(inp,msh,f,s)
            for i = s.u.s', p = s.u.p(i,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % >> Stencil cell/face indices, coordinates, etc.
                %  > Initialize...
                is_bnd = ~msh.f.logical(i);
                n      =  1;
                q      =  ceil(max(p)./2);
                sc     =  cell(1,q);
                sf     =  cell(1,q);
                wVdf   =  cell(1,2);
                
                %  > Loop through levels...
                while n <= ceil(q)
                    if n == 1
                        %  > Face "k"'s vertices.
                        kv    = msh.f.iv (i,:);
                        %  > Cells to which "kv" belong to.
                        sc{n} = RunLength(sort(cat(1,msh.v.ic{kv})));
                        %  > Boundary faces to which "kv" belong to.
                        fv    = RunLength(sort(cat(1,msh.v.if{kv})));
                        sf{n} = fv(~msh.f.logical(fv));
                        %  > Remove faces that belong to the same cell...
                        if is_bnd
                            ck    = msh.f.ic{i};
                            sf{n} = sf{n}(cat(1,msh.f.ic{sf{n}}) ~= ck | sf{n} == i);
                        end
                    else
                        if ~inp.m.nb
                            [sc{n},sf{n}] = B1_2D.fn(msh,sc{n-1},st); %  > [0]: Face   neighbours.
                        else
                            [sc{n},sf{n}] = B1_2D.vn(msh,sc{n-1},st); %  > [1]: Vertex neighbours.
                        end
                    end
                    st.c = cat(1,sc{:});
                    st.f = cat(1,sf{:});
                    
                    %  > Extend...
                    if n ~= 1 && ~isempty(st.f), break_extension = false;
                        while 1
                            e = B1_2D.Extend_1(msh,2*n-1,st);
                            if all(~e.f)
                                break;
                            end
                            %  > Loop through x/y-directions...
                            for l = 1:numel(e.f)
                                if e.f(l)
                                    %  > Direction.
                                    d  = func.setdiff(1:numel(e.f),l);
                                    %  > Update stencil elements...
                                    sn = B1_2D.Extend_2(inp,msh,d,e.L,sc{n},st);
                                    if ~isempty(sn.c)
                                        sc{n} = [sc{n};sn.c]; st.c = cat(1,sc{:});
                                    else
                                        break_extension = true;
                                    end
                                    if ~isempty(sn.f)
                                        sf{n} = [sf{n};sn.f]; st.f = cat(1,sf{:});
                                    end
                                    %  > Update stencil length...
                                    if all(e.f)
                                        e.L = B1_2D.Compute_L(msh,st);
                                    end
                                end
                            end
                            if break_extension
                                break;
                            end
                        end
                    end
                    n = n+1;
                end
                %  > Assign logical indexing.
                [logical,st] = B1_2D.Assign_li(inp.m.cls,is_bnd,i,st);
                %  > Update field "s".
                if ~isempty(st.f)
                    s.i   {i} = [st.c;st.f];
                else
                    s.i   {i} = st.c;
                end
                s.logical {i} = logical;
                if inp.m.cls && is_bnd
                    sf    {1} = sf{1}(sf{1} ~= i);
                end
                s.sc      {i} = sc;
                s.sf      {i} = sf;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %  > xf, xv and xt.
                xf    = msh.f.xy.c(i,:);
                xv    = msh.f.xy.v{i};
                xt    = B1_2D.s_xt(msh,st);
                %  > D.
                W     = B1_2D.Assemble_W(inp,xf,xt);
                c_Df  = B1_2D.C(p,1);
                D     = B1_2D.Assemble_D(inp,msh,c_Df,f,logical,p,st.f,xf,xt,W);
                %  > 1D quadrature.
                Q     = A3_2D.Q_1D_2(q);
                xg    = f.qd.xu(Q.x,xv);
                dx_g  = xg-xf;
                %  > wVdf.
                for j = 1:size(s.x.s.wVdf,2)
                    switch j
                        case 1
                            %  > \phi([x,y]).
                            c_df  = c_Df;
                            x_df  = c_df.c.*dx_g(:,1).^c_df.e(1,:).*dx_g(:,2).^c_df.e(2,:);
                            for k = 1:size(inp.c,2)
                                wVdf{j}(k,:) = Q.w'./2*(inp.c{j,k}(xg).*x_df);
                            end
                        case 2
                            %  > \nabla\phi([x,y]).
                            c_df  = B1_2D.C(p,j);
                            for k = 1:size(inp.c,2)
                                wVdf{j}(k,:) = Q.w'./2*(inp.c{j,k}(xg).*c_df{k}.c.*dx_g(:,1).^c_df{k}.e(1,:).*dx_g(:,2).^c_df{k}.e(2,:));
                            end
                    end
                end
                %  > tfV.
                if inp.m.cls && is_bnd
                    % >> w/  constraint(s).
                    %  > #1.
                    cls_m = B1_2D.Assemble_cls_m(inp,msh,f,i,p,xf,xg);
                    cls_t = func.cls_t(cls_m.b,cls_m.C,D.DTWD,D.DTW);
                    %  > #2.
                    for j = 1:size(s.x.s.wVdf,2)
                        for k = 1:numel(cls_t)
                            tfV{j}{k} = wVdf{j}*cls_t{k};
                        end
                    end
                else
                    % >> w/o constraint(s).
                    %  > #1.
                    DTWD_DTW = func.backslash(D.DTWD,D.DTW);
                    %  > #2.
                    for j = 1:size(s.x.s.wVdf,2)
                        tfV{j}{1} = wVdf{j}*DTWD_DTW; %  > ...from Cf.
                        tfV{j}{2} = [];               %  > ...from kf.
                    end
                end
                %  > Update field "s".
                s.D       {i}   = D;
                s.x.s.tfV (i,:) = tfV;
                s.x.s.wVdf(i,:) = wVdf;
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Assign logical indexing.
        function [logical,st] = Assign_li(cls,is_bnd,f,st)
            if cls && is_bnd
                st.f = st.f(st.f ~= f);
            end
            nc = numel(st.c);
            nf = numel(st.f);
            if nf ~= 0
                logical = [true(nc,1);false(nf,1)];
            else
                logical = true(nc,1);
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > 1.2.2.1. -----------------------------------------------------
        %  > Face neighbours.
        function [sc,sf] = fn(msh,c,t)
            %  > Stencil cell(s).
            sc   = func.setdiff(RunLength(sort(cat(1,msh.c.c.nb.f{c}))),t.c);
            %  > Stencil face(s).
            fn_f = RunLength(sort(reshape(msh.c.f.if(c,:),[],1)));
            %  > Check faces to be added...
            vb_f = fn_f(~msh.f.logical(fn_f));
            if ~isempty(t.f)
                sf = func.setdiff(vb_f,t.f);
            else
                sf = vb_f;
            end
        end
        %  > 1.2.2.2. -----------------------------------------------------
        %  > Vertex neighbours.
        function [sc,sf] = vn(msh,c,t)
            %  > Stencil cell(s).
            sc   = func.setdiff(RunLength(sort(cat(1,msh.c.c.nb.v{c}))),t.c);
            %  > Stencil face(s).
            v    = RunLength(sort(reshape(msh.struct.ConnectivityList(c,:),[],1)));
            vn_f = RunLength(sort(cat(1,msh.v.if{v})));
            %  > Check faces to be added...
            vb_f = vn_f(~msh.f.logical(vn_f));
            if ~isempty(t.f)
                sf = func.setdiff(vb_f,t.f);
            else
                sf = vb_f;
            end
        end
        %  > 1.2.3. -------------------------------------------------------
        %  > 1.2.3.1. -----------------------------------------------------
        %  > Compute coordinates.
        function [xt] = s_xt(msh,st)
            xc = msh.c.c.xy.c(st.c,:);
            if ~isempty(st.f)
                xt = [xc;msh.f.xy.c(st.f,:)];
            else
                xt = xc;
            end
        end
        %  > 1.2.3.2. -----------------------------------------------------
        %  > Compute stencil limits.
        function [L] = Compute_L(msh,st)
            xt = B1_2D.s_xt(msh,st);
            for i = 1:size(xt,2)
                [L(1,i),L(2,i)] = MinMaxElem(xt(:,i),'finite');
            end
        end
        %  > 1.2.4. -------------------------------------------------------
        %  > 1.2.4.1. -----------------------------------------------------
        %  > Compute stencil length/aimensional parameters in the x/y-direction(s).
        function [s] = Compute_adp(msh,st)
            s.L = B1_2D.Compute_L(msh,st);
            s.n = ceil((s.L(2,:)-s.L(1,:))./func.mean(msh.c.h.xy(st.c,:),1));
        end
        %  > 1.2.4.2. -----------------------------------------------------
        %  > Return extension flag.
        function [e] = Extend_1(msh,p,st)
            %  > Adimensional parameters in the x/y-direction(s).
            s   = B1_2D.Compute_adp(msh,st);
            %  > Extend(?).
            e.f = s.n < p;
            e.L = s.L;
        end
        %  > 1.2.4.3. -----------------------------------------------------
        %  > Perform extension in the appropriate direction.
        function [sn] = Extend_2(inp,msh,d,L,c,t)
            %  > Compute "new layer"...
            if ~inp.m.nb
                [sn.c,sn.f] = B1_2D.fn(msh,c,t); %  > [0]: Face   neighbours.
            else
                [sn.c,sn.f] = B1_2D.vn(msh,c,t); %  > [1]: Vertex neighbours.
            end
            %  > Select direction...
            n = 10;
            L = round(L,n);
            %  > c.
            sn.c      = func.setdiff(sn.c,t.c);
            xc        = round(msh.c.c.xy.c(sn.c,:),n);
            logical_c = xc(:,d) >= L(1,d) & xc(:,d) <= L(2,d);
            sn.c      = sn.c(logical_c);
            %  > f.
            if ~isempty(sn.f)
                sn.f      = func.setdiff(sn.f,t.f);
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
        function [cls_m] = Assemble_cls_m(inp,msh,f,k,p,xf,x_ip)
            %  > dx.
            dx = x_ip-xf;
            %  > Check boundary type of face "k".
            bd = B1_2D.Check_bd_t(inp,msh,f,k);
            switch bd.t
                case "Dirichlet"
                    %  > b.
                    cls_b = f.fh.f.f(x_ip);
                    %  > C.
                    c     = B1_2D.C(p,1);
                    cls_C = c.c.*dx(:,1).^c.e(1,:).*dx(:,2).^c.e(2,:);
                case {"Neumann","Robin"}
                    switch bd.t
                        case "Neumann"
                            c = B1_2D.C(p,2);
                            for j = 1:size(x_ip,2)
                                %  > v.
                                v  (:,j) = f.fh.f.d{j}(x_ip);
                                %  > Sf*C.
                                C_Sf {j} = bd.Sf(j).*c{j}.c.*dx(:,1).^c{j}.e(1,:).*dx(:,2).^c{j}.e(2,:);
                            end
                        case "Robin"
                            c = B1_2D.C(p,1);
                            d = B1_2D.C(p,2);
                            for j = 1:size(x_ip,2)
                                %  > v.
                                gov (:,j) = inp.c{2,j}(x_ip)./inp.c{1,j}(x_ip);
                                v   (:,j) = f.fh.f.f(x_ip)+gov(j).*f.fh.f.d{j}(x_ip);
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
        function [W] = Assemble_W(inp,xf,xt)
            %  > dx.
            dx = xt-xf;
            %  > W.
            if ~inp.m.wls
                W  = eye(size(dx,1));
            else
                dd = sqrt(sum(dx.^2,2)); dd(dd == 0) = min(dd(dd ~= 0));
                W  = diag(inp.m.wf(dd));
            end
        end
        %  > 1.3.3.2. -----------------------------------------------------
        function [M] = Assemble_D(inp,msh,c,f,logical,p,sf_t,xf,xt,W)
            %  > dx.
            dx = xt-xf;
            n  = find(~logical);
            %  > Assemble matrices D and Dw.
            if ~isempty(sf_t)
                %  > Check boundary type of face(s) "sf_t".
                bd = B1_2D.Check_bd_t(inp,msh,f,sf_t);
                if any(bd.t == ["Neumann","Robin"], 'all')
                    d  = B1_2D.C(p,2);
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
            %  > Assign to structure "M".
            M.D    = D;
            M.DTW  = bsxfun(@times,diag(W),D)';    %  > D'*W.
            M.DTWD = D.'*bsxfun(@times,D,diag(W)); %  > D'*W*D.
            M.W    = W;
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Compute polynomial regression coefficients/exponents.
        function [t] = C(p,k)
            switch k
                case 1
                    % >> \phi: x^n*y^n.
                    %  > Auxiliary variables.
                    [n{2},n{1}]  = func.meshgrid(0:max(p));
                    %  > Assign to structure "t".
                    c            = ones(size(n{1}));
                    logical      = sum (cat(3,n{:}),3) <= max(p);
                    t.c    (1,:) = c   (logical); %  > c.
                    t.e    (1,:) = n{1}(logical); %  > x.
                    t.e    (2,:) = n{2}(logical); %  > y.
                case 2
                    % >> \nabla\phi([x,y]): (dx)^n*y^n or x^n*(dy)^n.
                    %  > Auxiliary variables.
                    i = 0:max(p)-1; i = [0,i];
                    j = 0:max(p);
                    [n{2,2},n{1,1}] = func.meshgrid(i);
                    [n{1,2},n{2,1}] = func.meshgrid(j); c{1} = n{1,2}'; c{2} = n{2,1}';
                    %  > Assign to structure "t".
                    for l = 1:size(n,1)
                        v{l} = c{l}; v{l}(c{l} ~= 0) = 1;
                        switch l
                            case 1, a = p(1)-1; b = p(2);
                            case 2, a = p(1);   b = p(2)-1;
                        end
                        logical{l}        = n{l,1} <= a & n{l,2} <= b & v{l}.*sum(cat(3,n{l,:}),3) <= max(p)-1;
                        t      {l}.c(1,:) = c{l}  (logical{l}); % > c.
                        t      {l}.e(1,:) = n{l,1}(logical{l}); % > x.
                        t      {l}.e(2,:) = n{l,2}(logical{l}); % > y.
                    end
                otherwise
                    return;
            end
        end
    end
end