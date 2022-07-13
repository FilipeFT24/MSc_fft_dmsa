classdef B1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "s" (stencil cell/face indices, coordinates, coefficents, etc.).
        function [s] = Initialize_s(inp,f,Nf,XYc)
            %  > Fiels: "i", "logical", "sc", "sf" and "tfV".
            s.D        = cell(Nf,1);
            s.i        = cell(Nf,1);
            s.logical  = cell(Nf,1);
            s.pf       = cell(Nf,1);
            s.sc       = cell(Nf,1);
            s.sf       = cell(Nf,1);
            %  > Field: "u".
            s.u.p      = repelem(inp.P,Nf,1);
            s.u.s      = (1:Nf)';
            %  > Field: "x".
            s.x.s.wVtf = cell(Nf,2); %  > Convection/diffusion.
            s.x.s.wVdf = cell(Nf,2); %  > Convection/diffusion.
            for j = ["a","x"]
                s.x.x.cf.(j) = cell(Nf,1);
                switch j
                    case "a", s.x.x.nv.(j) = f.fh.f.f(XYc);
                    case "x", s.x.x.nv.(j) = zeros   (size(XYc,1),1);
                end
                s.x.x.vf .(j)      = cell (Nf,1);
                s.x.x.xfT.(j){1,1} = zeros(Nf,3); %  >   \phi      (x).
                s.x.x.xfT.(j){1,2} = zeros(Nf,3); %  >   \phi      (y).
                s.x.x.xfT.(j){2,1} = zeros(Nf,3); %  >   \nabla\phi(x).
                s.x.x.xfT.(j){2,2} = zeros(Nf,3); %  >   \nabla\phi(y).
                s.x.x.xfV.(j){1}   = zeros(Nf,2); %  > V*\phi      ([x,y]).
                s.x.x.xfV.(j){2}   = zeros(Nf,2); %  > G*\nabla\phi([x,y]).
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil cell/face indices, coordinates, coefficents, etc.
        function [s] = Update_ss(inp,msh,f,s)
            %  > Initialize...
            K = repmat(inp.M.K,msh.f.Nf,1);
            l = 1;
            while l <= numel(s.u.s), i = s.u.s(l); p = s.u.p(i,:);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % >> Stencil cell/face indices, coordinates, etc.
                %  > Initialize...
                is_b  = ~msh.f.logical(i);
                n     =  1;
                q     =  ceil(max(p)./2);
                sc    =  cell(1,q);
                sf    =  cell(1,q);
                c_df  =  cell(1,size(s.u.p,2));
                wVtf  =  cell(1,size(s.u.p,2));
                wVdf  =  cell(1,size(s.u.p,2));

                %  > Cells to which "kv's" vertices belong to.
                kv    = msh.f.iv (i,:);
                sc{n} = RunLength(sort(cat(1,msh.v.ic{kv})));
                %  > Boundary faces to which "kv's" vertices belong to.
                fv    = RunLength(sort(cat(1,msh.v.if{kv})));
                sf{n} = fv(~msh.f.logical(fv));
                %  > Remove faces that belong to the same cell...
                if is_b
                    ck    = msh.f.ic{i};
                    sf{n} = sf{n}([msh.f.ic{sf{n}}]' ~= ck | sf{n} == i);
                end
                
                if any(p > 1)
                    %  > Initialize...
                    cp   = sc{n};
                    NCFs = B1_2D.NCFs_l(K(i),p);
                    nt   = numel(sc{n});
                    %  > Loop through levels...
                    while n < q
                        %  > Compute face neighbours.
                        [cn,fn] = B1_2D.d_fn(msh,cp,ceil(p./2) > n,cat(1,sc{:}),cat(1,sf{:}));
                        %  > Update...
                        nt      = nt+numel([cn;fn]);
                        sc{n+1} = cat(1,sc{n+1},cn); cp = sc{n+1};
                        sf{n+1} = cat(1,sf{n+1},fn);
                        %  > Check whether an extension is required...
                        if nt > NCFs(n)
                            n = n+1;
                        end
                    end
                end
                st.c = cat(1,sc{:});
                st.f = cat(1,sf{:});
                
                %  > Assign logical indexing.
                [logical,st] = B1_2D.Assign_li(inp.M.Cls,is_b,i,st);
                %  > Update field "s".
                if ~isempty(st.f)
                    s.i   {i} = [st.c;st.f];
                else
                    s.i   {i} = st.c;
                end
                s.logical {i} = logical;
                if inp.M.Cls && is_b
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
                %  > Check whether matrix D is ill-conditioned...
                if isIllConditioned(decomposition(D.D))
                    K(i) = K(i)+0.05; continue;
                else
                    l    = l+1;
                end
                %  > 1D quadrature.
                Q     = A3_2D.Q_1D_2(q);
                xg    = f.qd.xu(Q.x,xv);
                dx_g  = xg-xf;
                %  > c_df and wVdf.
                for j = 1:size(s.x.s.wVdf,2)
                    %  > c_df.
                    switch j
                        case 1, c_df{1}{1} = c_Df;
                                c_df{1}{2} = c_Df;
                        case 2, c_df{2}    = B1_2D.C(p,j);
                    end
                    %  > wVdf.
                    for k = 1:size(inp.c,2)
                        wVdf{j}(k,:) = Q.w'./2*(inp.c{j,k}(xg).*c_df{j}{k}.c.*dx_g(:,1).^c_df{j}{k}.e(1,:).*dx_g(:,2).^c_df{j}{k}.e(2,:));
                    end
                end
                %  > tfV.
                if inp.M.Cls && is_b
                    % >> w/  constraint(s).
                    %  > Apply constraints...
                    cls_m = B1_2D.Assemble_cls_m(inp,msh,f,i,p,xf,xg);
                    cls_t = func.cls_t(cls_m.b,cls_m.C,D.DTW,D.DTWD,D.DTWD_DTW);
                    for j = 1:size(s.x.s.wVdf,2)
                        wVtf{j}{1} = wVdf{j}*cls_t{1};   %  > Pf.
                        wVtf{j}{2} = wVdf{j}*cls_t{2};   %  > kf (add to the RHS).
                    end
                    %  > Assign matrices "cls_m{:}" and "cls_t{:}" to structure "D"...
                    D.Cf = cls_m.C;
                    D.bf = cls_m.b;
                    D.Pf = cls_t  {1};
                    D.kf = cls_t  {2};
                else
                    for j = 1:size(s.x.s.wVdf,2)
                        wVtf{j}{1} = wVdf{j}*D.DTWD_DTW; %  > Pf.
                        wVtf{j}{2} = [];
                    end
                end
                %  > Update field "s".
                s.D       {i}   = D;
                s.pf      {i}   = c_Df.t;
                s.x.s.wVtf(i,:) = wVtf;
                s.x.s.wVdf(i,:) = wVdf;
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Compute coordinates.
        function [xt] = s_xt(msh,st)
            xc = msh.c.c.xy.c(st.c,:);
            if ~isempty(st.f)
                xt = [xc;msh.f.xy.c(st.f,:)];
            else
                xt = xc;
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > Assign logical indexing.
        function [logical,st] = Assign_li(cls,is_b,f,st)
            if cls && is_b
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
        %  > 1.2.3. -------------------------------------------------------
        %  > 1.2.3.1. -----------------------------------------------------
        %  > Face neighbours (NOT USED).
        %  function [sc,sf] = fn(msh,c,t)
        %      %  > Stencil cell(s).
        %      sc   = func.setdiff(RunLength(sort(cat(1,msh.c.c.nb.f{c}))),t.c);
        %      %  > Stencil face(s).
        %      fn_f = RunLength(sort(reshape(msh.c.f.if(c,:),[],1)));
        %      %  > Check faces to be added...
        %      vb_f = fn_f(~msh.f.logical(fn_f));
        %      if ~isempty(t.f)
        %          sf = func.setdiff(vb_f,t.f);
        %      else
        %          sf = vb_f;
        %      end
        %  end
        %  > 1.2.3.2. -----------------------------------------------------
        %  > Vertex neighbours (NOT USED).
        %  function [sc,sf] = vn(msh,c,t)
        %      %  > Stencil cell(s).
        %      sc   = func.setdiff(RunLength(sort(cat(1,msh.c.c.nb.v{c}))),t.c);
        %      %  > Stencil face(s).
        %      v    = RunLength(sort(reshape(msh.struct.ConnectivityList(c,:),[],1)));
        %      vn_f = RunLength(sort(cat(1,msh.v.if{v})));
        %      %  > Check faces to be added...
        %      vb_f = vn_f(~msh.f.logical(vn_f));
        %      if ~isempty(t.f)
        %          sf = func.setdiff(vb_f,t.f);
        %      else
        %          sf = vb_f;
        %      end
        %  end
        %  > 1.2.3.3. -----------------------------------------------------
        %  > (Directional face) neighbours.
        function [sc,sf] = d_fn(msh,c,d,tc,tf)
            %  > Stencil cell(s).
            for i = 1:numel(c)
                nb_f{i} = msh.c.c.nb.f{c(i)};
                j   {i} = B1_2D.fd(msh.c.c.xy.c(c(i),:),msh.c.c.xy.c(msh.c.c.nb.f{c(i)},:));
                sc_j{i} = nb_f{i}(any(j{i}(:,d),2));
            end
            sc = func.setdiff(RunLength(sort(cat(1,sc_j{:}))),tc);
            %  > Stencil face(s).
            fn_f = RunLength(sort(reshape(msh.c.f.if(c,:),[],1)));
            vb_f = fn_f(~msh.f.logical(fn_f));
            %  > Check faces to be added...
            if ~isempty(tf)
                sf = func.setdiff(vb_f,tf);
            else
                sf = vb_f;
            end
        end
        %  > 1.2.3.3.1. ---------------------------------------------------
        %  > Auxiliary function.
        function [flag] = fd(a,b)
            %  > x-direction.
            flag(:,1)   =  ...
               b(:,2)  >= -1.*(b(:,1)-a(1))+a(2) & b(:,2) <= 1.*(b(:,1)-a(1))+a(2) | ...
               b(:,2)  <= -1.*(b(:,1)-a(1))+a(2) & b(:,2) >= 1.*(b(:,1)-a(1))+a(2);
            %  > y-direction.
            flag(:,2)   =  ...
               b(:,2)  >= -1.*(b(:,1)-a(1))+a(2) & b(:,2) >= 1.*(b(:,1)-a(1))+a(2) | ...
               b(:,2)  <= -1.*(b(:,1)-a(1))+a(2) & b(:,2) <= 1.*(b(:,1)-a(1))+a(2);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        %  > Check boundary type.
        function [bd] = Check_bd_t(inp,msh,f,k)
            for i = 1:numel(k)
                bd_k(i,1) = inp.B.T(f.bd.t(f.bd.i == k(i)));
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
            dd = sqrt(sum(dx.^2,2)); dd(dd == 0) = min(dd(dd ~= 0));
            W  = diag(inp.M.Wf(dd));
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
            M.D        = D;
            M.DTW      = bsxfun(@times,diag(W),D)';    %  > D'*W.
            M.DTWD     = D.'*bsxfun(@times,D,diag(W)); %  > D'*W*D.
            M.DTWD_DTW = func.backslash(M.DTWD,M.DTW); %  >(D'*W*D)\(D'W).
            M.W        = W;
        end
        % >> 1.4. ---------------------------------------------------------
        %  > 1.4.1. -------------------------------------------------------
        %  > Compute number of (#) coefficients/level.
        function [NCFs] = NCFs_l(K,p)
            [x{2},x{1}] = func.meshgrid(0:max(p)); y = (max(p)+1)./2-1;
            for i = 1:y
                z       = repmat(2.*i+1,1,2); z(z > p) = p(z > p);
                NCFs(i) = nnz(x{1} <= z(1) & x{2} <= z(2) & sum(cat(3,x{:}),3) <= max(z));
            end
            NCFs(y) = NCFs(y).*K;
        end
        %  > 1.4.2. -------------------------------------------------------
        %  > Compute polynomial regression coefficients/exponents.
        function [t] = C(p,k)
            switch k
                case 1
                    % >> \phi: x^n*y^n.
                    %  > Auxiliary variables.
                    [n{2},n{1}] = func.meshgrid(0:max(p));
                    %  > Assign to structure "t".
                    l        = n{1} <= p(1) & n{2} <= p(2) & sum(cat(3,n{:}),3) <= max(p);
                    c        = ones(max(p)+1);
                    t.c(1,:) = c   (l); %  > c.
                    t.e(1,:) = n{1}(l); %  > x.
                    t.e(2,:) = n{2}(l); %  > y.
                    m        = 1:size(t.c,2);
                    t.t{1}   = m(t.e(1,:) >= t.e(2,:) | (t.e(1,:) == 1 & t.e(2,:) == 0 | t.e(1,:) == 0 & t.e(2,:) == 1)); %  > x (w/ linear profile).
                    t.t{2}   = m(t.e(2,:) >= t.e(1,:) | (t.e(1,:) == 1 & t.e(2,:) == 0 | t.e(1,:) == 0 & t.e(2,:) == 1)); %  > y (w/ linear profile).
                case 2
                    % >> \nabla\phi([x,y]): (dx)^n*y^n or x^n*(dy)^n.
                    %  > Auxiliary variables.
                    i = 0:max(p)-1; i = [0,i];
                    j = 0:max(p);
                    [n{2,2},n{1,1}] = func.meshgrid(i);
                    [n{1,2},n{2,1}] = func.meshgrid(j); c{1} = n{1,2}'; c{2} = n{2,1}';
                    %  > Assign to structure "t".
                    for k = 1:size(n,1)
                        v{k} = c{k}; v{k}(c{k} ~= 0) = 1;
                        switch k
                            case 1, a = p(1)-1; b = p(2);
                            case 2, a = p(1);   b = p(2)-1;
                        end
                        l{k}        = n{k,1} <= a & n{k,2} <= b & v{k}.*sum(cat(3,n{k,:}),3) <= max(p)-1;
                        t{k}.c(1,:) = c{k}  (l{k}); %  > c.
                        t{k}.e(1,:) = n{k,1}(l{k}); %  > x.
                        t{k}.e(2,:) = n{k,2}(l{k}); %  > y.
                    end
                otherwise
                    return;
            end
        end
    end
end