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
            for i = 1:size(u.s,2)
                if isempty(u.s{i})
                    continue;
                else
                    for j = 1:size(u.s{i},2)
                        for k = u.s{i}{j}'
                            %  > Initialize...
                            p  = u.p{i}(k,j);
                            sc = [];
                            sf = [];
                            
                            % >> Level 1.
                            %  > Cells to which face "k"'s vertices belong to.
                            fv    = msh.f.iv(k,:);
                            sc{1} = unique(cat(1,msh.v.ic{fv}));
                            
                            %  > Are these boundary cells?
                            bc    = msh.c.logical(sc{1});
                            %  > Add boundary faces...
                            fc    = unique(reshape(msh.c.f.if(sc{1}(bc),:),1,[]));
                            sf{1} = fc(~msh.f.logical(fc))';
                            %  > Check whether the selected boundary faces contain any of "fv"'s vertices.
                            if ~isempty(sf{1})
                                %  > Initialize array "l".
                                m = false(numel(sf{1}),numel(fv));
                                %  > For each of the vertices in "fv"...
                                for l = 1:numel(fv)
                                    m(:,l) = any(bsxfun(@eq,msh.f.iv(sf{1},:),fv(l)),2);
                                end
                                sf{1} = sf{1}(any(m,2));
                                %  > Remove face that shares 1 vertex w/ face "j" if it belongs to the same cell.
                                if any(sf{1} == k)
                                    fv_c             = unique (cat(1,msh.v.if{fv}));
                                    b                = ismembc(sf{1},fv_c);
                                    b   (sf{1} == k) = false;
                                    %  > Check whether face "j" belongs to the same cell...
                                    c_b              = cat(1,msh.f.ic{sf{1}});
                                    c_j              = msh.f.ic{k};
                                    b   (c_b ~= c_j) = false;
                                    %  > Remove...
                                    if nnz(b) == 1
                                        sf{1}(b) = [];
                                    end
                                end
                            end
%                             % >> Level N.
%                             if p > 1
%                                 if ~inp.p.nb_t
%                                     %  > Face neighbours.
%                                     [sc{2},sf{2}] = B1_2D.fn(msh,sc{1},sf{1});
%                                 else
%                                     %  > Vertex neighbours.
%                                     [sc{2},sf{2}] = B1_2D.vn(msh,sc{1},sf{1});
%                                 end
%                             end
                            %  > ------------------------------------------
                            %  > Compute coordinates...
                            xc = B1_2D.xt_c(msh.c.c.xy.c,cat(1,sc{:}));
                            if ~isempty(sf{1})
                                %  > xt.
                                xf      = B1_2D.xt_f(msh.f.xy.c,cat(1,sf{:}));
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
                            if ~isempty(sf{1})
                                s.i  {k,i}{j} = cat(1,cat(1,sc{:}),cat(1,sf{:}));
                            else
                                s.i  {k,i}{j} = cat(1,sc{:});
                            end
                            s.logical{k,i}{j} = logical;
                            s.sc     {k,i}{j} = sc;
                            s.sf     {k,i}{j} = sf;
                            s.xt     {k,i}{j} = xt;
                        end
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
        %  > 2.1.3. -------------------------------------------------------
        %  > 2.1.3.1. -----------------------------------------------------
        %  > Vertex neighbours.
        function [sc,sf] = fn(msh,c,f)
            %  > Stencil cell(s).
            sc   = Tools_1.setdiff(unique(cat(1,msh.c.c.nb.v{c})),c);
            %  > Stencil face(s).
            v    = unique(reshape(msh.struct.ConnectivityList(c,:),[],1));
            vn_f = unique(cat(1,msh.v.if{v}));
            vb_f = vn_f  (~msh.f.logical(vn_f));
            is_f = false (numel(vb_f),1);
            %  > Check faces to be added...
            for i = 1:numel(vb_f)
                is_f(i) = nnz(ismembc(msh.f.iv(vb_f(i),:),v)) > 1;
            end
            if ~isempty(f)
                sf = Tools_1.setdiff(vb_f(is_f),f);
            else
                sf = vb_f(is_f);
            end
        end
        %  > 2.1.3.2. -----------------------------------------------------
        %  > Face neighbours.
        function [sc,sf] = vn(msh,c,f)
            %  > Stencil cell(s).
            sc   = Tools_1.setdiff(unique(cat(1,msh.c.c.nb.v{c})),c);
            %  > Stencil face(s).
            v    = unique(reshape(msh.struct.ConnectivityList(c,:),[],1));
            vn_f = unique(cat(1,msh.v.if{v}));
            vb_f = vn_f  (~msh.f.logical(vn_f));
            %  > Check faces to be added...
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
                x{i}.Df   = cell(Nf,nc(2));
                x{i}.gf   = cell(Nf,nc(2));
                x{i}.Pf   = cell(Nf,nc(2));
                x{i}.Tf_V = cell(Nf,nc(2));
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
            for i = 1:size(u.s,2)
                if isempty(u.s{i})
                    continue;
                else
                    for j = 1:size(u.s{i},2)
                        for k = u.s{i}{j}'
                            %  > p.
                            p       = u.p{i}(k,j);
                            %  > Pf.
                            xy_fc   = msh.f.xy.c(k,:);
                            xy_t    = s.xt      {k,i}{j};
                            [Df,Pf] = B1_2D.Assemble_DPf(inp,p,xy_fc,xy_t);
                            %  > Tf_V.
                            xy_fv = msh.f.xy.v{k};
                            switch i
                                case 1, t = Tools_2.Terms_1(p);
                                case 2, t = Tools_2.Terms_2(p,j);
                            end
                            gf    = B1_2D.Gauss_f(inp.c{i,j},f.qd{i,j},xy_fv);
                            Tf_V  = B1_2D.Assemble_Tf_V(gf,Pf,t,xy_fv);
                            %  > Update field x".
                            x.Df  {k,i}{j} = Df;
                            x.gf  {k,i}{j} = gf;
                            x.Pf  {k,i}{j} = Pf;
                            x.Tf_V{k,i}{j} = Tf_V;
                        end
                    end
                end
            end
        end
        %  > 2.3.1. -------------------------------------------------------
        %  > Check boundary type and assemble matrices Df and Pf.
        function [Df,Pf] = Assemble_DPf(inp,p,xy_fc,xy_t)
            % >> Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
            d_ft = xy_t-xy_fc;
            t    = Tools_2.Terms_1(p);
            Df   = t.c.*d_ft(:,1).^t.e(1,:).*d_ft(:,2).^t.e(2,:);
            
            % >> Pf = inv(Dwf_T*Df)*Dwf_T, where Dwf = W*Df.
            if ~inp.p.wls
                Dwf = Df;
            else
                %  > Compute distances to face centroid and check for nil d's.
                d   = sqrt(sum((xy_t-xy_fc).^2,2));
                d_n = d == 0;
                if any(d_n)
                    d (d_n) = min(d(~d_n),[],'all')./2;
                end
                Dwf = diag(inp.p.wf(d,max(d)))*Df;
            end
            Pf = inv(transpose(Dwf)*Df)*transpose(Dwf);
        end
        %  > 2.3.2. -------------------------------------------------------
        %  > 2.3.2.1. -----------------------------------------------------
        function [g] = Gauss_f(inp_c,qd,xy_fv)
            %  > Auxiliary variables.
            p = qd.pw.Points;
            w = qd.pw.Weights;
            
            %  > (xg,yg).
            for i = 1:size(xy_fv,2)
                g.xy(:,i) = qd.xu(p,xy_fv(:,i));
            end
            %  > V.
            g.V = 1./2.*w.*inp_c(g.xy);
        end
        %  > 2.3.2.2. -----------------------------------------------------
        %  > Assemble matrix Tf_V.
        function [Tf_V] = Assemble_Tf_V(g,Pf,t,xy_fv)
            %  > df_V.
            d_cg = mean(xy_fv,1)-g.xy;
            df   = t.c.*d_cg(:,1).^t.e(1,:).*d_cg(:,2).^t.e(2,:);
            df_V = df.*g.V;
            %  > Tf_V.
            Tf_V = sum(df_V,1)*Pf;
        end
    end
end