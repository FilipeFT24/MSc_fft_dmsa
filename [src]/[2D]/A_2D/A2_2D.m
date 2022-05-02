classdef A2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up 'msh' structure.
        function [msh] = Set_msh(h)
            %  > Auxiliary variables.
            inp = A1_2D.Set_inp_1(h);
            u   = inp.m.Uniform;
            h   = inp.m.h;
            
            if u
                switch inp.m.p{1}
                    case "s"
                        %  > w/ squares.
                        switch inp.m.p{2}
                            case 1
                                %  > Uniform w/ domain: (x,y)={[0,1],[0,1]}.
                                struct = A2_2D.msh_s_1(h);
                            otherwise
                                return;
                        end
                    case "v"
                        %  > Uniform w/ triangles.
                        struct = A2_2D.msh_v_1(h,XLim,YLim);
                    otherwise
                        return;
                end
            else
            end
            %  > ----------------------------------------------------------
            % >> struct.
            msh.struct          = struct;
            CL_c                = struct.ConnectivityList;
            %  > ----------------------------------------------------------
            % >> c.
            msh.c.Nc            = size(CL_c,1);
            %  > c.c.
            %  #1: Cell "i": face/vertex neighbours + auxiliary array "V".
            [msh.c.c.nb,V]      = A2_2D.cc_nb   (CL_c);
            %  #2: Cell "i": centroid/vertex coordinates (for each cell).
            msh.c.c.xy          = A2_2D.cc_xy   (struct);
            %  #3: Cell "i": volume.
            msh.c.Volume        = A2_2D.c_Volume(msh.c.c.xy.v);
            %  #4: Cell "i": reference length.
            msh.c.h             = A2_2D.c_h     (msh.c.c.xy.v,msh.c.Volume);
            %  > c.f.
            %  #1: Identify all/boundary/bulk faces.
            F{1}                = A2_2D.cf_F1   (CL_c,msh.c.c.nb.f);
            %  #2: Cell "i": centroid/vertex coordinates (for each face).
            msh.c.f.xy          = A2_2D.cf_xy   (F{1}.c.v,struct);
            %  #3: Cell "i": outer face normals (Sf).
            msh.c.f.Sf          = A2_2D.cf_Sf_1 (msh.c.c.xy.c,msh.c.f.xy);
            %  > ----------------------------------------------------------
            % >> f.
            %  #1: List faces.
            F{2}                = A2_2D.f_F2(F{1});
            %  -1) Number of (unique) faces.
            msh.f.Nf            = size(F{2}.ic,1);
            %  -2) Face "i": cell indices.
            msh.f.ic            = F{2}.ic;
            %  -3) Face "i": vertex indices.
            msh.f.iv            = F{2}.iv;
            %  -4) Face "i": is it a boundary or bulk face ?
            msh.f.logical       = F{2}.logical;
            %  -5) Boundary cell indices.
            ic_b                = unique([msh.f.ic{~F{2}.logical}]');
            msh.c.logical       = false (msh.c.Nc,1);
            msh.c.logical(ic_b) = true;
            %  #2: Face "i": centroid/vertex coordinates.
            msh.f.xy            = A2_2D.f_xy (msh.f.iv,struct);
            %  #3: Cell "i": face indices.
            msh.c.f.if          = A2_2D.cf_if(CL_c,F);
            %  > ----------------------------------------------------------
            % >> v.
            %  #1: Vertex "i": cell indices.
            msh.v.ic            = A2_2D.v_c      (CL_c,V);
            %  #2: Vertex "i": face indices.
            msh.v.if            = A2_2D.v_f      (F{2});
            %  #3: Identify boundary/bulk vertices.
            msh.v.logical       = A2_2D.v_logical(F{2});
            %  > ----------------------------------------------------------
            % >> d.
            %  #1: Domain reference length.
            msh.d.h             = Tools_1.mean(msh.c.h.h,1);
            %  > Sort fields...
            msh                 = Tools_1.Sort_msh_2D(msh);
            %  > ----------------------------------------------------------
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [Nd,xy_v] = msh_grid(h,XLim,YLim)
            %  > X/Y_Lim.
            XLim    = sort(XLim);
            YLim    = sort(YLim);
            %  > Nd.
            Nd  (1) = round(1./h.*(XLim(2)-XLim(1)));
            Nd  (2) = round(1./h.*(YLim(2)-YLim(1)));
            %  > Vd.
            Vd  {1} = linspace(XLim(1),XLim(2),Nd(1)+1);
            Vd  {2} = linspace(YLim(1),YLim(2),Nd(2)+1);
            %  > X/Yd.
            [Xd,Yd] = meshgrid(Vd{1},Vd{2});
            %  > X/Yv.
            xy_v    = cat(2,reshape(Xd,[],1),reshape(Yd,[],1));
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Example #1 (w/ squares): uniform w/ domain: (x,y)={[0,1],[0,1]}.
        function [struct] = msh_s_1(h)
            %  > Grid limits.
            XLim      = [0,1];
            YLim      = [0,1];
            %  > xy_v.
            [Nd,xy_v] = A2_2D.msh_grid(h,XLim,YLim);
            %  > 'struct'.
            struct    = A2_2D.Set_s(Nd,xy_v);
        end
        % >> 3.2. ---------------------------------------------------------
        %  > Auxiliary function (square grid 'struct').
        function [struct] = Set_s(Nc,xy_v)
            %  > Number of cells/vertices.
            Nv = Nc+1;
            NC = Nc(1).*Nc(2);
            NV = Nv(1).*Nv(2);
            %  > Connectivity list.
            CList                        = reshape(1:NV,Nv(2),Nv(1));
            struct.ConnectivityList(:,1) = reshape(CList(1:Nv(2)-1,1:Nv(1)-1),NC,1); % > SW.
            struct.ConnectivityList(:,2) = reshape(CList(1:Nv(2)-1,2:Nv(1)  ),NC,1); % > SE.
            struct.ConnectivityList(:,3) = reshape(CList(2:Nv(2)  ,2:Nv(1)  ),NC,1); % > NE.
            struct.ConnectivityList(:,4) = reshape(CList(2:Nv(2)  ,1:Nv(1)-1),NC,1); % > NW.
            %  > Points.
            struct.Points                = xy_v;
        end
        
        %% > 4. -----------------------------------------------------------
        % >> 4.1. ---------------------------------------------------------
        function [struct] = msh_v_1(h,XLim,YLim)
            %  > xy_v.
            [~,xy_v] = A2_2D.msh_grid(h,XLim,YLim);
            %  > 'struct'.
            struct   = delaunayTriangulation(xy_v(:,1),xy_v(:,2));
        end
        
        %% > 5. -----------------------------------------------------------
        % >> 5.1. ---------------------------------------------------------
        %  > Field: "c".
        %  > 5.1.1. -------------------------------------------------------
        %  > #1: Cell "i": face/vertex neighbours + auxiliary array "V".
        function [nb,V] = cc_nb(CL_c)
            %  > Auxiliary variables.
            CL_s    = sort (CL_c,2);
            sz      = size (CL_c);
            sz  (3) = max  (CL_c,[],'all');
            V   {1} = false(sz(1),sz(3));
            V   {2} = false(sz(1),sz(2));
            
            %  > Vertex neighbours.
            for i = 1:sz(1)
                V{1}(i,CL_c(i,:)) = true;
            end
            V{1} = sparse(V{1});
            for i = 1:sz(1)
                for j = 1:sz(2)
                    V{2}(i,V{1}(:,CL_c(i,j)) ~= 0) = true;
                end
                V{2}     (i,i) = false;
                nb.v{i,1}(:,1) = find(V{2}(i,:));
            end
            %  > Face neighbours.
            for i = 1:sz(1)
                j = numel(nb.v{i});
                for k = 1:j
                    f_ij{i,1}(k,1) = nnz(ismembc(CL_s(i,:),CL_s(nb.v{i}(k),:))) > 1;
                end
                nb.f{i,1}(:,1) = nb.v{i,1}(f_ij{i,1});
            end
        end
        %  > 5.1.2. -------------------------------------------------------
        %  > #2: Cell "i": centroid/vertex coordinates (for each cell).
        function [xy] = cc_xy(struct)
            %  > Auxiliary variables.
            sz_CL = size(struct.ConnectivityList);
            sz_Pt = size(struct.Points);
            
            for i = 1:sz_CL(1)
                j              = 1:sz_Pt(2);
                xy.v{i,1}(:,j) = struct.Points(struct.ConnectivityList(i,:),j);
                xy.c     (i,j) = Tools_1.mean(xy.v{i,1}(:,j),1);
            end
        end
        %  > 5.1.3. -------------------------------------------------------
        %  > #3: Cell "i": volume.
        function [Volume] = c_Volume(xy_v)
            %  > Auxiliary variables.
            sz = size(xy_v);
            
            for i = 1:sz(1)
                Volume(i,1) = polyarea(xy_v{i}(:,1),xy_v{i}(:,2));
            end
        end
        %  > 5.1.4. -------------------------------------------------------
        %  > #4: Cell "i": reference length.
        function [h] = c_h(xy_v,Volume)
            %  > Auxiliary variables.
            sz = size(Volume);
            
            for i = 1:sz(1)
                j = size(xy_v{i});
                %  > h.
                k = [1:j(1);circshift(1:j,j(1)-1)]';
                for l = 1:j(1)
                    Length{i,1}(l,1) = Tools_1.dist(xy_v{i}(k(l,:),:));
                end
                h.h (i,1) = 4.*Volume(i)./sum(Length{i});
                %  > h(x,y).
                for l = 1:j(2)
                    [L(1,l),L(2,l)] = MinMaxElem(xy_v{i}(:,l),'finite');
                end
                h.xy(i,:) = L(2,:)-L(1,:);
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  > 5.1.5. -------------------------------------------------------
        %  > Identify all/boundary/bulk faces (only face neighbours are evaluated).
        function [F] = cf_F1(CL_c,nb_f)
            %  > Auxiliary variables.
            CL_s  = sort(CL_c,2);
            sz    = size(CL_c);
            vj    = [1:sz(2);circshift(1:sz(2),sz(2)-1)]';
            [~,n] = size(vj);
            
            for i = 1:sz(1)
                %  > Deal all faces (through connectivity matrix).
                for j = 1:sz(2)
                    f_v{i,1}(j,1:n) = sort(CL_c(i,vj(j,:)));
                    f_v{i,1}(j,n+1) = i;
                end
                f_i{i,1} = false(sz(2),1);
                %  > Check and increment accordingly (if bulk face)...
                l = 1;
                for j = 1:numel(nb_f{i})
                    k = ismembc(CL_s(i,:),CL_s(nb_f{i}(j),:));
                    %  > If the evaluated cells have more than 1 vertex in common...
                    if nnz(k) > 1
                        f_blk{i,1}(l) = find(all(bsxfun(@eq,f_v{i,1}(:,1:n),CL_s(i,k)),2)); l = l+1;
                    end
                end
                %  > f_i: 0-bnd.
                %         1-blk.
                f_i{i,1}(f_blk{i,1},1) = true;
            end
            %  > Assign to structure "F".
            F.c.i = f_i;
            F.c.v = f_v;
            F.f.i = cat(1,F.c.i{:});
            F.f.v = cat(1,F.c.v{:});
        end
        %  > 5.1.6. -------------------------------------------------------
        %  > #2: Cell "i": centroid/vertex coordinates (for each face).
        function [xy] = cf_xy(F1,struct)
            %  > Auxiliary variables.
            sz = size(F1,1);
            
            for i = 1:sz
                j = size(F1{i});
                for k = 1:j(1)
                    l                   = 1:j(2)-1;
                    xy.v{i,1}{k,1}(l,:) = struct.Points(F1{i}(k,l),:);         %  > Vertex #k(x,y).
                    xy.c{i,1}     (k,l) = Tools_1.mean(xy.v{i,1}{k,1}(:,l),1); %  > Face   #j(x,y).
                end
            end
        end
        %  > 5.1.7. -------------------------------------------------------
        %  > #3: Cell "i": outer face normals (Sf).
        function [Sf] = cf_Sf_1(xy_m,cf_xy)
            %  > Auxiliary variables.
            sz = size(cf_xy.v);
            
            for i = 1:sz(1)
                j = size(cf_xy.v{i});
                for k = 1:j(1)
                    Sf{i,1}(k,:) = A2_2D.cf_Sf_2(xy_m(i,:),cf_xy.c{i}(k,:),cf_xy.v{i}{k});
                end
            end
        end
        %  > 5.1.7.1. -----------------------------------------------------
        %  > Auxiliary function.
        function [Sf] = cf_Sf_2(xy_cm,xy_fm,xy_vv)
            %  > \vec{FC}.
            FC    = xy_cm-xy_fm;
            FC    = bsxfun(@rdivide,FC,sqrt(sum(FC.^2)));
            %  > \vec{Nf}.
            Sf(2) = xy_vv(1,1)-xy_vv(2,1);
            Sf(1) = xy_vv(2,2)-xy_vv(1,2);
            Nf    = bsxfun(@rdivide,Sf,sqrt(sum(Sf.^2)));
            %  > Check...
            if dot(FC,Nf) > 0
                Sf = -Sf;
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % >> 5.2. ---------------------------------------------------------
        %  > Field: "f".
        %  > 5.2.1. -------------------------------------------------------
        %  > Auxiliary function (identify vertex indices of bulk faces).
        function [blk] = blk_f(F1)
            %  > Auxiliary variables.
            f_blk = F1.f.v(F1.f.i,:);
            [~,n] = size  (f_blk);
            
            %  > Face "i" belongs to cells "j(1)" and "j(2)".
            [~,~,a] = unique(f_blk(:,1:n-1),'rows');
            [~,b,~] = unique(a);
            for i = 1:numel (b)
                blk.if(i,:) = find(a == a(b(i)),2)';
            end
            blk.if = sortrows(blk.if,1);
            %  > Get indices...
            for i = 1:size(blk.if,1)
                blk.ic(i,:) = f_blk(blk.if(i,:),n);
            end
        end
        %  > 5.2.2. -------------------------------------------------------
        %  > List faces.
        function [F2] = f_F2(F1)
            %  > Auxiliary variables.
            f_blk = A2_2D.blk_f(F1);
            l     = 1;
            n     = size(f_blk.if,1);
            
            %  > List of unique faces.
            j     = find(F1.f.i);
            k     = 1;
            [a,b] = sortrows(cat(1,j(f_blk.if(:,k)),find(~F1.f.i)),1);
            %  > Assign to structure "F2".
            F2.iv                = F1.f.v(a(:,1),1:end-1);
            F2.logical(b >  n,1) = false;
            F2.logical(b <= n,1) = true;
            for i = 1:size(a,1)
                if ~F2.logical(i)
                    F2.ic{i,1} = F1.f.v(a(i,1),end);
                else
                    F2.ic{i,1} = f_blk.ic(l,:); l = l+1;
                end
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  > 5.2.3. -------------------------------------------------------
        %  > #2: Face "i": centroid/vertex coordinates.
        function [xy] = f_xy(f,struct)
            %  > Auxiliary variables.
            sz = size(f);
            
            for i = 1:sz(1)
                j              = 1:sz(2);
                xy.v{i,1}(j,:) = struct.Points(f(i,j),:);        %  > Vertex #j(x,y).
                xy.c     (i,j) = Tools_1.mean(xy.v{i,1}(:,j),1); %  > Face   #j(x,y).
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  > 5.2.4. -------------------------------------------------------
        %  #3: Cell "i": face indices.
        function [ic_f] = cf_if(CL_c,F)
            %  > Auxiliary variables.
            sz_c = size (F{1}.c.i,1);
            sz_e = size (CL_c,2);
            sz_f = size (F{2}.ic ,1);
            U    = false(sz_c,sz_f);
            
            for i = 1:sz_f
                U(F{2}.ic{i},i) = true;
            end
            for i = 1:sz_c
                ic_f(i,:) = find(U(i,:),sz_e,'first');
                c         = A2_2D.i_AB(F{1}.c.v{i}(:,1:end-1),F{2}.iv(ic_f(i,:),:));
                ic_f(i,c) = ic_f(i,:);
            end
        end
        %  > 5.2.4.1. -----------------------------------------------------
        %  > Auxiliary function (faster version of built-in function "intersect" w/ 'rows' option).
        function [is_match] = i_AB(A,B)
            %  > Sort/keep track of indices of A and [A;B].
            [A,is_A] = sortrows(A);
            [~,is_T] = sortrows([A;B]);
            %  > Match w/ indices of B.
            is_B     = is_T(is_T > size(A,1))-size(A,1);
            is_match = is_A(is_B);
        end
        % >> 5.3. ---------------------------------------------------------
        %  > Field: "v".
        %  > 5.3.1. -------------------------------------------------------
        %  > #1: Vertex "i" cell indices (similar to routines in "A2_2D.nb_cc").
        function [iv_c] = v_c(CL_c,V)
            %  > Auxiliary variables.
            sz = max(CL_c,[],'all');
            U  = transpose(V{1});
            
            for i = 1:sz
                iv_c{i,1}(:,1) = find(U(i,:));
            end
        end
        %  > 5.3.2. -------------------------------------------------------
        %  > #2: Vertex "i" face indices (similar to routines in "A2_2D.nb_cc").
        function [iv_f] = v_f(F2)
            %  > Auxiliary variables.
            sz_f = size (F2.iv,1);
            sz_v = max  (F2.iv,[],'all');
            U    = false(sz_v,sz_f);
            
            for i = 1:sz_f
                U(F2.iv(i,:),i) = true;
            end
            for i = 1:sz_v
                iv_f{i,1}(:,1) = find(U(i,:),2);
            end
        end
        %  > 5.3.3. -------------------------------------------------------
        %  > #3: Identify boundary/bulk vertices.
        function [logical] = v_logical(F2)
            sz             = max   (F2.iv,[],'all');
            bf_v           = unique(F2.iv(~F2.logical,:));
            logical        = false (sz,1);
            logical(bf_v)  = true;
        end
    end
end