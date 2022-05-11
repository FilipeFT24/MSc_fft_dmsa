classdef A2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up "msh" structure.
        function [msh] = Set_msh(h)
            %  > Auxiliary variables.
            inp_m = A1_2D.Set_inp_1(h);
            
            %  > ----------------------------------------------------------
            % >> struct.
            struct              = Tools_1.Set_struct(inp_m);
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
            %  #2: Cell "i": outer face normals.
            msh.c.f.Sf          = A2_2D.cf_Sf   (F{1},struct,msh.c.c.xy.c);
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
            msh.v.ic            = A2_2D.v_c      (CL_c,V{1});
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
        %  > Field: "c".
        %  > 2.1.1. -------------------------------------------------------
        %  > #1: Cell "i": face/vertex neighbours + auxiliary array "V".
        function [nb,V] = cc_nb(CL_c)
            %  > Auxiliary variables.
            CL_s    = sort (CL_c,2);
            sz      = size (CL_c);
            sz  (3) = max  (CL_c,[],'all');
            V   {1} = false(sz(1),sz(3));
            V   {2} = false(sz(1),sz(1));
            l       = 1:sz(1);
            
            %  > Assemble sparse matrix...
            a     = repelem(1:sz(1),sz(2))';
            b     = reshape(CL_c',[],1);
            V{1}  = sparse (a,b,true(numel(a),1));
            %  > Vertex neighbours.
            for i = 1:sz(1)
                for j = 1:sz(2)
                    V{2}(i,V{1}(:,CL_c(i,j))) = true;
                end
                V{2}     (i,i) = false;
                nb.v{i,1}(:,1) = l(V{2}(i,:));
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
        %  > 2.1.2. -------------------------------------------------------
        %  > #2: Cell "i": centroid/vertex coordinates (for each cell).
        function [xy] = cc_xy(struct)
            %  > Auxiliary variables.
            sz_CL = size(struct.ConnectivityList);
            sz_Pt = size(struct.Points);
            
            for i = 1:sz_CL(1)
                j              = 1:sz_Pt(2);
                xy.v{i,1}(:,j) = struct.Points(struct.ConnectivityList(i,:),j);
                xy.c     (i,j) = Tools_1.mean (xy.v{i,1}(:,j),1);
            end
        end
        %  > 2.1.3. -------------------------------------------------------
        %  > #3: Cell "i": volume.
        function [Volume] = c_Volume(xy_v)
            %  > Auxiliary variables.
            sz = size(xy_v);
            
            for i = 1:sz(1)
                Volume(i,1) = polyarea(xy_v{i}(:,1),xy_v{i}(:,2));
            end
        end
        %  > 2.1.4. -------------------------------------------------------
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
        %  > 2.1.5. -------------------------------------------------------
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
                        f_blk{i,1}(l) = find(all(bsxfun(@eq,f_v{i,1}(:,1:end-1),CL_s(i,k)),2)); l = l+1;
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
        %  > 2.1.6. -------------------------------------------------------
        %  > #2: Cell "i": outer face normals (Sf).
        function [Sf] = cf_Sf(F1,struct,xy_cm)
            for i = 1:size(F1.c.v,1)
                %  > xy_v.
                for k = 1:size(F1.c.v{i},1)
                    xy_v{i,1}{k,1} = struct.Points(F1.c.v{i}(k,1:end-1),:);
                end
                %  > Sf.
                for j  = 1:size(xy_v{i},1)
                    %  > \vec{FC}.
                    vec_FC      = xy_cm(i,:)-Tools_1.mean(xy_v{i}{j},1);
                    %  > \vec{Sf}.
                    vec_Sf(j,2) = xy_v{i}{j}(1,1)-xy_v{i}{j}(2,1);
                    vec_Sf(j,1) = xy_v{i}{j}(2,2)-xy_v{i}{j}(1,2);
                    %  > Check...
                    if dot(vec_FC,vec_Sf(j,:)) > 0
                        vec_Sf(j,:) = -vec_Sf(j,:);
                    end
                end
                Sf{i,1} = vec_Sf;
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % >> 2.2. ---------------------------------------------------------
        %  > Field: "f".
        %  > 2.2.1. -------------------------------------------------------
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
        %  > 2.2.2. -------------------------------------------------------
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
        %  > 2.2.3. -------------------------------------------------------
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
        %  > 2.2.4. -------------------------------------------------------
        %  > #3: Cell "i": face indices.
        function [ic_f] = cf_if(CL_c,F)
            %  > Auxiliary variables.
            sz_c = size (F{1}.c.i,1);
            sz_e = size (CL_c,2);
            sz_f = size (F{2}.ic ,1);
            U    = false(sz_c,sz_f);
            j    = 1:sz_c;
            k    = 1:sz_f;
            
            for i = k
                U(F{2}.ic{i},i) = true;
            end
            for i = j
                ic_f(i,:) = k(U(i,:));
                c         = A2_2D.AB(F{2}.iv(ic_f(i,:),:),F{1}.c.v{i}(:,1:end-1));
                ic_f(i,c) = ic_f(i,:);
            end
        end
        %  > 2.2.4.1. -----------------------------------------------------
        %  > Auxiliary function (faster version of built-in function "intersect" w/ 'rows' option).
        function [l] = AB(A,B)
            [~,i] = sortrows(A);
            [~,j] = sortrows(B);
            [~,k] = sortrows(i);
            l     = j(k);
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Field: "v".
        %  > 2.3.1. -------------------------------------------------------
        %  > #1: Vertex "i" cell indices (similar to routines in "A2_2D.nb_cc").
        function [iv_c] = v_c(CL_c,V)
            %  > Auxiliary variables.
            U  = transpose(V);
            sz = size(U);
            j  = 1:sz(2);
            
            for i = 1:sz(1)
                iv_c{i,1}(:,1) = j(U(i,:));
            end
        end
        %  > 2.3.2. -------------------------------------------------------
        %  > #2: Vertex "i" face indices (similar to routines in "A2_2D.nb_cc").
        function [iv_f] = v_f(F2)
            %  > Auxiliary variables.
            sz_f = size (F2.iv,1);
            sz_v = max  (F2.iv,[],'all');
            U    = false(sz_v,sz_f);
            j    = 1:sz_f;
            k    = 1:sz_v;
            
            for i = j
                U(F2.iv(i,:),i) = true;
            end
            for i = k
                iv_f{i,1}(:,1) = j(U(i,:));
            end
        end
        %  > 2.3.3. -------------------------------------------------------
        %  > #3: Identify boundary/bulk vertices.
        function [logical] = v_logical(F2)
            sz             = max   (F2.iv,[],'all');
            bf_v           = unique(F2.iv(~F2.logical,:));
            logical        = false (sz,1);
            logical(bf_v)  = true;
        end
    end
end