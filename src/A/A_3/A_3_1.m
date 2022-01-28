classdef A_3_1
    methods (Static)
        %% > Wrap-up A_3_1.
        % >> --------------------------------------------------------------
        % >> 1.   Set domain faces.
        % >> 2.   Set grid properties.
        %  > 2.1. Set cell neighbouring cells.
        %  > 2.2. Select bulk faces.
        %  > 2.3. Select boundary faces.
        %  > 2.4. Match bulk face cells.
        % >> 3.   Compute face normals/identify to which boundary do they belong to (based on its direction).
        %  > 3.1. Compute face normals.
        %  > 3.2. Face normal tools.
        % >> 4.   Compute reference length.
        % >> --------------------------------------------------------------
        function [msh] = WrapUp_A_3_1(struct,msh)
            % >> ----------------------------------------------------------
            % >> 1. Wrap up cell (cell components).
            % >> 2. Wrap up face (face components).
            %       Remark: 'WrapUp_Cell' is called first since some of the components of this field help setting up 'WrapUp_Face'.
            % >> ----------------------------------------------------------
            % >> 1.
            for i = 1:msh.c.NC
                %  > Cell centroid.
                msh.c.mean(1,i) = A_Tools.fft_mean(msh.c.xy_v{i}(:,1));
                msh.c.mean(2,i) = A_Tools.fft_mean(msh.c.xy_v{i}(:,2));
            end
            %  > Connectivity list.
            Cn_c = sort(struct.ConnectivityList,2);
            %  > Cell neighbouring cells.
            msh = A_3_1.Set_CellNeighbours(msh,Cn_c);
            
            % >> 2.
            %  > Face coordinates.
            msh = A_3_1.Set_DomainFaces(struct,msh,Cn_c);
            %  > Cell reference length.
            msh = A_3_1.Set_ReferenceLength(msh);
            
            for j = 1:size(msh.f.xy_v,2)
                %  > Face centroid.
                msh.f.mean(1,j) = A_Tools.fft_mean(msh.f.xy_v{j}(:,1));
                msh.f.mean(2,j) = A_Tools.fft_mean(msh.f.xy_v{j}(:,2));
            end
            %  > Face normals.
            msh = A_3_1.Set_FaceNormals(msh);
        end
                                      
        %% > 1. -----------------------------------------------------------
        function [msh] = Set_CellNeighbours(msh,Cn_c)
            %  > Initialize x_ij.
            x_ij = zeros(size(Cn_c,1),size(Cn_c,1));
            
            %  > If any(A) belongs to B, then A and B are neighbours.
            %    Remark: Some constraints (if clauses) were set to speed up evaluation process (VERY time consuming for large grids).
            for i = 1:size(Cn_c,1)
                for j = i+1:size(Cn_c,1)
                    x_ij(i,j) = any(ismembc(Cn_c(i,:),Cn_c(j,:)));
                    x_ij(j,i) = x_ij(i,j);
                end
                msh.c.c{i} = find(x_ij(i,:));
            end
        end
        
        %% > 2. -----------------------------------------------------------
        function [msh] = Set_DomainFaces(struct,msh,Cn_c)
            % >> Domain faces.
            %  > Bulk faces (contains duplicates).
            blk_ij = A_3_1.Pick_BlkFaces(msh,Cn_c);
            %  > Boundary faces.
            bnd_ij = A_3_1.Pick_BndFaces(msh,Cn_c,blk_ij);
            
            % >> Process faces...
            %  > blk.
            blk_f = reshape(blk_ij,[msh.c.NC,1]);
            blk_f = cat(1,blk_f{:});
            shr_f = A_3_1.Match_blkFaces(blk_f);
            %  > bnd.
            bnd_f = reshape(bnd_ij,[size(bnd_ij,2),1]);
            bnd_f = cat(1,bnd_f{:});
            %  > Concatenate arrays.
            for i = 1:size(blk_f,1)+size(bnd_f,1)
                if i <= size(blk_f,1)
                    dom_f(i,1:2) = blk_f(i,1:2);
                    dom_f(i,3)   = blk_f(i,3);
                    dom_f(i,4)   = 0;
                else
                    dom_f(i,1:2) = bnd_f(i-size(blk_f,1),1:2);
                    dom_f(i,3)   = bnd_f(i-size(blk_f,1),3);
                    dom_f(i,4)   = 1;
                end
            end
            %  > Eliminate 'repeated' (bulk) faces.
            uni_f = unique(sort(dom_f(:,1:2),2),'rows','stable');
            
            %  > Re-arranged cell array.
            jj = size(uni_f,1)-size(bnd_f,1);
            for i = 1:size(uni_f,1)
                fin_f{i,1}     = uni_f(i,1);
                fin_f{i,2}     = uni_f(i,2);
                if i <= jj
                    fin_f{i,3} = shr_f(i,:);
                    fin_f{i,4} = 0;
                else
                    fin_f{i,3} = bnd_f(-jj+i,3);
                    fin_f{i,4} = 1;
                end
            end
            
            %  > Number of faces.
            msh.f.NF = size(fin_f,1);
            
            % >> 'Unique' faces.
            %  > Initialize field.
            msh.c.f.f = cell(1,msh.c.NC);
            for i = 1:msh.f.NF
                %  > (Xv,Yv).
                msh.f.xy_v{i}(1,:) = [struct.Points(fin_f{i,1},1),struct.Points(fin_f{i,1},2)];
                msh.f.xy_v{i}(2,:) = [struct.Points(fin_f{i,2},1),struct.Points(fin_f{i,2},2)];
                %  > Track faces that belong to a given cell, e.g.: cell X is composed by faces A, B, C, ...
                for j = 1:length(fin_f{i,3})
                    msh.c.f.f{fin_f{i,3}(j)} = [msh.c.f.f{fin_f{i,3}(j)},i];
                end
            end
            %  > Cell faces.
            dom_f = sortrows(dom_f,3);
            for i = 1:msh.c.NC
                n(i) = mcount(dom_f(:,3),i,'==');
                if i == 1
                    for j = 1:n(i)
                        %  > (Xv,Yv).
                        msh.c.f.xy_v{i}{j}(1,:) = [struct.Points(dom_f(j,1),1),struct.Points(dom_f(j,1),2)];
                        msh.c.f.xy_v{i}{j}(2,:) = [struct.Points(dom_f(j,2),1),struct.Points(dom_f(j,2),2)];
                    end
                    %  > Vertex indices.
                    msh.c.v{i} = unique([dom_f(1:n(i),1);dom_f(1:n(i),2)])';
                else
                    for j = 1:n(i)
                        %  > (Xv,Yv).
                        msh.c.f.xy_v{i}{j}(1,:) = [struct.Points(dom_f(sum(n(1:i-1))+j,1),1),struct.Points(dom_f(sum(n(1:i-1))+j,1),2)];
                        msh.c.f.xy_v{i}{j}(2,:) = [struct.Points(dom_f(sum(n(1:i-1))+j,2),1),struct.Points(dom_f(sum(n(1:i-1))+j,2),2)];
                    end
                    %  > Vertex indices.
                    msh.c.v{i} = unique([dom_f(sum(n(1:i-1))+1:sum(n(1:i-1))+n(i),1);dom_f(sum(n(1:i-1))+1:sum(n(1:i-1))+n(i),2)])';
                end  
            end
            %  > Face vertex indices.
            for i = 1:size(fin_f,1)
                msh.f.v{i} = [fin_f{i,1:2}];
            end
            
            % >> Set neighbours.
            %  > Cells.
            %    Remark: Already done.
            %  > Faces.
            for i = 1:msh.f.NF
                msh.f.c{i} = fin_f{i,3};
            end
            
            % >> Set boundaries.
            %  > Cells.
            uni_c = A_Tools.fft_unique(bnd_f(:,3)');
            for i = 1:length(uni_c)
                %  > (Xv,Yv).
                msh.bnd.c{1,i} = msh.c.f.xy_v{uni_c(i)};
                %  > Cell index.
                msh.bnd.c{2,i} = uni_c(i);
            end
            %  > Faces.
            for i = 1:size(bnd_f)
                %  > (Xv,Yv).
                msh.bnd.f{1,i}(1,:) = [struct.Points(fin_f{jj+i,1},1),struct.Points(fin_f{jj+i,1},2)];
                msh.bnd.f{1,i}(2,:) = [struct.Points(fin_f{jj+i,2},1),struct.Points(fin_f{jj+i,2},2)];
                %  > Face index.
                msh.bnd.f{2,i}      = jj+i;
                %  > Cell index.
                msh.bnd.f{3,i}      = fin_f{jj+i,3};
            end
        end
        % >> 2.1. ---------------------------------------------------------
        function [blk_ij] = Pick_BlkFaces(msh,Cn_c)
            %  > Find bulk faces (i.e. faces w/ 0/1 boundary vertices).
            %    Remark: Only neighbouring cells are evaluated.
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.c{i},2)
                    blk_ij{i}(j,:) = A_Tools.fft_ismember_2(Cn_c(i,:),Cn_c(msh.c.c{i}(j),:));
                end
                %  > Remove '0' rows.
                blk_ij{i}(all(~blk_ij{i},2),:) = [];
                %  > Keep cell index.
                blk_ij{i}(:,3) = i;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [bnd_ij] = Pick_BndFaces(msh,Cn_c,blk_ij)
            for i = 1:msh.c.NC
                %  > Find boundary vertices.
                %    Remark: '0': #matching faces = 0 (incomplete) - Cell w/ 2 faces on boundaries.
                %            '1': #matching faces = 1 (incomplete) - Cell w/ 1 face  on boundaries.
                %            '2': #matching faces = 2 (  complete) - Bulk cell.
                for j = 1:size(Cn_c(i,:),2)
                    Cn_n(i,j) = mcount(blk_ij{i}(:,1:2),Cn_c(i,j),'==');
                end
                %  > Check for any '0s'/'1s'.
                if ~all(Cn_n(i,:) == 2)
                    if mcount(Cn_n(i,:),0,'==')
                        %  > Get indices of '0s'/'1s'.
                        i_0 = find(Cn_n(i,:) == 0);
                        i_1 = find(Cn_n(i,:) == 1);
                        %  > Fill f_ij.
                        bnd_jk{i}(1,:) = [Cn_c(i,i_0),Cn_c(i,i_1(1))];
                        bnd_jk{i}(2,:) = [Cn_c(i,i_0),Cn_c(i,i_1(2))];
                    else
                        %  > Get indices of '1s'.
                        i_1 = find(Cn_n(i,:) == 1);
                        %  > Fill f_ij.
                        bnd_jk{i}(1,:) = [Cn_c(i,i_1(1)),Cn_c(i,i_1(2))];
                    end
                end
            end
            %  > Eliminate empty cells (bulk faces).
            Empty  = cellfun(@isempty,bnd_jk) == 0;
            for i = 1:length(Empty)
                if Empty(i)
                    bnd_jk{i}(:,3) = i;
                end
            end
            bnd_ij = bnd_jk(Empty);
        end
        % >> 2.3. ---------------------------------------------------------
        function [Array_i] = Match_blkFaces(blk_f)
            %  > Find common faces.
            [~,~,id_f]   = unique(sort(blk_f(:,1:2),2),'rows');
            
            %  > Array.
            len          = length(id_f);
            Array_i(:,1) = id_f;
            Array_i(:,2) = 1:1:len;
            Array_o      = zeros(2,len./2);
            Array_v      = zeros(1,len./2);
            
            j = 1; %  > Growth of Array_o(1,:).
            k = 1; %  > Growth of Array_o(2,:).
            for i = 1:len
                if j <= len./2
                    if (i == 1 || ~A_Tools.fft_ismember_1(id_f(i),Array_v(1,1:j))) && j <= len./2
                        Array_o(1,j)    = i;
                        Array_v(1,j)    = id_f(i);
                        j               = j+1;
                    else
                        l(k)            = find(id_f(i) == Array_v(1,1:j));
                        Array_o(2,l(k)) = i;
                        k               = k+1;
                    end
                else
                    l(k)            = find(id_f(i) == Array_v(1,:));
                    Array_o(2,l(k)) = i;
                    k               = k+1;
                end
            end
            Array_o = Array_o';
            
            for i = 1:size(Array_o,1)
                for j = 1:size(Array_o,2)
                    Array_i(i,j) = blk_f(Array_o(i,j),3);
                end
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [msh] = Set_FaceNormals(msh)
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.xy_v{i},1)
                    %  > Face centroid.
                    msh.c.f.mean{i}(1,j) = A_Tools.fft_mean(msh.c.f.xy_v{i}{j}(:,1));
                    msh.c.f.mean{i}(2,j) = A_Tools.fft_mean(msh.c.f.xy_v{i}{j}(:,2));
                    %  > Outer (unit)/(face) normal.
                    [msh.c.f.Nf{i}(:,j),msh.c.f.Sf{i}(:,j)] = ...
                        A_3_1.Tools_FaceNormals(msh.c.f.xy_v{i}{j},msh.c.mean(:,i));
                end
            end
        end     
        % >> 3.2. ---------------------------------------------------------
        function [Nf,Sf] = Tools_FaceNormals(fv,mean_ic)
            % >> Centroids.
            %  > Face centroid.
            if_mean(1) = A_Tools.fft_mean(fv(:,1));
            if_mean(2) = A_Tools.fft_mean(fv(:,2));
            %  > Cell centroid.
            ic_mean(1) = mean_ic(1);
            ic_mean(2) = mean_ic(2);
            
            % >> Arrays.
            %  > \vec{FC}.
            FC(1) = ic_mean(1)-if_mean(1);
            FC(2) = ic_mean(2)-if_mean(2);
            %  > \vec{Nf}.
            Sf(2) = fv(1,1)-fv(2,1);
            Sf(1) = fv(2,2)-fv(1,2);
            %  > Normalize arrays.
            FC = bsxfun(@rdivide,FC,sqrt(sum(FC.^2)));
            Nf = bsxfun(@rdivide,Sf,sqrt(sum(Sf.^2)));
            %  > Check...
            if dot(FC,Nf) > 0
                Nf = -Nf;
                Sf = -Sf;
            end
        end
        
        %% > 4. -----------------------------------------------------------
        function [msh] = Set_ReferenceLength(msh)
            %  > Face length.
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.f.xy_v{i},2)
                    msh.c.f.len{i}(j) = A_Tools.fft_dist(msh.c.f.xy_v{i}{j});
                end
                %  > Face perimeter.
                p        (i) = sum(msh.c.f.len{i});
                %  > Face volume.
                msh.c.vol(i) = polyarea(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2));
            end
            %  > Hydraulic diameter (h) and reference length(h_ref).
            i           = 1:msh.c.NC;
            msh.c.h (i) = 4.*msh.c.vol(i)./p(i);
            msh.d.h_ref = sum(msh.c.h)./msh.c.NC;
        end
    end
end