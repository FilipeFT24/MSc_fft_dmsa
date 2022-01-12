classdef SubClass_1_3
    methods (Static)
        %% > Tools.
        %  > --------------------------------------------------------------
        % >> 1.   Set cells' neighbouring cells.
        % >> 2.   Set domain faces (Main).
        %  > 2.1. Pick_BlkFaces (Select     bulk faces).
        %  > 2.2. Pick_BndFaces (Select boundary faces).
        %  > 2.3. Match bulk faces' cells.
        % >> 3.   Set face normals/identify to which boundary does it belong to based on its direction.
        %  > 3.1. Compute face normals.
        %  > 3.2. Identify boundary -> 1) Nf = [ 1, 0] -> East (E) boundary.
        %                              2) Nf = [ 0, 1] -> North(N) boundary.
        %                              3) Nf = [-1, 0] -> West (W) boundary.
        %                              4) Nf = [ 0,-1] -> South(S) boundary.
        % >> 4.   Stencil setup.
        %  > 4.1. Set stencil neighbours.
        %  > 4.2. Set stencil limits. 
        % >> 5.   Set reference length.
        % >> 6.   Other tools.
        %  > 6.1. Modified version(#1) of ismember/ismembc: Compare array  B w/ sorted array  A.
        %  > 6.2. Modified version(#2) of ismember/ismembc: Compare array  B w/ sorted array  A and return nnz elements of B if nnz(Array) > 2.
        %  > 6.3. Modified version(#1) of isequal: Compare array   B w/ matrix A.
        %  > 6.4. Modified version(#3) of isequal: Compare matrix  B w/ matrix A.
        %  > 6.5. Modified version(#1) of setdiff.
        %  > --------------------------------------------------------------

        %% > 1.) ----------------------------------------------------------
        function [msh] = Set_CellNeighbours(struct,msh)
            %  > Connectivity.
            Cn_c = struct.ConnectivityList;
            %  > Initialize.
            x_ij = zeros(size(Cn_c,1),size(Cn_c,1));

            %  > If any(A) belongs to B, then A and B are neighbours.
            % -> Remark: Some constraints (if clauses) were set to speed up evaluation process (VERY time consuming for large grids).
            for i = 1:size(Cn_c,1)
                for j = i+1:size(Cn_c,1)
                    x_ij(i,j) = SubClass_1_3.fft_ismember_1(Cn_c(i,:),Cn_c(j,:));
                    x_ij(j,i) = x_ij(i,j);
                end
                msh.c.nb{i} = find(x_ij(i,:) ~= 0);
            end
        end
        
        %% > 2.) ----------------------------------------------------------
        function [msh] = Set_DomainFaces(struct,msh)
            %  > Connectivity.
            Cn_c = struct.ConnectivityList;

            % >> Domain faces.
            %  > Bulk faces (contains duplicates).
            blk_ij      = SubClass_1_3.Pick_BlkFaces(msh,Cn_c);
            %  > Boundary faces.
            [bnd_ij,CT] = SubClass_1_3.Pick_BndFaces(msh,Cn_c,blk_ij);
            
            % >> Process faces...
            %  > blk.
            blk_f = cell2mat(reshape(blk_ij,[msh.c.NC,1]));
            shr_f = SubClass_1_3.Match_blkFaces(blk_f);
            %  > bnd.
            bnd_f = cell2mat(reshape(bnd_ij,[size(bnd_ij,2),1]));
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
            %  > Eliminate "repeated" (bulk) faces.
            uni_f = unique(sort(dom_f(:,1:2),2),'rows','stable');
            
            %  > Re-arranged cell.
            jj = size(uni_f,1)-size(bnd_f,1);
            for i = 1:size(uni_f,1)
                fin_f{i,1} = uni_f(i,1);
                fin_f{i,2} = uni_f(i,2);
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
            %  > "Unique" faces.
            for i = 1:msh.f.NF
                %  > (Xv,Yv).
                msh.f.XY_v{i}(1,:) = [struct.Points(fin_f{i,1},1),struct.Points(fin_f{i,1},2)];
                msh.f.XY_v{i}(2,:) = [struct.Points(fin_f{i,2},1),struct.Points(fin_f{i,2},2)];
            end
            %  > Cell faces.
            dom_f = sortrows(dom_f,3);
            for i = 1:msh.c.NC
                n(i) = mcount(dom_f(:,3),i,'==');
                if i == 1
                    for j = 1:n(i)
                        %  > (Xv,Yv).
                        msh.c.f.XY_v{i}{j}(1,:) = [struct.Points(dom_f(j,1),1),struct.Points(dom_f(j,1),2)];
                        msh.c.f.XY_v{i}{j}(2,:) = [struct.Points(dom_f(j,2),1),struct.Points(dom_f(j,2),2)];
                    end
                else
                    for j = 1:n(i)
                        %  > (Xv,Yv).
                        msh.c.f.XY_v{i}{j}(1,:) = [struct.Points(dom_f(sum(n(1:i-1))+j,1),1),struct.Points(dom_f(sum(n(1:i-1))+j,1),2)];
                        msh.c.f.XY_v{i}{j}(2,:) = [struct.Points(dom_f(sum(n(1:i-1))+j,2),1),struct.Points(dom_f(sum(n(1:i-1))+j,2),2)];
                    end
                end
            end
                       
            % >> Set neighbours.
            %  > Cells. 
            % -> Remark: 1) Already done (see call of previous function on 3.2.1.) of SubClass_1_1.
            %            2) That routine "triggers" the remaining routines coded in this SubClass.
            %  > Faces.
            for i = 1:msh.f.NF
                msh.f.nb{i} = fin_f{i,3};
            end
            
            % >> Set boundaries.
            %  > Cells.
            uni_c = unique(bnd_f(:,3));
            for i = 1:length(uni_c)
                %  > (Xv,Yv).
                msh.bnd.c{1,i} = msh.c.f.XY_v{uni_c(i)};
                %  > Cell index.
                msh.bnd.c{2,i} = uni_c(i);
                %  > Type.
                msh.bnd.c{3,i} = CT(uni_c(i));
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
                %  > Identify boundary.
                [msh.bnd.f{4,i},...
                    msh.bnd.f{5,i}] = SubClass_1_3.Identify_bnd(msh,msh.bnd.f{1,i},msh.bnd.f{3,i});
            end
        end
        % >> 2.1.) --------------------------------------------------------
        function [blk_ij] = Pick_BlkFaces(msh,Cn_c)
            % >> Find bulk faces (i.e. faces w/ 0/1 boundary vertices).
            %  > Remark: Only neighbouring cells are evaluated.
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.nb{i},2)
                    blk_ij{i}(j,:) = SubClass_1_3.fft_ismember_2(Cn_c(i,:),Cn_c(msh.c.nb{i}(j),:));
                end
                %  > Remove '0' rows.
                blk_ij{i}(all(~blk_ij{i},2),:) = [];
                %  > Keep cell index.
                blk_ij{i}(:,3) = i;
            end
        end
        % >> 2.2.) --------------------------------------------------------
        function [bnd_ij,Cell_Type] = Pick_BndFaces(msh,Cn_c,blk_ij)
            %  > Initialize.
            Cell_Type = 2.*ones(1,msh.c.NC);
            
            % >> Find boundary faces.
            for i = 1:msh.c.NC
                %  > Find boundary vertices.
                %  > NOTE: '0' -> #matching faces = 0 (incomplete) -> Cell w/ 2 faces on boundaries.
                %          '1' -> #matching faces = 1 (incomplete) -> Cell w/ 1 face  on boundaries.
                %          '2' -> #matching faces = 2 (  complete) -> Bulk cell.
                for j = 1:size(Cn_c(i,:),2)
                    Cn_n(i,j) = mcount(blk_ij{i}(:,1:2),Cn_c(i,j),'==');
                end
                %  > Check for any '0s'/'1s'.
                if ~all(Cn_n(i,:) == 2)
                    if mcount(Cn_n(i,:),0,'==')
                        %  > Keep information.
                        Cell_Type(i) = 0;
                        
                        %  > Get indices of '0'.
                        i_0 = find(Cn_n(i,:) == 0);
                        %  > Get indices of '1s'.
                        i_1 = find(Cn_n(i,:) == 1);
                        %  > Fill f_ij.
                        bnd_jk{i}(1,:) = [Cn_c(i,i_0),Cn_c(i,i_1(1))];
                        bnd_jk{i}(2,:) = [Cn_c(i,i_0),Cn_c(i,i_1(2))];
                    else
                        %  > Keep information.
                        Cell_Type(i) = 1;
                        
                        %  > Get indices of '1s'.
                        i_1 = find(Cn_n(i,:) == 1);
                        %  > Fill f_ij.
                        bnd_jk{i}(1,:) = [Cn_c(i,i_1(1)),Cn_c(i,i_1(2))];
                    end
                end
            end
            %  > Eliminate empty cells (bulk faces)/keep cell index.
            Empty  = cellfun(@isempty,bnd_jk) == 0;
            for i = 1:length(Empty)
                if Empty(i)
                    bnd_jk{i}(:,3) = i;
                end
            end
            bnd_ij = bnd_jk(Empty); 
        end
        % >> 2.3.) --------------------------------------------------------
        function [Array_i] = Match_blkFaces(blk_f)
            %  > Find common faces.
            [~,~,id_f] = unique(sort(blk_f(:,1:2),2),'rows');
            
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
                    if (i == 1 || ~ismembc(id_f(i),sort(Array_v(1,1:j)))) && j <= len./2
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
        
        %% > 3.) ----------------------------------------------------------
        % >> 3.1.) --------------------------------------------------------
        function [msh] = Set_FaceNormals(msh)
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.XY_v{i},1)
                    %  > Face centroid.
                    msh.c.f.mean{i}(1,j) = mean(msh.c.f.XY_v{i}{j}(:,1));
                    msh.c.f.mean{i}(2,j) = mean(msh.c.f.XY_v{i}{j}(:,2));
                    %  > \vec{FC}.
                    FC{i}(1,j) = msh.c.mean(1,i)-msh.c.f.mean{i}(1,j);
                    FC{i}(2,j) = msh.c.mean(2,i)-msh.c.f.mean{i}(2,j);
                    %  > \vec{Nf}.
                    Nf{i}(2,j) = msh.c.f.XY_v{i}{j}(1,1)-msh.c.f.XY_v{i}{j}(2,1);
                    Nf{i}(1,j) = msh.c.f.XY_v{i}{j}(2,2)-msh.c.f.XY_v{i}{j}(1,2);
                    %  > Normalize arrays.
                    FC{i}(:,j) = bsxfun(@rdivide,FC{i}(:,j),sqrt(sum(FC{i}(:,j).^2)));
                    Nf{i}(:,j) = bsxfun(@rdivide,Nf{i}(:,j),sqrt(sum(Nf{i}(:,j).^2)));
                    %  > Check...
                    if dot(FC{i}(:,j),Nf{i}(:,j)) > 0
                        msh.c.f.Nf{i}(:,j) = -Nf{i}(:,j);
                    else
                        msh.c.f.Nf{i}(:,j) =  Nf{i}(:,j);
                    end
                end
            end          
        end
        % >> 3.2.) --------------------------------------------------------
        function [i_bnd,Nf_bnd] = Identify_bnd(msh,f_v,i_c)            
            % >> Centroids.
            %  > Face centroid.
            i_f_mean(1) = mean(f_v(:,1));
            i_f_mean(2) = mean(f_v(:,2));
            %  > Cell centroid.
            i_c_mean(1) = msh.c.mean(1,i_c);
            i_c_mean(2) = msh.c.mean(2,i_c);
            
            % >> Arrays.
            %  > \vec{FC}.
            FC(1) = i_c_mean(1)-i_f_mean(1);
            FC(2) = i_c_mean(2)-i_f_mean(2);
            %  > \vec{Nf}.
            Nf(2) = f_v(1,1)-f_v(2,1);
            Nf(1) = f_v(2,2)-f_v(1,2);
            %  > Normalize arrays.
            FC = bsxfun(@rdivide,FC,sqrt(sum(FC.^2)));
            Nf = bsxfun(@rdivide,Nf,sqrt(sum(Nf.^2)));
            %  > Check...
            if dot(FC,Nf) > 0
                Nf = -Nf;
            end
            
            % >> Set boundary.
            switch true
                case isequal(Nf,[ 1, 0])
                    i_bnd  = 'E';
                case isequal(Nf,[ 0, 1])
                    i_bnd = 'N';
                case isequal(Nf,[-1, 0])
                    i_bnd = 'W';
                case isequal(Nf,[ 0,-1])
                    i_bnd = 'S';  
            end
            Nf_bnd = Nf;
        end
                
        %% > 4.) ----------------------------------------------------------
        % >> 4.1.) --------------------------------------------------------
        function [msh] = Set_FaceNeighbours(msh,Type,Nlev)
            %  > Auxiliary arrays.
            for i = 1:size(msh.bnd.f,2)
                bnd_faces(i) = msh.bnd.f{2,i};
            end
            for j = 1:size(msh.bnd.c,2)
                bnd_cells(j) = msh.bnd.c{2,j};
                bnd_types(j) = msh.bnd.c{3,j};
            end
                       
            % >> Cells to be evaluated (stencil 1).
            %  > Evaluate stencil for neighbours of neighbouring cells of face i (purely for code efficiency).
            for i = 1:msh.f.NF
                j  = length(msh.f.nb{i});
                ki = length(msh.c.nb{msh.f.nb{i}(1)});
                if j == 1
                    nb{i}(1)       = msh.f.nb{i};
                    nb{i}(2:ki+1)  = msh.c.nb{msh.f.nb{i}};
                elseif j == 2
                    kj             = ki+length(msh.c.nb{msh.f.nb{i}(2)});
                    nb{i}(1:ki)    = msh.c.nb{msh.f.nb{i}(1)};
                    nb{i}(ki+1:kj) = msh.c.nb{msh.f.nb{i}(2)};
                end
                nb{i} = unique(nb{i});
            end
            
            % >> Stencils.
            for i = 1:msh.f.NF
                l = 1;
                % >> Level 1.
                for j = 1:length(nb{i})
                    %  > 'Direct' face neighbouring cells.
                    if SubClass_1_3.fft_ismember_1(nb{i}(j),msh.f.nb{i})
                        ngh {1,i}(l) = nb{i}(j);
                        x_ij{1,i}(j) = 1;
                        l            = l+1;
                    else
                        if SubClass_1_3.fft_isequal_1(msh.f.XY_v{i},msh.c.XY_v{nb{i}(j)})
                            ngh {1,i}(l) = nb{i}(j);
                            x_ij  {i}(j) = 1;
                            l            = l+1;
                        else
                            x_ij{i}(j) = 0;
                        end
                    end
                end
                msh.s.st{1,i} = ngh{1,i};
                
                % >> Level N.
                if Nlev > 1
                    for j = 2:Nlev
                        % >> Select stencil type.
                        %  > Vertex type: Requires (at least) 1 common vertex (i.e. common vertex/face).
                        %  > Cell   type: Requires 2 common vertices (i.e. common face).
                        if strcmpi(Type,'Vertex')
                            for k = 1:length(msh.s.st{j-1,i})
                                if k == 1
                                    ngh{j,i} = msh.c.nb{msh.s.st{j-1,i}(k)};
                                else
                                    ngh{j,i} = [ngh{j,i},msh.c.nb{msh.s.st{j-1,i}(k)}];
                                end
                            end
                        elseif strcmpi(Type,'Face')
%                             %  > For each (previous) stencil cell, select eligible neighbouring cells.
%                             for k = 1:length(msh.s.st{j-1,i})
%                                 %  > Check cell neighbours that do not belong to previous stencil.
%                                 if j == 2
%                                     elig_c{j-1,i}{k} = SubClass_1_3.fft_setdiff(msh.c.nb{msh.s.st{j-1,i}(k)},msh.s.st{1,i});
%                                 else
%                                     elig_c{j-1,i}{k} = SubClass_1_3.fft_setdiff(msh.c.nb{msh.s.st{j-1,i}(k)},exc{j-2,i});
%                                 end
%                             end
%                             v_ijk{j-1,i} = elig_c{j-1,i}(~cellfun('isempty',elig_c{j-1,i})) ;
%                             v_ijk{j-1,i} = [v_ijk{j-1,i}{:}];
%                             
%                             % >> Check if it is a boundary face. 
%                             if ismembc(i,bnd_faces)
%                                 for k = 1:length(v_ijk{j-1,i})
%                                     numb{j-1,i}(k) = mcount(v_ijk{j-1,i},v_ijk{j-1,i}(k),'==');
%                                 end
%                                 %  > "ik" cells only share 1 face.
%                                 ik{j-1,i} = find(numb{j-1,i} == 1);
%                                 
%                                 %  > Evaluate whether cells "ik" share common face with stencil (j-1) cells.
%                                 for k = 1:length(ik{j-1,i})
%                                     for l = 1:length(msh.s.st{j-1,i})
%                                         i_Flag{j-1,i}(k,l) = SubClass_1_3.fft_isequal_2(msh.c.XY_v{v_ijk{j-1,i}(ik{j-1,i}(k))},msh.c.XY_v{msh.s.st{j-1,i}(l)});
%                                     end
%                                     if any(i_Flag{j-1,i}(k,:) == 1)
%                                         v_ijk{j-1,i} = [v_ijk{j-1,i},v_ijk{j-1,i}(ik{j-1,i}(k))];
%                                     end
%                                 end
%                             end
%                             %  > Remove unique elements.
%                             ngh{j,i} = find(accumarray(v_ijk{j-1,i}.',ones(size(v_ijk{j-1,i}))) > 1)';
                        end
                        %  > Exclude repeated cells.
                        ngh{j,i} = unique(ngh{j,i});
                        
                        %  > Exclude lower-order stencil cells.
                        if j == 2
                            exc{j-1,i} = msh.s.st{1,i};
                        else
                            exc{j-1,i} = [msh.s.st{j-1,i},exc{j-2,i}];
                        end
                        msh.s.st{j,i} = SubClass_1_3.fft_setdiff(ngh{j,i},exc{j-1,i});
                    end
                end
            end
        end
        % >> 4.2.) --------------------------------------------------------
        function [msh] = Set_Limits(msh)
            for i = 1:msh.f.NF
                %  > (x,y)_min.
                msh.s.xy_min(1,i) = min(msh.c.mean(1,msh.s.st{i}));
                msh.s.xy_min(2,i) = min(msh.c.mean(2,msh.s.st{i}));
                %  > (x,y)max.
                msh.s.xy_max(1,i) = max(msh.c.mean(1,msh.s.st{i}));
                msh.s.xy_max(2,i) = max(msh.c.mean(2,msh.s.st{i}));
            end
        end
        
        %% > 5.) ----------------------------------------------------------
        function [msh] = Set_ReferenceLength(msh)
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.f.XY_v{i},2)
                    %  > Length.
                    msh.c.f.len{i}(j) = pdist([msh.c.f.XY_v{i}{j}(:,1),msh.c.f.XY_v{i}{j}(:,2)],'euclidean');
                end
                N      (i) = size(msh.c.XY_v{i},1);
                p      (i) = sum(msh.c.f.len{i});
                msh.c.H(i) = N(i).*msh.c.Vol(i)./p(i);
            end
            msh.c.H_ref = sum(msh.c.H)./msh.c.NC;
        end
        
        %% > 6.) ----------------------------------------------------------
        % >> 6.1.) --------------------------------------------------------
        function [Flag] = fft_ismember_1(A,B)
            Flag = double(any(ismembc(B,sort(A))));
        end
        % >> 6.2.) --------------------------------------------------------
        function [Flag] = fft_ismember_2(A,B)
            Array = find(double(ismembc(B,sort(A))));
            if size(Array) < 2
                Flag = zeros(1,2);
            else
                Flag = B(Array);
            end  
        end
        % >> 6.3.) --------------------------------------------------------
        function [i_Flag] = fft_isequal_1(A,B)
            i_Flag = 0;
            for i = 1:size(A,1)
                for j = 1:size(B,1)
                    Flag(i,j) = isequal(A(i,:),B(j,:));
                end
                if any(Flag(i,:))
                    i_Flag = 1;
                    break;
                end
            end
        end
        % >> 6.4.) --------------------------------------------------------
        function [i_Flag] = fft_isequal_2(A,B)
            i_Flag = 0;
            for i = 1:size(A,1)
                for j = 1:size(B,1)
                    Flag(i,j) = isequal(A(i,:),B(j,:));
                    continue;
                end
            end
            if nnz(Flag) > 1
                i_Flag = 1;
            end
        end
        % >> 6.5.) --------------------------------------------------------
        function [Z] = fft_setdiff(X,Y)
            if ~isempty(X) && ~isempty(Y)
                check    = false(1, max(max(X), max(Y)));
                check(X) = true;
                check(Y) = false;
                Z        = X(check(X));
            else
                Z = X;
            end
        end
    end
end