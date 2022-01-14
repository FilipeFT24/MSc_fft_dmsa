classdef SubClass_1_3_2
    methods (Static)
        %% > SubClass_1_3_2.
        % >> --------------------------------------------------------------
        % >> 1.     Stencil setup.
        %  > 1.1.   Find vertex/face cell neighbours.
        %  > 1.2.   Add boundary faces/compute coordinates of stencil elements.
        %  > 1.2.1. Add boundary faces #1 (check faces to be added).
        %  > 1.2.2. Add boundary faces #2 (reshape cell array).
        %  > 1.3.   Compute stencil limits.
        %  > 1.4.   Perform stencil extension(s)/re-compute stencil limits.
        % >> --------------------------------------------------------------
  
        %% > 1. -----------------------------------------------------------
        function [msh] = Stencil_Setup(msh,Type,Nlev)
            % >> 1.
            %  > 1.1.
            msh = SubClass_1_3_2.Set_Neighbours(msh,Type,Nlev);
            %  > 1.2.
            msh = SubClass_1_3_2.Compute_StencilCoord(msh);
            %  > 1.3.
            msh = SubClass_1_3_2.Compute_Limits(msh,msh.s.xy_v);
        end
        % >> 1.1. ---------------------------------------------------------
        function [msh] = Set_Neighbours(msh,Type,Nlev)
            %  > Auxiliary arrays.
            for i = 1:size(msh.bnd.f,2)
                bnd_faces(i) = msh.bnd.f{2,i};
            end
            for j = 1:size(msh.bnd.c,2)
                bnd_cells(j) = msh.bnd.c{2,j};
            end
                       
            % >> Cells to be evaluated (stencil 1).
            %  > Evaluate stencil for neighbours of neighbouring cells of face i (purely for code efficiency).
            for i = 1:msh.f.NF
                j  = length(msh.f.cells{i});
                ki = length(msh.c.nb{msh.f.cells{i}(1)});
                if j == 1
                    nb{i}(1)       = msh.f.cells{i};
                    nb{i}(2:ki+1)  = msh.c.nb{msh.f.cells{i}};
                elseif j == 2
                    kj             = ki+length(msh.c.nb{msh.f.cells{i}(2)});
                    nb{i}(1:ki)    = msh.c.nb{msh.f.cells{i}(1)};
                    nb{i}(ki+1:kj) = msh.c.nb{msh.f.cells{i}(2)};
                end
                nb{i} = unique(nb{i});
            end
            
            % >> Stencils.
            for i = 1:msh.f.NF
                l = 1;
                % >> Level 1.
                for j = 1:length(nb{i})
                    %  > 'Direct' face neighbouring cells.
                    if SubClass_1_3_3.fft_ismember_1(nb{i}(j),msh.f.cells{i})
                        ngh {1,i}(l) = nb{i}(j);
                        x_ij{1,i}(j) = 1;
                        l            = l+1;
                    else
                        if SubClass_1_3_3.fft_isequal_1(msh.f.xy_v{i},msh.c.xy_v{nb{i}(j)})
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
                %  > Vertex type: Requires (at least) 1 common vertex (i.e. common vertex/face).
                %  > Cell   type: Requires 2 common vertices (i.e. common face).
                if Nlev > 1
                    for j = 2:Nlev
                        %  > Vertex stencil elements.
                        for k = 1:length(msh.s.st{j-1,i})
                            if k == 1
                                ngh{j,i} = msh.c.nb{msh.s.st{j-1,i}(k)};
                            else
                                ngh{j,i} = [ngh{j,i},msh.c.nb{msh.s.st{j-1,i}(k)}];
                            end
                        end
                        %  > Exclude repeated cells.
                        ngh{j,i} = unique(ngh{j,i});
                        
                        if strcmpi(Type,'Vertex')
                            %  > Do nothing...
                        elseif strcmpi(Type,'Face')
                            %  > Loop through previous stencil cells and select faces' index.
                            for k = 1:length(msh.s.st{j-1,i})
                                if k == 1
                                    prev_f{j-1,i} = msh.c.f.faces{msh.s.st{j-1,i}(k)};
                                else
                                    prev_f{j-1,i} = [prev_f{j-1,i},msh.c.f.faces{msh.s.st{j-1,i}(k)}];
                                end
                            end
                            %  > (Outer) faces of previous stencil.
                            prev_f{j-1,i} = unique(prev_f{j-1,i});
                            %  > Loop through 'ngh' cells and check whether cell k contains any element of 'prev_f'.
                            l = 'T';
                            for k = 1:length(ngh{j,i})
                                if l == 'T' && any(ismembc(msh.c.f.faces{ngh{j,i}(k)},prev_f{j-1,i}))
                                    v_ijk{j-1,i} = ngh{j,i}(k);
                                    l            = 'F';
                                else
                                    if any(ismembc(msh.c.f.faces{ngh{j,i}(k)},prev_f{j-1,i}))
                                        v_ijk{j-1,i} = [v_ijk{j-1,i},ngh{j,i}(k)];
                                    end
                                end
                            end
                            %  > Overwrite 'ngh'.
                            ngh{j,i} = v_ijk{j-1,i};
                        end
                        
                        %  > Exclude lower-order stencil cells.
                        if j == 2
                            exc{j-1,i} = msh.s.st{1,i};
                        else
                            exc{j-1,i} = [msh.s.st{j-1,i},exc{j-2,i}];
                        end
                        msh.s.st{j,i} = SubClass_1_3_3.fft_setdiff(ngh{j,i},exc{j-1,i});
                    end
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        function [msh] = Compute_StencilCoord(msh)
            % >> iD's.
            %  > Boundary cells indices.
            for i = 1:size(msh.bnd.c,2)
                bnd_cc(i) = msh.bnd.c{2,i};
            end
            %  > Boundary faces' cell index.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            
            % >> (Xv,Yv).
            for i = 1:msh.f.NF
                for j = 1:size(msh.s.st,1)
                    %  > Stencil cell centroids.
                    len_c{i}(j) = length(msh.s.st{j,i});
                    for k = 1:len_c{i}(j)
                        st_v{j,i}(1,k) = msh.c.mean(1,msh.s.st{j,i}(k));
                        st_v{j,i}(2,k) = msh.c.mean(2,msh.s.st{j,i}(k));
                    end
                    %  > Check whether face i's stencil cells belong to the boundary. If so, add the respective face to the stencil.
                    Flag{j,i} = zeros(1,len_c{i}(j));
                    Flag{j,i} = ismembc(msh.s.st{j,i},bnd_cc);
                    if any(Flag{j,i} == 1)
                        st_v{j,i} = SubClass_1_3_2.Add_bnd_1(Flag{j,i},st_v{j,i},bnd_ff,bnd_fc,len_c{i}(j),msh.s.st{j,i},msh.f.mean);
                    end
                    msh.s.xy_st{j,i} = st_v{j,i};
                end
            end
            %  > Reshape cell array.
            msh.s.xy_v = SubClass_1_3_2.Add_bnd_2(msh);            
        end
        %  > 1.2.1. -------------------------------------------------------
        function [st_v] = Add_bnd_1(Flag,st_v,bnd_ff,bnd_fc,len_c,st,mean_f)
            %  > Add respective boundary faces.
            %    Remark: A given boundary cell may contain more than 1 boundary face (see 3rd row of msh.bnd.f)
            j = 1;
            for i = 1:length(Flag)
                if Flag(i)
                    cf{j} = find(bnd_fc == st(i));
                    j     = j+1;
                end
            end
            %  > Faces to be added to the stencil.
            f_add = cell2mat(cf);
            for j = 1:length(f_add)
                st_v(1,j+len_c) = mean_f(1,bnd_ff(f_add(j)));
                st_v(2,j+len_c) = mean_f(2,bnd_ff(f_add(j)));
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [st_xy] = Add_bnd_2(msh)
            for i = 1:size(msh.s.st,2)
                for j = 1:size(msh.s.st,1)
                    st_xy{i}{j} = msh.s.xy_st{j,i};
                end
                st_xy{i} = cell2mat(st_xy{i});
            end
        end
        % >> 1.3. ---------------------------------------------------------
        function [msh] = Compute_Limits(msh,st_xy)                       
            for i = 1:msh.f.NF
                % >> Stencil limits.
                %  > (x,y)_min.
                msh.s.lim.x_min(i) = min(st_xy{i}(1,:));
                msh.s.lim.y_min(i) = min(st_xy{i}(2,:));
                %  > (x,y)_max.
                msh.s.lim.x_max(i) = max(st_xy{i}(1,:));
                msh.s.lim.y_max(i) = max(st_xy{i}(2,:));

                % >> Stencil dimensionless length.
                for j = 1:size(msh.s.st,1)
                    %  > Initialize local variables.
                    hx_ij(j,i) = 0;
                    hy_ij(j,i) = 0;
                    for k = 1:length(msh.s.st{j,i})
                        hx_ij(j,i) = hx_ij(j,i)+(max(msh.c.xy_v{msh.s.st{j,i}(k)}(:,1))-min(msh.c.xy_v{msh.s.st{j,i}(k)}(:,1)));
                        hy_ij(j,i) = hy_ij(j,i)+(max(msh.c.xy_v{msh.s.st{j,i}(k)}(:,2))-min(msh.c.xy_v{msh.s.st{j,i}(k)}(:,2)));
                    end
                    ni    (j,i) = length(msh.s.st{j,i});
                end
                nt          (i) = sum(ni(:,i));
                hx_i        (i) = sum(hx_ij(:,i));
                hy_i        (i) = sum(hy_ij(:,i));
                Ls_x        (i) = max(st_xy{i}(1,:))-min(st_xy{i}(1,:));
                Ls_y        (i) = max(st_xy{i}(2,:))-min(st_xy{i}(2,:));
                msh.s.hx_ref(i) = hx_i(i)./nt(i);
                msh.s.hy_ref(i) = hy_i(i)./nt(i);
                msh.s.nx    (i) = Ls_x(i)./msh.s.hx_ref(i);
                msh.s.ny    (i) = Ls_y(i)./msh.s.hy_ref(i);
            end            
        end
        % >> 1.4. ---------------------------------------------------------
        function [] = Perform_StencilExtension()
        end
        
        
    

        
%         % >> 1.4. ---------------------------------------------------------
%         function [add_to] = StencilExt_1(msh)
%             % >> Stencil indices.
%             for i = 1:size(msh.s.st,2)
%                 for j = 1:size(msh.s.st,1)
%                     st_i{i}{j} = msh.s.st{j,i};
%                 end
%                 st_i{i} = cell2mat(st_i{i});
%             end
%             % >> Add elements to stencil...
%             %  > Initialize.
%             add_to = cell(1,msh.f.NF);
%             for i = 1:msh.f.NF
%                 k = 1;
%                 for j = 1:msh.c.NC
%                     if (msh.c.mean(1,j) >= msh.s.lim.x_min(i) && msh.c.mean(1,j) <= msh.s.lim.x_max(i)) && ...
%                             (msh.c.mean(2,j) >= msh.s.lim.y_min(i) && msh.c.mean(2,j) <= msh.s.lim.y_max(i)) && ~ismembc(j,sort(st_i{i}))
%                         add_to{i}(k) = j;
%                         k = k+1;
%                     end
%                 end
%             end
%         end        
    end
end