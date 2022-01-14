classdef SubClass_1_3_2
    methods (Static)
        %% > SubClass_1_3_2.
        % >> --------------------------------------------------------------
        % >> 1.   Stencil setup.
        %  > 1.1. Set stencil neighbours.
        %  > 1.2. Set stencil limits        (called on SubClass_2_2).
        %  > 1.3. Perform stencil extension (called on SubClass_2_2).
        % >> --------------------------------------------------------------
  
        %% > 1.) ----------------------------------------------------------
        % >> 1.1.) --------------------------------------------------------
        function [msh] = Set_FaceNeighbours(msh,Type,Nlev)
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
                msh.s.st_i{1,i} = ngh{1,i};
                
                % >> Level N.
                %  > Vertex type: Requires (at least) 1 common vertex (i.e. common vertex/face).
                %  > Cell   type: Requires 2 common vertices (i.e. common face).
                if Nlev > 1
                    for j = 2:Nlev
                        %  > Vertex stencil elements.
                        for k = 1:length(msh.s.st_i{j-1,i})
                            if k == 1
                                ngh{j,i} = msh.c.nb{msh.s.st_i{j-1,i}(k)};
                            else
                                ngh{j,i} = [ngh{j,i},msh.c.nb{msh.s.st_i{j-1,i}(k)}];
                            end
                        end
                        %  > Exclude repeated cells.
                        ngh{j,i} = unique(ngh{j,i});
                        
                        if strcmpi(Type,'Vertex')
                            %  > Do nothing...
                        elseif strcmpi(Type,'Face')
                            %  > Loop through previous stencil cells and select faces' index.
                            for k = 1:length(msh.s.st_i{j-1,i})
                                if k == 1
                                    prev_f{j-1,i} = msh.c.f.faces{msh.s.st_i{j-1,i}(k)};
                                else
                                    prev_f{j-1,i} = [prev_f{j-1,i},msh.c.f.faces{msh.s.st_i{j-1,i}(k)}];
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
                            exc{j-1,i} = msh.s.st_i{1,i};
                        else
                            exc{j-1,i} = [msh.s.st_i{j-1,i},exc{j-2,i}];
                        end
                        msh.s.st_i{j,i} = SubClass_1_3_3.fft_setdiff(ngh{j,i},exc{j-1,i});
                    end
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        
        % >> 1.3. ---------------------------------------------------------
        function [msh] = Set_Limits(msh,st_xy)                       
            for i = 1:msh.f.NF
                for j = 1:size(msh.s.st_i,1)
                    %  > Number of cells/per stencil layer.
                    ni(j,i) = length(msh.s.st_i{j,i});
                end
                %  > Number of cells/per stencil.
                nt(i) = sum(ni(:,i));
            end
            
            for i = 1:msh.f.NF
                % >> Stencil limits.
                %  > (x,y)_min.
                msh.s.lim.x_min(i) = min(st_xy{i}(1,:));
                msh.s.lim.y_min(i) = min(st_xy{i}(2,:));
                %  > (x,y)_max.
                msh.s.lim.x_max(i) = max(st_xy{i}(1,:));
                msh.s.lim.y_max(i) = max(st_xy{i}(2,:));

                % >> Stencil dimensionless length.
                %  > Initialize.
                for j = 1:size(msh.s.st_i,1)
                    hx_ij(j,i) = 0;
                    hy_ij(j,i) = 0;
                    for k = 1:length(msh.s.st_i{j,i})
                        hx_ij(j,i) = hx_ij(j,i)+(max(msh.c.xy_v{msh.s.st_i{j,i}(k)}(:,1))-min(msh.c.xy_v{msh.s.st_i{j,i}(k)}(:,1)));
                        hy_ij(j,i) = hy_ij(j,i)+(max(msh.c.xy_v{msh.s.st_i{j,i}(k)}(:,2))-min(msh.c.xy_v{msh.s.st_i{j,i}(k)}(:,2)));
                    end
                end
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
        function [add_to] = StencilExt_1(msh)
            % >> Stencil indices.
            for i = 1:size(msh.s.st_i,2)
                for j = 1:size(msh.s.st_i,1)
                    st_i{i}{j} = msh.s.st_i{j,i};
                end
                st_i{i} = cell2mat(st_i{i});
            end
            % >> Add elements to stencil...
            %  > Initialize.
            add_to = cell(1,msh.f.NF);
            for i = 1:msh.f.NF
                k = 1;
                for j = 1:msh.c.NC
                    if (msh.c.mean(1,j) >= msh.s.lim.x_min(i) && msh.c.mean(1,j) <= msh.s.lim.x_max(i)) && ...
                            (msh.c.mean(2,j) >= msh.s.lim.y_min(i) && msh.c.mean(2,j) <= msh.s.lim.y_max(i)) && ~ismembc(j,sort(st_i{i}))
                        add_to{i}(k) = j;
                        k = k+1;
                    end
                end
            end
        end        
    end
end