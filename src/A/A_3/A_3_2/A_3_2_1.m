classdef A_3_2_1
    methods (Static)
        %% > A_3_2_1.
        % >> --------------------------------------------------------------
        % >> 1.     Set stencil cell/face indices.
        % >> 2.     Tools.
        %  > 2.1.   Check face index (i.e. check whether a given stencil cell is a boundary cell).
        %  > 2.2.   Reshape element array.
        %  > 2.3.   Compute stencil coordinates.
        %  > 2.3.1. Compute stencil coordinates (cells).
        %  > 2.3.2. Compute stencil coordinates (faces).
        % >> --------------------------------------------------------------
        
        %% > 1. -----------------------------------------------------------
        function [msh] = Set_Stencil(msh,bnd_cc,Type,NLay)
            % >> Cells to be evaluated (stencil 1).
            %  > Evaluate stencil for neighbours of neighbouring cells of face i (purely for code efficiency).
            for i = 1:msh.f.NF
                j  = length(msh.f.c{i});
                ki = length(msh.c.c{msh.f.c{i}(1)});
                if j == 1
                    nb{i}(1)       = msh.f.c{i};
                    nb{i}(2:ki+1)  = msh.c.c{msh.f.c{i}};
                elseif j == 2
                    kj             = ki+length(msh.c.c{msh.f.c{i}(2)});
                    nb{i}(1:ki)    = msh.c.c{msh.f.c{i}(1)};
                    nb{i}(ki+1:kj) = msh.c.c{msh.f.c{i}(2)};
                end
                nb{i} = unique(nb{i});
            end
            
            for i = 1:msh.f.NF
                % >> Stencil cell indices.
                l = 1;
                % >> Level 1.
                for j = 1:length(nb{i})
                    %  > 'Direct' face neighbouring cells.
                    if A_Tools.fft_ismember_1(nb{i}(j),msh.f.c{i})
                        ngh {1,i}(l) = nb{i}(j);
                        x_ij{1,i}(j) = 1;
                        l            = l+1;
                    else
                        if A_Tools.fft_isequal_1(msh.f.xy_v{i},msh.c.xy_v{nb{i}(j)})
                            ngh {1,i}(l) = nb{i}(j);
                            x_ij  {i}(j) = 1;
                            l            = l+1;
                        else
                            x_ij{i}(j) = 0;
                        end
                    end
                end
                msh.s.c{1,i} = ngh{1,i};
                
                % >> Level N.
                %  > Vertex type: Requires (at least) 1 common vertex (i.e. common vertex/face).
                %  > Cell   type: Requires 2 common vertices (i.e. common face).
                if strcmpi(Type,'Face')
                    %  > Initialize.
                    prev_f = cell(NLay-1,msh.f.NF);
                end
                if NLay > 1
                    for j = 2:NLay
                        %  > Vertex stencil elements.
                        for k = 1:length(msh.s.c{j-1,i})
                            if k == 1
                                ngh{j,i} = msh.c.c{msh.s.c{j-1,i}(k)};
                            else
                                ngh{j,i} = [ngh{j,i},msh.c.c{msh.s.c{j-1,i}(k)}];
                            end
                        end
                        %  > Exclude repeated cells.
                        ngh{j,i} = unique(ngh{j,i});
                        
                        if strcmpi(Type,'Vertex')
                            %  > Do nothing...
                        elseif strcmpi(Type,'Face')
                            %  > Loop through previous stencil cells and select faces' index.
                            for k = 1:length(msh.s.c{j-1,i})
                                prev_f{j-1,i} = [prev_f{j-1,i},msh.c.f.f{msh.s.c{j-1,i}(k)}];
                            end
                            %  > (Outer) faces of previous stencil.
                            prev_f{j-1,i} = unique(prev_f{j-1,i});
                            
                            %  > Loop through 'ngh' cells and check whether cell k contains any element of 'prev_f'.
                            m = true;
                            for k = 1:length(ngh{j,i})
                                if m && any(ismembc(msh.c.f.f{ngh{j,i}(k)},prev_f{j-1,i}))
                                    v_ijk{j-1,i} = ngh{j,i}(k);
                                    m            = false;
                                else
                                    if any(ismembc(msh.c.f.f{ngh{j,i}(k)},prev_f{j-1,i}))
                                        v_ijk{j-1,i} = [v_ijk{j-1,i},ngh{j,i}(k)];
                                    end
                                end
                            end
                            %  > Overwrite 'ngh'.
                            ngh{j,i} = v_ijk{j-1,i};
                        end
                        
                        %  > Exclude lower-order stencil cells.
                        if j == 2
                            exc{j-1,i} = msh.s.c{1,i};
                        else
                            exc{j-1,i} = [msh.s.c{j-1,i},exc{j-2,i}];
                        end
                        msh.s.c{j,i} = A_Tools.fft_setdiff(ngh{j,i},exc{j-1,i});
                    end
                end
                % >> Stencil face indices.
                for j = 1:size(msh.s.c,1)
                    %  > Initialize.
                    Flag {j,i}   = false(1,length(msh.s.c{j,i}));
                    msh.s.f(j,i) = cell (1,1);
                    
                    %  > Check whether face i's stencil cells belong to the boundary. If so, add the respective face to the stencil.
                    Flag {j,i}  = ismembc(msh.s.c{j,i},bnd_cc);
                    if any(Flag{j,i})
                        msh.s.f{j,i} = A_3_2_1.Add_Face(Flag{j,i},msh,msh.s.c{j,i});
                    end
                end
                
                % >> Stencil coordinates.
                %  > Cells.
                msh.s.xy_v_c{i} = A_3_2_1.Compute_Coordinates_c(msh.s.c(:,i),msh.c.mean);
                %  > Faces.
                if any(~cellfun(@isempty,msh.s.f(:,i)))
                    msh.s.xy_v_f{i} = A_3_2_1.Compute_Coordinates_f(msh.s.f(:,i),msh.f.mean);
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [add_f] = Add_Face(Flag,msh,c)
            %  > Auxiliary arrays.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            
            %  > Add respective boundary faces.
            %    Remark: A given boundary cell may contain more than 1 boundary face (see 3rd row of msh.bnd.f)
            j = 1;
            for i = 1:length(Flag)
                if Flag(i)
                    cf{j} = find(bnd_fc == c(i));
                    j     = j+1;
                end
            end
            %  > Faces to be added to the stencil.
            add_f = bnd_ff(cell2mat(cf));
        end
        % >> 2.2. ---------------------------------------------------------
        function [arr] = Deal_StencilElem(st)
            %  > Initialize.
            ijk = 1;
            
            %  > Cell(s).
            for i = 1:size(st,1)
                if isempty(st{i})
                    continue;
                else
                    if ijk == 1
                        arr = st{i};
                        ijk = ijk+1;
                    else
                        arr = [arr,st{i}];
                        ijk = ijk+1;
                    end
                end
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > 2.3.1. -------------------------------------------------------
        function [xy_v_c] = Compute_Coordinates_c(st_c,mean_c)
            %  > Deal elements.
            arr_c = A_3_2_1.Deal_StencilElem(st_c);
            
            len_c = length(arr_c);
            for i = 1:len_c
                xy_v_c(1,i) = mean_c(1,arr_c(i));
                xy_v_c(2,i) = mean_c(2,arr_c(i));
            end
        end
        %  > 2.3.2. -------------------------------------------------------
        function [xy_v_f] = Compute_Coordinates_f(st_f,mean_f)
            %  > Deal elements.
            arr_f = A_3_2_1.Deal_StencilElem(st_f);
            
            len_f = length(arr_f);
            for i = 1:len_f
                xy_v_f(1,i) = mean_f(1,arr_f(i));
                xy_v_f(2,i) = mean_f(2,arr_f(i));
            end
        end
        %  > 2.3.3. -------------------------------------------------------
        function [xy_v_t] = Compute_Coordinates_t(st_c,st_f)
            xy_v_t = cat(2,st_c,st_f);
        end
    end
end