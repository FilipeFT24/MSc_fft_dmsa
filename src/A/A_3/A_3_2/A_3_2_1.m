classdef A_3_2_1
    methods (Static)
        %% > A_3_2_1.
        % >> --------------------------------------------------------------
        % >> 1.     Set stencil cell/face indices.
        % >> 2.     Tools.  
        %  > 2.1.   Check face index (i.e. check whether a given stencil cell is a boundary cell).
        %  > 2.2.   Reshape element array.
        %  > 2.3.   Compute stencil coordinates.
        %  > 2.3.1. Compute stencil (x,y) coordinates (cells/faces).
        %  > 2.3.2. Compute stencil (x,y) coordinates (total).
        % >> --------------------------------------------------------------
        
        %% > 1. -----------------------------------------------------------
        function [msh] = Set_Stencil(msh,bnd_cc,bnd_ff,bnd_fc,nt,nl)
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
                nb{i} = A_Tools.fft_unique(sort(nb{i}));
            end
                     
            for i = 1:msh.f.NF
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
                if nl > 1
                    for j = 2:nl
                        % >> Vertex stencil...
                        %  > ...cell neighbours.
                        ngh{j,i} = A_Tools.fft_unique(sort([msh.c.c{msh.s.c{j-1,i}}]));
                        %  > ...outer layer cell neighbours.
                        if j == 2
                            ngh{j,i} = A_Tools.fft_setdiff(ngh{j,i},msh.s.c{1,i});
                        else
                            ngh{j,i} = A_Tools.fft_setdiff(ngh{j,i},cat(2,msh.s.c{:,i})); 
                        end
                        
                        % >> Face stencil...
                        if nt
                            %  > ...outer cell layer's faces.
                            out_f{j-1,i} = A_Tools.fft_unique(sort([msh.c.f.f{ngh{j-1,i}}]));
                            
                            l = 1;
                            for k = 1:length(ngh{j,i})
                                if any(ismembc(msh.c.f.f{ngh{j,i}(k)},out_f{j-1,i}))
                                    %  >...add cell to the stencil.
                                    out_c{j-1,i}(l) = ngh{j,i}(k);
                                    l               = l+1;
                                end
                            end
                            %  > Overwrite...
                            ngh{j,i} = out_c{j-1,i};
                        end
                        %  > Deal elements...
                        msh.s.c{j,i} = ngh{j,i};
                    end
                end
                % >> Stencil face indices.
                for j = 1:size(msh.s.c,1)
                    %  > Initialize.
                    Flag   {j,i} = false(1,length(msh.s.c{j,i}));
                    msh.s.f(j,i) = cell (1,1);
                    
                    %  > Check whether face i's stencil cells belong to the boundary. If so, add the respective face to the stencil.
                    Flag{j,i} = ismembc(msh.s.c{j,i},bnd_cc);
                    if any(Flag{j,i})
                        msh.s.f{j,i} = A_3_2_1.Add_Face(Flag{j,i},msh.s.c{j,i},bnd_ff,bnd_fc);
                    end
                end
                
                % >> Stencil coordinates.
                %  > Cells.
                msh.s.xy_v_c{i} = A_3_2_1.Compute_Coordinates_cf(msh.s.c(:,i),msh.c.mean);
                %  > Faces.
                if any(~cellfun(@isempty,msh.s.f(:,i)))
                    msh.s.xy_v_f{i} = A_3_2_1.Compute_Coordinates_cf(msh.s.f(:,i),msh.f.mean);
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------   
        % >> 2.1. ---------------------------------------------------------
        function [add_f] = Add_Face(Flag,c,bnd_ff,bnd_fc)
            %  > Add respective boundary faces.
            %    Remark: A given boundary cell may contain more than 1 boundary face (see 3rd row of msh.bnd.f)
            j = 1;
            for i = 1:length(Flag)
                if Flag(i)
                    cf{j} = find(bnd_fc == c(i));
                    j     = j+1;
                end
            end
            add_f = bnd_ff([cf{:}]);
        end
        % >> 2.2. ---------------------------------------------------------
        function [arr] = Deal_StencilElem(st)
            %  > Initialize.
            ijk = 1;

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
        function [xy] = Compute_Coordinates_cf(st,mean_cf)
            arr_c     = A_3_2_1.Deal_StencilElem(st);            
            i         = 1:length(arr_c);
            xy(:,i)   = mean_cf(:,arr_c(i));
        end
        %  > 2.3.2. -------------------------------------------------------
        function [xy] = Compute_Coordinates_tt(st_c,st_f)
            xy = [st_c,st_f];
        end
    end
end