classdef SubClass_1_3_2
    methods (Static)
        %% > SubClass_1_3_2.
        % >> --------------------------------------------------------------
        % >> 1.     Stencil setup.
        %  > 1.1.   Set stencil cell/face indices.
        %  > 1.1.1. Check face index (i.e. check whether a given stencil cell is a boundary cell).
        %  > 1.1.2. Compute stencil elements' coordinates.
        % >> 2.     Compute stencil limits and check if and extension is required.
        %  > 2.1.   Compute limits.
        %  > 2.2.   Perform extension.
        % >> --------------------------------------------------------------
        
        %% > 1. -----------------------------------------------------------
        function [msh] = Stencil_Setup(msh,Type,p,NLay)
            %  > Auxiliary arrays.
            for i = 1:size(msh.bnd.c,2)
                bnd_cc(i) = msh.bnd.c{2,i};
            end
            
            % >> 1.
            msh = SubClass_1_3_2.Set_Stencil(msh,bnd_cc,Type,NLay);
            % >> 2.
            msh = SubClass_1_3_2.Extend_Stencil(msh,bnd_cc,p);
        end
        % >> 1.1. ---------------------------------------------------------
        function [msh] = Set_Stencil(msh,bnd_cc,Type,NLay)
            % >> Cells to be evaluated (stencil 1).
            %  > Evaluate stencil for neighbours of neighbouring cells of face i (purely for code efficiency).
            for i = 1:msh.f.NF
                j  = length(msh.f.c{i});
                ki = length(msh.c.nb{msh.f.c{i}(1)});
                if j == 1
                    nb{i}(1)       = msh.f.c{i};
                    nb{i}(2:ki+1)  = msh.c.nb{msh.f.c{i}};
                elseif j == 2
                    kj             = ki+length(msh.c.nb{msh.f.c{i}(2)});
                    nb{i}(1:ki)    = msh.c.nb{msh.f.c{i}(1)};
                    nb{i}(ki+1:kj) = msh.c.nb{msh.f.c{i}(2)};
                end
                nb{i} = unique(nb{i});
            end
            
            for i = 1:msh.f.NF
                % >> Stencil cell indices.
                l = 1;
                % >> Level 1.
                for j = 1:length(nb{i})
                    %  > 'Direct' face neighbouring cells.
                    if SubClass_1_3_3.fft_ismember_1(nb{i}(j),msh.f.c{i})
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
                                ngh{j,i} = msh.c.nb{msh.s.c{j-1,i}(k)};
                            else
                                ngh{j,i} = [ngh{j,i},msh.c.nb{msh.s.c{j-1,i}(k)}];
                            end
                        end
                        %  > Exclude repeated cells.
                        ngh{j,i} = unique(ngh{j,i});
                        
                        if strcmpi(Type,'Vertex')
                            %  > Do nothing...
                        elseif strcmpi(Type,'Face')
                            %  > Loop through previous stencil cells and select faces' index.
                            for k = 1:length(msh.s.c{j-1,i})
                                prev_f{j-1,i} = [prev_f{j-1,i},msh.c.f.faces{msh.s.c{j-1,i}(k)}];
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
                            exc{j-1,i} = msh.s.c{1,i};
                        else
                            exc{j-1,i} = [msh.s.c{j-1,i},exc{j-2,i}];
                        end
                        msh.s.c{j,i} = SubClass_1_3_3.fft_setdiff(ngh{j,i},exc{j-1,i});
                    end
                end
                % >> Stencil face indices.
                for j = 1:size(msh.s.c,1)
                    %  > Initialize.
                    Flag {j,i}  = zeros(1,length(msh.s.c{j,i}));
                    %  > Check whether face i's stencil cells belong to the boundary. If so, add the respective face to the stencil.
                    Flag {j,i}  = ismembc(msh.s.c{j,i},bnd_cc);
                    if any(Flag{j,i} == 1)
                        msh.s.f{j,i} = SubClass_1_3_2.Add_Face(Flag{j,i},msh,msh.s.c{j,i});
                    end
                end
                % >> Stencil coordinates.
                msh.s.xy_v_t{i} = SubClass_1_3_2.Compute_Coordinates(msh,msh.s.c(:,i),msh.s.f(:,i));
            end
        end
        %  > 1.1.1. -------------------------------------------------------
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
        %  > 1.1.2. --------------------------------------------------------
        function [arr_c,arr_f] = Deal_StencilElem(st_c,st_f)
            %  > Initialize.
            ijk_c = 1;
            ijk_f = 1;
            %  > Cell(s).
            for i = 1:size(st_c,1)
                if isempty(st_c{i})
                    continue;
                else
                    if ijk_c == 1
                        arr_c = st_c{i};
                        ijk_c = ijk_c+1;
                    else
                        arr_c = [arr_c,st_c{i}];
                        ijk_c = ijk_c+1;
                    end
                end
            end
            %  > Face(s).
            for i = 1:size(st_f,1)
                if isempty(st_f{i})
                    continue;
                else
                    if ijk_f == 1
                        arr_f = st_f{i};
                        ijk_f = ijk_f+1;
                    else
                        arr_f = [arr_f,st_f{i}];
                        ijk_f = ijk_c+1;
                    end
                end
            end
        end
        %  > 1.1.3. -------------------------------------------------------
        function [xy_v] = Compute_Coordinates(msh,st_c,st_f)
            %  > Deal elements.
            [arr_c,arr_f] = SubClass_1_3_2.Deal_StencilElem(st_c,st_f);
            
            %  > Initialize.
            len_c = length(arr_c);
            len_f = length(arr_f);
            %  > Cell(s).
            for i = 1:len_c
                xy_v(1,i) = msh.c.mean(1,arr_c(i));
                xy_v(2,i) = msh.c.mean(2,arr_c(i));
            end
            %  > Face(s).
            for i = len_c+1:len_c+len_f
                xy_v(1,i) = msh.f.mean(1,arr_f(i-len_c));
                xy_v(2,i) = msh.f.mean(2,arr_f(i-len_c));
            end
        end
        
        %% > 2. -----------------------------------------------------------
        function [msh] = Extend_Stencil(msh,bnd_cc,p)
            % >> 2.1.
            par = SubClass_1_3_2.Compute_Parameters(msh,msh.s.c,msh.s.f);
            % >> 2.2.
            %  > Initialize.
            msh.s.par.n_e = zeros(2,msh.f.NF);
            msh.s.c_e     = cell (1,msh.f.NF);
            msh.s.c_f     = cell (1,msh.f.NF);
            %  > Extend stencil until...
            Cont = 1;
            while (any(par.n_x) < p-1./2 || any(par.n_y) < p-1./2) && Cont
                [msh,Cont] = SubClass_1_3_2.Perform_Extension (msh,bnd_cc,par,p);
                par        = SubClass_1_3_2.Compute_Parameters(msh,msh.s.c,msh.s.f);
            end
            %msh.s.par = par;
        end
        % >> 2.1. ---------------------------------------------------------
        function [par] = Compute_Parameters(msh,st_c,st_f)
            for i = 1:msh.f.NF
                % >> Deal elements.
                [arr_c{i},arr_f{i}] = SubClass_1_3_2.Deal_StencilElem(st_c(:,i),st_f(:,i));
                
                % >> Compute/re-compute 'par' fields.
                %  > (x,y)_min.
                par.l_x(1,i) = min(msh.s.xy_v_t{i}(1,:));
                par.l_y(1,i) = min(msh.s.xy_v_t{i}(2,:));
                %  > (x,y)_max.
                par.l_x(2,i) = max(msh.s.xy_v_t{i}(1,:));
                par.l_y(2,i) = max(msh.s.xy_v_t{i}(2,:));
                %  > (h_x,h_y).
                for j = 1:length(arr_c{i})
                    h_x(i,j) = max(msh.c.xy_v{arr_c{i}(j)}(:,1))-min(msh.c.xy_v{arr_c{i}(j)}(:,1));
                    h_y(i,j) = max(msh.c.xy_v{arr_c{i}(j)}(:,2))-min(msh.c.xy_v{arr_c{i}(j)}(:,2));
                end
                par.h_x(i)   = sum(h_x(i,:))./length(arr_c{i});
                par.h_y(i)   = sum(h_y(i,:))./length(arr_c{i});
                %  > (n_x,n_y).
                par.n_x(i) = (par.l_x(2,i)-par.l_x(1,i))./par.h_x(i);
                par.n_y(i) = (par.l_y(2,i)-par.l_y(1,i))./par.h_y(i);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [msh,Cont] = Perform_Extension(msh,bnd_cc,par,p)
            %  > Initialize.
            Fl_id = zeros(msh.f.NF,2);
            Cont  = 0;
            k     = 0;
            len   = size(msh.s.c,1);
                        
            for i = 1:msh.f.NF
                %  > Stencil elements.
                if ~isempty(msh.s.c{len,i})
                    if par.n_x(i) < p-1./2 || par.n_y(i) < p-1./2
                        k   = k+1;
                        Lst = msh.s.c{len,i};
                        for j = 1:len
                            st{k,j} = msh.s.c{j,i};
                        end
                        st{k} = cell2mat(st(k,:));
                    end
                    if par.n_x(i) < p-1./2
                        % >> Add cell indices.
                        %  > Stencil limits.
                        y_min = par.l_y(1,i);
                        y_max = par.l_y(2,i);
                        Ext_X = SubClass_1_3_2.Extension_1(msh,Lst,st{k},'x',y_min,y_max);
                        %  > if none added, continue...
                        if isempty(Ext_X)
                            continue;
                        end
                        
                        % >> Update/add...
                        %  > Number of extensions (y-direction).
                        msh.s.par.n_e(1,i) = msh.s.par.n_e(1,i)+1;
                        %  > Cell indices.
                        msh.s.c{len+1,i} = Ext_X;
                        msh.s.c_e    {i} = [msh.s.c_e{i},Ext_X];
                        %  > Face indices.
                        Flag = ismembc(Ext_X,bnd_cc) == 1;
                        if any(Flag)
                            msh.s.c_f{i} = SubClass_1_3_2.Add_Face(Flag,msh,Ext_X);
                        end
                        Fl_id(i,1) = length(Ext_X);
                    end
                    if par.n_y(i) < p-1./2
                        % >> Add cell indices.
                        % > Stencil limits.
                        x_min = par.l_x(1,i);
                        x_max = par.l_x(2,i);
                        Ext_Y = SubClass_1_3_2.Extension_1(msh,Lst,st{k},'y',x_min,x_max);
                        %  > if none added, continue...
                        if isempty(Ext_Y)
                            continue;
                        end
                        
                        % >> Update/add...
                        %  > Number of extensions (y-direction).
                        msh.s.par.n_e(2,i) = msh.s.par.n_e(2,i)+1;
                        %  > Cell indices.
                        msh.s.c{len+1,i} = Ext_Y;
                        msh.s.c_e    {i} = [msh.s.c_e{i},Ext_Y];
                        %  > Face indices.
                        Flag = ismembc(Ext_Y,bnd_cc) == 1;
                        if any(Flag)
                            msh.s.c_f{i} = SubClass_1_3_2.Add_Face(Flag,msh,Ext_Y);
                        end
                        Fl_id(i,2) = length(Ext_Y);
                    end
                end
                msh.s.xy_v_t{i} = SubClass_1_3_2.Compute_Coordinates(msh,msh.s.c(:,i),msh.s.f(:,i));
            end           
            if any(any(Fl_id))
                Cont = 1;
            end
        end
        %  > 1.4.1. -------------------------------------------------------
        function [add_to] = Extension_1(msh,Lst,st,Dir,v_min,v_max)
            % >> Outer cells.
            %  > Stencil (outer) cells.
            outer_c = Lst;
            %  > ...neighbours.
            for i = 1:length(outer_c)
                nb_c_out{i} = msh.c.nb{outer_c(i)};
            end
            %  > Outer cells (NOT in the stencil).
            nb_diff = setdiff(unique(cell2mat(nb_c_out)),st);
            
            %  > Select direction.
            if strcmpi(Dir,'x')
                k = 2;
            elseif strcmpi(Dir,'y')
                k = 1;
            end
            %  > Select cells within stencil limits.
            j = 1;
            for i = 1:length(nb_diff)
                if msh.c.mean(k,nb_diff(i)) >= v_min && msh.c.mean(k,nb_diff(i)) <= v_max
                    add_to(j) = nb_diff(i);
                    j         = j+1;
                end
            end
            if j == 1
                add_to = [];
            end
        end
        %  > 1.4.2. -------------------------------------------------------
        function [msh] = Extension_2(msh)
        end
    end
end