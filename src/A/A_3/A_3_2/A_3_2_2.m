classdef A_3_2_2
    methods (Static)
        %% > A_3_2_2.
        % >> --------------------------------------------------------------
        % >> 1.     Compute stencil limits and check if and extension is required.
        %  > 1.1.   Compute parameters.
        %  > 1.2.   Perform extension.
        %  > 1.2.1. Extension #1.
        %  > 1.2.2. Extension #2.
        % >> --------------------------------------------------------------
        function [msh] = Extend_Stencil(msh,bnd_cc,Type,p,et)
            %  > Initialize.
            msh.s.par.n_e = zeros(2,msh.f.NF);
            msh.s.c_e     = cell (1,msh.f.NF);
            msh.s.f_e     = cell (1,msh.f.NF);
            
            %  > hg_x,hg_y.
            [hg_x,hg_y] = A_3_2_2.Compute_hgx_hgy(msh);
            
            % >> 1.
            %  > 1.1.
            for i = 1:msh.f.NF
                [par.n_x(i),par.n_y(i),par.ng_x(i),par.ng_y(i),par.l_x(:,i),par.l_y(:,i)] = ...
                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),hg_x,hg_y);
            end
            %  > 1.2.
            %  > Extend stencil until it is incomplete and domain limits have not been reached...
            for i = 1:msh.f.NF
                %  > Initialize.
                Continue = true;
                
                %  > Check...
                if (~et && single(par.n_x(i)) >= p-1./2 && single(par.n_y(i)) >= p-1./2) || ...
                        (et && single(par.n_x(i)) >= p-1./2 && single(par.n_y(i)) >= p-1./2 && single(par.ng_x(i)) >= p-1./2 && single(par.ng_y(i)) >= p-1./2)
                    continue;
                else
                    while ((~et && (single(par.n_x(i)) < p-1./2 || single(par.n_y(i)) < p-1./2)) || ...
                            (et && (single(par.n_x(i)) < p-1./2 || single(par.n_y(i)) < p-1./2 || single(par.ng_x(i)) < p-1./2 || single(par.ng_y(i)) < p-1./2))) && Continue
                        %% > x-direction.
                        Continue_X = false;
                        if (~et && single(par.n_x(i)) < p-1./2) || ...
                                (et && (single(par.n_x(i)) < p-1./2 || single(par.ng_x(i)) < p-1./2))
                            %  > Add/update...
                            [msh,Continue_X] = ...
                                A_3_2_2.Perform_Extension(i,msh,bnd_cc,'x',par.l_y(1,i),par.l_y(2,i),Type);
                            [par.n_x(i),par.n_y(i),par.ng_x(i),par.ng_y(i),par.l_x(:,i),par.l_y(:,i)] = ...
                                A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),hg_x,hg_y);
                            %  > Number of extensions (x-direction).
                            msh.s.par.n_e(1,i) = msh.s.par.n_e(1,i)+1;
                        end
                        %% > y-direction.
                        Continue_Y = false;
                        if (~et && single(par.n_y(i)) < p-1./2) || ...
                                (et && (single(par.n_y(i)) < p-1./2 || single(par.ng_y(i)) < p-1./2))
                            %  > Add/update...
                            [msh,Continue_Y] = ...
                                A_3_2_2.Perform_Extension(i,msh,bnd_cc,'y',par.l_x(1,i),par.l_x(2,i),Type);
                            [par.n_x(i),par.n_y(i),par.ng_x(i),par.ng_y(i),par.l_x(:,i),par.l_y(:,i)] = ...
                                A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),hg_x,hg_y);
                            %  > Number of extensions (y-direction).
                            msh.s.par.n_e(2,i) = msh.s.par.n_e(2,i)+1;
                        end
                        Continue = any([Continue_X,Continue_Y]);
                    end
                end
            end
            %  > Deal fields...
            msh.s.par.n_x  = par.n_x;
            msh.s.par.n_y  = par.n_y;
            msh.s.par.ng_x = par.ng_x;
            msh.s.par.ng_y = par.ng_y;
            msh.s.par.l_x  = par.l_x;
            msh.s.par.l_y  = par.l_y;
            %  > Total stencil points.
            for i = 1:msh.f.NF
                msh.s.xy_v_t{i} = A_3_2_1.Compute_Coordinates_t(msh.s.xy_v_c{i},msh.s.xy_v_f{i});
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [n_x,n_y,ng_x,ng_y,l_x,l_y] = Compute_Parameters(i,msh,st_c,hg_x,hg_y)
            % >> Deal elements.
            arr_c = A_3_2_1.Deal_StencilElem(st_c);
            
            % >> Compute/re-compute 'par' fields.
            %  > (x,y)_min, (x,y)_max.
            [lc_x(1),lc_x(2)] = MinMaxElem(msh.s.xy_v_c{i}(1,:));
            [lc_y(1),lc_y(2)] = MinMaxElem(msh.s.xy_v_c{i}(2,:));
            %  > .. for face w/o boundary faces.
            if ~isempty(msh.s.xy_v_f{i})
                [lf_x(1),lf_x(2)] = MinMaxElem(msh.s.xy_v_f{i}(1,:));
                [lf_y(1),lf_y(2)] = MinMaxElem(msh.s.xy_v_f{i}(2,:));
                [ l_x(1),~      ] = MinMaxElem([lc_x(1),lf_x(1)]);
                [      ~, l_x(2)] = MinMaxElem([lc_x(2),lf_x(2)]);
                [ l_y(1),~      ] = MinMaxElem([lc_y(1),lf_y(1)]);
                [      ~, l_y(2)] = MinMaxElem([lc_y(2),lf_y(2)]);
            else
                l_x(1) = lc_x(1);
                l_y(1) = lc_y(1);
                l_x(2) = lc_x(2);
                l_y(2) = lc_y(2);
            end
            %  > (n_x,n_y).
            len_c = length(arr_c);
            for j = 1:len_c
                [A_x(j),B_x(j)] = MinMaxElem(msh.c.xy_v{arr_c(j)}(:,1));
                [A_y(j),B_y(j)] = MinMaxElem(msh.c.xy_v{arr_c(j)}(:,2));
                h_x(j)          = B_x(j)-A_x(j);
                h_y(j)          = B_y(j)-A_y(j);
            end
            h_x = sum(h_x)./len_c;
            h_y = sum(h_y)./len_c;
            n_x = (l_x(2)-l_x(1))./h_x;
            n_y = (l_y(2)-l_y(1))./h_y;
            %  > (ng_x,ng_y).
            ng_x = (l_x(2)-l_x(1))./hg_x;
            ng_y = (l_y(2)-l_y(1))./hg_y;
        end
        % >> 1.1.1. -------------------------------------------------------
        function [hg_x,hg_y] = Compute_hgx_hgy(msh)
            for j = 1:msh.c.NC
                [C_x(j),D_x(j)] = MinMaxElem(msh.c.xy_v{j}(:,1));
                [C_y(j),D_y(j)] = MinMaxElem(msh.c.xy_v{j}(:,2));
                hg_x(j)         = D_x(j)-C_x(j);
                hg_y(j)         = D_y(j)-C_y(j);
            end
            hg_x = sum(hg_x)./msh.c.NC;
            hg_y = sum(hg_y)./msh.c.NC;
        end
        
        % >> 1.2. ---------------------------------------------------------
        function [msh,Flag] = Perform_Extension(i,msh,bnd_cc,Dir,v_min,v_max,Type)
            %  > Initialize.
            Flag = false;
            len  = nnz(~cellfun(@isempty,msh.s.c(:,i)));
            
            %  > Previous layers' stencil elements.
            for j = 1:len
                st_el{j} = msh.s.c{j,i};
            end
            st_el = cell2mat(st_el);
            %  > Select direction.
            if strcmpi(Dir,'x')
                Add = A_3_2_2.Extension_1(msh,st_el,Dir,v_min,v_max,Type);
            elseif strcmpi(Dir,'y')
                Add = A_3_2_2.Extension_1(msh,st_el,Dir,v_min,v_max,Type);
            end
            
            % >> Update/add...
            %  > Cell indices.
            msh.s.c{len+1,i} = Add;
            msh.s.c_e    {i} = [msh.s.c_e{i},Add];
            %  > Face indices.
            i_face = ismembc(Add,bnd_cc) == true;
            if any(i_face)
                msh.s.f_e{i}       = A_3_2_1.Add_Face(i_face,msh,Add);
                msh.s.f  {len+1,i} = msh.s.f_e{i};
            end
            Flag = ~isempty(Add);
            
            if Flag
                %  > Cells.
                msh.s.xy_v_c{i} = A_3_2_1.Compute_Coordinates_c(msh.s.c(:,i),msh.c.mean);
                %  > Faces.
                if any(~cellfun(@isempty,msh.s.f(:,i)))
                    msh.s.xy_v_f{i} = A_3_2_1.Compute_Coordinates_f(msh.s.f(:,i),msh.f.mean);
                end
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        function [add_to] = Extension_1(msh,st_el,Dir,v_min,v_max,Type)
            % >> Stencil cell neighbours.
            for i = 1:length(st_el)
                nb_c{i} = msh.c.c{st_el(i)};
            end
            nb_c    = cell2mat(nb_c);
            nb_c_un = unique(nb_c);
            %  > Outer cells (NOT in the stencil).
            nb_diff = setdiff(nb_c_un,st_el);
            
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
                    add_to_v(j) = nb_diff(i);
                    j         = j+1;
                end
            end
            if strcmpi(Type,'Vertex')
                %  > Do nothing...
            elseif strcmpi(Type,'Face')
                if j == 1
                    %  > Continue...
                else
                    %  > Outer cell layer's faces.
                    for i = 1:length(st_el)
                        out_f{i} = msh.c.f.f{st_el(i)};
                    end
                    out_f = unique(cell2mat(out_f));
                    
                    %  > Add cell to extended stencil if it has a common face with the previous stencil layer.
                    k = 1;
                    for i = 1:length(add_to_v)
                        if any(ismembc(msh.c.f.f{add_to_v(i)},out_f))
                            add_to(k) = add_to_v(i);
                            k         = k+1;
                        end
                    end
                end
            end
                        
            %  > No added elements...
            if j == 1 || (k == 1 && strcmpi(Type,'Face'))
                add_to = [];
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [msh] = Extension_2(msh)
        end
    end
end