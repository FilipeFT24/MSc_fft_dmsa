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
        function [msh] = Extend_Stencil(msh,bnd_cc,p)
            %  > Initialize.
            msh.s.par.n_e = zeros(2,msh.f.NF);
            msh.s.c_e     = cell (1,msh.f.NF);
            msh.s.c_f     = cell (1,msh.f.NF);
            
            % >> 1.
            %  > 1.1.
            for i = 1:msh.f.NF
                [par.n_x(i),par.n_y(i),par.l_x(:,i),par.l_y(:,i)] = ...
                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i));
            end
            %  > 1.2.
            %  > Extend stencil until it is incomplete and domain limits have not been reached...
            for i = 1:msh.f.NF
                %  > Initialize.
                Continue = true;
                
                %  > Check...
                if par.n_x(i) >= p-1./2 && par.n_y(i) >= p-1./2
                    continue;
                else
                    while (par.n_x(i) < p-1./2 || par.n_y(i) < p-1./2) && Continue
                        %% > x-direction.
                        Continue_X = false;
                        if par.n_x(i) < p-1./2
                            %  > Add/update...
                            [msh,Continue_X] = ...
                                A_3_2_2.Perform_Extension(i,msh,bnd_cc,'x',par.l_y(1,i),par.l_y(2,i));
                            [par.n_x(i),par.n_y(i),par.l_x(:,i),par.l_y(:,i)] = ...
                                A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i));
                            %  > Number of extensions (x-direction).
                            msh.s.par.n_e(1,i) = msh.s.par.n_e(1,i)+1;
                        end
                        %% > y-direction.
                        Continue_Y = false;
                        if par.n_y(i) < p-1./2
                            %  > Add/update...
                            [msh,Continue_Y] = ...
                                A_3_2_2.Perform_Extension(i,msh,bnd_cc,'y',par.l_x(1,i),par.l_x(2,i));
                            [par.n_x(i),par.n_y(i),par.l_x(:,i),par.l_y(:,i)] = ...
                                A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i));
                            %  > Number of extensions (y-direction).
                            msh.s.par.n_e(2,i) = msh.s.par.n_e(2,i)+1;
                        end
                        Continue = any([Continue_X,Continue_Y]);
                    end
                end
            end
            %  > Deal fields...
            msh.s.par.l_x = par.l_x;
            msh.s.par.l_y = par.l_y;
            msh.s.par.n_x = par.n_x;
            msh.s.par.n_y = par.n_y;
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [n_x,n_y,l_x,l_y] = Compute_Parameters(i,msh,st_c)
            % >> Deal elements.
            arr_c = A_3_2_1.Deal_StencilElem(st_c);
            
            % >> Compute/re-compute 'par' fields.
            %  > (x,y)_min.
            l_x(1) = min(msh.s.xy_v_t{i}(1,:));
            l_y(1) = min(msh.s.xy_v_t{i}(2,:));
            %  > (x,y)_max.
            l_x(2) = max(msh.s.xy_v_t{i}(1,:));
            l_y(2) = max(msh.s.xy_v_t{i}(2,:));
            %  > (n_x,n_y).
            len_c = length(arr_c);
            for j = 1:len_c
                h_x(j) = max(msh.c.xy_v{arr_c(j)}(:,1))-min(msh.c.xy_v{arr_c(j)}(:,1));
                h_y(j) = max(msh.c.xy_v{arr_c(j)}(:,2))-min(msh.c.xy_v{arr_c(j)}(:,2));
            end
            h_x = sum(h_x)./len_c;
            h_y = sum(h_y)./len_c;
            n_x = (l_x(2)-l_x(1))./h_x;
            n_y = (l_y(2)-l_y(1))./h_y;
        end
        % >> 2.2. ---------------------------------------------------------
        function [msh,Flag] = Perform_Extension(i,msh,bnd_cc,Dir,v_min,v_max)
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
                Add = A_3_2_2.Extension_1(msh,st_el,Dir,v_min,v_max);
            elseif strcmpi(Dir,'y')
                Add = A_3_2_2.Extension_1(msh,st_el,Dir,v_min,v_max);
            end
            
            % >> Update/add...
            %  > Cell indices.
            msh.s.c{len+1,i} = Add;
            msh.s.c_e    {i} = [msh.s.c_e{i},Add];
            %  > Face indices.
            i_face = ismembc(Add,bnd_cc) == true;
            if any(i_face)
                msh.s.c_f{i}       = A_3_2_1.Add_Face(i_face,msh,Add);
                msh.s.f  {len+1,i} = msh.s.c_f{i};
            end
            Flag = ~isempty(Add);
            
            if Flag
                msh.s.xy_v_t{i} = A_3_2_1.Compute_Coordinates(msh,msh.s.c(:,i),msh.s.f(:,i));
            end
        end
        %  > 2.2.1. -------------------------------------------------------
        function [add_to] = Extension_1(msh,st_el,Dir,v_min,v_max)
            % >> Stencil cell neighbours.
            for i = 1:length(st_el)
                nb_c{i} = msh.c.nb{st_el(i)};
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
                    add_to(j) = nb_diff(i);
                    j         = j+1;
                end
            end
            if j == 1
                add_to = [];
            end
        end
        %  > 2.2.2. -------------------------------------------------------
        function [msh] = Extension_2(msh)
        end
    end
end