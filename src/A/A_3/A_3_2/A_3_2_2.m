classdef A_3_2_2
    methods (Static)
        %% > A_3_2_2.
        % >> --------------------------------------------------------------
        % >> 1.     Compute stencil limits and check if and extension is required.
        % >> 2.     Tools.
        %  > 2.1.   Compute parameters.
        %  > 2.2.   Perform extension.
        %  > 2.2.1. Extension #1.
        %  > 2.2.2. Extension #2.
        % >> --------------------------------------------------------------
                
        %% > 1. -----------------------------------------------------------
        function [msh] = Extend_Stencil(msh,bnd_cc,p)
            % >> 2.1.
            for i = 1:msh.f.NF
                [par.n_x(i),par.n_y(i),par.l_x(:,i),par.l_y(:,i)] = ...
                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),msh.s.f(:,i));
            end
            % >> 2.2.
            %  > Initialize.
            msh.s.par.n_e = zeros(2,msh.f.NF);
            msh.s.c_e     = cell (1,msh.f.NF);
            msh.s.c_f     = cell (1,msh.f.NF);
            %  > Extend stencil until...
            for i = 1:msh.f.NF
                if par.n_x(i) >= p-1./2 && par.n_y(i) >= p-1./2
                    continue;
                else
                    %  > Do while stencil is incomplete and domain limits have not been reached...
                    Continue_X = 1;
                    %while par.n_x(i) < p-1./2 || Continue_X 
                        [msh,Continue_X] = ...
                            A_3_2_2.Perform_Extension(i,msh,bnd_cc,par,p);
                        [par.n_x(i),par.l_x(:,i)] = ...
                            A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),msh.s.f(:,i));
                    %end
%                     Continue_Y = 1;
%                     while par.n_y(i) < p-1./2 || Continue_Y 
%                     end
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
        function [n_x,n_y,l_x,l_y] = Compute_Parameters(i,msh,st_c,st_f)
            % >> Deal elements.
            [arr_c,~] = A_3_2_1.Deal_StencilElem(st_c,st_f);
            
            % >> Compute/re-compute 'par' fields.
            %  > (x,y)_min.
            l_x(1) = min(msh.s.xy_v_t{i}(1,:));
            l_y(1) = min(msh.s.xy_v_t{i}(2,:));
            %  > (x,y)_max.
            l_x(2) = max(msh.s.xy_v_t{i}(1,:));
            l_y(2) = max(msh.s.xy_v_t{i}(2,:));
            %  > (n_x,n_y).
            for j = 1:length(arr_c)
                h_x(j) = max(msh.c.xy_v{arr_c(j)}(:,1))-min(msh.c.xy_v{arr_c(j)}(:,1));
                h_y(j) = max(msh.c.xy_v{arr_c(j)}(:,2))-min(msh.c.xy_v{arr_c(j)}(:,2));
            end
            h_x = sum(h_x)./length(arr_c);
            h_y = sum(h_y)./length(arr_c);
            n_x = (l_x(2)-l_x(1))./h_x;
            n_y = (l_y(2)-l_y(1))./h_y;
        end
        % >> 2.2. ---------------------------------------------------------
        function [msh,Continue] = Perform_Extension(i,msh,bnd_cc,par,p)
            %  > Initialize.
            Ext_Flag = 0;
            Continue = 0;
            k        = 1;
            len      = nnz(~cellfun(@isempty,msh.s.c(:,i)));
            
            %  > Previous layers' stencil elements.
            Lst = msh.s.c{len,i};
            for j = 1:len
                st{j} = msh.s.c{j,i};
            end
            st = cell2mat(st);
                        
            
            
            
            if par.n_x(i) < p-1./2
                % >> Add cell indices.
                %  > Stencil limits.
                y_min = par.l_y(1,i);
                y_max = par.l_y(2,i);
                Ext_X = A_3_2_2.Extension_1(msh,Lst,st,'x',y_min,y_max);
                
                % >> Update/add...
                %  > Number of extensions (y-direction).
                msh.s.par.n_e(1,i) = msh.s.par.n_e(1,i)+1;
                %  > Cell indices.
                msh.s.c{len+1,i} = Ext_X;
                msh.s.c_e    {i} = [msh.s.c_e{i},Ext_X];
                %  > Face indices.
                Flag = ismembc(Ext_X,bnd_cc) == 1;
                if any(Flag)
                    msh.s.c_f{i}       = A_3_2_1.Add_Face(Flag,msh,Ext_X);
                    msh.s.f  {len+1,i} = msh.s.c_f{i};
                end
                Ext_Flag(1) = length(Ext_X);
            end
            if par.n_y(i) < p-1./2
                % >> Add cell indices.
                % > Stencil limits.
                x_min = par.l_x(1,i);
                x_max = par.l_x(2,i);
                Ext_Y = A_3_2_2.Extension_1(msh,Lst,st,'y',x_min,x_max);
                
                % >> Update/add...
                %  > Number of extensions (y-direction).
                msh.s.par.n_e(2,i) = msh.s.par.n_e(2,i)+1;
                %  > Cell indices.
                msh.s.c{len+1,i} = Ext_Y;
                msh.s.c_e    {i} = [msh.s.c_e{i},Ext_Y];
                %  > Face indices.
                Flag = ismembc(Ext_Y,bnd_cc) == 1;
                if any(Flag)
                    msh.s.c_f{i}       = A_3_2_1.Add_Face(Flag,msh,Ext_Y);
                    msh.s.f  {len+1,i} = msh.s.c_f{i};
                end
                Ext_Flag(2) = length(Ext_Y);
            end
            Continue = any(Ext_Flag);
            if Continue
                msh.s.xy_v_t{i} = A_3_2_1.Compute_Coordinates(msh,msh.s.c(:,i),msh.s.f(:,i));
            end
        end
        %  > 2.2.1. -------------------------------------------------------
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
        %  > 2.2.2. -------------------------------------------------------
        function [msh] = Extension_2(msh)
        end
    end
end