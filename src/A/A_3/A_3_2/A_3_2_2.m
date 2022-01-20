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
        function [msh] = Extend_Stencil(msh,bnd_cc,nt,p,et)
            % >> Auxiliary arrays.
            %  > bnd_ff,bnd_fc.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            %  >  L_(x,y).
            L = A_3_2_2.Compute_Lx_Ly(msh);
            %  > hg_(x,y).
            if ~et
                %  > Initialize.
                hg = zeros(2,1);
            else
                hg = A_3_2_2.Compute_hgx_hgy(msh);
            end
            
            % >> 1.
            %  > 1.1.
            for i = 1:msh.f.NF
                if ~et
                    %  > w/o extension #1.
                    [par.nx(i),~,par.lx(:,i)] = A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(1,:),hg(1),1,et);
                    [par.ny(i),~,par.ly(:,i)] = A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(2,:),hg(2),2,et);
                else
                    %  > w/  extension #1.
                    [par.nx(i),par.ngx(i),par.lx(:,i)] = A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(1,:),hg(1),1,et);
                    [par.ny(i),par.ngy(i),par.ly(:,i)] = A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(2,:),hg(2),2,et);
                end
            end
            %  > 1.2.
            %  > Initialize.
            msh.s.par.ne = zeros(2,msh.f.NF);
            msh.s.c_e    = cell (1,msh.f.NF);
            msh.s.f_e    = cell (1,msh.f.NF);
            %  > Extend stencil until it is incomplete and domain limits have not been reached...
            for i = 1:msh.f.NF
                %  > Initialize.
                Continue = true;
                %  > Check...
                if (~et && round(par.nx(i)) >= p && round(par.ny(i)) >= p) || ...
                    (et && round(par.nx(i)) >= p && round(par.ny(i)) >= p && round(par.ngx(i)) >= p && round(par.ngy(i)) >= p)
                    continue;
                else
                    while ((~et && (round(par.nx(i)) < p || round(par.ny(i)) < p)) || ...
                            (et && (round(par.nx(i)) < p || round(par.ny(i)) < p   || round(par.ngx(i)) < p || round(par.ngy(i)) < p))) && Continue
                        %% > x-direction.
                        Continue_X = false;
                        if (~et &&  round(par.nx(i)) < p) || ...
                            (et && (round(par.nx(i)) < p  || round(par.ngx(i)) < p))
                            %  > Add/update...
                            [msh,Continue_X] = A_3_2_2.Perform_Extension(i,2,msh,bnd_cc,par.ly(1,i),par.ly(2,i),nt,bnd_ff,bnd_fc);
                            %  > Re-compute parameters.
                            if ~et
                                [par.nx(i),~,par.lx(:,i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(1,:),hg(1),1,et);
                            else
                                [par.nx(i),par.ngx(i),par.lx(:,i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(1,:),hg(1),1,et);
                            end
                            %  > Number of extensions (x-direction).
                            msh.s.par.ne(1,i) = msh.s.par.ne(1,i)+1;
                        end
                        %% > y-direction.
                        Continue_Y = false;
                        if (~et &&  round(par.ny(i)) < p) || ...
                            (et && (round(par.ny(i)) < p  || round(par.ngy(i)) < p))
                            %  > Add/update...
                            [msh,Continue_Y] = A_3_2_2.Perform_Extension(i,1,msh,bnd_cc,par.lx(1,i),par.lx(2,i),nt,bnd_ff,bnd_fc);
                            %  > Re-compute parameters.
                            if ~et
                                [par.ny(i),~,par.ly(:,i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(2,:),hg(2),2,et);
                            else
                                [par.ny(i),par.ngy(i),par.ly(:,i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L(2,:),hg(2),2,et);
                            end
                            %  > Number of extensions (y-direction).
                            msh.s.par.ne(2,i) = msh.s.par.ne(2,i)+1;
                        end
                        Continue = any([Continue_X,Continue_Y]);
                    end
                end
            end
            %  > Deal fields...
            msh.s.par.nx = par.nx;
            msh.s.par.ny = par.ny;
            if et
                msh.s.par.ngx = par.ngx;
                msh.s.par.ngy = par.ngy;
            end
            msh.s.par.lx = par.lx;
            msh.s.par.ly = par.ly;
            %  > Total stencil points.
            for i = 1:msh.f.NF
                msh.s.xy_v_t{i} = A_3_2_1.Compute_Coordinates_tt(msh.s.xy_v_c{i},msh.s.xy_v_f{i});
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [n,ng,l] = Compute_Parameters(i,msh,st_c,L,hg,k,et)
            % >> Deal elements.
            arr_c = A_3_2_1.Deal_StencilElem(st_c);
            
            % >> Compute/re-compute 'par' fields.
            %  > (x,y)_min,max.
            [lc(1),lc(2)] = MinMaxElem(msh.s.xy_v_c{i}(k,:),'finite');
            %  > .. for face w/o boundary faces.
            if ~isempty(msh.s.xy_v_f{i})
                [lf(1),lf(2)] = MinMaxElem(msh.s.xy_v_f{i}(k,:),'finite');
                [ l(1), ~   ] = MinMaxElem([lc(1),lf(1)],'finite');
                [    ~, l(2)] = MinMaxElem([lc(2),lf(2)],'finite');
            else
                l(1) = lc(1);
                l(2) = lc(2);
            end
            %  > n_(x,y).
            len_c = length(arr_c);
            j     = 1:len_c;
            L     = sum(L(arr_c(j)))./len_c;
            n     = (l(2)-l(1))./L;
            %  > ng_(x,y).
            if ~et
                %  > Initialize.
                ng = 0;
            else
                ng = (l(2)-l(1))./hg;
            end
        end
        % >> 1.1.1. -------------------------------------------------------
        function [Lt] = Compute_Lx_Ly(msh)
            for j = 1:msh.c.NC
                %  > (Lx,Ly)_min,max.
                [Lx(1,j),Lx(2,j)] = MinMaxElem(msh.c.xy_v{j}(:,1),'finite');
                [Ly(1,j),Ly(2,j)] = MinMaxElem(msh.c.xy_v{j}(:,2),'finite');
            end
            %  > (Lt)_x,y.
            [Lt(1,:),Lt(2,:)] = deal(Lx(2,:)-Lx(1,:),Ly(2,:)-Ly(1,:));
        end
        % >> 1.1.2. -------------------------------------------------------
        function [hgt] = Compute_hgx_hgy(msh)
            for j = 1:msh.c.NC
                %  > (Xv,Yv)_min,max.
                [Xv(1,j),Xv(2,j)] = MinMaxElem(msh.c.xy_v{j}(:,1),'finite');
                [Yv(1,j),Yv(2,j)] = MinMaxElem(msh.c.xy_v{j}(:,2),'finite');
                %  > (Xv)_(max-min)-(Yv)_(max-min).
                [hx(j)  ,hy(j)  ] = deal(Xv(2,j)-Xv(1,j),Yv(2,j)-Yv(1,j));
            end
            %  > (hgt)_x,y.
            [hgt(1,:),hgt(2,:)] = deal(sum(hx)./msh.c.NC,sum(hy)./msh.c.NC);
        end
        % >> 1.2. ---------------------------------------------------------
        function [msh,Flag] = Perform_Extension(i,k,msh,bnd_cc,v_min,v_max,Type,bnd_ff,bnd_fc)
            %  > Initialize.
            Flag = false;
            len  = nnz(~cellfun(@isempty,msh.s.c(:,i)));
            
            %  > Previous layers' stencil elements.
            for j = 1:len
                st_el{j} = msh.s.c{j,i};
            end
            st_el = cat(2,msh.s.c{:,i});
            Add   = A_3_2_2.Extension_1(msh,st_el,k,v_min,v_max,Type);
            
            % >> Update/add...
            %  > Cell indices.
            msh.s.c{len+1,i} = Add;
            msh.s.c_e    {i} = [msh.s.c_e{i},Add];
            %  > Face indices.
            i_face = ismembc(Add,bnd_cc) == true;
            if any(i_face)
                msh.s.f_e{i}       = A_3_2_1.Add_Face(i_face,Add,bnd_ff,bnd_fc);
                msh.s.f  {len+1,i} = msh.s.f_e{i};
            end
            Flag = ~isempty(Add);
            
            if Flag
                %  > Cells.
                msh.s.xy_v_c{i} = A_3_2_1.Compute_Coordinates_cf(msh.s.c(:,i),msh.c.mean);
                %  > Faces.
                if any(~cellfun(@isempty,msh.s.f(:,i)))
                    msh.s.xy_v_f{i} = A_3_2_1.Compute_Coordinates_cf(msh.s.f(:,i),msh.f.mean);
                end
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        function [add_to] = Extension_1(msh,st_el,k,v_min,v_max,nt)
            % >> Stencil cell neighbours.
            for i = 1:length(st_el)
                nb_c{i} = msh.c.c{st_el(i)};
            end
            nb_c    = cat(2,nb_c{:});
            nb_c_un = unique(nb_c);
            %  > Outer cells (NOT in the stencil).
            nb_diff = A_Tools.fft_setdiff(nb_c_un,st_el);
            
            %  > Select cells within stencil limits.
            j = 1;
            for i = 1:length(nb_diff)
                if msh.c.mean(k,nb_diff(i)) >= v_min && msh.c.mean(k,nb_diff(i)) <= v_max
                    add_to_v(j) = nb_diff(i);
                    j           = j+1;
                end
            end
            
            l = 1;
            if ~nt
                %  > Do nothing...
                add_to = add_to_v;
            else
                if j ~= 1
                    %  > Outer cell layer's faces.
                    for i = 1:length(st_el)
                        out_f{i} = msh.c.f.f{st_el(i)};
                    end
                    out_f = unique(cat(2,out_f{:}));
                    
                    %  > Add cell to extended stencil if it has a common face with the previous stencil layer.
                    for i = 1:length(add_to_v)
                        if any(ismembc(msh.c.f.f{add_to_v(i)},out_f))
                            add_to(l) = add_to_v(i);
                            l         = l+1;
                        end
                    end
                end
            end
            
            %  > No added elements...
            if j == 1 || (l == 1 && nt)
                add_to = [];
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [msh] = Extension_2(msh)
        end
    end
end