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
        function [msh] = Extend_Stencil(msh,bnd_cc,bnd_ff,bnd_fc,nt,p,et_1,et_2)
            %  >  L_(x,y): Face length on the x,y-directions.
            L = A_3_2_2.Compute_Lx_Ly(msh);
            %  > hg_(x,y): Reference length of each domain cell.
            if et_1
                hg = A_3_2_2.Compute_hgx_hgy(msh);
            end
            
            % >> 1.
            %  > 1.1.
            for i = 1:msh.f.NF
                if ~et_1
                    %  > w/o extension #1.
                    [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i)] = ...
                        A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L);
                else
                    %  > w/  extension #1.
                    [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i),par.ngx(i),par.ngy(i)] = ...
                        A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L,hg);
                end
            end
            %  > 1.2.
            %  > Initialize.
            msh.s.par.ne = zeros(2,msh.f.NF);
            msh.s.c_e    = cell (1,msh.f.NF);
            msh.s.f_e    = cell (1,msh.f.NF);
            
            for i = 1:msh.f.NF
                %% > Extension #1 (...until it is incomplete and domain limits have not been reached.).
                %  > Initialize.
                Continue = true;
                
                % >> Check...
                if (~et_1 && round(par.nx(i)) >= p && round(par.ny(i)) >= p) || ...
                   ( et_1 && round(par.nx(i)) >= p && round(par.ny(i)) >= p && round(par.ngx(i)) >= p && round(par.ngy(i)) >= p)
                    if ~et_2
                        continue;
                    end
                else
                    while ((~et_1 && (round(par.nx(i)) < p || round(par.ny(i)) < p)) || ...
                           ( et_1 && (round(par.nx(i)) < p || round(par.ny(i)) < p   || round(par.ngx(i)) < p || round(par.ngy(i)) < p))) && Continue
                        %% > x-direction.
                        Continue_X = false;
                        if (~et_1 &&  round(par.nx(i)) < p) || ...
                           ( et_1 && (round(par.nx(i)) < p  || round(par.ngx(i)) < p))
                            % >> Add/update...
                            [msh,Continue_X] = A_3_2_2.Perform_Extension(i,2,msh,bnd_cc,bnd_ff,bnd_fc,par.ly(1,i),par.ly(2,i),nt);
                            
                            % >> Re-compute parameters.
                            if ~et_1
                                %  > w/o extension #1.
                                [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L);
                            else
                                %  > w/  extension #1.
                                [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i),par.ngx(i),par.ngy(i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L,hg);
                            end
                            % >> Increment number of extensions (x-direction).
                            msh.s.par.ne(1,i) = msh.s.par.ne(1,i)+1;
                        end
                        %% > y-direction.
                        Continue_Y = false;
                        if (~et_1 &&  round(par.ny(i)) < p) || ...
                           ( et_1 && (round(par.ny(i)) < p  || round(par.ngy(i)) < p))
                            % >> Add/update...
                            [msh,Continue_Y] = A_3_2_2.Perform_Extension(i,1,msh,bnd_cc,bnd_ff,bnd_fc,par.lx(1,i),par.lx(2,i),nt);
                            
                            % >> Re-compute parameters.
                            if ~et_1
                                %  > w/o extension #1.
                                [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L);
                            else
                                %  > w/  extension #1.
                                [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i),par.ngx(i),par.ngy(i)] = ...
                                    A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L,hg);
                            end
                            % >> Increment number of extensions (y-direction).
                            msh.s.par.ne(2,i) = msh.s.par.ne(2,i)+1;
                        end
                        Continue = any([Continue_X,Continue_Y]);
                    end
                end
                %% > Extension #2 (...until all cell centroids lie within the stencil limits).
                if et_2
                    % >> Add cell(s) and respective face(s)...
                    msh = A_3_2_2.Extension_2(i,msh,par.lx(:,i),par.ly(:,i),[msh.s.c{:,i}],bnd_cc,bnd_ff,bnd_fc);
                    
                    % >> Re-compute parameters (w/o incrementing the number of extensions).
                    if ~et_1
                        %  > w/o extension #1.
                        [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i)] = ...
                            A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L);
                    else
                        %  > w/  extension #1.
                        [par.lx(:,i),par.ly(:,i),par.nx(i),par.ny(i),par.ngx(i),par.ngy(i)] = ...
                            A_3_2_2.Compute_Parameters(i,msh,msh.s.c(:,i),L,hg);
                    end
                end
            end
            %  > Deal fields...
            msh.s.par.lx = par.lx;
            msh.s.par.ly = par.ly;
            msh.s.par.nx = par.nx;
            msh.s.par.ny = par.ny;
            if et_1
                msh.s.par.ngx = par.ngx;
                msh.s.par.ngy = par.ngy;
            end
            %  > Total stencil points.
            for i = 1:msh.f.NF
                msh.s.xy_v_t{i} = A_3_2_1.Compute_Coordinates_tt(msh.s.xy_v_c{i},msh.s.xy_v_f{i});
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [lx,ly,nx,ny,ngx,ngy] = Compute_Parameters(i,msh,stl_c,L,hg)
            % >> Deal elements.
            arr_c = A_3_2_1.Deal_StencilElem(stl_c);
            
            % >> (x,y)_min,max.
            [lcx(1),lcx(2)] = MinMaxElem(msh.s.xy_v_c{i}(1,:),'finite');
            [lcy(1),lcy(2)] = MinMaxElem(msh.s.xy_v_c{i}(2,:),'finite');
            %  > .. for face w/o boundary faces.
            if ~isempty(msh.s.xy_v_f{i})
                [lfx(1),lfx(2)] = MinMaxElem(msh.s.xy_v_f{i}(1,:),'finite');
                [lfy(1),lfy(2)] = MinMaxElem(msh.s.xy_v_f{i}(2,:),'finite');
                [ lx(1),~     ] = MinMaxElem([lcx(1),lfx(1)],'finite');
                [    ~,  lx(2)] = MinMaxElem([lcx(2),lfx(2)],'finite');
                [ ly(1),~     ] = MinMaxElem([lcy(1),lfy(1)],'finite');
                [    ~,  ly(2)] = MinMaxElem([lcy(2),lfy(2)],'finite');
            else
                lx = lcx;
                ly = lcy;
            end
            
            % >> n_(x,y).
            %  > n_j.
            len = length(arr_c);
            j   = 1:len;
            Lxj = L(1,arr_c(j));
            Lyj = L(2,arr_c(j));
            %  > n.
            Lx  = sum(Lxj)./len;
            Ly  = sum(Lyj)./len;
            nx  = (lx(2)-lx(1))./Lx;
            ny  = (ly(2)-ly(1))./Ly;
            
            % >> ng_(x,y).
            if nargin == 5
                ngx = (lx(2)-lx(1))./hg(1);
                ngy = (ly(2)-ly(1))./hg(2);
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
        function [msh,Flag] = Perform_Extension(i,k,msh,bnd_cc,bnd_ff,bnd_fc,v_min,v_max,nt)
            % >> Previous stencils'...
            %  > ...length.
            len_c  = nnz(~cellfun(@isempty,msh.s.c(:,i)));
            %  > ...cell(s).
            stl_c  = [msh.s.c{:,i}];
            %  > ...cell neighbour(s).
            stl_cn = [msh.c.c{[msh.s.c{:,i}]}];
            stl_cn = A_Tools.fft_unique(sort(stl_cn));
            Add_c  = A_3_2_2.Extension_1(stl_c,stl_cn,msh.c.mean,k,v_min,v_max,msh.c.f.f,nt);
            
            % >> Update/add...
            [msh,Flag] = A_3_2_2.Add_Update(i,msh,Add_c,len_c,bnd_cc,bnd_ff,bnd_fc);
        end
        %  > 1.2.1. -------------------------------------------------------
        function [msh,Flag] = Add_Update(i,msh,Add_c,len_c,bnd_cc,bnd_ff,bnd_fc)
            %  > Initialize.
            Flag_c = false;
            
            % >> Update/add...
            %  > ...cell indices.
            msh.s.c{len_c+1,i} = Add_c;
            msh.s.c_e      {i} = [msh.s.c_e{i},Add_c];
            %  > ...face indices.
            i_face = ismembc(Add_c,bnd_cc) == true;
            if any(i_face)
                Add_f              = A_3_2_1.Add_Face(i_face,Add_c,bnd_ff,bnd_fc);
                msh.s.f{len_c+1,i} = Add_f;
                msh.s.f_e      {i} = [msh.s.f_e{i},Add_f];
            end
            %  > Check whether any additions were made.
            Flag = ~isempty(Add_c);
            
            % >> Compute coordinates of...
            if Flag
                %  > ...cell(s).
                msh.s.xy_v_c{i} = A_3_2_1.Compute_Coordinates_cf(msh.s.c(:,i),msh.c.mean);
                %  > ...face(s) if there is any to be added.
                if any(~cellfun(@isempty,msh.s.f(:,i)))
                    msh.s.xy_v_f{i} = A_3_2_1.Compute_Coordinates_cf(msh.s.f(:,i),msh.f.mean);
                end
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [add_to] = Extension_1(stl_c,stl_cn,mean_c,k,v_min,v_max,cell_f,nt)
            % >> Select...
            %  > ...outer cells (NOT in the current stencil) within the stencil limits.
            out_c = A_Tools.fft_setdiff(stl_cn,stl_c);
            j     = 1;
            for i = 1:length(out_c)
                if mean_c(k,out_c(i)) >= v_min && mean_c(k,out_c(i)) <= v_max
                    add_to_v(j) = out_c(i);
                    j           = j+1;
                end
            end
            
            %  > ...stencil type.
            l = 1;
            if ~nt
                % >> Vertex stencil.
                add_to = add_to_v;
            else
                % >> Face stencil.
                if j ~= 1
                    %  > Outer cell layer's faces.
                    for i = 1:length(stl_c)
                        out_f{i} = [cell_f{stl_c(i)}];
                    end
                    out_f = A_Tools.fft_unique(sort([out_f{:}]));
                    
                    for i = 1:length(add_to_v)
                        if any(ismembc(cell_f{add_to_v(i)},out_f))
                            add_to(l) = add_to_v(i);
                            l         = l+1;
                        end
                    end
                end
            end
            
            %  > If no added elements...
            if j == 1 || (l == 1 && nt)
                add_to = [];
            end
        end
        %  > 1.2.3. -------------------------------------------------------
        function [msh] = Extension_2(i,msh,x_lim,y_lim,stl_c,bnd_cc,bnd_ff,bnd_fc)
            %  > Select cells within stencil x,y-limits.
            k = 1;
            for j = 1:msh.c.NC
                if (msh.c.mean(1,j) >= x_lim(1) && msh.c.mean(1,j) <= x_lim(2)) && (msh.c.mean(2,j) >= y_lim(1) && msh.c.mean(2,j) <= y_lim(2))
                    bulk(k) = j;
                    k       = k+1;
                end
            end
            %  > Select cells in 'bulk' that don't belong to the stencil.
            Add_c = A_Tools.fft_setdiff(bulk,stl_c);
            
            %  > Update/add...
            if ~isempty(Add_c)
                len_c = nnz(~cellfun(@isempty,msh.s.c(:,i)));
                msh   = A_3_2_2.Add_Update(i,msh,Add_c,len_c,bnd_cc,bnd_ff,bnd_fc);
            end
        end
    end
end