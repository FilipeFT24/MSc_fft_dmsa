classdef B_2_1D
    methods (Static)
        %% > Wrap-up B_2 (1D).
        function [pde] = WrapUp_B_2_1D(msh,pde,ft,st,ng,np,v,g,bnd_w,bnd_e)
            % >> Set...
            %  > ...A and B.
            [msh,f,gradf,A,B] = B_2_1D.Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e);
            %  > ...PDE solution.
            switch ft
                case 'Implicit'
                    %  > Flux reconstruction: Implicit.
                    pde = B_2_1D.PDE_Implicit(msh,pde,A,B);
                    pde = B_2_1D.PDE_fv      (msh,pde,f,gradf);
                case 'Explicit'
                    %  > Flux reconstruction: Explicit.
                    pde = B_2_1D.PDE_Explicit(msh,pde,A,B);
                otherwise
                    return;
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [msh,f,gradf,A,B] = Assemble_AB(msh,pde,v,g,st,np,ng,bnd_w,bnd_e)
            %% > Face 'i'...
            for i = 1:msh.f.NF
                %  > ...stencil indices.
                [msh.s.c{i},msh.s.f{i}] = ...
                    A_2_1D.SetUp_Stencil_f(i,msh,np);
                %  > ...stencil coordinates.
                [msh.s.x_v_c{i},msh.s.x_v_f{i},msh.s.x_v_t{i}] = ...
                    A_2_1D.Compute_Coordinates_cft(msh,msh.s.c{i},msh.s.f{i});
                %  > ...Phi_f & GradPhi_f.
                [f{i},gradf{i},xf.v{i},xf.g{i}] = ...
                    A_2_1D.x_f(v,g,np,msh.s.x_v_t{i},msh.f.Xv(i));
            end
            
            %% > Assemble matrix A.
            % >> A (Cell dependent coefficients).
            %  > ... for cell #i: Phi(C) = Phi_f(e)-Phi_f(w).
            %  > Initialize.
            A = zeros(msh.c.NC);
            
            for i = 1:msh.c.NC
                %  > West face (w).
                nf_cw        {i}  = length(msh.s.c{i});
                kf_cw        {i}  = msh.s.c{i};
                xf_vw        {i}  = xf.v{i};
                xf_gw        {i}  = xf.g{i};
                A    (i,kf_cw{i}) = A(i,kf_cw{i})-( xf_vw{i}(1:nf_cw{i}));
                A    (i,kf_cw{i}) = A(i,kf_cw{i})-(-xf_gw{i}(1:nf_cw{i}));
                %  > East face (e).
                nf_ce        {i}  = length(msh.s.c{i+1});
                kf_ce        {i}  = msh.s.c{i+1};
                xf_ve        {i}  = xf.v{i+1};
                xf_ge        {i}  = xf.g{i+1};
                A    (i,kf_ce{i}) = A(i,kf_ce{i})+( xf_ve{i}(1:nf_ce{i}));
                A    (i,kf_ce{i}) = A(i,kf_ce{i})+(-xf_ge{i}(1:nf_ce{i}));
            end
            
            %% > Assemble matrix B.
            % >> B (Face dependent coefficients).
            %  > Initialize.
            B     = zeros(msh.c.NC,1);
            np_2  = np./2;
            %  > Auxiliary arrays.
            f_v   = pde.an.f_v;
            df_v  = pde.an.df_v;
            
            % >> Boundary contribution(s).
            %  > West boundary(w).
            kc_iw = 1;
            kc_fw = np_2;
            for i = kc_iw:kc_fw
                %  > A.
                Aw      = B_2_1D.BC_A('w',bnd_w,msh.c.NC,msh.s.f{i},nf_cw{i},msh.s.f{i+1},nf_ce{i},xf_vw{i},xf_gw{i},xf_ve{i},xf_ge{i},xf.g{1},xf.g{end},v,g);
                A (i,:) = A(i,:)+Aw;
                %  > B.
                Bw      = B_2_1D.BC_B('w',bnd_w,msh.s.f{i},nf_cw{i},msh.s.f{i+1},nf_ce{i},xf_vw{i},xf_gw{i},xf_ve{i},xf_ge{i},f_v,df_v,xf.v{1},xf.g{1},xf.v{end},xf.g{end},v,g);
                B (i,1) = B(i,1)-Bw;
            end
            %  > East boundary(e).
            kc_ie = msh.c.NC-np_2+1;
            kc_fe = msh.c.NC;
            for i = kc_ie:kc_fe
                %  > A.
                Ae      = B_2_1D.BC_A('e',bnd_e,msh.c.NC,msh.s.f{i},nf_cw{i},msh.s.f{i+1},nf_ce{i},xf_vw{i},xf_gw{i},xf_ve{i},xf_ge{i},xf.g{1},xf.g{end},v,g);
                A (i,:) = A(i,:)+Ae;
                %  > B.
                Be      = B_2_1D.BC_B('e',bnd_e,msh.s.f{i},nf_cw{i},msh.s.f{i+1},nf_ce{i},xf_vw{i},xf_gw{i},xf_ve{i},xf_ge{i},f_v,df_v,xf.v{1},xf.g{1},xf.v{end},xf.g{end},v,g);
                B (i,1) = B(i,1)-Be;
            end
            
            % >> Source term contribution(s).
            i     = 1:msh.c.NC;
            F_Vol = B_1_2_1D.Compute_SourceTerm(msh,pde,st,ng);
            B     = B+F_Vol(i);
            
            %  > Transform to sparse notation...
            A = sparse(A);
            B = sparse(B);
        end
        %  > 1.2. ---------------------------------------------------------
        %  > 1.2.1. -------------------------------------------------------
        function [Add] = BC_A(bnd,bnd_type,NC,f_w,nc_w,f_e,nc_e,xf_vw,xf_gw,xf_ve,xf_ge,xf_gi,xf_gf,v,g)
            %  > Initialize.
            Add = zeros(1,NC);
            
            switch bnd
                %  > Select boundary...
                case 'w'
                    %  > West boundary(w).
                    n_fw = length(f_w)+nc_w;
                    if ~isempty(f_e)
                        n_fe = length(f_e)+nc_e;
                    end
                    li      =  length(xf_gi)-1;
                    F(1:li) = -xf_gi(1:li)./xf_gi(li+1);

                    switch bnd_type
                        %  > Select boundary type...
                        case 'Dirichlet'
                            %  > Do nothing.
                        case 'Neumann'
                            Add(1:nc_w) = -F(1:nc_w).*(xf_vw(n_fw)-xf_gw(n_fw));
                            if ~isempty(f_e)
                                Add(1:nc_e) = Add(1:nc_w)+F(1:nc_e).*(xf_ve(n_fe)-xf_ge(n_fe));
                            end
                        case 'Robin'
                        otherwise
                            return;
                    end
                case 'e'
                    %  > East boundary(e).
                    n_fe = length(f_e)+nc_e;
                    if ~isempty(f_w)
                        n_fw = length(f_w)+nc_w;
                    end
                    len_f      =  length(xf_gf)-1;
                    F(1:len_f) = -xf_gf(1:len_f)./xf_gf(len_f+1);
                    
                    switch bnd_type
                        %  > Select boundary type...
                        case 'Dirichlet'
                            %  > Do nothing.
                        case 'Neumann'
                            Add(1:nc_e) = F(1:nc_e).*(xf_ve(n_fe)-xf_ge(n_fe));
                            if ~isempty(f_w)
                                Add(1:nc_w) = Add(1:nc_e)-F(1:nc_w).*(xf_vw(n_fw)-xf_gw(n_fw));
                            end
                            Add = circshift(Add,nnz(~Add));
                        case 'Robin'
                        otherwise
                            return;
                    end
                otherwise
                    return;
            end 
        end
        %  > 1.2.2. -------------------------------------------------------
        function [Add_t] = BC_B(bnd,bnd_type,f_w,nc_w,f_e,nc_e,xf_vw,xf_gw,xf_ve,xf_ge,f_v,df_v,xf_vi,xf_gi,xf_vf,xf_gf,v,g)
            switch bnd
                %  > Select boundary...
                case 'w'
                    %  > West boundary(w).
                    n_fw = length(f_w)+nc_w;
                    if ~isempty(f_e)
                        n_fe = length(f_e)+nc_e;
                    end
                    switch bnd_type
                        %  > Select boundary type...
                        case 'Dirichlet'
                            f_0 = f_v(1);
                        case 'Neumann'
                            f_0 = df_v(1)./(xf_vi(n_fw)-xf_gi(n_fw));
                        case 'Robin'
                            f_0 = (v.*f_v(1)-g.*df_v(1))./(v-g.*xf_gi(n_fw));
                        otherwise
                            return;
                    end
                    %  > Add contribution(s)...
                    Add_v = -f_0.*xf_vw(n_fw);
                    Add_g =  f_0.*xf_gw(n_fw);
                    if ~isempty(f_e)
                        n_fe  = length(f_e)+nc_e;
                        Add_v = Add_v+f_0.*xf_ve(n_fe);
                        Add_g = Add_g-f_0.*xf_ge(n_fe);
                    end
                case 'e'
                    %  > East boundary(e).
                    n_fe = length(f_e)+nc_e;
                    if ~isempty(f_w)
                        n_fw = length(f_w)+nc_w;
                    end
                    switch bnd_type
                        %  > Select boundary type...
                        case 'Dirichlet'
                            f_0 = f_v(end);
                        case 'Neumann'
                            f_0 = df_v(end)./xf_gf(n_fe);
                        case 'Robin'
                            f_0 = (v.*f_v(end)-g.*df_v(end))./(v-g.*xf_gf(n_fe));
                        otherwise
                            return;
                    end
                    %  > Add contribution(s)...
                    Add_v =  f_0.*xf_ve(n_fe);
                    Add_g = -f_0.*xf_ge(n_fe);
                    if ~isempty(f_w)
                        n_fw  = length(f_w)+nc_w;
                        Add_v = Add_v-f_0.*xf_vw(n_fw);
                        Add_g = Add_g+f_0.*xf_gw(n_fw);
                    end
                otherwise
                    return;
            end
            Add_t = Add_v+Add_g;
        end
        % >> 1.3. ---------------------------------------------------------
        function [X] = SetUp_Solver(A,B,str)
            if strcmpi(str,'backslash')
                X       = A\B; % Equivalent to: V*((U'*B)./s).
            elseif strcmpi(str,'Tikhonov')
                lmbd    = 1E-06;
                [U,S,V] = svd_lapack(A);
                s       = diag(S);
                X       = V*(s.*(U'*B)./(s.^2+lmbd));
            else
                return;
            end
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Implicit flux reconstruction.
        %  > 1.4.1. -------------------------------------------------------
        function [pde] = PDE_Implicit(msh,pde,A,B)
            %  > Phi(C).
            pde.Phi = B_2_1D.SetUp_Solver(A,B,'backslash');
            %  > Error norms.
            i       = 1:msh.c.NC;
            Xc(i)   = abs(pde.an.c_v(i)-pde.Phi(i));
            pde.Ec  = B_2_1D.Compute_ErrorNorms(msh,Xc);
        end
        %  > 1.4.2. -------------------------------------------------------
        function [pde] = PDE_fv(msh,pde,f,gradf)
            for i = 1:msh.f.NF
                % >> Implicit flux reconstruction values...
                %  > ...cell.
                n_c(i) = length(msh.s.c{i});
                k_c{i} = msh.s.c{i};
                v_c{i} = pde.Phi(k_c{i});
                v_t{i} = v_c{i};
                %  > ...face.
                if ~isempty(msh.s.f{i})
                    n_f(i) = length(msh.s.f{i})+n_c(i);
                    k_f{i} = msh.s.f{i};
                    v_f{i} = pde.an.f_v(k_f{i});
                    v_t{i} = [v_t{i};[v_f{i}]];
                end
                % >> PDE face values...
                %  > ...Phi_f.
                Phi_f(i,1)     = f{i}*v_t{i};
                %  > ...gradPhi_f.
                GradPhi_f(i,1) = gradf{i}*v_t{i};
            end
            i              = 1:msh.f.NF;
            pde.Ef.EA(i,1) = abs(Phi_f    (i,1)-pde.an.f_v (i,1));
            pde.Ef.EA(i,2) = abs(GradPhi_f(i,1)-pde.an.df_v(i,1));
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Explicit flux reconstruction.
        function [pde] = PDE_Explicit(msh,pde,A,B)
            %  > Phi(C).
            Phic   = pde.an.blk;
            %  > Error norms.
            Xc      = abs(A*Phic-B)';
            pde.Ec  = B_2_1D.Compute_ErrorNorms(msh,Xc);
        end
        
        %% > 2. -----------------------------------------------------------
        function [Ec]   = Compute_ErrorNorms(msh,X)
            i           = 1:msh.c.NC;
            Ec.EA(i,1)  = X(i);
            Ec_iX{1}(i) = X(i).*msh.c.Vol(i);
            Ec_iX{2}(i) = X(i).^2.*msh.c.Vol(i).^2;
            Ec.EN(1)    = sum(Ec_iX{1})./sum(msh.c.Vol);
            Ec.EN(2)    = sum(sqrt(Ec_iX{2}))./sum(sqrt(msh.c.Vol.^2));
            Ec.EN(3)    = max(X);
        end
    end
end