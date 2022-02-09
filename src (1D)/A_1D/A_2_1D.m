classdef A_2_1D
    methods (Static)
        %% > Wrap-up A_2 (1D).
        function [inp,msh] = WrapUp_A_2_1D(inp)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            h    = inp.msh.h;
            eg   = inp.msh.eg;
            
            switch eg
                case "1"
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_Uniform(Xv_i,Xv_f,h);
                case "2"
                    Nf_X        = inp.msh.s_nu.Nf_X;
                    Ks_X        = inp.msh.s_nu.Ks_X;
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_1(Xv_i,Xv_f,h,Nf_X,Ks_X);
                case "3"
                    Ks_X        = inp.msh.s_nu.Ks_X;
                    Lc_X        = inp.msh.s_nu.Lc_X;
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_2(Xv_i,Xv_f,h,Ks_X,Lc_X);
                otherwise
                    return;
            end
            % >> Grid properties.
            %  > Number of cells/faces.
            msh.c.NC     = NX_c;
            msh.f.NF     = NX_c+1;
            %  > Cell vertex coordinates.
            msh.f.Xv     = Xd_x;
            %  > Xc/Cell volume.
            i            = 1:msh.c.NC;
            msh.c.Xc (i) = 1./2.*(msh.f.Xv(i)+msh.f.Xv(i+1));
            msh.c.Vol(i) = msh.f.Xv(i+1)-msh.f.Xv(i);
            %  > Reference length.
            msh.d.H_ref  = (msh.f.Xv(msh.f.NF)-msh.f.Xv(1))./msh.c.NC;
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [NX_c,Xd_x] = SquareMesh_Uniform(Xv_i,Xv_f,h)
            % >> Vertex coordinates.
            %  > Xd.
            NX_c = round(1./h.*(Xv_f-Xv_i));
            Xd_x = linspace(Xv_i,Xv_f,NX_c+1);
        end
        % >> 1.2. ---------------------------------------------------------
        function [NX_c,Xd_x] = SquareMesh_NonUniform_1(Xv_i,Xv_f,h,Nf_X,Ks_X)
            % >> Transformation parameters.
            NX_v     = round(1./h.*(Xv_f-Xv_i))+1;
            Nf_Unf_X = linspace(0,1,NX_v);
            Pt_X     = (Nf_X-0)./(1-0);
            B_X      = 1./(2.*Ks_X).*log((1+(exp(Ks_X)-1).*(Pt_X))./(1+(exp(-Ks_X)-1).*(Pt_X)));
            i        = 1:NX_v;
            NF_X     = (Nf_X-0).*(1+sinh(Ks_X.*(Nf_Unf_X(i)-B_X))./sinh(Ks_X.*B_X))+0;
            
            % >> Vertex coordinates.
            %  > Xd.
            Xd_x = NF_X(1:NX_v).*(Xv_f-Xv_i);
            NX_c = NX_v-1;
        end
        % >> 1.3. ---------------------------------------------------------
        function [NX_c,Xd_x] = SquareMesh_NonUniform_2(Xv_i,Xv_f,h,Ks_X,bnd)
            % >> Transformation parameters.
            NX_v   = round(1./h.*(Xv_f-Xv_i))+1;
            Nf_Unf = linspace(0,1,NX_v);
            Bp     = Ks_X+1;
            Bm     = Ks_X-1;
            i      = 1:NX_v;
            X  (i) = (Bp./Bm).^(1-Nf_Unf(i));
            Num(i) = Bp-Bm.*X(i);
            Den(i) = X(i)+1;
            x_h(i) = Num(i)./Den(i);
            
            % >> Vertex coordinates.
            %  > Xd.
            Xd_x = Xv_i+(Xv_f-Xv_i).*x_h;
            NX_c = NX_v-1;
            
            switch bnd
                case "w"
                    %  > Do nothing.
                case "e"
                    %  > Invert mesh.
                    fl_x = flip(Xd_x);
                    Xd_x = ones(1,length(fl_x)).*Xv_i+Xv_f-fl_x;
                otherwise
                    return;
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [s] = Problem_SetUp(msh,stl_p,stl_s,stl_t,sn,bnd,v,g)
            %  > Auxiliary variables.
            Xc     = msh.c.Xc;
            NC     = msh.c.NC;
            Xv     = msh.f.Xv;
            NF     = msh.f.NF;
            f(:,1) = sn.f(:,1);
            f(:,2) = sn.f(:,2);

            % >> Loop through selected faces...
            for i = stl_s
                switch stl_t(i)
                    case "CDS"
                        %  > Auxiliary variables.
                        P1    = stl_p(i)+1;
                        P2    = P1./2;
                        LHS   = P2;
                        RHS   = NF+1-P2;
                        add_f = false;
                                                
                        % >> Add...
                        %  > Check whether face i's stencil contains a boundary face...
                        if i <= LHS
                            %  > ...face(s).
                            add_f = true;
                            bnd_f = 1;
                            stl_f = bnd_f;
                            bnd_i = bnd(1);
                            %  > ...cell(s).
                            stl_c = 1:P1-1;
                        elseif i >= RHS
                            %  > ...face(s).
                            add_f = true;
                            bnd_f = NF;
                            stl_f = bnd_f;
                            bnd_i = bnd(2);
                            %  > ...cell(s).
                            stl_c = NF+1-P1:NC;
                        else
                            %  > ...face(s).
                            stl_f = [];
                            bnd_i = string([]);
                            %  > ...cell(s).
                            stl_c = i-P2:i+P2-1;
                        end
                    case "UPW"
                        %  > Auxiliary variables.
                        P1    = stl_p(i)+1;
                        LHS   = P1;
                        add_f = false;
                        
                        % >> Add...
                        %  > Check whether face i's stencil contains a boundary face...
                        if i <= LHS
                            %  > ...face(s).
                            add_f = true;
                            bnd_f = 1;
                            stl_f = bnd_f;
                            bnd_i = bnd(1);
                            %  > ...cell(s).
                            stl_c = 1:P1-1;
                        else
                            %  > ...cell(s).
                            stl_c = i-P1:i-1;
                            %  > ...face(s).
                            stl_f = [];
                            bnd_i = string([]);
                        end
                    case "DOW"
                    otherwise
                        return;
                end
                %  > Compute stencil coordinates.
                xc = Xc(stl_c);
                if ~add_f
                    xt = xc;
                else
                    xt = [xc,Xv(bnd_f)];
                end
                %  > Stencil length.
                Ls  = (max(xt)-min(xt))./length(stl_c);

                % >> Tf.  
                %  > Auxiliary variables.
                len = length(xt); 
                j   = 1:len; 
                fx  = Xv(i);                
                %  > Df.
                switch bnd_i
                    case "Dirichlet"
                        %  > Boundary value.
                        bnd_v   = f(bnd_f,1);
                        %  > Df.
                        Df(j,j) = (xt(j)-fx)'.^(j-1);
                    case "Neumann"
                        %  > Boundary value.
                        bnd_v   = f(bnd_f,2);
                        %  > Df.
                        k       = 1:len-1;
                        Df(k,j) = (xt(k)-fx)'.^(j-1);
                        l       = len;
                        m       = k+1;
                        Df(l,1) = 0;
                        Df(l,m) = k.*(xt(l)-fx).^(k-1);
                    case "Robin"
                        %  > Boundary value.
                        g_v     = g./v;
                        bnd_v   = f(bnd_f,1)+g_v.*f(bnd_f,2);
                        %  > Df.
                        k       = 1:len-1;
                        Df(k,j) = (xt(k)-fx)'.^(j-1);
                        l       = len;
                        lv  (j) = (xt(l)-fx)'.^(j-1);
                        m       = k+1;
                        lg  (1) = 0;
                        lg  (m) = k.*(xt(l)-fx).^(k-1);
                        Df(l,j) = lv(j)+g./v.*lg(j);
                    otherwise
                        %  > Boundary value.
                        bnd_v   = [];
                        %  > Df.
                        Df(j,j) = (xt(j)-fx)'.^(j-1);
                end
                df       = zeros(2,len);
                df(1,1)  = 1;
                df(2,2)  = df(1,1);
                Inv      = inv(Df);
                Tf       = df*Inv;
                xf_v     = v.*Tf(1,:);
                xf_g     = g.*Tf(2,:);
                xf       = xf_v-xf_g;
                
                %  > Update 'msh' structure...
                s.c  {i} = stl_c;
                s.f  {i} = stl_f;
                s.xt {i} = xt;
                s.Ls (i) = Ls;
                s.bnd{i} = bnd_v;
                s.Tf {i} = Tf;
                s.xf {i} = xf;
            end
        end
    end
end