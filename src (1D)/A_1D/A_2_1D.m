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
                case '1'
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_Uniform(Xv_i,Xv_f,h);
                case '2'
                    Nf_X        = inp.msh.s_nu.Nf_X;
                    Ks_X        = inp.msh.s_nu.Ks_X;
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_1(Xv_i,Xv_f,h,Nf_X,Ks_X);
                case '3'
                    Ks_X        = inp.msh.s_nu.Ks_X;
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_2(Xv_i,Xv_f,h,Ks_X,'e');
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
                case 'w'
                    %  > Do nothing.
                case 'e'
                    %  > Invert mesh.
                    fl_x = flip(Xd_x);
                    Xd_x = ones(1,length(fl_x)).*Xv_i+Xv_f-fl_x;
                otherwise
                    return;
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [sc,sf,bnd] = SetUp_Stencil_f(msh,i,np,bnd_w,bnd_e)
            %  > Auxiliary variables.
            Xc   = msh.c.Xc;
            NC   = msh.c.NC;
            Xf   = msh.f.Xv;
            NF   = msh.f.NF;
            np_2 = np./2;
            
            % >> Add...
            %  > ...face(s).
            if i <= np_2
                sf  = 1;
                bnd = bnd_w;
            elseif i >= NF-np_2+1
                sf  = NF;
                bnd = bnd_e;
            else
                sf  = [];
                bnd = [];
            end
            %  > ...cell(s).
            df     = Xf(i)-Xc;
            [~,ic] = min(abs(df));
            vc     = Xc(ic);
            if i  <= np_2
                sc = 1:1:np-1;
            elseif i >= NF-np_2+1
                sc = NC-np+2:1:NC;
            else
                sc = i-np_2:1:i+np_2-1;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [xc,xf,xt] = Compute_Coordinates_cft(Xc,sc,Xv,sf)
            xc = Xc(sc);
            if isempty(sf)
                xf = [];
                xt = xc;
            else
                xf = Xv(sf);
                xt = [xc,xf];
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > 2.3.1.--------------------------------------------------------
        function [df] = Compute_df(np)
            %  > φ(x) = C1 + C2(x-xf)^1 +  C3(x-xf)^2 +  C4(x-xf)^3 +   C5(x-xf)^4 +    C6(x-xf)^5.
            %  > φ(x) = 0  + C2         + 2C3(x-xf)^1 + 3C4(x-xf)^2 +  4C5(x-xf)^3 +   5C6(x-xf)^4.
            %  > φ(x) = 0  + 0          + 2C3         + 6C4(x-xf)^1 + 12C5(x-xf)^2 +  20C6(x-xf)^3.
            %  > φ(x) = 0  + 0          + 0           + 6C4         + 24C5(x-xf)^1 +  60C6(x-xf)^2.
            %  > φ(x) = 0  + 0          + 0           + 0           + 24C5         + 120C6(x-xf)^1.
            %  > φ(x) = 0  + 0          + 0           + 0           +  0           + 120C6        .
            %  > ...
            %  > Initialize.
            df = zeros(np);
            
            for i = 1:np
                if i == 1
                    df(i,i) = 1;
                else
                    df(i,:) = circshift(df(i-1,:),1).*(i-1);
                end
            end
        end
        %  > 2.3.2.--------------------------------------------------------
        function [Tf,x_f,bnd_v] = x_f(fv,dfv,np,stlx,fx,bnd_i,bnd_t,v,g)
            df = A_2_1D.Compute_df(np);
            j  = 1:np;
            
            switch char(bnd_t)
                case 'Dirichlet'
                    %  > Boundary value.
                    bnd_v   = fv(bnd_i);
                    %  > Df.
                    Df(j,j) = (stlx(j)-fx)'.^(j-1);
                case 'Neumann'
                    %  > Boundary value.
                    bnd_v   = dfv(bnd_i);
                    %  > Df.
                    k       = 1:np-1;
                    Df(k,j) = (stlx(k)-fx)'.^(j-1);
                    l       = np;
                    m       = k+1;
                    Df(l,1) = 0;
                    Df(l,m) = k.*(stlx(l)-fx).^(k-1);
                case 'Robin'
                    %  > Boundary value.
                    g_v     = g./v;
                    bnd_v   = fv(bnd_i)-g_v.*dfv(bnd_i);
                    %  > Df.
                    k       = 1:np-1;
                    Df(k,j) = (stlx(k)-fx)'.^(j-1);
                    l       = np;
                    lv  (j) = (stlx(l)-fx)'.^(j-1);
                    m       = k+1;
                    lg  (1) = 0;
                    lg  (m) = k.*(stlx(l)-fx).^(k-1);
                    Df(l,j) = lv(j)-g./v.*lg(j);
                otherwise
                    %  > Boundary value.
                    bnd_v   = [];
                    %  > Df.
                    Df(j,j) = (stlx(j)-fx)'.^(j-1);
            end
            Inv  = inv(Df);
            Tf   = df*Inv;
            x_fv = v.*Tf(1,:);
            x_fg = g.*Tf(2,:);
            x_f  = x_fv-x_fg;
        end
    end
end