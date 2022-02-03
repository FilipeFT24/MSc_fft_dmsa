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
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_2(Xv_i,Xv_f,h,Ks_X,'East');
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
                case 'West'
                    %  > Do nothing.
                case 'East'
                    %  > Invert mesh.
                    fl_x = flip(Xd_x);
                    Xd_x = ones(1,length(fl_x)).*Xv_i+Xv_f-fl_x;
                otherwise
                    return;
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [sc,sf] = SetUp_Stencil_f(i,msh,np)
            %  > Auxiliary variables.
            Xc = msh.c.Xc;
            NC = msh.c.NC;
            Xf = msh.f.Xv;
            NF = msh.f.NF;
            
            %  > Add faces...
            np_2 =  np./2;
            if i <= np_2
                sf = 1;
            elseif i >= NF-np_2+1
                sf = NF;
            else
                sf = [];
            end
            %  > Add cells...
            df     = Xf(i)-Xc;
            [~,ic] = min(abs(df));
            vc     = Xc(ic);
            
            if i <= np_2
                sc = 1:1:np-1;
            elseif i >= NF-np_2+1
                sc = NC-np+2:1:NC;
            else
                sc = i-np_2:1:i+np_2-1;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [xc,xf,xt] = Compute_Coordinates_cft(msh,sc,sf)
            xc = msh.c.Xc(sc);
            xf = msh.f.Xv(sf);
            xt = [xc,xf];
        end
        % >> 2.3. ---------------------------------------------------------
        function [Phi_f,GradPhi_f,x_f_v,x_f_g] = x_f(v,g,np,stl_x,xf)
            df         = zeros(2,np);
            df  (1,1)  = 1; % > Phi_f.
            df  (2,2)  = 1; % > gradPhi_f.
            j          = 1:np;
            k          = 1:np;
            Df  (j,k)  = (stl_x(j)-xf)'.^(k-1);
            Phi_f      = df(1,:)*inv(Df);
            GradPhi_f  = df(2,:)*inv(Df);
            x_f_v      = v.*Phi_f;
            x_f_g      = g.*GradPhi_f;
            x_f        = x_f_v-x_f_g;
        end
    end
end