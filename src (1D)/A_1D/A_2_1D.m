classdef A_2_1D
    methods (Static)
        %% > Wrap-up A_2 (1D).
        function [inp,msh] = WrapUp_A_2_1D(inp)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            h    = inp.msh.h;
            t    = inp.msh.t;
            
            switch t
                case "1"
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_Uniform(Xv_i,Xv_f,h);
                case "2"
                    Nf_X        = inp.msh.Nf_X;
                    Ks_X        = inp.msh.Ks_X;
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_1(Xv_i,Xv_f,h,Nf_X,Ks_X);
                case "3"
                    Ks_X        = inp.msh.Ks_X;
                    Lc_X        = inp.msh.Lc_X;
                    [NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_2(Xv_i,Xv_f,h,Ks_X,Lc_X);
                otherwise
                    return;
            end
            % >> Grid properties.
            %  > Number of cells/faces.
            msh.c.NC       = NX_c;
            msh.f.NF       = NX_c+1;
            %  > Cell vertex coordinates.
            i              = 1:msh.f.NF;
            msh.f.Xv(i,1)  = Xd_x;
            %  > Xc/Cell volume.
            j              = 1:msh.c.NC;
            msh.c.Xc (j,1) = 1./2.*(msh.f.Xv(j)+msh.f.Xv(j+1));
            msh.c.Vol(j,1) = msh.f.Xv(j+1)-msh.f.Xv(j);
            %  > Reference length.
            msh.d.h_ref    = (msh.f.Xv(msh.f.NF)-msh.f.Xv(1))./msh.c.NC;
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [NX_c,Xd_x] = SquareMesh_Uniform(Xv_i,Xv_f,h)
            % >> Vertex coordinates.
            NX_c = round(1./h.*(Xv_f-Xv_i));
            Xd_x = linspace(Xv_i,Xv_f,NX_c+1);
        end
        % >> 1.2. ---------------------------------------------------------
        function [NX_c,Xd_x] = SquareMesh_NonUniform_1(Xv_i,Xv_f,h,Nf_X,Ks_X)
            % >> Transformation parameters/vertex coordinates.
            NX_v     = round(1./h.*(Xv_f-Xv_i))+1;
            Nf_Unf_X = linspace(0,1,NX_v);
            Pt_X     = (Nf_X-0)./(1-0);
            B_X      = 1./(2.*Ks_X).*log((1+(exp(Ks_X)-1).*(Pt_X))./(1+(exp(-Ks_X)-1).*(Pt_X)));
            i        = 1:NX_v;
            NF_X     = (Nf_X-0).*(1+sinh(Ks_X.*(Nf_Unf_X(i)-B_X))./sinh(Ks_X.*B_X))+0;
            Xd_x     = NF_X(1:NX_v).*(Xv_f-Xv_i);
            NX_c     = NX_v-1;
        end
        % >> 1.3. ---------------------------------------------------------
        function [NX_c,Xd_x] = SquareMesh_NonUniform_2(Xv_i,Xv_f,h,Ks_X,bnd)
            % >> Transformation parameters/vertex coordinates.
            NX_v   = round(1./h.*(Xv_f-Xv_i))+1;
            Nf_Unf = linspace(0,1,NX_v);
            Bp     = Ks_X+1;
            Bm     = Ks_X-1;
            i      = 1:NX_v;
            X  (i) = (Bp./Bm).^(1-Nf_Unf(i));
            Num(i) = Bp-Bm.*X(i);
            Den(i) = X(i)+1;
            x_h(i) = Num(i)./Den(i);
            Xd_x   = Xv_i+(Xv_f-Xv_i).*x_h;
            NX_c   = NX_v-1;
            
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
        %  > 2.1.1. -------------------------------------------------------
        %  > Initialize 'stl' structure.
        function [stl] = Initialize_stl(msh,p,s)
            for i = 1:length(p)
                stl.p   (:,i) = repelem(p(i),msh.f.NF);
                stl.s{i}(:,1) = 1:msh.f.NF;
                stl.t   (:,i) = repelem(s(i),msh.f.NF);
            end
        end
        %  > 2.1.2. -------------------------------------------------------
        %  > Compute method's order.
        function [p] = Compute_p(stl_p,stl_t)
            for i = 1:length(stl_p)
                switch stl_t(i)
                    case "c"
                        p(i) = 2.*stl_p(i);
                    otherwise
                        p(i) = 2.*stl_p(i)-1;
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [s] = Assemble_stl(inp,msh,pde,s,stl,stl_s)
            %  > Auxiliary variables.
            Xc      = msh.c.Xc';
            NC      = msh.c.NC;
            Xv      = msh.f.Xv;
            NF      = msh.f.NF;
            f (:,1) = pde.av.f(:,1);
            f (:,2) = pde.av.f(:,2);
            vg  (1) = inp.ps.v  (1);
            vg  (2) = inp.ps.v  (2);
            
            % >> Loop through selected faces...
            for n = 1:size(stl_s,2)
                if isempty(stl_s{n})
                    continue;
                else
                    for q = 1:size(stl_s{n},1)
                        %  > Auxiliary variables.
                        o     = stl_s{n}(q);
                        p     = stl.p(o,n);
                        add_f = false;
                        switch stl.t(o,n)
                            case "c"
                                LHS = p;
                                RHS = NF+1-p;
                            case "u"
                                LHS = p;
                                RHS = NF+1-(p-1);
                            case "d"
                                LHS = p-1;
                                RHS = NF+1-p;
                            otherwise
                                return;
                        end
                        %  > Add...
                        if o <= LHS
                            %  > ...face(s).
                            add_f = true;
                            bnd_f = 1;
                            stl_f = bnd_f;
                            bnd_i = inp.ps.b(1);
                            %  > ...cell(s).
                            switch stl.t(o,n)
                                case "c"
                                    stl_c = 1:2.*p-1;
                                otherwise
                                    stl_c = 1:2.*(p-1);
                            end
                        elseif o >= RHS
                            %  > ...face(s).
                            add_f = true;
                            bnd_f = NF;
                            stl_f = bnd_f;
                            bnd_i = inp.ps.b(2);
                            %  > ...cell(s).
                            switch stl.t(o,n)
                                case "c"
                                    stl_c = NF+1-2.*p:NC;
                                otherwise
                                    stl_c = NF-2.*(p-1):NC;
                            end
                        else
                            %  > ...face(s).
                            stl_f = [];
                            bnd_i = string([]);
                            %  > ...cell(s).
                            switch stl.t(o,n)
                                case "c"
                                    stl_c = o-p:o+p-1;
                                case "u"
                                    stl_c = o-p:o-1+(p-1);
                                case "d"
                                    stl_c = o-(p-1):o+(p-1);
                                otherwise
                                    return;
                            end
                        end
                        %  > Compute stencil coordinates.
                        xc = Xc(stl_c);
                        if ~add_f
                            xt = xc;
                        else
                            xt = [xc,Xv(bnd_f)];
                        end
                        
                        % >> Tf.
                        %  > Auxiliary variables.
                        len = length(xt);
                        j   = 1:len;
                        fx  = Xv(o);
                        %  > Df.
                        Df  = zeros(len,len);
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
                                g_v     = inp.g./inp.v;
                                bnd_v   = f(bnd_f,1)+g_v.*f(bnd_f,2);
                                %  > Df.
                                k       = 1:len-1;
                                Df(k,j) = (xt(k)-fx)'.^(j-1);
                                l       = len;
                                lv  (j) = (xt(l)-fx)'.^(j-1);
                                m       = k+1;
                                lg  (1) = 0;
                                lg  (m) = k.*(xt(l)-fx).^(k-1);
                                Df(l,j) = lv(j)+g_v*lg(j);
                            otherwise
                                %  > Boundary value.
                                bnd_v   = [];
                                %  > Df.
                                Df(j,j) = (xt(j)-fx)'.^(j-1);
                        end
                        if len == 1
                            %  > 1st order UDS/DDS cannot be used to discretize the diffusive term.
                            if n == 2
                                warndlg('1st order UDS/DDS cannot be used to discretize the diffusive term.');
                                break;
                            else
                                df   = 1;
                                Inv  = inv(Df);
                                Tf   = df*Inv;
                                xf   = Tf(n,:);
                            end
                        else
                            df       = zeros(1,len);
                            df(1,n)  = 1;
                            Inv      = inv(Df);
                            xf       = df*Inv;
                        end
                        %  > Update 'msh' structure...
                        s.c    {n,o} = stl_c;
                        s.f    {n,o} = stl_f;
                        s.bnd_i{n,o} = bnd_i;
                        s.bnd_v{n,o} = bnd_v;
                        s.xf   {n,o} = xf;
                        s.xt   {n,o} = xt;
                        s.Inv  {n,o} = Inv;
                    end
                end
            end
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Evaluate truncated error terms' magnitude (re-assemble Df).
        function [tm] = Compute_tm(inp,msh,s,stl,df,p,nt) 
            for i = 1:size(stl.s,2)
                if isempty(stl.s{i})
                    continue;
                else
                    for j = 1:size(stl.s{i},1)
                        %  > Auxiliary variables.
                        k  = stl.s   {i}(j);
                        xt = s.xt    {i,k};
                        fx = msh.f.Xv(k);
                        
                        %  > Df_T.
                        bnd_i = s.bnd_i{i,j};
                        switch bnd_i
                            case "Neumann"
                                l         = 1:nt(i);
                                m         = 1:length(xt);
                                n         = 1:length(xt)-1;
                                o         = stl.o(i)-1+l;
                                q         = n-1;
                                r         = length(xt);
                                t         = o-1;
                                Df  (n,l) = (xt(n)-fx)'.^o;
                                Df  (r,l) = (xt(r)-fx)'.^t.*o;
                                Df_T(l,m) = transpose(Df(m,l));
                            case "Robin"
                                g_v       = inp.g./inp.v;
                                l         = 1:nt(i);
                                m         = 1:length(xt);
                                n         = 1:length(xt)-1;
                                o         = p(i)-1+l;
                                q         = n-1;
                                r         = length(xt);
                                t         = o-1;
                                Df  (m,l) = (xt(m)-fx)'.^o;
                                Df  (r,l) = Df(r,l)+g_v.*(xt(r)-fx)'.^t.*o;
                                Df  (l,m) = transpose(Df(m,l));
                            otherwise
                                l         = 1:nt(i);
                                m         = 1:length(xt);
                                n         = p(i)-1+l;
                                Df  (m,l) = (xt(m)-fx)'.^n;
                                Df_T(l,m) = transpose(Df(m,l));
                        end
                        %  > Truncated terms.
                        tm{i}(k,l) = transpose(Df_T(l,m)*s.xf{i,k}').*df{i}(k,l);
                    end
                end
                tm{i} = abs(tm{i});
            end
        end
    end
end