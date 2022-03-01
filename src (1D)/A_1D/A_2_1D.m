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
            msh.c.NC       = NX_c;
            msh.f.NF       = NX_c+1;
            %  > Cell vertex coordinates.
            msh.f.Xv       = Xd_x;
            %  > Xc/Cell volume.
            i              = 1:msh.c.NC;
            msh.c.Xc (i)   = 1./2.*(msh.f.Xv(i)+msh.f.Xv(i+1));
            msh.c.Vol(i,1) = msh.f.Xv(i+1)-msh.f.Xv(i);
            %  > Reference length.
            msh.d.H_ref    = (msh.f.Xv(msh.f.NF)-msh.f.Xv(1))./msh.c.NC;
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
        %  > 2.1.1. -------------------------------------------------------
        %  > Initialize 'stl' structure.
        function [p,s,t] = Initialize_stl(msh,t1,t2)
            p = repelem(t2,msh.f.NF);
            s = 1:msh.f.NF;
            t = repelem(t1,msh.f.NF);
        end
        %  > 2.1.2. -------------------------------------------------------
        %  > Compute method's order.
        function [p] = Compute_p(stl_p,stl_t)
            switch stl_t
                case "CDS"
                    p = 2.*stl_p;
                otherwise
                    p = 2.*stl_p-1;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [s] = Assemble_stl(obj,msh,pde,s,stl,stl_s)
            %  > Auxiliary variables.
            Xc      = msh.c.Xc;
            NC      = msh.c.NC;
            Xv      = msh.f.Xv;
            NF      = msh.f.NF;
            f (:,1) = pde.a.f(:,1);
            f (:,2) = pde.a.f(:,2);
            vg  (1) = obj.v;
            vg  (2) = obj.g;
            
            % >> Loop through selected faces...
            for n = 1:size(stl_s,2)
                if isempty(stl_s{n})
                    continue;
                else
                    for q = 1:size(stl_s{n},1)
                        %  > Auxiliary variables.
                        o     = stl_s{n}(q);
                        p     = stl.p{n}(o);
                        add_f = false;
                        switch stl.t{n}(o)
                            case "CDS"
                                LHS = p;
                                RHS = NF+1-p;
                            case "UDS"
                                LHS = p;
                                RHS = NF+1-(p-1);
                            case "DDS"
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
                            bnd_i = obj.bnd(1);
                            %  > ...cell(s).
                            switch stl.t{n}(o)
                                case "CDS"
                                    stl_c = 1:2.*p-1;
                                otherwise
                                    stl_c = 1:2.*(p-1);
                            end
                        elseif o >= RHS
                            %  > ...face(s).
                            add_f = true;
                            bnd_f = NF;
                            stl_f = bnd_f;
                            bnd_i = obj.bnd(2);
                            %  > ...cell(s).
                            switch stl.t{n}(o)
                                case "CDS"
                                    stl_c = NF+1-2.*p:NC;
                                otherwise
                                    stl_c = NF-2.*(p-1):NC;
                            end
                        else
                            %  > ...face(s).
                            stl_f = [];
                            bnd_i = string([]);
                            %  > ...cell(s).
                            switch stl.t{n}(o)
                                case "CDS"
                                    stl_c = o-p:o+p-1;
                                case "UDS"
                                    stl_c = o-p:o-1+(p-1);
                                case "DDS"
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
                                g_v     = obj.g./obj.v;
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
                                return;
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
        function [ttm] = Compute_ttm(obj,msh,s,stl_s,p,nt,dfn)
            for i = 1:size(stl_s,2)
                if isempty(stl_s{i})
                    continue;
                else
                    for j = 1:size(stl_s{i},1)
                        %  > Auxiliary variables.
                        k   = stl_s{i}(j);
                        xt  = s.xt{i,k};
                        fx  = msh.f.Xv(k);
                        
                        %  > Df_T.
                        bnd_i = s.bnd_i{i,j};
                        switch bnd_i
                            case "Neumann"
                                l         = 1:nt(i);
                                m         = 1:length(xt);
                                n         = 1:length(xt)-1;
                                o         = p(i)-1+l;
                                q         = n-1;
                                r         = length(xt);
                                t         = o-1;
                                Df  (n,l) = (xt(n)-fx)'.^o;
                                Df  (r,l) = (xt(r)-fx)'.^t.*o;
                                Df_T(l,m) = transpose(Df(m,l));
                            case "Robin"
                                g_v       = obj.g./obj.v;
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
                        ttm{i}(k,l) = transpose(Df_T(l,m)*s.xf{i,k}').*dfn{i}(k,l);
                    end
                end
                ttm{i} = abs(ttm{i});
            end
        end
    end
end