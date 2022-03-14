classdef A_2_1D
    methods (Static)
        function [inp,msh] = Set_A(h)
            %  > Set 'inp' structure.
            inp  = A_1_1D.Set_inp(h);
            %  > Auxiliary variables.
            u    = inp.msh.Uniform;
            XLim = inp.msh.XLim;
            h    = inp.msh.h;
            NC   = round(1./h.*(XLim(2)-XLim(1)));
            
            switch u
                case true
                    Xv = A_2_1D.msh_1(XLim,NC);
                case false
                    Xv = A_2_1D.msh_2(XLim,NC,inp.msh.A,inp.msh.c);
                otherwise
                    return;
            end
            msh = A_2_1D.Set_Grid(XLim,Xv);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Uniformly spaced grid.
        function [Xv] = msh_1(XLim,NC)
            Xv = linspace(XLim(1),XLim(2),NC+1);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Smooth non-uniform grid (transformation taken from reference in PDF).
        function [Xv] = msh_2(XLim,NC,A,c)
            T  = (0:NC)./NC;
            B  = 1./(2.*A).*log((1+(exp(A)-1).*c)./(1-(1-exp(-A)).*c));
            Xv = c.*(1+sinh(A.*(T-B))./sinh(A.*B));
            Xv = XLim(1)+diff(XLim).*Xv;
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Set remaining grid properties.
        function [msh] = Set_Grid(XLim,Xv)
            %  > Field: 'c'.
            NC            = length(Xv)-1;
            msh.c.NC      = NC;
            i             = 1:msh.c.NC;
            msh.c.Xc(i,1) = 1./2.*(Xv(i)+Xv(i+1));
            msh.c.Vc(i,1) = Xv(i+1)-Xv(i);
            %  > Field: 'd'.
            msh.d.href    = 1./NC.*(XLim(2)-XLim(1));
            %  > Field: 'f'.
            msh.f.NF      = NC+1;
            j             = 1:msh.f.NF;
            msh.f.Xv(j,1) = Xv;
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
        function [p] = Compute_p(stl_p,stl_s)
            [m,n] = size(stl_p);
            for i = 1:m
                for j = 1:n
                    switch stl_s(i,j)
                        case "c"
                            p(i,j) = 2.*stl_p(i,j);
                        otherwise
                            p(i,j) = 2.*stl_p(i,j)-1;    
                    end
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
            vg  (1) = inp.pv.v  (1);
            vg  (2) = inp.pv.v  (2);
            
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
                            bnd_i = inp.pv.b(1);
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
                            bnd_i = inp.pv.b(2);
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
        %  > Evaluate truncated error terms' magnitude w/ analytic solution (re-assemble Df).
        function [TTM] = Compute_TTM(inp,msh,s,stl,df)
            %  > Analytic derivatives.
            nt = inp.pl.nt;

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
                                o         = s.p(i)-1+l;
                                q         = n-1;
                                r         = length(xt);
                                t         = o-1;
                                Df  (m,l) = (xt(m)-fx)'.^o;
                                Df  (r,l) = Df(r,l)+g_v.*(xt(r)-fx)'.^t.*o;
                                Df  (l,m) = transpose(Df(m,l));
                            otherwise
                                l         = 1:nt(i);
                                m         = 1:length(xt);
                                n         = s.p(i)-1+l;
                                Df  (m,l) = (xt(m)-fx)'.^n;
                                Df_T(l,m) = transpose(Df(m,l));
                        end
                        %  > Term's magnitude.
                        TTM{i}(k,l) = transpose(Df_T(l,m)*s.xf{i,k}').*df{i}(k,l);
                    end
                end
                TTM{i} = abs(TTM{i});
            end
        end
    end
end