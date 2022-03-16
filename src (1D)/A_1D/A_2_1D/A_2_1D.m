classdef A_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Uniformly spaced grid.
        function [Xv] = msh_1(XLim,NC)
            Xv = linspace(XLim(1),XLim(2),NC+1);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Smooth non-uniform grid (transformation taken from reference).
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
        %  > Initialize 'stl' structure.
        function [stl] = Initialize_stl(msh,p,t)
            %  > Auxiliary variables.
            l  = length(p);
            NF = msh.f.NF;
                       
            for i = 1:l
                j = i*l-1;
                k = i*l;
                if ~rem(p(i),2)
                    stl.p(:,k) = ones (NF,1).*t(i);
                else
                    stl.p(:,k) = zeros(NF,1);
                end
                stl.p    (:,j) = repelem(p(i),NF);
                stl.s{i} (:,1) = 1:NF;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        %  > Set 's' structure.
        function [s] = Set_s(inp,msh,pde,s,stl)
            %  > Auxiliary variables.
            Xc = msh.c.Xc;
            Xv = msh.f.Xv;
            NF = size(Xv,1);
            av = pde.av;
            nv = size(stl.s,2);

            % >> Loop through selected faces...
            for i = 1:nv
                j = i.*nv-1;
                k = i.*nv;
                if isempty(stl.s{i})
                    continue;
                else
                    for l = 1:size(stl.s{i},1)
                        m = stl.s{i}(l);
                        %  > Check polynomial order (odd/even).
                        nb =  ceil(stl.p(m,j)./2);
                        nl = -nb;
                        nr =  nb; 
                        if ~rem(stl.p(m,i),2)
                            if stl.p(m,k) < 0
                                %  > UDS.
                                nl = nl+stl.p(m,k);
                            else
                                %  > DDS.
                                nr = nr+stl.p(m,k);
                            end
                        end
                        ns = m+nl:m+nr-1;
                        %  > Check whether stencil reaches any of the boundaries and shift it accordingly.
                        %  > WB.
                        if any(ns < 0)
                            nw = mcount(ns,0,'<');
                            ns = ns+nw;
                        end
                        %  > EB.
                        if any(ns > NF)
                            ne = mcount(ns,NF,'>');
                            ns = ns-ne;
                        end
                        %  > Set up stencil/assemble matrix Df.
                        stl_f = A_2_1D.Set_stl_f  (inp,ns,Xc,Xv);
                        stl_f = A_2_1D.Assemble_Df(stl_f,Xv(m),av.f);
                        
                        %  > Update structure 's'.
                        Df = stl_f.Df;
                        lt = length(stl_f.xt);
                        if lt == 1
                            if i == 2
                                warndlg('1st order UDS/DDS cannot be used to discretize the diffusive term.');
                                break;
                            else
                                df   = 1;
                                Inv  = inv(Df);
                                Tf   = df*Inv;
                                xf   = Tf(n,:);
                            end
                        else
                            df       = zeros(1,lt);
                            df(1,i)  = 1;
                            Inv      = inv(Df);
                            xf       = df*Inv;
                        end
                        s.c    {i,m} = stl_f.sc;
                        s.f    {i,m} = stl_f.sf;
                        s.bnd_i{i,m} = stl_f.vi;
                        s.bnd_v{i,m} = stl_f.vf;
                        s.xf   {i,m} = xf;
                        s.xt   {i,m} = stl_f.xt;
                        s.Inv  {i,m} = Inv;
                    end
                end
            end
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Set 'stl_f' structure.
        function [stl_f] = Set_stl_f(inp,ns,Xc,Xv)
            %  > Auxiliary variables.
            NF   = size(Xv,1);
            is_f = ismembc(ns,[0,NF]);
            
            if any(is_f)
                if ns(is_f) == 0
                    sc = ns(ns ~= 0);
                    sf = 1;
                    vi = inp.pv.b(1);
                else
                    sc = ns(ns ~= NF);
                    sf = NF;
                    vi = inp.pv.b(2);
                end
                xt = [Xc(sc);Xv(sf)]';
            else
                sc = ns;
                sf = [];
                vi = [];
                xt = Xc(sc)';
            end
            stl_f.sc = sc;
            stl_f.sf = sf;
            stl_f.vi = vi;
            stl_f.xt = xt;
        end
        %  > 2.2.3. -------------------------------------------------------
        %  > Update matrix 'stl_f' (assemble matrix 'Df').
        function [stl_f] = Assemble_Df(stl_f,fx,func)
            %  > Auxiliary variables.
            vi = stl_f.vi;
            xt = stl_f.xt;
            lt = length(xt);
            j  = 1:lt;
            Df = zeros(lt);
            
            switch string(vi)
                case "Dirichlet"
                    Df(j,j) = (xt(j)-fx)'.^(j-1);
                    vf      = func(stl_f.sf,1);
                case "Neumann"
                    k       = 1:lt-1;
                    Df(k,j) = (xt(k)-fx)'.^(j-1);
                    l       = lt;
                    m       = k+1;
                    Df(l,1) = 0;
                    Df(l,m) = k.*(xt(l)-fx).^(k-1);
                    vf      = func(stl_f.sf,2);
                case "Robin"
                    g_v     = inp.g./inp.v;
                    k       = 1:len-1;
                    Df(k,j) = (xt(k)-fx)'.^(j-1);
                    l       = len;
                    lv  (j) = (xt(l)-fx)'.^(j-1);
                    m       = k+1;
                    lg  (1) = 0;
                    lg  (m) = k.*(xt(l)-fx).^(k-1);
                    Df(l,j) = lv(j)+g_v*lg(j);
                    vf      = f(stl_f.sf,1)+g_v.*f(stl_f.sf,2);
                otherwise
                    Df(j,j) = (xt(j)-fx)'.^(j-1);
                    vf      = [];
            end
            stl_f.Df = Df;
            stl_f.vf = vf;
        end
        %  > 2.2.4. -------------------------------------------------------
        %  > Evaluate truncated error terms' magnitude w/ analytic solution (re-assemble Df).
        function [TTM] = Compute_TTM(inp,msh,s,stl,df)
            %  > Analytic derivatives.
            nt = inp.pl.nt;
            ls = size(stl.s,2);

            for i = 1:ls
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
                        switch string(bnd_i)
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
                                n         = stl.p(i*ls-1)+l;
                                Df  (m,l) = (xt(m)-fx)'.^n;
                                Df_T(l,m) = transpose(Df(m,l));
                        end
                        TTM{i}(k,l) = transpose(Df_T(l,m)*s.xf{i,k}').*df{i}(k,l);
                    end
                end
                TTM{i} = abs(TTM{i});
            end
        end
    end
end