classdef A_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Uniformly spaced grid.
        function [Xv] = msh_1(XLim,Nc)
            Xv = linspace(XLim(1),XLim(2),Nc+1);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Smooth non-uniform grid (transformation taken from reference).
        function [Xv] = msh_2(XLim,Nc,A,c)
            T  = (0:Nc)./Nc;
            B  = 1./(2.*A).*log((1+(exp(A)-1).*c)./(1-(1-exp(-A)).*c));
            Xv = c.*(1+sinh(A.*(T-B))./sinh(A.*B));
            Xv = XLim(1)+diff(XLim).*Xv;
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Set remaining grid properties.
        function [msh] = Set_Grid(XLim,Xv)
            %  > Field: 'c'.
            Nc       = length(Xv)-1;
            msh.c.Nc = Nc;
            for i = 1:msh.c.Nc
                msh.c.f.f {i,1} = [i,i+1];
                msh.c.f.Nf{i,1} = [-1,1];
            end
            j             = 1:msh.c.Nc;
            msh.c.Xc(j,1) = 1./2.*(Xv(j)+Xv(j+1));
            msh.c.Vc(j,1) = Xv(j+1)-Xv(j);
            %  > Field: 'd'.
            msh.d.h       = 1./Nc.*(XLim(2)-XLim(1));
            %  > Field: 'f'.
            Nf            = Nc+1;
            msh.f.Nf      = Nf;
            for i = 1:Nf
                msh.f.c{i,1} = [i-1,i];
                if i == 1 || i == Nf
                    msh.f.c{i,1} = setdiff(msh.f.c{i,1},[0,Nf]);
                end
            end            
            j             = 1:msh.f.Nf;
            msh.f.Xv(j,1) = Xv;
        end

        %% > 2. -----------------------------------------------------------
        %  > Initialize 'upd' structure.
        function [upd] = Initialize_upd(msh,p)         
            for i = 1:length(p)
                upd.p    (:,i) = repelem(p(i),msh.f.Nf);
                upd.s{i} (:,1) = 1:msh.f.Nf;
            end
        end
    end
end