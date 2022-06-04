classdef A2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up 'msh' structure.
        function [msh] = Set_msh(h)
            %  > Auxiliary variables.
            inp  = A1_1D.Set_inp_1(h);
            u    = inp.m.Uniform;
            h    = inp.m.h;
            XLim = inp.m.XLim;
            Nc   = round(1./h.*(XLim(2)-XLim(1)));
            
            if u
                Xv  = A2_1D.msh_1(XLim,Nc);
            else
                Xv  = A2_1D.msh_2(XLim,Nc,inp.m.A,inp.m.c);
            end
            msh     = A2_1D.Set_Grid(XLim,Xv);
            msh.d.h = h;
        end
  
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Uniformly spaced grid.
        function [Xv] = msh_1(XLim,Nc)
            Xv = linspace(XLim(1),XLim(2),Nc+1);
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Smooth non-uniform grid (transformation taken from reference).
        function [Xv] = msh_2(XLim,Nc,A,c)
            T  = (0:Nc)./Nc;
            B  = 1./(2.*A).*log((1+(exp(A)-1).*c)./(1-(1-exp(-A)).*c));
            Xv = c.*(1+sinh(A.*(T-B))./sinh(A.*B));
            Xv = XLim(1)+diff(XLim).*Xv;
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Highly distorted grid (transformation taken from reference).
        function [Xv] = msh_3(XLim,Nc)
            q     = 1.1;
            Xv    = (0:Nc)./Nc;
            i     = setdiff(1:Nc,[1,Nc+1]);
            Xv(i) = Xv(i)+(rand(1,numel(i))-0.5).*Nc.^(-q);
        end
        
        % >> 2.4. ---------------------------------------------------------
        %  > Set remaining grid properties.
        function [msh] = Set_Grid(XLim,Xv)
            %  > Field: 'c'.
            Nc       = length(Xv)-1;
            msh.c.Nc = Nc;
            for i = 1:msh.c.Nc
                msh.c.f.if(i,:) = [i,i+1];
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
                msh.f.ic{i,1} = [i-1,i];
                if i == 1 || i == Nf
                    msh.f.ic{i,1} = setdiff(msh.f.ic{i,1},[0,Nf]);
                end
            end            
            j             = 1:msh.f.Nf;
            msh.f.Xv(j,1) = Xv;
        end
    end
end