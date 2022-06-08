classdef A2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up structure "msh".
        function [msh] = Set_msh(h,t)
            %  > Auxiliary variables.
            inp = A1_1D.Set_msh(h,t);
            
            %  > "msh".
            Xv  = A2_1D.Set_msh_1(inp);
            msh = A2_1D.Set_msh_2(inp,Xv);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set grid vertices.
        %  > switch inp.t
        %  >    case 1, uniformly spaced grid.
        %  >    case 2, smooth non-uniform grid.
        %  >    case 3, highly distorted non-uniform grid.
        %  > end
        function [Xv] = Set_msh_1(inp)
            %  > Auxiliary variables.
            Xi = inp.Lim(1)+(0:diff(inp.Lim)./inp.h).*inp.h;
            
            %  > Select grid type...
            switch inp.t
                case 1
                    Xv      = Xi;
                case 2
                    i (1)   = 1-(1-exp( inp.x(1))).*inp.x(2);
                    i (2)   = 1-(1-exp(-inp.x(1))).*inp.x(2);
                    B       = 1./(2.*inp.x(1)).*log(i(1)./i(2));
                    j (1,:) = inp.x(1).*(Xi-B);
                    j (2,:) = inp.x(1).*B;
                    Xv      = inp.x(2).*(1+sinh(j(1,:))./sinh(j(2,:)));
                case 3
                    rng default;
                    Nv      = numel(Xi);
                    i       = [1,Nv];
                    j       = src_Tools.setdiff(1:Nv,i);
                    Xv(i)   = Xi(i);
                    Xv(j)   = Xi(j)+(rand(1,Nv-2)-0.5).*inp.h.^inp.x(1);
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Set remaining grid properties.
        function [msh] = Set_msh_2(inp,Xv)
            %  > Field: "c".
            Nc       = numel(Xv)-1;
            msh.c.Nc = Nc;
            for i = 1:Nc
                msh.c.f.if  (i,:) = [i,i+1];
                msh.c.f.Sf  (i,:) = [-1,1];
                msh.c.Xc    (i,1) = 1./2.*(Xv(i)+Xv(i+1));
                msh.c.Volume(i,1) = Xv(i+1)-Xv(i);
            end
            %  > Field: "d".
            msh.d.h = diff(inp.Lim)./Nc;
            %  > Field: "f".
            Nf       = Nc+1;
            msh.f.Nf = Nf;
            for i = 1:Nf
                msh.f.ic{i,1} = [i-1,i];
                if i == 1 || i == Nf
                    msh.f.ic{i,1} = setdiff(msh.f.ic{i,1},[0,Nf]);
                end
                msh.f.Xv(i,1) = Xv(i);
            end
            %  > Sort...
            msh = src_Tools.Sort_1D_msh(msh);
        end
    end
end