classdef B1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "s" (stencil cell/face indices, coordinates, coefficents, etc.).
        function [s] = Initialize_s(inp,f,nc,ns,Nc,Nf,Xc,Xv)
            A = [0,2];
            for i = 1:ns
                %  > Fiels: "i", "logical" and "xt".
                s{i}.i       = cell(Nf,nc);
                s{i}.logical = cell(Nf,nc);
                s{i}.xt      = cell(Nf,nc);
                %  > Field: "u".
                for j = 1:numel(inp.p.p)
                    s{i}.u.p   (:,j) = repelem(inp.p.p(j)+A(i),Nf,1);
                    s{i}.u.s{j}(:,1) = 1:Nf;
                end
                %  > Field: "x".
                for j = ["a","x"]
                    s{i}.x.vf. (j) = cell (Nf,nc);
                    s{i}.x.xfV.(j) = zeros(Nf,nc); %  > \phi_f*V and \nabla\phi_f*V.
                    if j == "a"
                        s{i}.x.nv.(j) = f.fh.f.f(Xc);
                        for k = 1:nc
                            s{i}.x.xfV.(j)(:,k) = f.fh.f.f(Xv).*inp.c{k}(Xv);
                        end
                    else
                        s{i}.x.nv.(j) = zeros(Nc,1);
                    end
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil cell/face indices, coordinates, etc.
        function [s] = Update_ss(inp,msh,f,s,ext)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            Xc = msh.c.Xc;
            Xv = msh.f.Xv;
            
            for i = 1:numel(s.u.s)
                if isempty(s.u.s{i})
                    continue;
                else
                    for j  = s.u.s{i}', p = s.u.p(j,i);
                        %  > Stencil cell/face indices.
                        nb = ceil(p./2);
                        ns = j-nb:j+nb-1;
                        %  > Check whether the stencil reaches any of the boundaries...
                        if any(ns <= 0)
                            % >> WB.
                            %  > Add cell(s) to the right(?).
                            if ~ext || ~ B1_1D.Change_bd(i,j,f.bd.t(1),Nf)
                                sc = 1:p;
                            else
                                sc = 1:p+1;
                            end
                            si      = [1;sc'];
                            logical = [false;true(numel(sc),1)];
                            xt      = [Xv(1),Xc(sc)'];
                        elseif any(ns >= Nf)
                            % >> EB.
                            %  > Add cell(s) to the left(?).
                            if ~ext || ~ B1_1D.Change_bd(i,j,f.bd.t(2),Nf)
                                sc = Nc-(p-1):Nc;
                            else
                                sc = Nc-p:Nc;
                            end
                            si      = [sc';Nf];
                            logical = [true(numel(sc),1);false];
                            xt      = [Xc(sc)',Xv(Nf)];
                        else
                            si      = ns';
                            logical = true(numel(ns),1);
                            xt      = Xc(ns)';
                        end
                        %  > "tfV".
                        xf       = msh.f.Xv(j);
                        Df       = B1_1D.Assemble_Df(inp,f.bd,xf,xt);
                        df       = zeros(1,numel(xt));
                        df (1,i) = 1;
                        tfV      = inp.c{i}(xf).*df*inv(Df);
                        %  > Update field "s".
                        s.i      {j,i} = si;
                        s.logical{j,i} = logical;
                        s.xt     {j,i} = xt;
                        s.tfV    {j,i} = tfV;
                    end
                end
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Auxiliary function #1: change boundary formulation.
        function [flag] = Change_bd(i,j,bd_t,Nf)
            if ((i == 1 && bd_t == "Dirichlet") || (i == 2 && bd_t == "Neumann")) && (j == 1 || j == Nf)
                flag = false;
            else
                flag = true;
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > Auxiliary function #2: check boundary type and assemble matrix Df accordingly.
        function [Df] = Assemble_Df(inp,f_bd,xf,xt)
            %  > Auxiliary variables.
            bd = ismembc(xt,f_bd.x);
            p  = numel  (xt);
            Df = zeros  (p);
            
            if any(bd) && f_bd.t(arrayfun(@(x) find(f_bd.x == x),xt(bd))) ~= "Dirichlet"
                switch f_bd.t(f_bd.x == xt(bd))
                    case "Neumann"
                        Df(~bd,:) = (xt(~bd)-xf)'.^((1:p)-1);
                        Df( bd,:) = (xt( bd)-xf)'.^((1:p)-2).*[0:p-1]; Df(isinf(Df) | isnan(Df)) = 0;
                    case "Robin"
                        Df        = (xt     -xf)'.^((1:p)-1);
                        DG        = (xt( bd)-xf)'.^((1:p)-2).*[0:p-1]; DG(isinf(DG) | isnan(DG)) = 0;
                        Df( bd,:) = Df(bd,:)+inp.c{2}(xf)./inp.c{1}(xf).*DG;
                    otherwise
                        return;
                end
            else
                Df = (xt-xf)'.^((1:p)-1);
            end
        end
    end
end