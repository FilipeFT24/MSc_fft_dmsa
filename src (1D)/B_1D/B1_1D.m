classdef B1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update stencil coordinates,etc.
        function [s] = Update_1(msh,f,s,u,ch)
            [s]      = B1_1D.Update_sc(msh,f,s,u,ch);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update stencil coefficients,etc.
        function [x] = Update_2(inp,msh,f,s,u,x)
            x        = B1_1D.Update_sx(inp,msh,f,s,u,x);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Update stencil coordinates.
        function [s] = Update_sc(msh,f,s,u,ch)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            Xc = msh.c.Xc;
            Xv = msh.f.Xv;
            
            for i = 1:size(u.s,2)
                if isempty(u.s{i})
                    continue;
                else
                    for j  = u.s{i}'
                        p  = u.p(j,i);
                        
                        %  > Stencil cell/face indices.
                        nb = ceil(p./2);
                        ns = j-nb:j+nb-1;
                        %  > Check whether stencil reaches any of the boundaries and shift it accordingly.
                        if any(ns <= 0)
                            % >> WB.
                            %  > Add cell(s) to the right(?).
                            if ~ch || ~ B1_1D.Change_bd(i,j,f.bd.t(1),Nf)
                                sc = 1:p;
                            else
                                sc = 1:p+1;
                            end
                            xt = [Xc(sc)',f.bd.x(1)];
                        elseif any(ns >= Nf)
                            % >> EB.
                            %  > Add cell(s) to the left(?).
                            if ~ch || ~ B1_1D.Change_bd(i,j,f.bd.t(2),Nf)
                                sc = Nc-(p-1):Nc;
                            else
                                sc = Nc-p:Nc;
                            end
                            xt = [Xc(sc)',f.bd.x(2)];
                        else
                            sc = ns;
                            xt = Xc(sc)';
                        end
                        %  > Update 's' field.
                        s.c{j,i} = sc;
                        s.t{j,i} = sort(xt);
                    end
                end
            end
        end
        %  > 2.1.1. -------------------------------------------------------
        %  > Auxiliary function: change boundary formulation (flag).
        function [flag] = Change_bd(i,j,bd_t,Nf)
            switch i
                case 1
                    if (j == 1 || j == Nf) && bd_t == "Dirichlet"
                        flag = 0;
                    else
                        flag = 1;
                    end
                case 2
                    if (j == 1 || j == Nf) && bd_t == "Neumann"
                        flag = 0;
                    else
                        flag = 1;
                    end
                otherwise
                    return;
            end
        end       
        % >> 2.2. ---------------------------------------------------------
        %  > Update stencil coefficients.
        function [x] = Update_sx(inp,msh,f,s,u,x)
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j = u.s{i}'
                        %  > Df.
                        xt        = s.t{j,i};
                        p         = 1:length(xt);
                        gv        = inp.c(2)./inp.c(1);
                        Df        = B1_1D.Assemble_Df(f.bd,p,xt,msh.f.Xv(j),gv);
                        %  > Tf.
                        df        = zeros(1,length(xt));
                        df  (1,i) = 1;
                        Inv       = inv(Df);
                        Tf        = df*Inv;
                        %  > Update 'x' field.
                        x.if{j,i} = Inv;
                        x.Tf{j,i} = Tf;
                    end
                end
            end
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Check boundary type and assemble matrix Df accordingly.
        function [Df] = Assemble_Df(f_bd,p,xt,x,gv)
            %  > Auxiliary variables.
            a  = length(p);
            b  = length(xt);
            Df = zeros (b,a);
            
            %  > Assemble Df...
            bd = ismembc(xt,f_bd.x);
            if any(bd) && f_bd.t(f_bd.x == xt(bd)) ~= "Dirichlet"
                switch f_bd.t(f_bd.x == xt(bd))
                    case "Neumann"
                        i          = p;
                        j          = p(1:end-1);
                        k          = xt == f_bd.x(f_bd.x == xt(bd));
                        Df(~k,1:a) = (xt(~k)-x)'.^(i-1);
                        Df( k,1:a) = [0,j.*(xt(k)-x)'.^(j-1)];
                    case "Robin"
                        i          = p;
                        j          = p(1:end-1);
                        k          = xt == f_bd.x(f_bd.x == xt(bd));
                        Df( i,1:a) = (xt(i)-x)'.^(i-1);
                        Df( k,1:a) = Df( k,i)-[0,j.*(xt(k)-x)'.^(j-1)].*gv;
                    otherwise
                        return;
                end
            else
                i         = p;
                j         = 1:b;
                Df(j,1:a) = (xt(j)-x)'.^(i-1);
            end
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Update truncated terms' magnitude (w/ analytic field).
        function [em] = Update_t_terms(inp,msh,em,f,n,s,u,x)
            for i = 1:size(u.s,2)
                if ~isempty(u.s{i})
                    for j  = u.s{i}'
                        %  > Df(extended).
                        xt = s.t{j,i};
                        p  = length(xt)+1:length(xt)+n(i);
                        gv = inp.c(2)./inp.c(1);
                        Df = B1_1D.Assemble_Df(f.bd,p,xt,msh.f.Xv(j),gv);
                        %  > Derivatives up to...
                        for k = 1:length(p)
                            Df(:,k) = Df(:,k).*f.fh.f.d2{p(k)-1}(msh.f.Xv(j));
                        end
                        %  > Truncated terms' magnitude (approximated/finite truncation error).
                        em.f    {i}(j,p) = x.Tf{j,i}*Df;
                        em.f_abs{i}      = abs(em.f{i});
                    end
                end
            end
        end
    end
end