classdef A1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp] = Set_msh(h)
            inp.h               = h;            %  > Grid size.
            inp.Lim(1,:)        = [0,1];        %  > Grid limits (x-direction).
            inp.Lim(2,:)        = [0,1];        %  > Grid limits (y-direction).
            inp.p               = "s";          %  > Cell polyhedral type.
            inp.t               = 0;            %  > #Example.
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp(t,v)
            %  > ----------------------------------------------------------
            %  > Boundary conditions.
            %  > NOTE: hard coded for square domain (boundaries are identified by outher face normals (Sf)).
            inp.b.t(1)      = "Dirichlet";      %  > East (E).
            inp.b.t(2)      = "Dirichlet";      %  > North(N).
            inp.b.t(3)      = "Dirichlet";      %  > West (W).
            inp.b.t(4)      = "Dirichlet";      %  > South(S).
            if ~all(ismember(inp.b.t,["Dirichlet","Neumann","Robin"]))
                return;
            end
            %  > ----------------------------------------------------------
            %  > Coefficients/analytic function handle(s).
            fh              = A1_2D.fh_cf(t,v);
            inp.c           = fh.c;             %  > c.
            inp.f           = fh.f;             %  > f.
            %  > ----------------------------------------------------------
            %  > Method(s).
            inp.m.cls       = 1;
            inp.m.wf        = A1_2D.fh_wf(2);   %  > Weight function.
            %  > NOTE: if 0: Neighbour type: vertex.
            %          if 1: Neighbour type: face.
            inp.m.nb        = 1;                %  > Neighbouring type.
            %  > ----------------------------------------------------------
            %  > Polynomial fit.
            inp.p.p{1}(1,:) = [1,1];            %  > Convection(x): [x,y].
            inp.p.p{1}(2,:) = [1,1];            %  > Convection(y): [x,y].
            inp.p.p{2}(1,:) = [1,1];            %  > Diffusion (x): [x,y].
            inp.p.p{2}(2,:) = [1,1];            %  > Diffusion (y): [x,y].
            %  > ----------------------------------------------------------
            %  > P-Adaptation.
            inp.p.t         = 0;
            inp.p.n         = 2;                %  > Maximum number of cycles.
            inp.p.e         = 1.0E-10;          %  > Minimum global discretization/truncation error.
            inp.p.trsh      = 0.95;             %  > Treshold for face selection based on maximum face truncation error (%).
            if ~(inp.p.trsh <= 1)
                return;
            end
            %  > ----------------------------------------------------------
            %  > Plot.
            inp.plot      = [1,1];
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        function [fh] = fh_cf(t,v)
            %  > Auxiliary variables.
            c(1,:) =  [1,1];
            c(2,:) = -[1,1];
            
            %  > ch.
            for i = 1:size(c,1)
                for j = 1:size(c,2)
                    switch t(1)
                        case 1, fh.c{i,j} = @(x) repmat(c(i,j),size(x,1),1);
                        otherwise
                            return;
                    end
                end
            end
            %  > fh.
            switch t(2)
                case 1
                    fh.f = @(x) exp(-v(1).*((x(:,1)-v(2)).^2+(x(:,2)-v(3)).^2));
                case 2
                    if any(c(:,1) == 0)
                        return;
                    else
                        A = i./(2.*pi.*c(2,1));
                        B = c(1,1)./(2.*i);
                    end
                    fh.f = @(x) A.*bessely(0,A.*sqrt((x(:,1)-v(2)).^2+(x(:,2)-v(3)).^2)).*exp(B.*x(:,1));
                otherwise
                    return;
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        function [wf] = fh_wf(p)
            %  > Auxiliary variables.
            e  = 1;
            k  = 1./2;
            c  = exp(-(1./k).^2);
            
            %  > wf.
            a  = @(d) exp(-(d./(k.*(1+e).*max(d))).^2)-c;
            b  = 1-c;
            g  = @(d) a(d)./b;
            wf = @(d) g(d)./d.^p;
        end
    end
end