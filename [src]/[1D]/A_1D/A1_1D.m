classdef A1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp] = Set_msh(h,t)
            inp.h   = h;                            %  > Grid size.
            inp.Lim = [0,1];                        %  > Grid limits(x).
            inp.t   = t;                            %  > #Example.
            if inp.t ~= 1
                switch inp.t
                    case 2, inp.x(1) = 5.00;        %  > Stretching (α).
                            inp.x(2) = 0.50;        %  > Clustered location (x0).
                    case 3, inp.x(1) = 1.25;        %  > Distortion (q).
                    otherwise
                        return;
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp(t,v)
            %  > ----------------------------------------------------------
            %  > Boundary conditions.
            inp.b.t(1)         = "Dirichlet";       %  > West(W).
            inp.b.t(2)         = "Dirichlet";       %  > East(E).
            inp.b.extension(1) = 1;                 %  > Extension layer for boundary face fit.
            inp.b.extension(2) = 0;                 %  > Extension layer for boundary face fit(error estimator).
            if ~all(ismember(inp.b.t,["Dirichlet","Neumann","Robin"]))
                return;
            end
            %  > ----------------------------------------------------------
            %  > Coefficients/analytic function handle(s).
            fh                 = A1_1D.fh_cf(t,v);
            inp.c              = fh.c;              %  > c.
            inp.f              = fh.f;              %  > f.
            %  > ----------------------------------------------------------
            %  > Polynomial fit.
            inp.p.p(1)         = 1;                 %  > Convection(x).
            inp.p.p(2)         = 1;                 %  > Diffusion (x).
            %  > ----------------------------------------------------------
            %  > P-Adaptation.
            inp.p.t            = 0;
            inp.p.n            = 2;                 %  > Maximum number of cycles.
            inp.p.e            = 1.0E-10;           %  > Minimum global discretization/truncation error.
            inp.p.trsh         = 0.95;              %  > Treshold for face selection based on maximum face truncation error (%).
            if ~(inp.p.trsh <= 1)
                return;
            end
            %  > ----------------------------------------------------------
            %  > Plot.
            inp.plot           = [0,1,0];
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        function [fh] = fh_cf(t,v)
            %  > Auxiliary variables.
            c = [1,-1];
            
            %  > ch.
            for i = 1:numel(c)
                switch t(1)
                    case 1, fh.c{i} = @(x) repmat(c(i),size(x));
                    otherwise
                        return;
                end
            end
            %  > fh.
            switch t(2)
                case 1, fh.f = @(x) sin(pi.*(v(1).*x+v(2)));
                case 2, fh.f = @(x) exp(-v(1).*(x-v(2)).^2);
                otherwise
                    return;
            end
        end
    end
end