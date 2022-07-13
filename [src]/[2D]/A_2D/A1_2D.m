classdef A1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp] = Set_msh(h)
            inp.h              = h;                   %  > Grid size.
            inp.Lim(1,:)       = [0,1];               %  > Grid limits (x-direction).
            inp.Lim(2,:)       = [0,1];               %  > Grid limits (y-direction).
            inp.p              = "s";                 %  > Cell polyhedral type.
            if inp.p == "s"
                inp.t          = 0;                   %  > #Example.
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp(t,v)
            %  > ----------------------------------------------------------
            %  > Boundary conditions.
            %  > NOTE: hard coded for square domain (boundaries are identified by outher face normals (Sf)).
            inp.B.T(1)         = "Neumann";           %  > East (E).
            inp.B.T(2)         = "Dirichlet";         %  > North(N).
            inp.B.T(3)         = "Dirichlet";         %  > West (W).
            inp.B.T(4)         = "Dirichlet";         %  > South(S).
            if ~all(ismember(inp.B.T,["Dirichlet","Neumann","Robin"]))
                return;
            end
            %  > ----------------------------------------------------------
            %  > Coefficients/analytic function handle(s).
            fh                 = A1_2D.fh_cf(t{1},v);
            inp.c              = fh.c;                %  > c.
            inp.f              = fh.f;                %  > f.
            %  > ----------------------------------------------------------
            %  > Method(s).
            inp.M.Cls          = 1;                   %  > Constrain polynomial at the boundary.
            inp.M.K            = 1.05;                %  > Set stencil w/ +(%) of elements.
            inp.M.Wf           = A1_2D.fh_wf(t{2});   %  > Weight function.
            %  > ----------------------------------------------------------
            %  > Polynomial fit.
            inp.P              = [1,1];               %  > [x,y].
            %  > ----------------------------------------------------------
            %  > #Test.
            inp.T              = 2;
            %  > ----------------------------------------------------------
            %  > P-Adaptation.
            %  > Selection criteria.
            inp.P_Ad.Config    = 3;                   %  > Neighbouring configuration type.
            inp.P_Ad.Isotropic = 1;                   %  > Isotropic coarsening/refinement.
            if inp.T == 2 && (inp.P_Ad.Isotropic && range(inp.P) ~= 0)
                return;
            end
            inp.P_Ad.Qrt       = [0.000,0.950];       %  > Treshold for face selection (%of faces): coarsening/refinement.
            inp.P_Ad.Qrt_PT    = 0.50;                %  > Treshold for face selection (%of faces): premature termination.
            if any(inp.P_Ad.Qrt < 0 | inp.P_Ad.Qrt > 1) || any(inp.P_Ad.Qrt_PT < 0 | inp.P_Ad.Qrt_PT > 1)
                return;
            end
            %  > Stopping criteria.
            inp.P_Ad.ec        = 1.00E-07;            %  > Minimum global discretization error.
            inp.P_Ad.Nc        = 1;                  %  > Maximum number of cycles.
            inp.P_Ad.NNZ       = 50E3;                %  > Maximum nnz(A).
            %  > ----------------------------------------------------------
            %  > Plot.
            inp.Plot{1}        = [0,0];
            inp.Plot{2}        = [1,1];
            inp.Plot{3}        = 1;
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        function [fh] = fh_cf(t,v)
            %  > Auxiliary variables.
            c(1,:) =  [5,0];
            c(2,:) = -[1,1];
            
            % >> ch.
            for i = 1:size(c,1)
                for j = 1:size(c,2)
                    switch t(1)
                        case 1, fh.c{i,j} = @(x) repmat(c(i,j),size(x,1),1);
                        otherwise
                            return;
                    end
                end
            end
            % >> fh.
            if ~t(2)
                %  > ...case   isotropic.
                for j = 1:size(v{2},2)
                    rhs{j} = @(x) (x(:,1)-v{2}(1,j)).^2+(x(:,2)-v{2}(2,j)).^2;
                end
                switch t(3)
                    case 1, fh.f = @(x) exp (-v{1}(1).*sum(cellfun(@(f) f(x),rhs)));
                    case 2, fh.f = @(x) sqrt( v{1}(1).*sum(cellfun(@(f) f(x),rhs)));
                    otherwise
                        return;
                end
            else
                %  > ...case anisotropic.
                switch t(3)
                    case 1
                        for j = 1:size(v{2},2)
                            rhs.a{j} = @(x) (x(:,1)-v{2}(1,j)).^2;
                            rhs.b{j} = @(x) (x(:,2)-v{2}(2,j));
                        end
                        a    = -v{1}(1);
                        b    =  v{1}(2);
                        fh.f =  @(x) exp(a.*sum(cellfun(@(f) f(x),rhs.a))).*sin(b.*sum(cellfun(@(f) f(x),rhs.b)));
                    case 2
                        for j = 1:size(v{2},2)
                            rhs.a{j} = @(x) (x(:,1)-v{2}(1,j));
                            rhs.b{j} = @(x) (x(:,2)-v{2}(2,j)).^2;
                        end
                        a    =  v{1}(2);
                        b    = -v{1}(1);
                        fh.f =  @(x) sin(a.*sum(cellfun(@(f) f(x),rhs.a))).*exp(b.*sum(cellfun(@(f) f(x),rhs.b)));
                    otherwise
                        return;
                end
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        function [wf] = fh_wf(t)
            switch t
                case 1, wf = @(d) (min(d)./d).^2;
                otherwise
                    return;
            end
        end
    end
end