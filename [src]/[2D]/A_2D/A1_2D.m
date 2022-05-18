classdef A1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp_m] = Set_msh(h)
            inp_m.h             = h;                                 %  > Grid size.
            inp_m.Lim(1,:)      = [0,1];                             %  > Grid limits (x-direction).
            inp_m.Lim(2,:)      = [0,1];                             %  > Grid limits (y-direction).
            inp_m.p             = "s";                               %  > Cell polyhedral (type).
            inp_m.t             = 0;                                 %  > #Example.
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp(t,v)
            %  > ----------------------------------------------------------
            %  > Analytic function/coefficients.
            h                   = A1_2D.cf_h(t,v);
            inp.c               = h.c;                               %  > c.
            inp.f               = h.f;                               %  > f.
            %  > Boundary conditions.
            %  > NOTE: hard coded for square domain (boundaries are identified by outher face normals (Sf)).
            inp.b.t(1)          = "Dirichlet";                       %  > East  (E).
            inp.b.t(2)          = "Dirichlet";                       %  > North (N).
            inp.b.t(3)          = "Dirichlet";                       %  > West  (W).
            inp.b.t(4)          = "Dirichlet";                       %  > South (S).
            if ~all(ismember(inp.b.t,["Dirichlet","Neumann","Robin"]))
                return;
            end
            %  > Polynomial fit.
            inp.p.p(1,:)        = [5,5];                             %  > p-convection(X/Y).
            inp.p.p(2,:)        = [5,5];                             %  > p-diffusion (X/Y).
            inp.p.nb_t          = 0;
            inp.p.wls           = 1;
            if inp.p.wls
                p               = 2;                                 %  > p.
                e               = 1;                                 %  > \epsilon.
                k               = 1./2;                              %  > k.
                c(1)            = k.*(1+e);
                c(2)            = exp(-(1./k).^2);
                inp.p.wf        = @(d) ((exp(-(d(1,:)./(c(1).*max(d))).^2)-c(2))./(1-c(2)))./d(1,:).^p;
            end
            %  > ----------------------------------------------------------
            %  > P-adaptation.
            %  > #1.
            inp.p_adapt.allow   = 0;                                 %  > Allow p-adaptation(?).
            inp.p_adapt.nc      = 50;                                %  > Maximum number of cycles.
            inp.p_adapt.ec_m    = 1.0E-10;                           %  > Minimum (global) discretization error.
            inp.p_adapt.lambda  = 0.75;                              %  > Treshold for face selection based on maximum face truncation error (%).
            if ~(inp.p_adapt.lambda < 1)
                return;
            end
            %  > #2.
            inp.p_adapt.opt(1)  = 0;                                 %  > Use higher-order solution(?).
            inp.p_adapt.opt(2)  = 0;                                 %  > Add lower-order (predicted) cell truncation error as source term(?).
            if all(inp.p_adapt.opt)                                  %  > Only allow w/ lower-order solution...
                return;
            end
            %  > ----------------------------------------------------------
            %  > Plot...
            inp.plot            = [0,1];
            %  > ----------------------------------------------------------
        end
        % >> 1.3. ---------------------------------------------------------
        function [h] = cf_h(t,v)
            %  > Auxiliary variables.
            c (1,:) =  [0,0];
            c (2,:) = -[1,1];
            i       =  v(1);
            xc      =  v(2);
            yc      =  v(3);
            
            %  > ch.
            switch t(1)
                case 1
                    h.c{1,1} = @(x) repelem(c(1,1),size(x,1),1); %  > V(x).
                    h.c{1,2} = @(x) repelem(c(1,2),size(x,1),1); %  > V(y).
                    h.c{2,1} = @(x) repelem(c(2,1),size(x,1),1); %  > G(x).
                    h.c{2,2} = @(x) repelem(c(2,2),size(x,1),1); %  > G(x).
                otherwise
                    return;
            end
            %  > fh.
            switch t(2)
                case 1
                    h.f = @(x) exp(-i.*((x(:,1)-xc).^2+(x(:,2)-yc).^2));
                case 2
                    if any(c(:,1) == 0)
                        return;
                    else
                        A = i./(2.*pi.*c(2,1));
                        B = c(1,1)./(2.*i);
                    end
                    h.f = @(x) A.*bessely(0,A.*sqrt((x(:,1)-xc).^2+(x(:,2)-yc).^2)).*exp(B.*x(:,1));
                otherwise
                    return;
            end
        end
    end
end