classdef A1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp_m] = Set_inp_1(h)
            inp_m.h             = h;                                 %  > Grid size.
            inp_m.Lim(1,:)      = [0,1];                             %  > Grid limits (x-direction).
            inp_m.Lim(2,:)      = [0,1];                             %  > Grid limits (y-direction).
            inp_m.p             = "s";                               %  > Cell polyhedral (type).
            inp_m.t             = 0;                                 %  > #Example.
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp_2(t,v)
            %  > ----------------------------------------------------------
            %  > Analytic function/coefficients.
            inp.c               = Tools_1.c   (t,v);                 %  > c.
            inp.f               = Tools_1.func(t,v);                 %  > f.
            %  > ----------------------------------------------------------
            %  > Boundary conditions.
            %  > NOTE: hard coded for square domain (boundaries are identified by outher face normals (Sf)).
            inp.b.t(1)          = "Dirichlet";                       %  > East  (E).
            inp.b.t(2)          = "Dirichlet";                       %  > North (N).
            inp.b.t(3)          = "Dirichlet";                       %  > West  (W).
            inp.b.t(4)          = "Dirichlet";                       %  > South (S).
            if ~all(ismember(inp.b.t,["Dirichlet","Neumann","Robin"]))
                return;
            end
            %  > ----------------------------------------------------------
            %  > Polynomial fit.
            inp.p.p(1,:)        = [3,1];                             %  > p-convection(X/Y).
            inp.p.p(2,:)        = [3,1];                             %  > p-diffusion (X/Y).
            if any(any(rem(inp.p.p,2) == 0)) || ...                  %  > Allow only p=1,3,5,7,9,etc.
                    ~all(inp.p.p(1,:) == inp.p.p(2,:))               %  > Treat in a unified manner...
                return;
            end
            inp.p.nb_t          = 1;
            inp.p.wls           = 1;
            if inp.p.wls
                syms d1 d2;
                p               = 2;                                 %  > p.
                e               = 1;                                 %  > \epsilon.
                k               = 1./2;                              %  > k.
                c(1)            = k.*(1+e);
                c(2)            = exp(-(1./k).^2);
                g               = ((exp(-(d1./(c(1).*d2)).^2)-c(2))./(1-c(2)));
                wf              = g./d1.^p;
                inp.p.wf        = matlabFunction(wf,'Vars',{d1,d2}); %  > wf (weight function) handle.
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
            inp.plot            = [1,1];
        end
    end
end