classdef A1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp] = Set_inp_1(h)
            %  > ----------------------------------------------------------
            %  > msh.
            %    ├─ v
            %        ├─ Example 1: Uniform.
            %                    └─ fft_Distmesh2D (https://github.com/ionhandshaker/distmesh/blob/master/distmesh.m).
            %        └─ Example 2: Non-uniform.
            %                    └─ Random.
            %    └─ s
            %        ├─ Example 1: Uniform.
            %        └─ Example 2: Non-uniform.
            %                    ├─ Random.
            %                    └─ Smooth non-uniform.
            %                        ├─  A: Stretching parameter.
            %                        └─  c: Clustered location.
            %  > ----------------------------------------------------------
            inp.m.Uniform       =  1;                                 %  > Set uniform grid(?).
            if ~inp.m.Uniform                                         %  > If non-uniform...
                inp.m.is_random =  0;
                if ~inp.m.is_random                                   %  > If not random...
                    inp.m.A     =  [5,5];                             %  > Stretching parameter(X/Y).
                    inp.m.c(1)  =  0.5;                               %  > Clustered location(X).
                    inp.m.c(2)  =  0.5;                               %  > Clustered location(Y).
                end
            end
            inp.m.p{1}          = "s";                                %  > Cell polyhedral.
            inp.m.p{2}          =  1;                                 %  > Example# (check A2_2D.m).
            inp.m.h             =  h;                                 %  > Grid size.
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp_2(v)
            %  > ----------------------------------------------------------
            %  > Analytic function.
            syms x y;
            xc                  =  v(1);
            yc                  =  v(2);
            i                   =  v(3);  
            f                   =  exp(-i.*((x-xc).^2+(y-yc).^2));    %  > f.
            inp.f               =  matlabFunction(f,'Vars',{x,y});    %  > f handle.
            %  > ----------------------------------------------------------
            %  > Boundary treatment.
            inp.b.change(1)     =  1;                                 %  > Add n extra points for boundary face fit.
            inp.b.change(2)     =  0;                                 %  > Add n extra points for boundary face fit (error estimator).
            %  > ----------------------------------------------------------
            %  > Coefficients.
            inp.c(1,:)          =  [0,0];                             %  > Convection(X/Y).
            inp.c(2,:)          = -[1,1];                             %  > Diffusion (X/Y).
            %  > ----------------------------------------------------------
            %  > Polynomial fit.
            inp.p.p(1)          =  1;                                 %  > p(X).
            inp.p.p(2)          =  1;                                 %  > p(Y).
            if any(rem(inp.p.p,2) == 0)                               %  > Allow only p=1,3,5,7,9,etc.
                return;
            end
            inp.p.nb            = "V";                                %  > F/V.
            inp.p.WLS           =  1;                                 %  > D/LSQ.
            if inp.p.WLS
                syms d1 d2;
                p               =  2;                                 %  > p.
                e               =  1;                                 %  > \epsilon.
                k               =  1./2;                              %  > k.
                c(1)            =  k.*(1+e);
                c(2)            =  exp(-(1./k).^2);
                g               =  ((exp(-(d1./(c(1).*d2)).^2)-c(2))./(1-c(2)));
                wf              =  g./d1.^p;                       
                inp.wf          =  matlabFunction(wf,'Vars',{d1,d2}); %  > w (weight function) handle.
            end
            %  > ----------------------------------------------------------
            %  > (Volumetric) source term.
            inp.p.st            =  0;                                 %  > w/ analytic integrand(0)/w/ 2D quadrature(1).
            %  > ----------------------------------------------------------
            %  > P-adaptation.
            %  > #1.
            inp.p_adapt.allow   =  0;                                 %  > Allow p-adaptation(?).
            inp.p_adapt.nc      =  50;                                %  > Maximum number of cycles.
            inp.p_adapt.ec_m    =  1.0E-10;                           %  > Minimum (global) discretization error.
            inp.p_adapt.lambda  =  0.75;                              %  > Treshold for face selection based on maximum face truncation error (%).
            if ~(inp.p_adapt.lambda < 1)
                return;
            end
            %  > #2. 
            inp.p_adapt.opt(1)  =  0;                                 %  > Use higher-order solution(?).
            inp.p_adapt.opt(2)  =  0;                                 %  > Add lower-order (predicted) cell truncation error as source term(?).   
            if all(inp.p_adapt.opt)                                   %  > Only allow w/ lower-order solution...
                return;
            end
        end         
    end
end