classdef B_1_2
    methods (Static)
        %% > Wrap-up B_1_2.
        % >> --------------------------------------------------------------
        % >> 1.   1D quadrature.
        %  > 1.1. Face mapping.
        %  > 1.2  Gauss points/weights per face (NOT used here).
        % >> 2.   2D quadrature.
        %  > 2.1. Cell mapping (computational domain): Triangle.
        %  > 2.2. Cell mapping (computational domain): Square.
        %  > 2.3. Isoparameteric mapping.
        %  > 2.4. Compute (individual) cell polygon integral (2D).
        %  > 2.5. Compute source term.
        % >> --------------------------------------------------------------
        function [pde] = WrapUp_B_1_2(msh,pde,ng)
            % >> 1.
            %  > 1.1.
            [pde.f.xy_fg,pde.f.j_fg,pde.f.Q_1D] = ...
                B_1_2.CD_1D(ng);
            % >> 2.
            [pde.c.Qp,pde.c.F_Vol] = ...
                B_1_2.Compute_SourceTerm(msh,pde.fn.func,ng);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1 ----------------------------------------------------------
        function [x,j,Q_1D] = CD_1D(ng)
            %  > Symbolic variables.
            syms a b csi;
            
            %  > x(csi) = a*(1-csi)/2+b*(1+csi)/2.
            x = a.*(1-csi)./2+b.*(1+csi)./2;
            x = matlabFunction(x);
            %  > j(csi) = d(x)/d(csi) = (b-a)/2.
            j = (b-a)./2;
            j = matlabFunction(j);
            %  > Quadrature abcissas/weights.
            Q_1D = quadGaussLegendre(ng);
        end
        % >> 1.2. ---------------------------------------------------------
        function [fg] = GaussFace_Points(Q_1D,xy_fg,j_fg,xy_f)
            %  > Points.
            i              = 1:length(Q_1D.Points);
            fg.Points(i,1) = xy_fg(xy_f(1,1),xy_f(2,1),Q_1D.Points(i));
            fg.Points(i,2) = xy_fg(xy_f(1,2),xy_f(2,2),Q_1D.Points(i));
            %  > Weights.
            fg.Weights     = 1./2.*Q_1D.Weights;
            fg.jac         = [j_fg(xy_f(1,1),xy_f(1,2)),j_fg(xy_f(2,1),xy_f(2,2))];
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [Qc,N] = CD_Triangle(n)
            %  > Symbolic variables.
            syms csi eta;
            
            %  > Quadrature abcissas/weights.
            Qc = quadtriangle(n,'Domain',[0,0;1,0;0,1]);
            %  > Vertex shape functions.
            N{1} = 1-csi-eta;
            N{2} = csi;
            N{3} = eta;
        end
        % >> 2.2. ---------------------------------------------------------
        function [Qc,N] = CD_Quadrilateral(n)
            %  > Symbolic variables.
            syms csi eta;
            
            %  > Quadrature abcissas/weights.
            Qc = quadsquare(n);
            %  > Vertex shape functions.
            N{1} = 1./4.*(1-csi).*(1-eta);
            N{2} = 1./4.*(1+csi).*(1-eta);
            N{3} = 1./4.*(1+csi).*(1+eta);
            N{4} = 1./4.*(1-csi).*(1+eta);
        end
        % >> 2.3. ---------------------------------------------------------
        function [x,y,det_J] = IsoMapping(N)
            %  > Symbolic variables.
            syms csi eta;
            
            %  > Isoparameteric mapping.
            x = 0;
            y = 0;
            a = sym('a',[1,size(N,2)]); %  > Xv(i).
            b = sym('b',[1,size(N,2)]); %  > Yv(i).
            for i = 1:size(N,2)
                x = x+N{i}.*a(i);
                y = y+N{i}.*b(i);
            end
            J     = jacobian([x,y],[csi,eta]);
            det_J = det(J);
            
            %  > Create function handles.
            x     = matlabFunction(x);
            y     = matlabFunction(y);
            det_J = matlabFunction(det_J);
        end
        % >> 2.4. ---------------------------------------------------------
        function [Qp,I] = Compute_Integral_2D(Qc,xx,yy,det_J,func)
            %  > Quadrature abcissas/weights.
            Qp.Points  = [xx,yy];
            Qp.Weights = Qc.Weights;
            %  > Integral.
            I = 0;
            for i = 1:size(Qp.Points,1)
                I = I+det_J(i).*Qp.Weights(i).*func(Qp.Points(i,1),Qp.Points(i,2));
            end
        end
        % >> 2.5. ---------------------------------------------------------
        function [Qp,F_Vol] = Compute_SourceTerm(msh,func,ng)
            %  > Cell mapping (computational domain).
            for i = 1:msh.c.NC
                l(i) = length(msh.c.f.f{i});
            end
            %  > Evaluate...
            if any(l == 3)
                [Qc_T,N_T]          = B_1_2.CD_Triangle     (ng);
                [xx_T,yy_T,det_J_T] = B_1_2.IsoMapping      (N_T);
            elseif any(l == 4)
                [Qc_S,N_S]          = B_1_2.CD_Quadrilateral(ng);
                [xx_S,yy_S,det_J_S] = B_1_2.IsoMapping      (N_S);
            end

            %  > Compute source term (based on cell polygon).
            for i = 1:msh.c.NC
                if size(msh.c.xy_v{i},1) == 3
                    j           = 1:size(N_T,2);
                    X    {i}(j) = msh.c.xy_v{i}(j,1);
                    Y    {i}(j) = msh.c.xy_v{i}(j,2);
                    Qc   {i}    = Qc_T;
                    xx   {i}    = xx_T   (X{i}(1),X{i}(2),X{i}(3),Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    yy   {i}    = yy_T   (Y{i}(1),Y{i}(2),Y{i}(3),Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    det_J{i}    = det_J_T(X{i}(1),X{i}(2),X{i}(3),Y{i}(1),Y{i}(2),Y{i}(3));
                    det_J{i}    = repelem(det_J{i},size(xx{i},1))';
                elseif size(msh.c.xy_v{i},1) == 4
                    j           = 1:size(N_S,2);
                    X    {i}(j) = msh.c.xy_v{i}(j,1);
                    Y    {i}(j) = msh.c.xy_v{i}(j,2);
                    Qc   {i}    = Qc_S;
                    xx   {i}    = xx_S   (X{i}(1),X{i}(2),X{i}(3),X{i}(4),Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    yy   {i}    = yy_S   (Y{i}(1),Y{i}(2),Y{i}(3),Y{i}(4),Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    det_J{i}    = det_J_S(X{i}(1),X{i}(2),X{i}(3),X{i}(4),Y{i}(1),Y{i}(2),Y{i}(3),Y{i}(4),Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                end
                [Qp{i},F_Vol(i)] = B_1_2.Compute_Integral_2D(Qc{i},xx{i},yy{i},det_J{i},func);
            end
        end
    end
end