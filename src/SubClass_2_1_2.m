classdef SubClass_2_1_2
    methods (Static)
        %% > SubClass_2_1.
        % >> --------------------------------------------------------------
        % >> 1.   1D quadrature.
        %  > 1.1. Face mapping.
        %  > 1.2  Gauss points/weights per face.
        % >> 2.   2D quadrature.
        %  > 2.1. Cell mapping (computational domain): Triangle.
        %  > 2.2. Cell mapping (computational domain): Square.
        %  > 2.3. Isoparameteric mapping.
        %  > 2.4. Compute (individual) cell polygon integral (2D).
        %  > 2.5. Compute source term.
        % >> --------------------------------------------------------------
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1 ----------------------------------------------------------
        function [x,j] = CD_1D()
            % >> Symbolic variables.
            syms a b csi;
            
            %  > x(csi) = a*(1-csi)/2+b*(1+csi)/2.
            x = a.*(1-csi)./2+b.*(1+csi)./2;
            x = matlabFunction(x);
            %  > j(csi) = d(x)/d(csi) = (b-a)/2.
            j = (b-a)./2;
            j = matlabFunction(j);
        end
        % >> 1.2. ---------------------------------------------------------
        function [fg] = GaussFace_Points(xy_fg,j_fg,ng,xy_f)
            %  > 1D Quadrature abcissas/weights.
            Q_1D = quadGaussLegendre(ng);
            
            for i = 1:length(Q_1D.Points)
                fg.Points(i,1) = xy_fg(xy_f(1,1),xy_f(2,1),Q_1D.Points(i));
                fg.Points(i,2) = xy_fg(xy_f(1,2),xy_f(2,2),Q_1D.Points(i));
            end
            fg.Weights = 1./2.*Q_1D.Weights;
            fg.jac     = [j_fg(xy_f(1,1),xy_f(1,2)),j_fg(xy_f(2,1),xy_f(2,2))];
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [Qc,N] = CD_Triangle(n)
            % >> Symbolic variables.
            syms csi eta;
            
            % >> Computational domain.
            %  > Quadrature abcissas/weights.
            Qc = quadtriangle(n,...
                'Type','product','Symmetry','allowAsymmetric','Weights','allowNegative','Points','allowOutside');
            %  > Vertex shape functions.
            N{1} = 1-csi-eta;
            N{2} = csi;
            N{3} = eta;
        end
        % >> 2.2. ---------------------------------------------------------
        function [Qc,N] = CD_Quadrilateral(n)
            % >> Symbolic variables.
            syms csi eta;
            
            %  > Quadrature abcissas/weights.
            Qc = quadsquare(n,...
                'Type','productLegendre','Symmetry','allowAsymmetric','Weights','allowNegative','Points','allowOutside');
            %  > Vertex shape functions.
            N{1} = 1./4.*(1-csi).*(1-eta);
            N{2} = 1./4.*(1+csi).*(1-eta);
            N{3} = 1./4.*(1+csi).*(1+eta);
            N{4} = 1./4.*(1-csi).*(1+eta);
        end
        % >> 2.3. ---------------------------------------------------------
        function [xx,yy,det_j] = IsoMapping(N)
            % >> Symbolic variables.
            syms csi eta;
            
            %  > Isoparameteric mapping.
            xx   = 0;
            yy   = 0;
            a    = sym('a',[1,size(N,2)]); % -> Xv(i).
            b    = sym('b',[1,size(N,2)]); % -> Yv(i).
            for i = 1:size(N,2)
                [xx,yy] = ...
                    deal(xx+N{i}.*a(i),yy+N{i}.*b(i));
            end
            j     = jacobian([xx,yy],[csi,eta]);
            det_j = det(j);
            
            %  > Create function handles.
            xx    = matlabFunction(xx);
            yy    = matlabFunction(yy);
            det_j = matlabFunction(det_j);
        end
        % >> 2.4. ---------------------------------------------------------
        function [Qp,I] = Compute_Integral_2D(Qc,xx,yy,det_J,func)
            %  > 2D Quadrature abcissas/weights.
            Qp.Points  = [xx,yy];
            Qp.Weights = Qc.Weights;
            %  > Func(phi).
            func_p     = func(Qp.Points(:,1),Qp.Points(:,2));
            %  > 2D integral.
            I = 0;
            for i = 1:size(Qp.Points,1)
                I = I+det_J(i).*Qp.Weights(i).*func_p(i);
            end
        end
        % >> 2.5. ---------------------------------------------------------
        function [Qp,F_Vol] = Compute_SourceTerm(n,msh,func)
            %  > Cell mapping (computational domain).
            [Qc_T,N_T]          = SubClass_2_1.CD_Triangle     (n);
            [xx_T,yy_T,det_j_T] = SubClass_2_1.IsoMapping      (N_T);
            [Qc_S,N_S]          = SubClass_2_1.CD_Quadrilateral(n);
            [xx_S,yy_S,det_j_S] = SubClass_2_1.IsoMapping      (N_S);
            
            %  > Compute source term (based on cell polygon).
            for i = 1:msh.c.NC
                if size(msh.c.XY_v{i},1) == 3
                    % >> (xx,yy,det_J,Qc).
                    %  > (X1,X2,X3):
                    for j = 1:size(N_T,2)
                        X{i}{j} = msh.c.XY_v{i}(j,1);
                        Y{i}{j} = msh.c.XY_v{i}(j,2);
                    end
                    Qc   {i} = Qc_T;
                    xx   {i} = xx_T   (X{i}{1},X{i}{2},X{i}{3},Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    yy   {i} = yy_T   (Y{i}{1},Y{i}{2},Y{i}{3},Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    det_J{i} = det_j_T(X{i}{1},X{i}{2},X{i}{3},Y{i}{1},Y{i}{2},Y{i}{3});
                elseif size(msh.c.XY_v{i},1) == 4
                    % >> (xx,yy,det_J,Qc).
                    %  > (X1,X2,X3):
                    for j = 1:size(N_S,2)
                        X{i}{j} = msh.c.XY_v{i}(j,1);
                        Y{i}{j} = msh.c.XY_v{i}(j,2);
                    end
                    Qc   {i} = Qc_S;
                    xx   {i} = xx_S   (X{i}{1},X{i}{2},X{i}{3},X{i}{4},Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    yy   {i} = yy_S   (Y{i}{1},Y{i}{2},Y{i}{3},Y{i}{4},Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                    det_J{i} = det_j_S(X{i}{1},X{i}{2},X{i}{3},X{i}{4},Y{i}{1},Y{i}{2},Y{i}{3},Y{i}{4},Qc{i}.Points(:,1),Qc{i}.Points(:,2));
                end
                %  > 2D Quadrature/cell source term.
                [Qp{i},F_Vol(i)] = SubClass_2_1.Compute_Integral_2D(Qc{i},xx{i},yy{i},det_J{i},func);
            end
        end
    end
end