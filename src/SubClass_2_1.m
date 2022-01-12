classdef SubClass_2_1
    methods (Static)
        %% > Wrap up SubClass_2_1.
        function [fn,bnd,blk] = WrapUp_2_1(inp,msh)     
            %  ------------------------------------------------------------
            % >> 1.   Compute analytic expressions of f,df and d2f.
            % >> 2.   Compute analytic solution.
            % >> 3.   Compute source term.
            %  > 3.1. Cell mapping (computational domain): Triangle.
            %  > 3.2. Cell mapping (computational domain): Square.
            %  > 3.3. Isoparameteric mapping.
            %  > 3.4. Compute (individual) cell polygon integral (2D).
            %  > 3.5. Compute source term.
            %  > ----------------------------------------------------------
            % >> Local variables.
            Order = inp.fr.n;
            
            % >> 1.
            fn = SubClass_2_1.Set_fn(inp);
            % >> 2.
            [bnd,blk] = SubClass_2_1.Compute_f_df_d2f(fn,msh); 
            % >> 3.
            % -> NOT compute here.
        end
        
        %% > Tools.
        %% > 1.) ----------------------------------------------------------
        function [fn] = Set_fn(inp)
            % >> Local variables.
            [vx,vy,gx,gy] = ...
                deal(inp.pr.vx,inp.pr.vy,inp.pr.gx,inp.pr.gy);
            % >> Symbolic variables.
            syms x y;
            
            %  > Phi.
            fn.f     = sin(3.*pi.*x).*sin(3.*pi.*y);
            %  > gradPhi.
            fn.df_x  = diff(fn.f,x);
            fn.df_y  = diff(fn.f,y);
            %  > lapPhi.
            fn.d2f_x = diff(fn.df_x,x);
            fn.d2f_y = diff(fn.df_y,y);
            %  > func.
            fn.func = (vx.*fn.df_x+vy.*fn.df_y)-(gx.*fn.d2f_x+gy.*fn.d2f_y);
            
            %  > Function handles.
            [fn.f,fn.df_x,fn.df_y,fn.d2f_x,fn.d2f_y,fn.func] = ...
                deal(matlabFunction(fn.f),matlabFunction(fn.df_x),matlabFunction(fn.df_y),matlabFunction(fn.d2f_x),matlabFunction(fn.d2f_y),matlabFunction(fn.func));
        end
        
        %% > 2.) ----------------------------------------------------------
        function [bnd,blk] = Compute_f_df_d2f(fn,msh)
            %  > Boundary faces.
            for i = 1:size(msh.bnd.f,2)               
                bnd.f    (i) = fn.f    (msh.f.mean(1,msh.bnd.f{2,i}),msh.f.mean(2,msh.bnd.f{2,i}));
                bnd.df_x (i) = fn.df_x (msh.f.mean(1,msh.bnd.f{2,i}),msh.f.mean(2,msh.bnd.f{2,i}));
                bnd.df_y (i) = fn.df_y (msh.f.mean(1,msh.bnd.f{2,i}),msh.f.mean(2,msh.bnd.f{2,i}));
                bnd.d2f_x(i) = fn.d2f_x(msh.f.mean(1,msh.bnd.f{2,i}),msh.f.mean(2,msh.bnd.f{2,i}));
                bnd.d2f_y(i) = fn.d2f_y(msh.f.mean(1,msh.bnd.f{2,i}),msh.f.mean(2,msh.bnd.f{2,i}));
            end
            %  > Domain cells.
            for i = 1:msh.c.NC
                blk.f    (i) = fn.f    (msh.c.mean(1,i),msh.c.mean(2,i));
                blk.df_x (i) = fn.df_x (msh.c.mean(1,i),msh.c.mean(2,i));
                blk.df_y (i) = fn.df_y (msh.c.mean(1,i),msh.c.mean(2,i));
                blk.d2f_x(i) = fn.d2f_x(msh.c.mean(1,i),msh.c.mean(2,i));
                blk.d2f_y(i) = fn.d2f_y(msh.c.mean(1,i),msh.c.mean(2,i));
            end
        end
        
        %% > 3.) ----------------------------------------------------------
        % >> 3.1.) --------------------------------------------------------
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
        % >> 3.2.) --------------------------------------------------------
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
        % >> 3.3.) --------------------------------------------------------
        function [xx,yy,det_j] = IsoMapping(N)
            % >> Symbolic variables.
            syms csi eta;
            
            %  > Isoparameteric mapping.
            xx   = 0;
            yy   = 0;
            a    = sym('a',[1,size(N,2)]); % -> Xv(i).
            b    = sym('b',[1,size(N,2)]); % -> Yv(i).
            for i = 1:size(N,2)
                xx = xx+N{i}.*a(i);
                yy = yy+N{i}.*b(i);
            end
            j     = jacobian([xx,yy],[csi,eta]);
            det_j = det(j);
            
            %  > Create function handles.
            xx      = matlabFunction(xx);
            yy      = matlabFunction(yy);
            det_j   = matlabFunction(det_j);           
        end
        % >> 3.4.) --------------------------------------------------------
        function [Qp,I] = Compute_Integral_2D(Qc,xx,yy,det_J,func)
            %  > Quadrature abcissas/weights.
            Qp.Points  = [xx,yy];
            Qp.Weights = Qc.Weights;
            %  > Func(phi).
            func_p     = func(Qp.Points(:,1),Qp.Points(:,2));
            
            I = 0;
            for i = 1:size(Qp.Points,1)
                I = I+det_J(i).*Qp.Weights(i).*func_p(i);
            end
        end
                
        % >> #2: Compute source term (2D cell integral).
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
                %  > Quadrature/cell source term.
                [Qp{i},F_Vol(i)] = SubClass_2_1.Compute_Integral_2D(Qc{i},xx{i},yy{i},det_J{i},func);
            end
        end
    end
end