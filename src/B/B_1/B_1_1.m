classdef B_1_1
    methods (Static)
        %% > Wrap-up B_1_1.
        % >> --------------------------------------------------------------
        % >> 1.   Compute analytic solution.
        %  > 1.1. Compute analytic expressions of f,df and d2f.
        %  > 1.2. Compute values.
        % >> --------------------------------------------------------------
        function [pde] = WrapUp_B_1_1(msh,vx,vy,gx,gy)
            % >> 1.
            %  > 1.1.
            pde.fn = B_1_1.Set_fn(vx,vy,gx,gy);
            %  > 1.2.
            [pde.bnd,pde.blk] = B_1_1.Compute_f_df_d2f(msh,pde.fn);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fn] = Set_fn(vx,vy,gx,gy)
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
        % >> 1.2. ---------------------------------------------------------
        function [bnd,blk] = Compute_f_df_d2f(msh,fn)
            %  > Boundary faces.
            i            = 1:size(msh.bnd.f,2);
            bnd.f    (i) = fn.f    (msh.f.mean(1,[msh.bnd.f{2,i}]),msh.f.mean(2,[msh.bnd.f{2,i}]));
            bnd.df_x (i) = fn.df_x (msh.f.mean(1,[msh.bnd.f{2,i}]),msh.f.mean(2,[msh.bnd.f{2,i}]));
            bnd.df_y (i) = fn.df_y (msh.f.mean(1,[msh.bnd.f{2,i}]),msh.f.mean(2,[msh.bnd.f{2,i}]));
            bnd.d2f_x(i) = fn.d2f_x(msh.f.mean(1,[msh.bnd.f{2,i}]),msh.f.mean(2,[msh.bnd.f{2,i}]));
            bnd.d2f_y(i) = fn.d2f_y(msh.f.mean(1,[msh.bnd.f{2,i}]),msh.f.mean(2,[msh.bnd.f{2,i}]));
            
            %  > Domain cells.
            j            = 1:msh.c.NC;
            blk.f    (j) = fn.f    (msh.c.mean(1,j),msh.c.mean(2,j));
            blk.df_x (j) = fn.df_x (msh.c.mean(1,j),msh.c.mean(2,j));
            blk.df_y (j) = fn.df_y (msh.c.mean(1,j),msh.c.mean(2,j));
            blk.d2f_x(j) = fn.d2f_x(msh.c.mean(1,j),msh.c.mean(2,j));
            blk.d2f_y(j) = fn.d2f_y(msh.c.mean(1,j),msh.c.mean(2,j));           
        end
    end
end