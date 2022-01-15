classdef SubClass_2_1_1
    methods (Static)
        %% > Wrap-up SubClass_2_1_1.
        % >> --------------------------------------------------------------
        % >> 1.   Compute analytic solution.
        %  > 1.1. Compute analytic expressions of f,df and d2f.
        %  > 1.2. Compute values.
        % >>
        % >> --------------------------------------------------------------
        function [bnd,blk] = WrapUp_2_1_1(inp,msh)
            % >> 1.
            %  > 1.1.
            fn = SubClass_2_1_1.Set_fn(inp);
            % >> 1.2.
            [bnd,blk] = SubClass_2_1_1.Compute_f_df_d2f(fn,msh);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
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
        % >> 1.2. ---------------------------------------------------------
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
    end
end