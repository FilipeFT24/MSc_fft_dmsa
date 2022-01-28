classdef B_1_1
    methods (Static)
        %% > Wrap-up B_1_1.
        % >> --------------------------------------------------------------
        % >> 1.   Compute analytic solution.
        %  > 1.1. Compute analytic expressions of f,df and d2f.
        %  > 1.2. Compute values.
        % >> --------------------------------------------------------------
        function [pde] = WrapUp_B_1_1(msh,Xv_i,Xv_f,Yv_i,Yv_f,vx,vy,gx,gy)
            % >> 1.
            %  > 1.1.
            [pde.fn] = B_1_1.Set_fn(Xv_i,Xv_f,Yv_i,Yv_f,vx,vy,gx,gy,'1');
            %  > 1.2.
            [pde.bnd,pde.blk] = B_1_1.Compute_f_df_d2f(msh,pde.fn);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fn] = Set_fn(Xv_i,Xv_f,Yv_i,Yv_f,vx,vy,gx,gy,as)
            % >> Symbolic variables.
            syms x y;
            
            %  > Phi.
            if as == '1'
                fn.f = sin(3.*pi.*x).*sin(3.*pi.*y);
            elseif as == '2'
                %  > Location.
                xc   = 1./2.*(Xv_f-Xv_i);
                yc   = 1./2.*(Yv_f-Yv_i);
                %  > Intensity.
                i    = 50;
                fn.f = exp(-i.*((x-xc).^2+(y-yc).^2));
            end
            
            %  > gradPhi.
            df_x  = diff(fn.f,x);
            df_y  = diff(fn.f,y);
            %  > lapPhi.
            d2f_x = diff(df_x,x);
            d2f_y = diff(df_y,y);
            %  > func.
            fn.func = (vx.*df_x+vy.*df_y)-(gx.*d2f_x+gy.*d2f_y);
            
            %  > Function handles.
            [fn.f,fn.func] = deal(matlabFunction(fn.f),matlabFunction(fn.func));
        end
        % >> 1.2. ---------------------------------------------------------
        function [bnd,blk] = Compute_f_df_d2f(msh,fn)
            %  > Boundary faces.
            i          = 1:size(msh.bnd.f,2);
            bnd.f(i,1) = fn.f(msh.f.mean(1,[msh.bnd.f{2,i}]),msh.f.mean(2,[msh.bnd.f{2,i}]));
            %  > Domain cells.
            j          = 1:msh.c.NC;
            blk.f(j,1) = fn.f(msh.c.mean(1,j),msh.c.mean(2,j));           
        end
    end
end