classdef B_1_1_1D
    methods (Static)
        %% > Wrap-up B_1_1 (1D).
        function [pde] = WrapUp_B_1_1_1D(msh,Xv_i,Xv_f,v,g,as)
            %  > Analytic functions.
            [pde.an.fn] = ...
                B_1_1_1D.Set_fn(Xv_i,Xv_f,v,g,as);
            %  > Analytic values.
            [pde.an.f_v,pde.an.df_v,pde.an.c_v] = ...
                B_1_1_1D.Compute_blkf_blkc(msh,pde.an.fn.f,pde.an.fn.df);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fn] = Set_fn(Xv_i,Xv_f,v,g,as)
            % >> Symbolic variable.
            syms x;
            
            %  > Phi.
            if as == 'sin'
                %  > Intensity.
                i    = 9;
                fn.f = sin(i.*pi.*x);
            elseif as == 'exp'
                %  > Location.
                xc   = 1./2.*(Xv_f-Xv_i);
                %  > Intensity.
                i    = 50;
                fn.f = exp(-i.*((x-xc).^2));
            end
            %  > gradPhi/lapPhi.
            fn.df    = diff(fn.f,x);
            d2f      = diff(diff(fn.f,x),x);
            %  > func/int(func).
            fn.func  = v.*fn.df-g.*d2f;
            fn.int   = int(fn.func);
            
            %  > Function handles.
            [fn.f,fn.df,fn.func,fn.int] = ...
                deal(matlabFunction(fn.f),matlabFunction(fn.df),matlabFunction(fn.func),matlabFunction(fn.int));
        end
        % >> 1.2. ---------------------------------------------------------
        function [f_v,df_v,c_v] = Compute_blkf_blkc(msh,f,df)
            %  > Face values.
            i         = 1:msh.f.NF;
            f_v (i,1) = f (msh.f.Xv(i));
            df_v(i,1) = df(msh.f.Xv(i));
            %  > Cell values.
            j         = 1:msh.c.NC;
            c_v (j,1) = f(msh.c.Xc(j));
        end
    end
end