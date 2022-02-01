classdef B_1_1_1D
    methods (Static)
        %% > Wrap-up B_1_1 (1D).
        function [pde] = WrapUp_B_1_1_1D(msh,Xv_i,Xv_f,v,g,as)
            %  > Analytic functions.
            [pde.an.fn] = ...
                B_1_1_1D.Set_fn(Xv_i,Xv_f,v,g,as);
            %  > Analytic values.
            [pde.an.bnd,pde.an.blk] = ...
                B_1_1_1D.Compute_blkf_blkc(msh,pde.an.fn.f);
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
            df       = diff(fn.f,x);
            d2f      = diff(diff(fn.f,x),x);
            %  > func/int(func).
            fn.func  = v.*df-g.*d2f;
            fn.int   = int(fn.func);
            
            %  > Function handles.
            [fn.f,fn.func,fn.int] = ...
                deal(matlabFunction(fn.f),matlabFunction(fn.func),matlabFunction(fn.int));
        end
        % >> 1.2. ---------------------------------------------------------
        function [blk_f,blk_c] = Compute_blkf_blkc(msh,f)
            %  > Face values.
            i          = 1:msh.f.NF;
            blk_f(i,1) = f(msh.f.Xv(i));
            %  > Cell values.
            j          = 1:msh.c.NC;
            blk_c(j,1) = f(msh.c.Xc(j));
        end
    end
end