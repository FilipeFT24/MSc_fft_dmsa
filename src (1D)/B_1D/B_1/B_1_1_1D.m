classdef B_1_1_1D
    methods (Static)
        %% > Wrap-up B_1_1 (1D).
        function [pde] = WrapUp_B_1_1_1D(msh,Xv_i,Xv_f,v,g,as)
            % >> Compute...
            %  > ...analytic functions.
            pde.fn = B_1_1_1D.Set_fn(Xv_i,Xv_f,v,g,as);
            %  > ...analytic values.
            pde.sn = B_1_1_1D.Compute_blkf_blkc(msh,pde.fn);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fn] = Set_fn(Xv_i,Xv_f,v,g,as)
            %  > Symbolic variable.
            syms x;
            
            switch char(as)
                case 'sin'
                    i    = 9;
                    f{1} = sin(i.*pi.*x);
                case 'exp'
                    c    = 1./2.*(Xv_f-Xv_i);
                    i    = 50;
                    f{1} = exp(-i.*((x-c).^2));
                otherwise
                    return;
            end
            f{2} = diff(f{1},x);
            func = v.*f{2}-g.*diff(f{2},x);
            intf = int(func);
            
            % >> Set...
            %  > ...symbolic functions.
            [fn.sym{1},fn.sym{2}] = deal(f{1},f{2});
            %  > ...function handles.
            [fn.f{1},fn.f{2},fn.func,fn.int] = deal(matlabFunction(f{1}),matlabFunction(f{2}),matlabFunction(func),matlabFunction(intf));
        end
        % >> 1.2. ---------------------------------------------------------
        function [s] = Compute_blkf_blkc(msh,fn)
            %  > Cell values.
            i        = 1:msh.c.NC;
            s.c(i,1) = fn.f{1}(msh.c.Xc(i));
            
            %  > Face values.
            j        = 1:msh.f.NF;
            s.f(j,1) = fn.f{1}(msh.f.Xv(j));
            s.f(j,2) = fn.f{2}(msh.f.Xv(j));
        end
    end
end