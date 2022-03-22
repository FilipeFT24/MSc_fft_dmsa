classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > 'Standard' run.
        function [obj,msh] = p_standard(inp,obj,msh,pde)
            % >> Update fields 'm', 's' and 'x'.
            obj.s        = B_2_1_1D.Update_1   (inp,msh,pde,obj.s,obj.u);
            obj.x        = B_2_1_1D.Update_2   (inp,msh,obj.s,obj.u,obj.x);
            obj.m        = B_2_1_1D.Update_3   (msh,pde,obj.m,obj.s,obj.u,obj.x);
            % >> Update...
            %  > ...solution.
            obj.x.nv.a   = pde.av;
            obj.x.nv.x.c = B_2_1_1D.Update_xc  (obj.m.At,obj.m.Bt);
            %  > ...field 'x'.
            obj.x        = B_2_1_1D.Update_4   (obj.s,obj.u,obj.x);
            % >> Update cell/face truncation and cell global discretization error distribution/norms.
            obj.e        = B_2_1_1D.Update_et_a(msh,obj.e,obj.s,obj.u,obj.x);
            obj.e        = B_2_1_1D.Update_ec_a(msh,obj.e,obj.x);
            obj.e        = B_2_1_1D.Update_et_p(inp,msh,pde,obj.e,obj.s,obj.u,obj.x);
            obj.e        = B_2_1_1D.Update_ec_p(msh,obj.e,obj.m);
            % >> Update structures.
            [obj,msh]    = B_2_1_1D.Set_struct (obj,msh);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 'Adaptative' run.
        function [obj,msh] = p_adaptative(inp,obj,msh,pde)
            i = 0;
            while 1
                % >> Update cycle count...
                i = i+1;
                % >> Update fields 'm', 's' and 'x'.
                obj.s        = B_2_1_1D.Update_1       (inp,msh,pde,obj.s,obj.u);
                obj.x        = B_2_1_1D.Update_2       (inp,msh,obj.s,obj.u,obj.x);
                obj.m        = B_2_1_1D.Update_3       (msh,pde,obj.m,obj.s,obj.u,obj.x);
                % >> Update...
                %  > ...solution(?).
                flag_1           = B_2_2_1D.Solve_AX_B (i);
                if flag_1
                    obj.x.nv.a   = pde.av;
                    obj.x.nv.x.c = B_2_1_1D.Update_xc  (obj.m.At,obj.m.Bt);
                end
                %  > ...field 'x'.
                obj.x            = B_2_1_1D.Update_4   (obj.s,obj.u,obj.x);
                % >> Update cell/face truncation and cell global discretization error distribution/norms.
                obj.e            = B_2_1_1D.Update_et_a(msh,obj.e,obj.s,obj.u,obj.x);
                obj.e            = B_2_1_1D.Update_ec_a(msh,obj.e,obj.x);
                obj.e            = B_2_1_1D.Update_et_p(inp,msh,pde,obj.e,obj.s,obj.u,obj.x);
                obj.e            = B_2_1_1D.Update_ec_p(msh,obj.e,obj.m);
                % >> Set structures.
                obj_m(i,:)       = obj.m;
                obj_e(i,:)       = obj.e.p;
                %  >> Stop adaptation(?).
                flag_2           = B_2_2_1D.Stop       (i,obj_e{i,1}.c.n_abs);
                if ~flag_2
                    obj.u        = B_2_2_1D.Update_u   (obj_e{i,1},obj.u);
                else
                    obj.m        = obj_m;
                    obj.e        = obj_e;
                    break;
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        %  > Select faces for p-coarsening/refinement.
        %  > nf    : Number of selected faces.
        %  > lambda: Parameter used to establish a treshold for the maximum face truncation error (group selection).        
        function [u_s] = Select_f(e)
            %  > Auxiliary variables.
            nf     = 2;
            lambda = 0.95;
            nv     = size(e.t.f_abs,2);
            
            [~,iM] = maxk(e.t.f_abs(:,1:nv-1),nf);
            for i  = 1:nv-1
                u_s{i}(:,1) = sort(iM(:,i));
            end
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Decrease/increase method's order (update field 'u').
        function [u] = Update_u(e,u)
            %  > Auxiliary variables.
            u_s = B_2_2_1D.Select_f(e);
            m   = size(u_s,2);
            A   = 2;
                       
            for i = 1:m
                if ~isempty(u_s{i})
                    for j = 1:size(u_s{i},1)
                        %  > Face/column index.
                        k        = u_s{i}(j);
                        l        = i*m-1;
                        %  > Increase polynomial order by "A".
                        u.p(k,l) = u.p(k,l)+A;
                    end
                end
                u.s{i} = u_s{i};
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Solve AX=B(?) criterion/criteria.
        function [flag] = Solve_AX_B(i)
            if i ~= 1
                flag = 0;
            else
                flag = 1;
            end
        end
        % >> 2.4. ---------------------------------------------------------
        %  > Set stopping criterion/criteria.
        %  > nc_M   : Maximum number of cycles.
        %  > ec_m_L1: Minimum cell global discretization error (L1 norm).
        %  > ec_M_L3: Maximum cell global discretization error (L_infinity norm).
        function [flag] = Stop(c,ec_abs)
            %  > Auxiliary variables.
            nc_M = 3;
            ec_m = 1e-10;
            
            if c > nc_M || ec_abs(1) <= ec_m
                flag = 1;
            else
                flag = 0;
            end
        end
    end
end