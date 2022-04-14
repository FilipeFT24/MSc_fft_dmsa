classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 'p-standard' run.
        function [obj] = p_standard(inp,msh,obj)
            %  > Auxiliary variables.
            add_b      = inp.pv.add;
            fs         = 1;
            obj.x.nv.a = obj.f.av;
            
            %  > Update fields 'e', 'f', 's', 'u' and 'x'.
            [obj.m,obj.s,obj.x] = B_2_1_1D.Update_all(inp,msh,obj.f,obj.m,obj.s,obj.u,obj.x,add_b,fs);
            [obj.e,obj.x]       = B_2_1_1D.Update_e  (inp,msh,obj.e,obj.f,obj.m,obj.s,obj.u,obj.x,add_b);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 'p-adaptative' run.
        % >> 2.1. ---------------------------------------------------------
        function [obj] = p_adaptative(inp,obj,msh,fld)
            %  > Auxiliary variables.
            add_b      = inp.ps.add;
            fld_u      = ["a","x"];
            obj.x.nv.a = obj.f.av;
            i          = 0;
            while 1
                %  > Update cycle count...
                i = i+1;
                %  > Update fields 'm', 's' and 'x'.
                [obj.s,obj.x,obj.m] = ...
                    B_2_1_1D.Update_123(inp,msh,obj.f,obj.m,obj.s,obj.u,obj.x,add_b);
                %  > Update field 'x'.
                flag_1           = B_2_2_1D.Solve_AX_B(i);
                if flag_1
                    obj.x.nv.x.c = B_2_1_1D.Update_xc (obj.m.At,obj.m.Bt);
                end
                obj.x            = B_2_1_1D.Update_4  (obj.f,obj.s,obj.u,obj.x,fld_u);
                %  > Update field 'e'.               
                obj.e            = B_2_1_1D.Update_e  (inp,msh,obj.e,obj.f,obj.m,obj.s,obj.u,obj.x,add_b);
                %  > Assign structures.
                obj_m(i,:)       = obj.m;
                obj_e(i,:)       = obj.e;
                %   > Stop adaptation(?).
                flag_2           = B_2_2_1D.Stop      (i,obj.e.(fld){1}.c.n_abs);
                if ~flag_2
                    obj.u        = B_2_2_1D.Update_u  (obj.e.(fld){1},obj.u);
                else
                    obj.m        = obj_m;
                    obj.e        = obj_e;
                    break;
                end
                fprintf("Loop: #%3d\n",i);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        %  > Select faces for p-coarsening/refinement.
        %  > nf    : Number of selected faces.
        %  > lambda: Parameter used to establish a treshold for the maximum face truncation error (group selection).        
        function [u_s] = Select_f(e)
            %  > Auxiliary variables.
            nf     = 1;
            lambda = 0.95;
            nv     = size(e.t.f_abs,2);
            u_s    = cell(1,nv-1);
            sf     = find(e.t.f_abs(:,nv) > lambda.*max(e.t.f_abs(:,nv)));
            
            %if ismembc(1,sf)
            %    sf = [1,2,3,4,5,6]';
            %    disp(3);
            %end
            
            for i  = 1:nv-1
                u_s{i} = sf;
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
            nc_M = 30;
            ec_m = 1.0E-10;
            
            if c > nc_M || ec_abs(1) <= ec_m
                flag = 1;
                if c > nc_M
                    fprintf("Stopping criterion: max. number of cycles.\n");
                elseif ec_abs(1) <= ec_m
                    fprintf("Stopping criterion: min. error treshold (L1 norm).\n");
                end
            else
                flag = 0;
            end
        end
    end
end