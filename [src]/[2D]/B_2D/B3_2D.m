classdef B3_2D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize structure "obj": {1} - present.        (".e",".f",".s",".m").
        %                                {2} - error estimator (     ".f",".s",".m").
        function [obj] = Initialize(inp,msh)
            %  > Auxiliary variables.
            nc = size(inp.c);
            ns = 2;
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            
            %  > Field: "e" (error).
            obj.e = B2_2D.Initialize_e(nc,Nc,Nf);
            %  > Field: "f" (analytic function(s), 1D/2D quadrature, etc.).
            obj.f = A3_2D.Initialize_f(inp,msh);
            %  > Field: "m" (matrices).
            obj.m = B2_2D.Initialize_m(nc,ns,Nc);
            %  > Field: "s" (stencil cell/face indices, coordinates, coefficents, etc.).
            obj.s = B1_2D.Initialize_s(inp,obj.f,nc,ns,Nc,Nf,msh.c.c.xy.c);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update all fields from structure "obj".
        function [m,s] = Update_all(inp,msh,f,m,s)
            s          = B1_2D.Update_ss(inp,msh,f,s);
            m          = B2_2D.Update_m (msh,f,m,s);
            s          = B2_2D.Update_sx(f,m,s);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set up "P-Standard" and "P-Adaptative" runs.
        function [obj] = Run_p(inp,msh)
            if ~inp.p.t
                %  > "P-Standard" run.
                obj = B3_2D.Initialize  (inp,msh);
                obj = B3_2D.P_Standard  (inp,msh,obj);
            else
                %  > "P-Adaptative" run.
                obj = B3_2D.Initialize  (inp,msh);
                obj = B3_2D.P_Adaptative(inp,msh,obj);
            end
        end
        %  > 2.1.1. -------------------------------------------------------
        function [obj] = P_Standard(inp,msh,obj)
            %  > Update fields "m" and "s".
            j                            = 1;
            [obj.m{j},obj.s{j}] = ...
                B3_2D.Update_all(inp,msh,obj.f,obj.m{j},obj.s{j});
            %  > Update field  "e".
            [obj.e] = ...
                B2_2D.Update_e  (inp,msh,obj.e,obj.m{j},obj.s{j});
            %  > Plot...
            Fig_2D_1.Plot(inp,msh,obj);
            Fig_2D_2.Plot(inp,msh,obj.e);
        end
        %  > 2.1.2. -------------------------------------------------------
        function [obj_i] = P_Adaptative(inp,msh,obj)
            %  > Initialize cycle count.
            fc = true;
            i  = 0;
            while 1
                %  > Update cycle count...
                i = i+1;
                fprintf("Loop: #%3d\n",i);
                
                %  > Update fields "m", "s" and "x".
                j                            = 1;
                [obj.m{j},obj.s{j},obj.x{j}] = ...
                    B3_2D.Update_all(inp,msh,obj.f,obj.m{j},obj.s{j},obj.u{j},obj.x{j},fc);
                %  > Update field  "e".
                [obj.e] = ...
                    B2_2D.Update_e  (inp,msh,obj.e,obj.m,obj.s,obj.x);
                %  > Assign structures (auxiliary variables).
                obj_i.e    (i,:) = obj.e;
                obj_i.m.nnz(i,:) = obj.m{j}.nnz;
                obj_i.u    (i,:) = obj.u{j};
                
                %   > Stop adaptation(?).
                if ~B2_2D.Stop(inp,i,obj.e.a.n_abs.c)
                    obj.u = B2_2D.Update_u(inp,obj.e.a,obj.u);
                else
                    break;
                end
            end
        end
    end
end