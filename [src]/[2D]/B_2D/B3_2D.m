classdef B3_2D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize "obj" structure: {1} - present.        (".s",".m",".x").
        %                                {2} - error estimator (".s",".m",".x").
        function [obj] = Initialize(inp,msh)
            %  > Auxiliary variables.
            A  = [0,2];
            nc = size(inp.c);
            ns = 2;
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            
            %  > Field: "e" (error).
            obj.e = B2_2D.Initialize_e (nc,Nc,Nf);
            %  > Field: "f" (analytic function, 1D/2D quadrature, etc.).
            obj.f = A3_2D.Initialize_f (inp,msh);
            %  > Field: "s" (stencil cell/face indices, coordinates, etc.).
            obj.s = B1_2D.Initialize_1 (nc,ns,Nf);
            %  > Field: "m" (matrices).
            obj.m = B2_2D.Initialize_3 (nc,ns,Nc);
            %  > Field: "u" (update flag).
            for i = 1:ns
                for j = 1:nc(1)
                    for k = 1:nc(2)
                        obj.u{i}.p{j}   (:,k) = repelem(inp.p.p(j,k)+A(i),Nf); %  > (x,y).
                        obj.u{i}.s{j}{k}(:,1) = 1:Nf;
                    end
                end
            end
            %  > Field: "x" (stencil coefficients/nodal solution,etc.).
            obj.x = B1_2D.Initialize_24(obj.f,obj.u,nc,ns,Nc,Nf);                  
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update all fields ("obj").
        function [m,s,x] = Update_all(inp,msh,f,m,s,u,x,f_xc)
            s            = B1_2D.Update_1(inp,msh,s,u);
            x            = B1_2D.Update_2(inp,msh,f,s,u,x);
            m            = B2_2D.Update_3(msh,f,m,s,u,x);
            x            = B2_2D.Update_4(f,m,s,u,x,f_xc);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set up 'P-standard' and 'P-adaptative' runs.
        function [obj] = Run_p(inp,msh)
            if ~inp.p_adapt.allow
                %  > 'P-standard' run.
                obj = B3_2D.Initialize  (inp,msh);
                obj = B3_2D.P_standard  (inp,msh,obj);
                %  > Plot...
                Fig_2D_1.Plot(inp,msh,obj);
                Fig_2D_2.Plot(inp,msh,obj.e);
            else
                %  > 'P-adaptative' run.
                obj = B3_2D.Initialize  (inp,msh);
                obj = B3_2D.P_adaptative(inp,msh,obj);
                %  > Plot...
                Fig_2D_2.Plot(inp,msh,obj.e(end,:));
            end
        end
        % >> 2.2. ---------------------------------------------------------
        % >> 'p-standard' run.
        function [obj] = P_standard(inp,msh,obj)
            %  > Update fields "m", "s" and "x".
            fc                           = true;
            j                            = 1;
            [obj.m{j},obj.s{j},obj.x{j}] = ...
                B3_2D.Update_all(inp,msh,obj.f,obj.m{j},obj.s{j},obj.u{j},obj.x{j},fc);
            %  > Update field  "e".
            [obj.e] = ...
                B2_2D.Update_e  (inp,msh,obj.e,obj.m,obj.x);
        end
        % >> 2.3. ---------------------------------------------------------
        % >> 'P-adaptative' run.
        function [obj_i] = P_adaptative(inp,msh,obj)
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
                    B2_2D.Update_e  (inp,msh,obj.e,obj.m,obj.x);
                %  > Assign structures (auxiliary variables).
                obj_i.e    (i,:) = obj.e;
                obj_i.m.nnz(i,:) = obj.m{j}.nnz;
                
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