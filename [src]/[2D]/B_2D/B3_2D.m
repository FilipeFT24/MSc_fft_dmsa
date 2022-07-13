classdef B3_2D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize structure "obj".
        function [obj] = Initialize_1(inp,msh)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            
            %  > Field: "e" (error).
            obj.e = B2_2D.Initialize_e(Nc,Nf);
            %  > Field: "f" (analytic function(s), 1D/2D quadrature, etc.).
            obj.f = A3_2D.Initialize_f(inp,msh);
            %  > Field: "m" (matrices).
            obj.m = B2_2D.Initialize_m(Nc);
            %  > Field: "s" (stencil cell/face indices, coordinates, coefficents, etc.).
            obj.s = B1_2D.Initialize_s(inp,obj.f,Nf,msh.c.c.xy.c);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update all fields from structure "obj".
        function [e,f,m,s] = Update_all(inp,msh,e,f,m,s)
            s              = B1_2D.Update_ss(inp,msh,f,s);
            m              = B2_2D.Update_m (msh,f,m,s);
            s              = B2_2D.Update_sx(inp,msh,f,m,s);
            e              = B2_2D.Update_e (inp,msh,e,m,s);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set up "P-Standard" and "P-Adaptative" runs.
        function [obj] = Run_p(inp,msh)
            switch inp.T
                %  > "P-Standard" run.
                case 1, obj = B3_2D.Initialize_1(inp,msh);
                    obj = B3_2D.P_Standard  (inp,msh,obj);
                    %  > "P-Adaptative" run.
                case 2, obj = B3_2D.Initialize_1(inp,msh);
                    obj = B3_2D.P_Adaptative(inp,msh,obj);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [obj] = P_Standard(inp,msh,obj)
            %  > Update fields "e", "m" and "s".
            [obj.e,obj.f,obj.m,obj.s] = ...
                B3_2D.Update_all(inp,msh,obj.e,obj.f,obj.m,obj.s);
            %  > Plot...
            Plot_2D_1.Plot(inp,msh,obj);
        end
        % >> 2.3. ---------------------------------------------------------
        function [obj] = P_Adaptative(inp,msh,obj)
            % >> Initialize...
            i = 1;
            for j = ["c","r"]
                block.(j) = false(msh.f.Nf,1);
            end
            %  > Fields "e", "f", "m" and "s".
            [obj.e,obj.f,obj.m,obj.s] = ...
                B3_2D.Update_all(inp,msh,obj.e,obj.f,obj.m,obj.s);
            
            % >> Until any stopping criterion has not been met...
            while 1
                %  > Re-do loop (if necessary)...
                if i > 1
                    block = B2_2D.rst(inp,msh,block,[obj(:).e],v(i-1));
                    if any(block.c | block.r)
                        i = i-1;
                    end
                end
                %  > Check stopping criteria.
                if ~(any(block.c | block.r)) && B2_2D.Stop(inp,obj)
                    break;
                else
                    %  > Select faces for coarsening/refinement...
                    [obj(i+1).s,v(i),stop] = ...
                        B2_2D.Update_u(inp,msh,obj(i).e.a.t.f_abs,obj(i).s,block);
                    if ~stop
                        %  > Print to terminal...
                        if ~(any(block.c | block.r))
                            fprintf("Loop #%2d\n",i);
                        end
                        %  > Update cycle count...
                        i = i+1;
                        %  > Update fields "e", "f", "m" and "s".
                        [obj(i).e,obj(i).f,obj(i).m,obj(i).s] = ...
                            B3_2D.Update_all(inp,msh,obj(i-1).e,obj(i-1).f,obj(i-1).m,obj(i).s);
                    else
                        if ~(any(block.c | block.r))
                            obj = obj(1:i);
                        else
                            obj = obj(1:i-1);
                        end
                        break;
                    end
                end
            end
            %  > Plot...
            Plot_2D_1.Plot(inp,msh,obj);
            Plot_2D_2.Plot(inp,msh,obj,v);
        end
    end
end