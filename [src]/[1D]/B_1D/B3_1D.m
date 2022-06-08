classdef B3_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize structure "obj": {1} - present.        (".e",".f",".s",".m").
        %                                {2} - error estimator (     ".f",".s",".m").
        function [obj] = Initialize(inp,msh)
            %  > Auxiliary variables.
            nc = numel(inp.c);
            ns = 2;
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;

            %  > Field: "e" (error).
            obj.e = B2_1D.Initialize_e(nc,Nc,Nf);
            %  > Field: "f" (analytic function(s), boundary values, etc.).
            obj.f = A3_1D.Initialize_f(inp,msh);
            %  > Field: "m" (matrices).
            obj.m = B2_1D.Initialize_m(nc,ns,Nc);
            %  > Field: "s" (stencil cell/face indices, coordinates, coefficents, etc.).
            obj.s = B1_1D.Initialize_s(inp,obj.f,nc,ns,Nc,Nf,msh.c.Xc,msh.f.Xv);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update all fields from structure "obj".
        function [m,s] = Update_all(inp,msh,f,m,s)
            s          = B1_1D.Update_ss(inp,msh,f,s,0);
            m          = B2_1D.Update_m (msh,f,m,s);
            s          = B2_1D.Update_sx(f,m,s);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set up "P-Standard" and "P-Adaptative" runs.
        function [obj] = Run_p(inp,msh)
            if ~inp.p.t
                %  > "P-Standard" run.
                obj = B3_1D.Initialize  (inp,msh);
                obj = B3_1D.P_Standard  (inp,msh,obj);
            else
                %  > "P-Adaptative" run.
                obj = B3_1D.Initialize  (inp,msh);
                obj = B3_1D.P_Adaptative(inp,msh,obj);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [obj] = P_Standard(inp,msh,obj)
            %  > Update fields "m" and "s".
            j                   = 1;
            [obj.m{j},obj.s{j}] = ...
                B3_1D.Update_all(inp,msh,obj.f,obj.m{j},obj.s{j});
            %  > Update field  "e".
            [obj.e] = ...
                B2_1D.Update_e  (inp,msh,obj.e,obj.m{j},obj.s{j});
            %  > Plot...
            Fig_1_1D.Plot(inp.plot(1),msh);
            Fig_2_1D.Plot(inp.plot(2),msh,obj);
        end
        % >> 2.3. ---------------------------------------------------------
        function [obj] = P_Adaptative(inp,msh,obj)
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > switch ...
        %  >    case 1, grid types.
        %  > end
        function [] = Load_p(t)
            %  > h.
            m     = 3;
            n     = 5;
            h_lim = [1.0E-2,1.0E-3];
            h     = exp(1).^(linspace(log(h_lim(1)),log(h_lim(2)),n));
            %  > Select...
            switch t
                case 1
                    %  > "inp".
                    inp      = A1_1D.Set_inp([1,2],[10,0.5]);       %  > c/f.
                    inp.plot = [0,0,1];
                    %  > "msh" and "obj".
                    for i = 1:n
                        for j = 1:m
                            msh(i,j) = A2_1D.Set_msh(h(i),j);       %  > h.
                            obj(i,j) = B3_1D.Run_p  (inp,msh(i,j)); %  > inp/msh.
                        end
                    end
                otherwise
                    return;
            end
            %  > Plot...
            Fig_3_1D.Plot(inp.plot(3),msh,obj);
        end
    end
end