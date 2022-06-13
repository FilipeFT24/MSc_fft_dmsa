classdef B3_2D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize structure "obj": {1} - present.        (".e",".f",".s",".m").
        %                                {2} - error estimator (     ".f",".s",".m").
        function [obj] = Initialize_1(inp,msh)
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
            switch inp.p.t
                %  > "P-Standard" run.
                case 1, obj = B3_2D.Initialize_1(inp,msh);
                        obj = B3_2D.P_Standard  (inp,msh,obj);
                    %  > "P-Adaptative" run.
                case 2, obj = B3_2D.Initialize_1(inp,msh);
                        obj = B3_2D.P_Adaptative(inp,msh,obj);
                otherwise
                    %  > Other tests.
                    h.lim = [3.5E-2,2.5E-2];
                    h.n   = 5;
                    load  = [1,1];
                    V1    = B3_2D.Initialize_2(h,load);
                    switch inp.p.t
                        %  > #1.
                        case 3, V2 = B3_2D.Load_p_1(V1.inp,V1.msh,load);
                        otherwise
                            return;
                    end
                    obj = V2.obj;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [obj] = P_Standard(inp,msh,obj)
            %  > Update fields "m" and "s".
            j                   = 1;
            [obj.m{j},obj.s{j}] = ...
                B3_2D.Update_all(inp,msh,obj.f,obj.m{j},obj.s{j});
            %  > Update field  "e".
            [obj.e] = ...
                B2_2D.Update_e  (inp,msh,obj.e,obj.m{j},obj.s{j});
            %  > Plot...
            Plot_2D_1.Plot(inp,msh,obj);
        end
        % >> 2.3. ---------------------------------------------------------
        function [obj_i] = P_Adaptative(inp,msh,obj)
            %  > Initialize cycle count.
            i = 0;
            while 1
                %  > Update cycle count...
                i = i+1;
                fprintf("Loop: #%3d\n",i);
                
                %  > Update fields "m" and "s".
                j                   = 1;
                [obj.m{j},obj.s{j}] = ...
                    B3_2D.Update_all(inp,msh,obj.f,obj.m{j},obj.s{j});
                %  > Update field  "e".
                [obj.e] = ...
                    B2_2D.Update_e  (inp,msh,obj.e,obj.m{j},obj.s{j});
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
        % >> 2.4. ---------------------------------------------------------
        %  > Set up other tests.
        %  > switch ...
        %  >    case #1, (...).
        %  > end
        %  > 2.4.1. -------------------------------------------------------
        %  > Set up structure(s) "inp" and "msh".
        function [V] = Initialize_2(h,load_V)
            %  > Working directory.
            wd             = "B_2D/[.mat Files]/[.msh]/";
            %  > Set structure(s) "inp".
            inp            = A1_2D.Set_inp([1,1],[100,0.5,0.5]); 
            inp.plot{1}(:) = 0;
            inp.p.t        = 1;
            %  > H.
            H              = exp(1).^(linspace(log(h.lim(1)),log(h.lim(2)),h.n));
            %  > Set structure(s) "msh".
            if ~load_V(1)
                for i = 1:numel(H)
                    msh(i) = A2_2D.Set_msh(H(i));
                    fprintf("Cycle #%3d (msh)\n",i);
                end
                save(strjoin([wd,"msh.mat"],''),"msh");
            else
                load(strjoin([wd,"msh.mat"],''));
            end
            %  > Assign to structure "V".
            V.inp = inp;
            V.msh = msh;           
        end
        %  > 2.4.2. -------------------------------------------------------
        function [V] = Load_p_1(inp,msh,load_V)
            %  > Working directory.
            wd = "B_2D/[.mat Files]/[.V]/[1]/";
            %  > Execute and assign to structure "V".
            if ~load_V(2)
                for i = 1:size(msh,2)
                    obj  (i)   = B3_2D.Run_p(inp,msh(i));
                    V.msh(i).d = msh(i).d;
                    V.obj(i).e = obj(i).e;
                    for j = 1:numel(obj(i).m)
                        V.obj(i).m{j}.nnz = obj(i).m{j}.nnz;
                    end
                    fprintf("Cycle #%3d (obj)\n",i);
                end
                save(strjoin([wd,"V.mat"],''),"V");
            else
                load(strjoin([wd,"V.mat"],''));
            end
            %  > Plot...
            Plot_2D_2.Plot(inp,V.msh,V.obj);
        end
    end
end