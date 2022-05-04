classdef B3_2D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize 'obj' structure.
        function [obj] = Initialize(inp,msh)
            %  > Auxiliary variables.
            A  = [0,2];
            nc = length(inp.c);
            ns = 2;
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            
            %  >
            
            
           
            %  > Field: 'f' (analytic function).
            obj.f.fh = A3_2D.Update_fh(inp);
            obj.f.av = A3_2D.Update_av(msh,obj.f.fh.f);
            obj.f.bd = A3_2D.Update_bd(inp,msh,obj.f.fh.f);
            obj.f.st = A3_2D.Update_st(msh,obj.f.fh.func.f,inp.p.st);
            %  > Field: 'q' (1D/2D quadrature).
            obj.f.qd = A3_2D.Q_1D(2);
            
            %  > ----------------------------------------------------------
            %  > Initialize fields...
            %  > {1} - present.        ('.s','.m','.x').
            %  > {2} - error estimator ('.s','.m','.x').
            %  > ----------------------------------------------------------
            %  > Field: "e" (error).
            obj.e = B2_2D.Initialize_e (nc,ns,Nc,Nf);
            %  > Field: "s" (stencil cell/face indices, coordinates, etc.).
            obj.s = B1_2D.Initialize_1 (nc,ns,Nf);
            %  > Field: "m" (matrices).
            obj.m = B2_2D.Initialize_3 (nc,ns,Nc);
            %  > Field: "u" (update flag).
            for i = 1:ns
                obj.u{i}.f = ["a","x"];
                for j = 1:numel(inp.p.p)
                    obj.u{i}.p   (:,j) = repelem(inp.p.p(j)+A(i),Nf); %  > X/Y.
                    obj.u{i}.s{j}(:,1) = 1:Nf;
                end
            end
            %  > Field: "x" (stencil coefficients/nodal solution,etc.).
            obj.x = B1_2D.Initialize_24(obj.f,obj.u,nc,ns,Nc,Nf); 
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update all fields ('obj').
        function [m,s,x] = Update_all(inp,msh,f,m,s,u,x,add,f_uc)
            %  > Update stencil.
            s            = B1_2D.Update_1   (msh,s,u);
            x            = B1_2D.Update_2   (inp,msh,f,s,u,x);
            %  > Update matrices(?).
            if f_uc(1)
                m        = B2_2D.Update_3   (inp,msh,f,m,s,u,x);
            end
            %  > Update nodal solution(?).
            if f_uc(2)
                x.nv.x.c = Tools_2.Update_xc(m.At,m.Bt);
            end
            %  > Update face values(?).
            if f_uc(3)
                x        = B1_2D.Update_4   (f,s,u,x);
            end
            f.av.c-x.nv.x.c
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set up 'p-standard' and 'p-adaptative' runs.
        function [obj] = Run_p(inp,msh)
            switch inp.p_adapt.allow
                case false
                    %  > 'p-standard' run.
                    obj = B3_2D.Initialize(inp,msh);
                    obj = B3_2D.p_standard(inp,msh,obj);
                    %  > Plot...
                    %  Fig_V1_0_2D.Plot(1,msh,obj,"blk");
                    Fig_V1_1_2D.Plot(1,msh,obj);
                case true
                    %  > 'p-adaptative' run.
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        % >> 'p-standard' run.
        function [obj] = p_standard(inp,msh,obj)
            %  > Auxiliary variables.
            ch   = inp.b.change(1);
            f_uc = [1,1,1];
            f_xc = 1;
            
            %  > Update fields 'm', 's' and 'x'.
            ic                              = 1;
            [obj.m{ic},obj.s{ic},obj.x{ic}] = ...
                B3_2D.Update_all(inp,msh,obj.f,obj.m{ic},obj.s{ic},obj.u{ic},obj.x{ic},ch,f_uc);
            %  > Update fields 'e', 'm' and 'x'.
            [obj.e] = ...
                B2_2D.Update_e  (inp,msh,obj.e,obj.f,obj.m,obj.s,obj.u,obj.x);
        end
    end
end