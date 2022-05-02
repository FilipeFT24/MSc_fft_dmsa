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
            
            %  > Fields: 'e.a', 'e.d' and 'e.p'.
            for i = ["a","da","p","d"]
                switch i
                    case "a"
                        o = ns;
                    otherwise
                        o = ns-1;
                end
                for j = 1:o
                    obj.e.(i){j}.t.c       = zeros(Nc,1);
                    obj.e.(i){j}.t.c_abs   = zeros(Nc,1);
                    obj.e.(i){j}.t.f       = zeros(Nf,nc+1);
                    obj.e.(i){j}.t.f_abs   = zeros(Nf,nc+1);
                    obj.e.(i){j}.t.n.c     = zeros(1 ,3);
                    obj.e.(i){j}.t.n_abs.c = zeros(1 ,3);
                    obj.e.(i){j}.t.n.f     = zeros(3 ,nc+1);
                    obj.e.(i){j}.t.n_abs.f = zeros(3 ,nc+1);
                    obj.e.(i){j}.c.c       = zeros(Nc,1);
                    obj.e.(i){j}.c.c_abs   = zeros(Nc,1);
                    obj.e.(i){j}.c.n       = zeros(1 ,3);
                    obj.e.(i){j}.c.n_abs   = zeros(1 ,3);
                end
            end
            %  > Field: 'f' (analytic function).
            obj.f.fh = A3_2D.Update_fh(inp);
            obj.f.av = A3_2D.Update_av(msh,obj.f.fh.f);
            obj.f.bd = A3_2D.Update_bd(inp,msh,obj.f.fh.f);
            obj.f.st = A3_2D.Update_st(msh,obj.f.fh.func.f,inp.p.st);
            %  > Field: 'q' (1D/2D quadrature).
            obj.f.qd = A3_2D.Q_1D(2);
            
            %  > Other fields...
            %  > {1} - present.        ('.s','.m','.x').
            %  > {2} - error estimator ('.s','.m','.x').
            for i = 1:ns
                %  > Field: 'm' (matrices).
                for j = 1:nc
                    obj.m{i}.Ac{j} = zeros(Nc);
                    obj.m{i}.Bc{j} = zeros(Nc,1);
                end
                obj.m{i}.At = zeros(Nc);
                obj.m{i}.Bt = obj.f.st;
                %  > Field: 's' (stencil cell/face indices,coordinates,etc.).
                obj.s{i}.i       = cell(Nf,nc);
                obj.s{i}.logical = cell(Nf,nc);
                obj.s{i}.xt      = cell(Nf,nc);
                %  > Field: 'u' (update flag).
                obj.u{i}.f  = ["a","x"];
                for j = 1:length(inp.p.p)
                    obj.u{i}.p   (:,j) = repelem(inp.p.p(j)+A(i),Nf); %  > X/Y.
                    obj.u{i}.s{j}(:,1) = 1:Nf;
                end
                %  > Field: 'x' (nodal solution/stencil coefficients,etc.).
                obj.x{i}.Df     = cell (Nf,nc);
                obj.x{i}.Tf     = cell (Nf,nc);
                obj.x{i}.Tf_V   = cell (Nf,nc);
                obj.x{i}.nv.a.f = zeros(Nf,1);
                for j = obj.u{i}.f
                    obj.x{i}.cf.(j) = cell (Nf,nc);
                    obj.x{i}.vf.(j) = cell (Nf,nc);
                    obj.x{i}.xf.(j) = zeros(Nf,nc);
                    if j == "a"
                        obj.x{i}.nv.(j)   = obj.f.av;
                    else
                        obj.x{i}.nv.(j).c = zeros(Nc,1);
                    end
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update all fields ('obj').
        function [m,s,x] = Update_all(inp,msh,f,m,s,u,x,add,f_uc)
            %  > Update stencil.
            s            = B1_2D.Update_1 (msh,s,u);
            x            = B1_2D.Update_2 (inp,msh,f,s,u,x);
            %  > Update matrices(?).
            if f_uc(1)
                m        = B2_2D.Update_3 (inp,msh,f,m,s,u,x);
            end
            %  > Update nodal solution(?).
            if f_uc(2)
                x.nv.x.c = Tools_2.Update_xc(m.At,m.Bt);
            end
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
                    Fig_V1_0_2D.Plot(1,msh,obj,"blk");
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
        end
    end
end