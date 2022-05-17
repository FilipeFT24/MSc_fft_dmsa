classdef B3_1D
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
                    if i == "a" && inp.t_terms.allow
                        obj.e.(i){j}.m.f     = cell(1,nc);
                        obj.e.(i){j}.m.f_abs = cell(1,nc);
                    end
                    obj.e.(i){j}.t.c         = zeros(Nc,1);
                    obj.e.(i){j}.t.c_abs     = zeros(Nc,1);
                    obj.e.(i){j}.t.f         = zeros(Nf,nc+1);
                    obj.e.(i){j}.t.f_abs     = zeros(Nf,nc+1);
                    obj.e.(i){j}.t.n.c       = zeros(1,3);
                    obj.e.(i){j}.t.n_abs.c   = zeros(1,3);
                    obj.e.(i){j}.t.n.f       = zeros(3,nc+1);
                    obj.e.(i){j}.t.n_abs.f   = zeros(3,nc+1);
                    obj.e.(i){j}.c.c         = zeros(Nc,1);
                    obj.e.(i){j}.c.c_abs     = zeros(Nc,1);
                    obj.e.(i){j}.c.n         = zeros(1,3);
                    obj.e.(i){j}.c.n_abs     = zeros(1,3);
                end
            end
            %  > Field: 'f' (analytic function).
            obj.f.fh = A3_1D.Update_fh(inp);
            obj.f.av = A3_1D.Update_av(msh,obj.f.fh.f);
            obj.f.bd = A3_1D.Update_bd(inp,msh,obj.f.fh.f);
            obj.f.st = A3_1D.Update_st(msh,obj.f.fh.func.i);             
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
                %  > Field: 's' (stencil coordinates,etc.).
                obj.s{i}.c       = cell (Nf,nc);
                obj.s{i}.f       = cell (Nf,nc);
                obj.s{i}.i       = cell (Nf,nc);
                obj.s{i}.logical = cell (Nf,nc);
                obj.s{i}.t       = cell (Nf,nc);
                %  > Field: 'u' (update flag).
                obj.u{i}.f  = ["a","x"];
                for j = 1:nc
                    obj.u{i}.p   (:,j) = repelem(inp.p(j)+A(i),Nf);
                    obj.u{i}.s{j}(:,1) = 1:Nf;
                end
                %  > Field: 'x' (nodal solution/ stencil coefficients,etc.).
                obj.x{i}.if     = cell (Nf,nc);
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
            s            = B1_1D.Update_1 (msh,f,s,u,add);
            x            = B1_1D.Update_2 (inp,msh,f,s,u,x);
            %  > Update matrices(?).
            if f_uc(1)
                m        = B2_1D.Update_3 (inp,msh,f,m,s,u,x);
            end
            %  > Update nodal solution(?).
            if f_uc(2)
                x.nv.x.c = B2_1D.Update_xc(m.At,m.Bt);
            end
            %  > Update face values(?).
            if f_uc(3)
                x        = B2_1D.Update_4 (f,s,u,x);
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Set up 'p-standard' and 'p-adaptative' runs.
        function [obj] = Run_p(inp,msh)
            switch inp.p_adapt.allow
                case false
                    %  > 'p-standard' run.
                    obj = B3_1D.Initialize(inp,msh);
                    
                    %obj.u{1}.p = [1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;3,3;3,3;3,3;3,3;3,3;3,3;3,3;3,3;3,3;3,3;3,3;3,3;3,3;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1;1,1];
                    
                    
                    
                    obj = B3_1D.p_standard(inp,msh,obj);                    
                    %  > Plot...
                    Fig_V1_1_1D.Plot([0,0],inp,msh,obj);
                case true
                    %  > 'p-adaptative' run.
                    obj_init    = B3_1D.Initialize(inp,msh);
                    for i       = "a"
                        obj.(i) = B3_1D.p_adaptative(inp,msh,obj_init,i);
                        %  > Plot...
                        Fig_V1_2_1D.Plot(1,inp,msh,obj.(i));
                    end
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
                B3_1D.Update_all(inp,msh,obj.f,obj.m{ic},obj.s{ic},obj.u{ic},obj.x{ic},ch,f_uc);
            %  > Update fields 'e', 'm' and 'x'.
            [obj.e,obj.m,obj.s,obj.x]       = ...
                B2_1D.Update_e  (inp,msh,obj.e,obj.f,obj.m,obj.s,obj.u,obj.x,f_xc);
        end        
        % >> 2.3. ---------------------------------------------------------
        % >> 'p-adaptative' run.
        function [obj] = p_adaptative(inp,msh,obj,f)
            %  > Auxiliary variables.
            ch   = inp.b.change(1);
            f_uc = [1,1,1];
%             if f == "a"
%                 f = "da";
%             end

            %  > Initialize cycle count.
            i = 0;
            while 1
                %  > Update cycle count...
                i = i+1;
                %  > Update fields 'm', 's' and 'x'.
                ic                              = 1;
                [obj.m{ic},obj.s{ic},obj.x{ic}] = ...
                    B3_1D.Update_all(inp,msh,obj.f,obj.m{ic},obj.s{ic},obj.u{ic},obj.x{ic},ch,f_uc);
                %  > Update fields 'e', 'm' and 'x'.
                f_xc                      = B3_1D.f_sol(i);
                [obj.e,obj.m,obj.s,obj.x]       = ...
                    B2_1D.Update_e  (inp,msh,obj.e,obj.f,obj.m,obj.s,obj.u,obj.x,f_xc);
                
                %  > Assign structures (auxiliary variables).
                obj_e(i,:) = obj.e;
                obj_m(i,:) = obj.m;
                obj_x(i,:) = obj.x;
                %   > Stop adaptation(?).
                if ~B3_1D.Stop(inp,i,obj.e.(f){1}.c.n_abs)
                    ic     = 1;
                    obj.u  = B3_1D.Set_u(inp,obj.e.(f){ic}.t.f_abs,obj.u);
                else
                    obj.e  = obj_e;
                    obj.m  = obj_m;
                    obj.x  = obj_x;
                    break;
                end
                fprintf("Loop: #%3d\n",i);
            end
        end 
        %  > 2.3.1. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
        function [u] = Set_u(inp,v,u)
            A     = 2;
            trsh  = inp.p_adapt.lambda.*max(v(:,3));
            fs_n  = find(v(:,3) > trsh);
            for i = 1:size(u,2)
                for j = 1:size(u{i}.p,2)
                    u_s{i}{j}      = fs_n;
                    u{i}.p(fs_n,j) = u{i}.p(fs_n,j)+A;
                end
                u{i}.s = u_s{i};
            end          
        end
        %  > 2.3.2. -------------------------------------------------------
        %  > Solve AX=B(?) criterion/criteria.
        function [flag] = f_sol(i)
            %if i ~= 1
            %    flag = 0;
            %else
                flag = 1;
            %end
        end
        %  > 2.3.3. -------------------------------------------------------
        %  > Set stopping criterion/criteria.
        function [flag] = Stop(inp,c,ec_abs)
            flag = 0;
            if c > inp.p_adapt.nc
                flag = 1;
                fprintf("Stopping criterion: max. number of cycles.\n");
            end
            if ec_abs(1) <= inp.p_adapt.ec_m
                flag = 1;
                fprintf("Stopping criterion: min. error treshold (L1 norm).\n");
            end
        end
    end
end