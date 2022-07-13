classdef B2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "m" (matrices).
        function [m] = Initialize_m(Nc)
            m.A      = zeros(Nc);
            m.b      = zeros(Nc,1);
            m.nnz.A  = zeros(Nc,1);
            m.nnz.At = 0;
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update matrix A and vector b.
        function [m] = Update_m(msh,f,m,s)
            %  > Rows to be updated (r).
            r        = RunLength(sort([msh.f.ic{s.u.s}]));
            m.A(r,:) = 0;
            m.b(r,1) = f.st(r);
            
            %  > For each cell...
            for i = r
                %  > For each face...
                for j  = 1:numel(msh.c.f.if(i,:)), k = msh.c.f.if(i,j);
                    %  > Auxiliary variables.
                    o  = s.logical {k};
                    a  = s.i       {k}( o);
                    b  = s.i       {k}(~o);
                    Sf = msh.c.f.Sf{i}( j,:);
                    %  > For each term...
                    for l = 1:size(s.x.s.wVtf,2)
                        %  > A.
                        m.A(i,a) = m.A(i,a)+Sf*s.x.s.wVtf{k,l}{1}(:,o);
                        %  > b.
                        if any(~o)
                            m.b(i,1) = m.b(i,1)-Sf*s.x.s.wVtf{k,l}{1}(:,~o)*f.bd.v(func.find_c(f.bd.i,b));
                        end
                        %  > Add kf (scalar) to the RHS...
                        if ~isempty(s.x.s.wVtf{k,l}{2})
                            m.b(i,1) = m.b(i,1)-Sf*s.x.s.wVtf{k,l}{2};
                        end
                    end
                end
                m.nnz.A(i) = nnz(m.A(i,:));
            end
            m.nnz.At = nnz(m.A);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update nodal solution/face values (multiplied by V).
        function [s] = Update_sx(inp,msh,f,m,s)
            %  > Auxiliary variables.
            fn = convertCharsToStrings(fieldnames(s.x.x.vf));
            if any(~cellfun(@isempty,s.x.x.vf.a))
                c.a = s.u.s';
            else
                c.a = 1:size(s.x.x.vf.a,1);
            end
            c.x = 1:size(s.x.x.vf.x,1);
            
            %  > Update nodal solution.
            s.x.x.nv.x = m.A\m.b;
            %  > Update face coefficients/face values.
            for i = 1:numel(fn)
                % >> Update nodal values used to fit face "j"...
                for j = c.(fn(i))
                    %  > Auxiliary variables.
                    n = s.logical{j};
                    a = s.i      {j}( n);
                    b = s.i      {j}(~n);
                    %  > Cell value(s).
                    s.x.x.vf.(fn(i)){j}      = zeros(numel(n),1);
                    s.x.x.vf.(fn(i)){j}(n,1) = s.x.x.nv.(fn(i))(s.i{j}(n));
                    %  > Face value(s).
                    if any(~n)
                        s.x.x.vf.(fn(i)){j}(~n,1) = f.bd.v(func.find_c(f.bd.i,b));
                    end
                end
                % >> Update coefficients/face values...
                cc = c.(fn(i))( inp.M.Cls & ~msh.f.logical(c.(fn(i))));
                cu = c.(fn(i))(~inp.M.Cls | (inp.M.Cls & msh.f.logical(c.(fn(i)))));
                %  > cf.
                for j = cc
                    s.x.x.cf.(fn(i)){j} = s.D{j}.Pf*s.x.x.vf.(fn(i)){j}+s.D{j}.kf;
                end
                for j = cu
                    s.x.x.cf.(fn(i)){j} = s.D{j}.DTWD_DTW*s.x.x.vf.(fn(i)){j};
                end
                %  > xfV.
                for j = c.(fn(i))
                    for k = 1:size(s.x.s.wVdf,2)
                        s.x.x.xfV.(fn(i)){k}(j,:) = s.x.s.wVdf{j,k}   *s.x.x.cf.(fn(i)){j};
                        %                         = s.x.s.wVtf{j,k}{1}*s.x.x.vf.(fn(i)){j};                     > w/o constraint(s).
                        %                         = s.x.s.wVtf{j,k}{1}*s.x.x.vf.(fn(i)){j}+s.x.s.wVtf{j,k}{2};  > w/  constraint(s).
                    end
                end
                %  > xfT (check "a posteriori" contribution of the polynomial terms in the x and y-directions w/ linear profile).
                for j = cc
                    for k = 1:numel(s.pf{j}), p{k} = s.pf{j}{k};
                        t{k} = func.cls_t(s.D{j}.bf,s.D{j}.Cf(:,p{k}),s.D{j}.DTW(p{k},:),s.D{j}.DTWD(p{k},p{k}),s.D{j}.DTWD(p{k},p{k})\s.D{j}.DTW(p{k},:));
                    end
                    for k = 1:size(s.x.x.xfT.(fn(i)),1)
                        for l = 1:size(s.x.x.xfT.(fn(i)),2)
                            for o = 1:size(s.x.x.xfT.(fn(i)){k,l},2)-1
                                s.x.x.xfT.(fn(i)){k,l}(j,o) = ...
                                    s.x.s.wVdf{j,k}(l,p{o})*(t{o}{1}*s.x.x.vf.(fn(i)){j}+t{o}{2});
                            end
                            s.x.x.xfT.(fn(i)){k,l}(j,size(s.x.x.xfT.(fn(i)){k,l},2)) = ...
                                s.x.x.xfV.(fn(i)){k}(j,l);
                        end
                    end
                end
                for j = cu
                    for k = 1:numel(s.pf{j}), p{k} = s.pf{j}{k};
                        t{k} = func.backslash(s.D{j}.DTWD(p{k},p{k}),s.D{j}.DTW(p{k},:));
                    end
                    for k = 1:size(s.x.x.xfT.(fn(i)),1)
                        for l = 1:size(s.x.x.xfT.(fn(i)),2)
                            for o = 1:size(s.x.x.xfT.(fn(i)){k,l},2)-1
                                s.x.x.xfT.(fn(i)){k,l}(j,o) = ...
                                    s.x.s.wVdf{j,k}(l,p{o})*(t{o}*s.x.x.vf.(fn(i)){j});
                            end
                            s.x.x.xfT.(fn(i)){k,l}(j,size(s.x.x.xfT.(fn(i)){k,l},2)) = ...
                                s.x.x.xfV.(fn(i)){k}(j,l);
                        end
                    end
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize field "e" (error): fields "e.a", "e.f" and "e.p".
        function [e] = Initialize_e(Nc,Nf)
            %  > Error distribution.
            for i = ["a","f","p"]
                if i ~= "f"
                    e.(i).c.c_abs = zeros(Nc,1);
                    e.(i).t.c     = zeros(Nc,1);
                    e.(i).t.c_abs = zeros(Nc,1);
                    e.(i).t.f_abs = zeros(Nf,3);
                else
                    e.(i) = cell(Nf,1);
                end
            end
            %  > Error norms.
            for i = ["a","p"]
                e.(i).n_abs.c   = zeros(3,1);
                e.(i).n_abs.t.c = zeros(3,1);
                e.(i).n_abs.t.f = zeros(3,3);
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Update field "e.(...)" (error).
        function [e] = Update_e(inp,msh,e,m,s)
            %  > #1: Update field "e.a".
            e.a = B2_2D.Update_ea(msh,e.a,m,s);
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Update "e.a(...)" field (analytic cell/face truncation error distribution/norms).
        function [e] = Update_ea(msh,e,m,s)
            % >> f.
            %  > \tau_f: {1}-\tau_f(\phi).
            %            {2}-\tau_f(\nabla\phi).
            %            {3}-\tau_f.
            k = 1;
            for i = 1:msh.f.Nf
                c       = msh.f.ic  {i}(k);
                Sf(i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == i,:);
            end
            [n,o] = size(s.x.x.xfT.a);
            for j = 1:o
                for i = 1:n
                    etf_ij{i,j} = s.x.x.xfT.a{i,j}-s.x.x.xfT.x{i,j};
                end
                etf_j {j} = sum(cat(3,etf_ij{1:n,j}),3);
                etf_jS{j} = etf_j{j}.*Sf(:,j);
            end
            e.t.f_abs = abs(sum(cat(3,etf_jS{:}),3));
            %  > Error norms.
            e.n_abs.t.f = func.Set_n(e.t.f_abs);
            
            % >> c.
            %  > \tau_c.
            for j = 1:o
                etf_jl(:,j) = etf_j{j}(:,o+1);
            end
            %  > Reset...
            if any(e.t.c)
                e.t.c = zeros(msh.c.Nc,1);
            end
            for i = 1:msh.c.Nc
                e.t.c(i,1) = sum(sum(msh.c.f.Sf{i}'.*etf_jl(msh.c.f.if(i,:),:)',1),2);
            end
            %  > \tau_c_abs.
            e.t.c_abs   = abs(e.t.c);
            %  > e_c_abs.
            e.c.c_abs   = abs(s.x.x.nv.a-s.x.x.nv.x);
            %           = abs(m.A\e.t.c);
            %  > Error norms.
            e.n_abs.c   = func.Set_n(e.c.c_abs,msh.c.Volume);
            e.n_abs.t.c = func.Set_n(e.t.c_abs,msh.c.Volume);
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Update "e.f(...)" field (fitting error distribution).
        %  > e = sqrt(W)*[\phi_s-D*b], where b = (DTWD\DTW)*\phi_s.
        function [e] = Update_ef(msh,e,s)
            for i = 1:size(s.D,1)
                b{i} = s.D{i}.DTWD\s.D{i}.DTW*s.x.x.vf.x{i};
                e{i} = abs(sqrt(s.D{i}.W)*(s.x.x.vf.x{i}-s.D{i}.D*b{i}));
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Update structure "u".
        function [s,v,flag] = Update_u(inp,msh,e,s,block)
            %  > Select variables.
            p = ceil(s.u.p./2);
            if inp.P_Ad.Isotropic
                x = e(:,3);
                y = p(:,1);
            else
                x = e(:,1:2);
                y = p(:,1:2);
            end
            %  > Set structures "f" and "v".
            for i = 1:msh.f.Nf
                F{i,1} = B2_2D.find_F(inp,msh,i);
                for j  = 1:size(x,2)
                    P(j).f(i,1) = y(i,j);    %  > Face "f".
                    P(j).F{i,1} = y(F{i},j); %  > Face "f" neighbours (depends on choosen configuration).
                end
            end
            [f,v] = B2_2D.Set_fv(inp,msh,block,F,P,x);
            
            %  > Update structure "u".
            if ~all(cellfun(@isempty,{f(:).c}) & cellfun(@isempty,{f(:).r}))
                flag = false;
                A    = 2;
                if inp.P_Ad.Isotropic
                    s.u.p(f.c,:) = s.u.p(f.c,:)-A;
                    s.u.p(f.r,:) = s.u.p(f.r,:)+A;
                else
                    for i = 1:size(f,2)
                        s.u.p(f(i).c,i) = s.u.p(f(i).c,i)-A;
                        s.u.p(f(i).r,i) = s.u.p(f(i).r,i)+A;
                    end
                end
                s.u.s = RunLength(sort([cat(1,f(:).r);cat(1,f(:).c)]));
            else
                flag = true;
                fprintf("Stopping criterion: couldn't coarsen/refine any further.\n");
            end
        end
        %  > 3.1.1. -------------------------------------------------------
        %  > Find face "f" neighbours (depends on choosen configuration).
        function [F] = find_F(inp,msh,f)
            if inp.P_Ad.Config == 1 || 3
                F1 = RunLength(sort(func.setdiff(cat(1,msh.v.if{msh.f.iv(f,:)})',f)));
            end
            if inp.P_Ad.Config == 2 || 3
                F2 = RunLength(sort(func.setdiff(reshape(msh.c.f.if(msh.f.ic{f},:),1,[]),f)));
            end
            switch inp.P_Ad.Config
                case 1, F = F1;
                case 2, F = F2;
                case 3, F = RunLength(sort(cat(2,F1,F2)));
                otherwise
                    return;
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        %  > Select faces for coarsening/refinement.
        function [f,v] = Set_fv(inp,msh,block,F,P,x)
            %  > Auxiliary variables.
            m = numel(F);
            n = size (P,2);
            
            % >> Set structure "list".
            for i = 1:n
                %  > "i" and "v".
                [list(i).v,list(i).i] = sort(x(:,i),'ascend');
                %  > "logical".
                for j = 1:m
                    list(i).logical.c(j,1) = B2_2D.is_elig(0,P(i).f(j,1),P(i).F{j,1});
                    list(i).logical.r(j,1) = B2_2D.is_elig(1,P(i).f(j,1),P(i).F{j,1});
                end
                list(i).logical.c = list(i).logical.c(list(i).i);
                list(i).logical.r = list(i).logical.r(list(i).i);
            end
            
            % >> Select faces and check whether they can be automatically...
            flag = B2_2D.f_select(block,list,P,inp.P_Ad.Qrt,inp.P_Ad.Qrt_PT);
            %  > ...refined.
            for i = 1:n
                if all(flag(i).r.logical)
                    %  > Identify valid faces for refinement...
                    f(i).r = flag(i).r.i;
                    %  > Keep track of flagged/selected faces (for refinement).
                    for j = 1:numel(f(i).r)
                        v(i).r.list{j}{1} = [f(i).r(j),true];
                    end
                else
                    %  > Identify invalid faces for refinement...
                    f_i(i).r      = flag(i).r.i(~flag(i).r.logical);
                    %  > Add neighbour(s)/append to structure "f".
                    f_v(i).r      = B2_2D.Add_nb_r(f_i(i).r,F,P(i));
                    f  (i).r      = RunLength(sort(cat(1,flag(i).r.i(flag(i).r.logical),f_v(i).r.i)));
                    %  > Keep track of flagged/selected faces (for refinement).
                    v  (i).r.list = f_v (i).r.list;
                end
                v(i).r.f = flag(i).r.i; %  > ...flagged  for refinement.
                v(i).r.s = f   (i).r;   %  > ...selected for refinement.
            end
            %  > ...coarsened.
            for i = 1:n
                if all(flag(i).c.logical)
                    f(i).c = flag(i).c.i;
                end
                v(i).c.f = flag(i).c.i; %  > ...flagged  for coarsening.
            end
        end
        %  > 3.2.2. -------------------------------------------------------
        %  > 3.2.2.1. -----------------------------------------------------
        %  > Check whether face "f" is eligible for coarsening/refinement...
        function [logical] = is_elig(k,Pf,PF)
            if ~k
                logical = all(Pf >= PF); %  > ...for coarsening.
            else
                logical = all(Pf <= PF); %  > ...for refinement.
            end
        end
        %  > 3.2.2.2. -----------------------------------------------------
        %  > Select faces for coarsening/refinement...
        function [flag] = f_select(block,list,P,Qrt,Qrt_PT)
            %  > Auxiliary variables.
            [~,j] = sort(cat(1,list.v),1,'ascend');
            for i = 1:size(list,2)
                list_i(:,i) = list(i).i;
                list_v(:,i) = list(i).v;
            end
            f(:,1) = list_i(j);
            f(:,2) = list_v(j);

            %  > Sort...
            [~,a]  = unique(f(:,1)); f = f(sort(a),:);
            %  > Exclude blocked faces...
            if all(~block.r)
                k = f(:,2);
            else
                [~,b] = func.find_c(f(:,1),find(block.r)); k = f(~b,2);
            end
            n = ceil(size(k,1).*Qrt);
            
            %  > Check...
            if n(2) < Qrt_PT.*size(f,1)
                % >> Terminate...
                s.c = false(size(list_v));
                s.r = false(size(list_v));
            else
                % >> Select faces...
                %  > ...for coarsening.
                if Qrt(1) == 0
                    v (1)  = 0;
                else
                    v (1)  = B2_2D.VT(0,3,k(n(1))); % > w/ tolerance.
                end
                s.c = list_v < v(1) & ~block.c(f(:,1),:);
                %  > ...for refinement.
                if Qrt(2) == 1
                    v (2)  = k(end);
                else
                    v (2)  = B2_2D.VT(1,3,k(n(2))); % > w/ tolerance.
                end
                s.r = list_v > v(2) & ~block.r(f(:,1),:);
            end
            %  > ...for each direction.
            for i = 1:size(list,2)
                %  > "i".
                flag(i).c.i      (:,1) = list(i).i(s.c(:,i) & P(i).f > 1); %  > Indices to coarsen.
                flag(i).r.i      (:,1) = list(i).i(s.r(:,i));              %  > Indices to refine.
                %  > "logical".
                flag(i).c.logical(:,1) = list(i).logical.c(s.c);           %  > Can it be automatically coarsened?
                flag(i).r.logical(:,1) = list(i).logical.r(s.r);           %  > Can it be automatically refined?
            end
        end
        %  > 3.2.2.3. -----------------------------------------------------
        %  > Add neighbours (depends on the selected configuration)...
        function [nb_f] = Add_nb_r(f,F,P)
            %  > Find neighbours and assign to "nb_f".
            nb     = B2_2D.Find_nb(1,f,F,P);
            nb_f.i = nb.i(nb.logical);
            %  > Keep track of flagged/selected faces (for refinement).
            o = 1;
            for i = 1:numel(f)
                nb_f.list{i}{o}   = [f(i),true];
                nb_f.list{i}{o+1} = [nb.list(i).i,nb.list(i).logical];
            end
            o = o+1;
            
            %  > Until all faces are valid or no more faces can be added...
            while ~all(nb.logical)
                %  > Find neighbours and assign to "nb_f".
                g      = nb.i(~nb.logical);
                nb     = B2_2D.Find_nb(1,g,F,P);
                nb_f.i = [nb_f.i;nb.i(nb.logical)];
                %  > Keep track of flagged/selected faces (for refinement).
                for i = 1:numel(f)
                    nb_f.list{i}{o+1} = [];
                end
                for i = 1:numel(g)
                    for j = 1:numel(f)
                        if ~isempty(nb_f.list{j}{o})
                            if ismembc(g(i),sort(nb_f.list{j}{o}(:,1)))
                                nb_f.list{j}{o+1} = cat(1,nb_f.list{j}{o+1},[nb.list(i).i,nb.list(i).logical]);
                            end
                        end
                    end
                end
                o = o+1;
            end
            for i = 1:numel(f)
                for j = 1:numel(nb_f.list{i})
                    nb_f.list{i}{j} = unique(nb_f.list{i}{j},'rows');
                end
            end
            nb_f.i = RunLength(sort(nb_f.i));
        end
        %  > 3.2.2.3.1. ---------------------------------------------------
        %  > Auxiliary function (to be used recursively in while loop).
        function [nb] = Find_nb(k,f,F,P)
            %  > Check what faces need to be coarsened/refined to coarsen/refine face "f"...
            if ~k
                for i = 1:numel(f)
                    nb_f(i).i = F{f(i)}(P.F{f(i)} > P.f(f(i))); %  > ...for coarsening.
                end
            else
                for i = 1:numel(f)
                    nb_f(i).i = F{f(i)}(P.F{f(i)} < P.f(f(i))); %  > ...for refinement.
                end
            end
            %  > Can its neighbour(s) be coarsened/refined(?).
            for i = 1:numel(f)
                for j = 1:numel(nb_f(i).i)
                    nb_f(i).logical(j) = B2_2D.is_elig(k,P.f(nb_f(i).i(j)),P.F{nb_f(i).i(j)});
                end
            end
            for i = ["i","logical"]
                %  > Sort.
                if i == "i"
                    [~,a] = sort([nb_f(:).i]);
                end
                nb.(i) = [nb_f(:).(i)]; nb.(i) = nb.(i)(a);
                %  > Exclude duplicates.
                if i == "i"
                    b = [true,nb.i(2:end)-nb.i(1:end-1)>0];
                end
                nb.(i) = nb.(i)(b)';
            end
            for i = 1:numel(f)
                for j = ["i","logical"]
                    nb.list(i).(j) = nb_f(i).(j)';
                end
            end
        end
        % >> 3.3. ---------------------------------------------------------
        %  > 3.3.1. -------------------------------------------------------
        %  > Check stopping criteria.
        function [flag] = Stop(inp,obj)
            %  > Auxiliary variables.
            flag = false;
            
            if inp.P_Ad.ec  >= obj(end).e.a.n_abs.c(1)
                flag = true;
                fprintf("Stopping criterion: min. error treshold (L1 norm).\n");
            end
            if inp.P_Ad.Nc  <  numel(obj)
                flag = true;
                fprintf("Stopping criterion: max. number of cycles.\n");
            end
            if inp.P_Ad.NNZ <= obj(end).m.nnz.At
                flag = true;
                fprintf("Stopping criterion: max. nnz(A).\n");
            end
        end
        %  > 3.3.2. -------------------------------------------------------
        %  > Monitor \tau_f.
        function [block] = rst(inp,msh,block,obj_e,v)
            %  > Auxiliary variables.
            x  = arrayfun(@(x) x.a.t.f_abs  (:,3),obj_e,'un',0);
            e  = x(end-1:end);
            y  = arrayfun(@(x) x.a.n_abs.t.f(1,3),obj_e);
            z  = arrayfun(@(x) x.a.n_abs.t.f(3,3),obj_e);
            nm = y(end-1:end);
            nM = z(end-1:end);
            
            %  > L1.
            f(1) = nm(1) < B2_2D.VT(1,3,nm(2));
            g{1} = v.r.s(e{2}(v.r.s) > e{1}(v.r.s) & e{1}(v.r.s) > 1E-10);
            if f(1) && isempty(g{1})
                g{1} = v.r.s;
            end
            %  > L3.
            f(2) = nM(1) < B2_2D.VT(1,3,nM(2));
            if f(2)
                g{2} = find(e{2} > B2_2D.VT(1,3,max(e{2})));
                if all(~ismembc(g{2},g{1}))
                    f(2) = false;
                end
            end
            %  > Select/reset...
            if all(~f)
                for j = ["c","r"]
                    block.(j)(:) = false;
                end
            else
                if f(1)
                    block = B2_2D.Block_f(inp,msh,block,g{1},v.r.list); %  > "f" (flagged).
                    block = B2_2D.Block_s(inp,msh,block,g{1});          %  > "s" (selected).
                end
                if f(2)
                    block = B2_2D.Block_s(inp,msh,block,g{2});          %  > "s" (selected).
                end
            end
        end
        %  > 3.3.2.1. -----------------------------------------------------
        function [block] = Block_f(inp,msh,block,k,list)
            for i = 1:numel(list)
                flag  = false;
                for j = 1:numel(list{i})
                    if ~isempty(list{i}{j})
                        if any(ismembc(k,list{i}{j}(:,1)))
                            flag = true;
                            break;
                        end
                    end
                end
                if flag
                    a          = cat   (1,list{i}{:});
                    b          = a     (:,1);
                    c          = unique(cat(1,b,B2_2D.find_F(inp,msh,list{i}{1}(1))'));
                    block.r(c) = true;
                end
            end
        end
        %  > 3.3.2.2. -----------------------------------------------------
        function [block] = Block_s(inp,msh,block,k)
            for i = 1:numel(k)
                block.r([k(i),B2_2D.find_F(inp,msh,k(i))]) = true;
            end
        end
        % >> 3.4. ---------------------------------------------------------
        %  > Add tolerance to value "v".
        function [v] = VT(k,l,v)
            if ~k
                v = v+10.^(ceil(log10(v)-1)-l); % > v* > v.
            else
                v = v-10.^(ceil(log10(v)-1)-l); % > v* < v.
            end
        end
    end
end