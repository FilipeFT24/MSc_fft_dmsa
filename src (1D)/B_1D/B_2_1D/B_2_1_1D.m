classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [pde,s,stl] = SetUp_p(obj,msh,pde,s,stl)
            %  > Update 's' and 'stl' structures and create 'x' field.
            [s,stl,x.c]   = B_2_1_1D.Update_s     (obj,msh,pde,s,stl);
            %  > Update 'x' field.
            [x.f.a,x.v.a] = B_2_1_1D.Update_pde_x (msh,s,pde.a);
            [x.f.f,x.v.f] = B_2_1_1D.Update_pde_x (msh,s,x);
            %  > Update 'e' field.
            [e.c]         = B_2_1_1D.Update_pde_ec(pde.a,x,msh.c.Vol);
            [e.f]         = B_2_1_1D.Update_pde_ef(pde.a,x);
            [e.t]         = B_2_1_1D.Update_pde_et(obj,msh,s,pde.a,x.f.a,pde.f.st); 
            %  > Set 'pde' fields.
            pde.x         = x;
            pde.e         = e;
            pde.e         = B_2_1_1D.Sort_pde_e(pde.e);
        end
               
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize matrix coefficients.
        function [u,A,B] = Initialize_coeffs(j,NC,stl_s,sc,xf,A,B,bnd_s,bnd_v)
            %  > Initialize flags.
            u.c.uA = false;
            u.w.uB = false;
            u.e.uB = false;
            
            if any(ismembc([j,j+1],stl_s))
                % >> A.
                %  > Update...
                A      = zeros(1,NC);
                u.c.uA = true;
                %  > West face contribution.
                u.w.kc = sc{j};
                l.w    = length(u.w.kc);
                u.w.Af = xf{j}(1:l.w);
                %  > East face contribution.
                u.e.kc = sc{j+1};
                l.e    = length(u.e.kc);
                u.e.Af = xf{j+1}(1:l.e);
                
                % >> B.
                if any(ismembc([j,j+1],bnd_s))
                    %  > Update...
                    B  = 0;
                    if ismembc(j,bnd_s)
                        u.w.uB = true;
                        u.w.fv = bnd_v(bnd_s == j);
                        u.w.Bf = xf{j}(l.w+1);
                    end
                    if ismembc(j+1,bnd_s)
                        u.e.uB = true;
                        u.e.fv = bnd_v(bnd_s == j+1);
                        u.e.Bf = xf{j+1}(l.e+1);
                    end
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Update matrix A (cell dependent coefficients).
        function [A] = Assemble_A(str_f,A,xf)
            switch str_f
                case "w"
                    %  > West face contribution.
                    A = A-xf;
                case "e"
                    %  > East face contribution.
                    A = A+xf;
                otherwise
                    return;
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Update matrix B (face dependent coefficients w/o source term contribution).
        function [B] = Assemble_B(str_f,B,f,xf)
            switch str_f
                case "w"
                    %  > West face contribution.
                    B = B+f.*xf;
                case "e"
                    %  > East face contribution.
                    B = B-f.*xf;
                otherwise
                    return;
            end
        end
        % >> 2.4. ---------------------------------------------------------
        %  > Update 's' structure.
        function [s,stl,nc] = Update_s(obj,msh,pde,s,stl)
            %  > Auxiliary variables.
            NC = msh.c.NC;
            k  = 1:size(stl.s,2);
            vg = [obj.v,-obj.g];
            Ac = zeros(NC);
            Bc = zeros(NC,1);
            
            %  > Update stencil and check for nil (i.e. below trsh=10e-6) coefficients in 's.xf'.
            [s]     = A_2_1D.Assemble_stl(obj,msh,pde,s,stl,stl.s);
            [s,stl] = B_2_1_1D.Check_nil_coeffs(obj,msh,pde,s,stl,10e-6);           

            %  > Update  xf.
            for i = 1:size(s.f,1)
                bnd_s{i} = find(~cellfun(@isempty,s.f(i,:)));
                bnd_v{i} = [s.bnd_v{i,bnd_s{i}}];
            end
            %  > Update A/B.
            for i = k
                if ~isempty(stl.s{i})
                    for j = 1:NC
                        %  > Initialize.
                        [u,A{i}(j,:),B{i}(j,1)] = ...
                            B_2_1_1D.Initialize_coeffs(j,NC,stl.s{i},s.c(i,:),s.xf(i,:),s.A{i}(j,:),s.B{i}(j,1),bnd_s{i},bnd_v{i});
                        
                        %  > A.
                        if u.c.uA
                            A{i}(j,u.w.kc) = B_2_1_1D.Assemble_A("w",A{i}(j,u.w.kc),u.w.Af);
                            A{i}(j,u.e.kc) = B_2_1_1D.Assemble_A("e",A{i}(j,u.e.kc),u.e.Af);
                        end
                        %  > B.
                        if u.w.uB
                            B{i}(j,1) = B_2_1_1D.Assemble_B("w",B{i}(j,1),u.w.fv,u.w.Bf);
                        end
                        if u.e.uB
                            B{i}(j,1) = B_2_1_1D.Assemble_B("e",B{i}(j,1),u.e.fv,u.e.Bf);
                        end
                    end
                else
                    A{i} = s.A{i};
                    B{i} = s.B{i};
                end
                %  > Add convective/diffusive contributions (cumulative matrices).
                Ac = Ac+vg(i).*A{i};
                Bc = Bc+vg(i).*B{i};
            end
            %  > Add source term.
            Bc   = Bc+pde.f.st;
            %  > Nodal solution.
            nc   = Ac\Bc;
            %  > Update fields.
            s.A  = A;
            s.B  = B;
            s.Ac = Ac;
            s.Bc = Bc;
        end
        
        
        
        
        
        
        
        
        
        
        
        
        % >> 2.5. ---------------------------------------------------------
        %  > Check for nil coefficients in 's.xf' and decrease/increase method's order accordingly.
        %    If we're trying to use an UDS/DDS for the diffusive term on a uniform grid (boundaries NOT included), use a CDS of higher-order instead.
        function [s,stl] = Check_nil_coeffs(obj,msh,pde,s,stl,trsh)
            %  > Initialize.
            [m,n]     = size(s.xf);
            flag(1:m) = false;
            
            %  > Select faces and update convective/diffusive term(s).
            for i = 1:m
                k = 0;
                for j = 1:n
                    if any(abs(s.xf{i,j}) < trsh)
                        k              = k+1;
                        ij_nil{i}(k,1) = j;
                        flag  (i)      = true;
                    end
                end
                %  > Overwrite 'stl.p' and/or 'stl.t'.
                if flag(i)
                    if i == 1
                        %  > Convective term.
                        
                        
                    else
                        %  > Diffusive  term.
                        stl.t{i}(ij_nil{i}) = "CDS";
                    end
                end
            end
            if any(flag)
                [s] = A_2_1D.Assemble_stl(obj,msh,pde,s,stl,ij_nil);
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Update analytic 'pde.x.(...)' field.
        function [rf,vf] = Update_pde_x(msh,s,n)
            %  > Auxiliary variables.
            l = size(s.xf,1);
            m = msh.f.NF;
            
            %  > Face (reconstructed) values.
            for i = 1:m
                for j = 1:l
                    if ~isempty(s.c{j,i})
                        %  > Nodal indices/values.
                        k = s.c{j,i};
                        v = n.c(k,1);
                        %  > Add boundary contribution?...
                        if ~isempty(s.f{j,i})
                            v =[v;s.bnd_v{j,i}];
                        end
                    else
                        v = s.bnd_v{j,i};
                    end
                    %  > Reconstructed face value.
                    vf{j,i} = v;
                    rf(i,j) = s.xf{j,i}*v;
                end
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        %  > Update 'pde.e.c' field (cell error).
        function [ec] = Update_pde_ec(a,x,Vol)
            %  > Error/absolute error distribution.
            ec.c    (:,1) = a.c(:,1)-x.c(:,1);
            ec.c_abs(:,1) = abs(ec.c);
            %  > Error/absolute error norms.
            ec_1    (:,1) = ec.c.*Vol;
            ec_1_abs(:,1) = ec.c_abs.*Vol;
            ec_2    (:,1) = ec_1.^2;
            ec_2_abs(:,1) = ec_1_abs.^2;
            ec.n    (1,1) = sum(ec_1)./sum(Vol);
            ec.n_abs(1,1) = sum(ec_1_abs)./sum(Vol);
            ec.n    (2,1) = sum(sqrt(ec_2))./sum(sqrt(Vol.^2));
            ec.n_abs(2,1) = sum(sqrt(ec_2_abs))./sum(sqrt(Vol.^2));
            ec.n    (3,1) = max(ec.c);
            ec.n_abs(3,1) = max(ec.c_abs);
        end
        %  > 3.2.2. -------------------------------------------------------
        %  > Update 'pde.e.f' field (face error).
        function [ef] = Update_pde_ef(a,x)
            %  > Convective/diffusive components.
            i = 1:2;
            %  > Error/absolute error distribution.
            ef.f    (:,i) = a.f(:,i)-x.f.f(:,i);
            ef.f_abs(:,i) = abs(ef.f(:,i));
            %  > Mean error/absolute error.
            ef.n    (1,i) = mean(ef.f(:,i));
            ef.n_abs(1,i) = mean(ef.f_abs(:,i));
        end
        %  > 3.2.3. -------------------------------------------------------
        %  > Update 'pde.e.t' field (truncation error).
        function [et] = Update_pde_et(obj,msh,s,a,x_fa,st)
            %  > Weighted (convective/diffusive) components.
            vg = [obj.v,obj.g];
            nc = length(vg);
            i  = 1:nc;
            j  = nc+1;
            k  = 1:j;
            
            % >> Face(s).
            %  > (Weighted) truncation/absolute truncation error distribution.
            et.f      (:,i) = vg(i).*(a.f(:,i)-x_fa(:,i));
            et.f      (:,j) = et.f(:,j-2)-et.f(:,j-1);
            et.f_abs  (:,k) = abs(et.f);
            %  > Mean/absolute mean truncation error.
            et.n.f    (:,k) = mean(et.f(:,k));
            et.n_abs.f(:,k) = mean(et.f_abs(:,k));
            
            % >> Cell(s).
            %  > Truncation/absolute truncation error distribution.
            %  > Equivalent formulations.
            %  l            = 1:msh.c.NC;
            %  m            = l+1;
            %  et.c   (l,1) = et.f(l,j)-et.f(m,j);
            et.c      (:,1) = s.Ac*a.c(:,1)-s.Bc(:,1);
            et.c_abs  (:,1) = abs(et.c);
           
            %  > Error/absolute error norms.
            Vol             = msh.c.Vol;
            ec_1            = et.c.*Vol;
            ec_1_abs        = et.c_abs.*Vol;
            ec_2            = ec_1.^2;
            ec_2_abs        = ec_1_abs.^2;
            et.n.c    (1,:) = sum(ec_1,1)./sum(Vol);
            et.n.c    (2,:) = sum(sqrt(ec_2),1)./sum(sqrt(Vol.^2));
            et.n.c    (3,:) = max(et.c);
            et.n_abs.c(1,:) = sum(ec_1_abs,1)./sum(Vol);
            et.n_abs.c(2,:) = sum(sqrt(ec_2_abs),1)./sum(sqrt(Vol.^2));
            et.n_abs.c(3,:) = max(et.c_abs);
        end
        % >> 3.3. ---------------------------------------------------------
        function [pde_e] = Sort_pde_e(pde_e)
            pde_e   = orderfields(pde_e  ,{'c','f','t'});
            pde_e.c = orderfields(pde_e.c,{'c','c_abs','n','n_abs'});
            pde_e.f = orderfields(pde_e.f,{'f','f_abs','n','n_abs'});
            pde_e.t = orderfields(pde_e.t,{'c','c_abs','f','f_abs','n','n_abs'});
        end
    end
end