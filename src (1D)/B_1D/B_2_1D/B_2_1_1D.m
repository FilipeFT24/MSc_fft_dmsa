classdef B_2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [pde,s] = SetUp_p(inp,msh,pde,s,stl)
            %  > Update 's' structure and 'x' field.
            s   = B_2_1_1D.Update_s (inp,msh,pde,s,stl);
            x.c = B_2_1_1D.Update_xc(s); 
            %  > Compute cell/face error.
            e.c = B_2_1_1D.Update_ec(pde.av.c,x.c,msh.c.Vc);
            e.f = B_2_1_1D.Update_ef(pde.av.f,x.c,s);
            
            %  > Compute analytic/estimated truncation error.  
            if ~inp.pa.ee
                %  > w/ analytic field (nodal values).
                n.c = pde.av.c;
                n.f = pde.av.f;
            else
                %  > w/ error estimators.
            end
            e.t   = B_2_1_1D.Update_et_f(inp,s,n);
            e.t   = B_2_1_1D.Update_et_c(inp,msh,e.t);
            pde   = B_2_1_1D.Set_pde    (pde,x,e);
            pde.e = B_2_1_1D.Sort_pde_e (pde.e);
        end
               
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize matrix coefficients.
        function [u,A,B] = Initialize_coeffs(j,NC,stl_s,sc,xf,A,B,bnd_s,bnd_v)
            %  > Initialize flags.
            u.c.uA = 0;
            u.w.uB = 0;
            u.e.uB = 0;
            
            if any(ismembc([j,j+1],stl_s))
                % >> A.
                %  > Update...
                A      = zeros(1,NC);
                u.c.uA = 1;
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
                        u.w.uB = 1;
                        u.w.fv = bnd_v(bnd_s == j);
                        u.w.Bf = xf{j}(l.w+1);
                    end
                    if ismembc(j+1,bnd_s)
                        u.e.uB = 1;
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
        function [s] = Update_s(inp,msh,pde,s,stl)
            %  > Auxiliary variables.
            NC = msh.c.NC;
            k  = 1:size(stl.s,2);
            vg = [inp.pv.v(1),-inp.pv.v(2)];
            Ac = zeros(NC);
            Bc = zeros(NC,1);
            
            %  > Update stencil and check for nil (i.e. below trsh=10e-6) coefficients in 's.xf'.
            s = A_2_1D.Assemble_stl(inp,msh,pde,s,stl,stl.s);

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
            Bc   = Bc+pde.fn.st;
            %  > Update fields.
            s.A  = A;
            s.B  = B;
            s.Ac = Ac;
            s.Bc = Bc;
        end
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        %  > Update 'pde.x.c' field (nodal solution).
        function [xc] = Update_xc(s)
            xc = s.Ac\s.Bc;
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > Update 'pde.x.v' and 'pde.x.f' fields (reconstructed face values).
        %  > xv.
        function [vf] = Update_xv(s,x_c)
            %  > Auxiliary variables.
            [m,n] = size(s.xf);

            for i = 1:m
                for j = 1:n
                    if ~isempty(s.c{i,j})
                        %  > Nodal indices/values.
                        k = s.c{i,j};
                        v = x_c(k,1);
                        %  > Add boundary contribution?...
                        if ~isempty(s.f{i,j})
                            v = [v;s.bnd_v{i,j}];
                        end
                    else
                        v = s.bnd_v{i,j};
                    end
                    vf{i,j} = v;
                end
            end
        end
        %  > xf.
        function [xf] = Update_xf(s,x_c)
            %  > Auxiliary variables.
            [m,n] = size(s.xf);
            
            %  > Nodal values used to reconstruct face.
            vf = B_2_1_1D.Update_xv(s,x_c);
            for i = 1:m
                for j = 1:n
                    %  > Reconstructed face value.
                    xf(j,i) = s.xf{i,j}*vf{i,j};
                end
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        %  > Update 'pde.e.c' field (cell error).
        function [ec] = Update_ec(a_c,x_c,Vc)
            %  > Error/absolute error distribution.
            ec.c    (:,1) = a_c(:,1)-x_c(:,1);
            ec.c_abs(:,1) = abs(ec.c);
            %  > Error/absolute error norms.
            ec_1    (:,1) = ec.c.*Vc;
            ec_1_abs(:,1) = ec.c_abs.*Vc;
            ec_2    (:,1) = ec_1.^2;
            ec_2_abs(:,1) = ec_1_abs.^2;
            ec.n    (1,1) = sum(ec_1)./sum(Vc);
            ec.n_abs(1,1) = sum(ec_1_abs)./sum(Vc);
            ec.n    (2,1) = sum(sqrt(ec_2))./sum(sqrt(Vc.^2));
            ec.n_abs(2,1) = sum(sqrt(ec_2_abs))./sum(sqrt(Vc.^2));
            ec.n    (3,1) = max(ec.c);
            ec.n_abs(3,1) = max(ec.c_abs);
        end
        %  > 3.2.2. -------------------------------------------------------
        %  > Update 'pde.e.f' field (face error).
        function [ef] = Update_ef(a_f,x_c,s)
            %  > Auxiliary variables.
            x_ff          = B_2_1_1D.Update_xf(s,x_c);
            %  > Convective/diffusive components.
            i             = 1:size(x_ff,2);
            %  > Error/absolute error distribution.
            ef.f    (:,i) = a_f(:,i)-x_ff(:,i);
            ef.f_abs(:,i) = abs(ef.f(:,i));
            %  > Mean error/absolute error.
            ef.n    (1,i) = mean(ef.f(:,i));
            ef.n_abs(1,i) = mean(ef.f_abs(:,i));
        end
        %  > 3.2.3. -------------------------------------------------------
        %  > Update 'pde.e.t.f' field (face truncation error).
        function [et] = Update_et_f(inp,s,x)
            %  > Weighted (convective/diffusive) components.
            vg = [inp.pv.v(1),inp.pv.v(2)];
            nc = length(vg);
            i  = 1:nc;
            j  = 1:nc+1;
            k  = j(end);

            %  > Reconstructed face values (w/ input nodal field).
            x.s             = B_2_1_1D.Update_xf(s,x.c);
            %  > (Weighted) truncation/absolute truncation error distribution.
            et.f      (:,i) = vg(i).*(x.f(:,i)-x.s(:,i));
            et.f      (:,k) = et.f(:,k-2)-et.f(:,k-1);
            et.f_abs  (:,j) = abs(et.f);
            %  > Mean/absolute mean truncation error.
            et.n.f    (:,j) = mean(et.f(:,j));
            et.n_abs.f(:,j) = mean(et.f_abs(:,j));
        end
        %  > 3.2.4. -------------------------------------------------------
        %  > Update 'pde.e.t.c' field (cell truncation error).
        function [et] = Update_et_c(inp,msh,et)
            %  > Weighted (convective/diffusive) components.
            vg = [inp.pv.v(1),inp.pv.v(2)];
            nc = length(vg);
            i  = 1:nc;
            j  = 1:nc+1;
            k  = nc+1;

            %  > Truncation/absolute truncation error distribution.
            %  > Equivalent formulations (need to import s).
            %  et.c   (:,1) = s.Ac*x.c(:,1)-s.Bc(:,1);
            l               = 1:msh.c.NC;
            m               = l+1;
            et.c      (l,1) = et.f(l,k)-et.f(m,k);
            et.c_abs  (:,1) = abs(et.c);
           
            %  > Error/absolute error norms.
            Vc              = msh.c.Vc;
            ec_1            = et.c.*Vc;
            ec_1_abs        = et.c_abs.*Vc;
            ec_2            = ec_1.^2;
            ec_2_abs        = ec_1_abs.^2;
            et.n.c    (1,:) = sum(ec_1,1)./sum(Vc);
            et.n.c    (2,:) = sum(sqrt(ec_2),1)./sum(sqrt(Vc.^2));
            et.n.c    (3,:) = max(et.c);
            et.n_abs.c(1,:) = sum(ec_1_abs,1)./sum(Vc);
            et.n_abs.c(2,:) = sum(sqrt(ec_2_abs),1)./sum(sqrt(Vc.^2));
            et.n_abs.c(3,:) = max(et.c_abs);
        end
        % >> 3.3. ---------------------------------------------------------
        %  > 3.3.1. -------------------------------------------------------
        %  function [et] = Set_et(ef,ec)
        %      n         = [fieldnames(ef.n)',fieldnames(ec.n)';struct2cell(ef.n)',struct2cell(ec.n)'];
        %      et.n      = struct(n{:});
        %      n_abs     = [fieldnames(ef.n_abs)',fieldnames(ec.n_abs)';struct2cell(ef.n_abs)',struct2cell(ec.n_abs)'];
        %      et.n_abs  = struct(n{:});
        %      et.f      = ef.f;
        %      et.f_abs  = ef.f_abs;
        %      et.c      = ec.c;
        %      et.c_abs  = ec.c_abs;
        %  end
        %  > 3.3.2. -------------------------------------------------------
        function [pde] = Set_pde(pde,x,e)
            pde.x = x;
            pde.e = e;
        end
        % >> 3.4. ---------------------------------------------------------
        function [pde_e] = Sort_pde_e(pde_e)
            pde_e   = orderfields(pde_e  ,{'c','f','t'});
            pde_e.c = orderfields(pde_e.c,{'c','c_abs','n','n_abs'});
            pde_e.f = orderfields(pde_e.f,{'f','f_abs','n','n_abs'});
            pde_e.t = orderfields(pde_e.t,{'c','c_abs','f','f_abs','n','n_abs'});
        end
    end
end