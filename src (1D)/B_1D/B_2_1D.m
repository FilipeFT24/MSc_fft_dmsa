classdef B_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [s,xn] = Update_stl(msh,stl_p,stl_s,stl_t,A,B,sn,bnd,v,g)
            % >> Compute/update...
            %  > ...stencil.
            s     = A_2_1D.Problem_SetUp(msh,stl_p,stl_s,stl_t,sn,bnd,v,g);
            xf    = s.xf;
            bnd_s = find(~cellfun(@isempty,s.f));
            bnd_v = cell2mat(s.bnd(bnd_s));
            
            %  > ...matrix A (cell dependent coefficients).
            %  > ...matrix B (face dependent coefficients w/o source term contribution).
            for i = 1:msh.c.NC
                %  > If cell i's west(w) face has been updated...
                if ismembc(i,stl_s)
                    %  > A.
                    kc_w       = s.c{i};
                    xf_w       = xf {i};
                    A(i,kc_w)  = A(i,kc_w)-xf_w(1:length(kc_w));
                    %  > B (boundary contribution(s)).
                    if ismembc(i,bnd_s)
                        vw     = bnd_v(bnd_s == i);
                        B(i,1) = B(i,1)+vw.*xf_w(length(kc_w)+1);
                    end
                end
                %  > If cell i's east(e) face has been updated...
                if ismembc(i+1,stl_s)
                    %  > A.
                    kc_e       = s.c{i+1};
                    xf_e       = xf {i+1};
                    A(i,kc_e)  = A(i,kc_e)+xf_e(1:length(kc_e));
                    %  > B (boundary contribution(s)).
                    if ismembc(i+1,bnd_s)
                        ve     = bnd_v(bnd_s == i+1);
                        B(i,1) = B(i,1)-ve.*xf_e(length(kc_e)+1);
                    end
                end
            end
            xn.c = A\B;
        end
        % >> 1.2. ---------------------------------------------------------
        function [e,x] = Update_pde(msh,x,sn,s)
            % >> Update pde...
            %  > ...face(s).
            for i  = 1:msh.f.NF
                kc = s.c{i};
                vc = x.c(kc);
                vt = vc;
                if ~isempty(s.f{i})
                    vf = s.bnd{i};
                    vt = [vt;vf];
                end
                j        = 1:2;
                x.f(i,j) = vt'*s.Tf{i}(j,:)';
            end
            %  > ...cell/face error/error norm(s).
            [e.c,e.f] = B_2_1D.Compute_Error(msh,sn,x,s);
        end
        % >> 1.3. ---------------------------------------------------------
        function [Ec,Ef] = Compute_Error(msh,sn,xn,s)
            % >> Cell(s).
            i         = 1:msh.c.NC;
            Vol(i,1)  = msh.c.Vol(i);
            %  > Absolute error.
            Ec.c(i,1) = abs(sn.c(i,1)-xn.c(i,1));
            %  > Error norms.
            Ec_1(i,1) = Ec.c(i,1).*Vol(i);
            Ec_2(i,1) = Ec_1(i,1).^2;
            Ec.n(1,1) = sum(Ec_1)./sum(Vol);
            Ec.n(2,1) = sum(sqrt(Ec_2))./sum(sqrt(Vol.^2));
            Ec.n(3,1) = max(Ec_1);
            
            % >> Face(s).
            i         = 1:msh.f.NF;
            Ls (i,1)  = s.Ls;
            %  > Absolute error: 1) Column #1: Absolute   error: Phi(f).
            %                    2) Column #2: Absolute   error: GradPhi(f).
            Ef.f(i,1) = abs(sn.f(i,1)-xn.f(i,1));
            Ef.f(i,2) = abs(sn.f(i,2)-xn.f(i,2));
            %  > Error norms.
            Ef_1(i,1) = Ef.f(i,1).*Ls(i,1);
            Ef_1(i,2) = Ef.f(i,2).*Ls(i,1);
            Ef_2(i,1) = Ef_1(:,1).^2;
            Ef_2(i,2) = Ef_1(:,2).^2;
            Ef.n(1,1) = sum(Ef_1(:,1))./sum(Ls);
            Ef.n(1,2) = sum(Ef_1(:,2))./sum(Ls);
            Ef.n(2,1) = sum(sqrt(Ef_2(:,1)))./sum(sqrt(Ls.^2));
            Ef.n(2,2) = sum(sqrt(Ef_2(:,2)))./sum(sqrt(Ls.^2));
            Ef.n(3,1) = max(Ef_1(:,1));
            Ef.n(3,2) = max(Ef_1(:,2));
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [msh,pde] = SetUp_Problem(msh,pde,v,g,bnd,p,p_adapt)
            %  > Initialize.
            An  = zeros(msh.c.NC);
            Bn  = zeros(msh.c.NC,1);
            Bn  = Bn+pde.FV;
            
            %  > Select...
            switch p_adapt
                case false
                    %  > Initialize.
                    stl_p  = repelem(p,msh.f.NF);
                    stl_s  = 1:msh.f.NF;
                    stl_t  = repelem("CDS",msh.f.NF);
                    %  > Solve PDE.
                    [s,xn] = B_2_1D.Update_stl(msh,stl_p,stl_s,stl_t,An,Bn,pde.sn,bnd,v,g);
                    [e,xn] = B_2_1D.Update_pde(msh,xn,pde.sn,s);
                    %  > Check error estimators.
                    B_2_1D.Check_EE("1","CDS",msh,stl_s,An,Bn,pde.sn,bnd,v,g);
                    %  > Update structure.
                    msh.s  = s;
                    pde.e  = e;
                    pde.x  = xn;   
                case true
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Check_EE(EE,DT,msh,stl_s,An,Bn,sn,bnd,v,g)
            switch EE
                %  > Higher-order face polynomial regression w/ lower-order values.
                case "1"
                    %  > Select discretization type...
                    switch DT
                        case "CDS"
                            %  > Auxiliary variables.
                            LO = 1;
                            HO = 7;
                            i  = LO:2:HO;
                            for j = 1:length(i)
                                ij         = i(j);
                                stl_p(j,:) = repelem(ij,msh.f.NF);
                                stl_t(j,:) = repelem(DT,msh.f.NF);
                            end
                            %  > Error estimation.
                            for j = 1:length(i)
                                %  > Exact error.
                                [SN1{j},XN{j}] = B_2_1D.Update_stl(msh,stl_p(j,:),stl_s,stl_t(j,:),An,Bn,sn,bnd,v,g);
                                [E1{j},XN{j}] = B_2_1D.Update_pde(msh,XN{j},sn,SN1{j});
                                if j == 1
                                    continue;
                                else
                                    %  > Estimated error.
                                    [E2{j-1},~] = B_2_1D.Update_pde(msh,XN{j-1},sn,SN1{j});
                                end
                            end
                            %  > Plot...
                            Fig_4_1D.WrapUp_Fig_4_1D(true,false,[4,5],msh.f.Xv,i,XN,E1,E2);
                        otherwise
                            return;
                    end
                case "2"
                    switch DT
                        case "CDS"
                            %  > Auxiliary variables.
                            LO = 1;
                            HO = 7;
                            i  = LO:2:HO;
                            for j = 1:length(i)
                                ij         = i(j);
                                stl_p(j,:) = repelem(ij,msh.f.NF);
                                stl_t(j,:) = repelem(DT,msh.f.NF);
                            end
                            %  > Error estimation.
                            for j = 1:length(i)
                                %  > Approximated values.
                                [SN1{j},XN{j}] = B_2_1D.Update_stl(msh,stl_p(j,:),stl_s,stl_t(j,:),An,Bn,sn,bnd,v,g);
                                [EN1{j},XN{j}] = B_2_1D.Update_pde(msh,XN{j},sn,SN1{j});
                                if j == 1
                                    continue;
                                else
                                    %  > Estimated error.
                                    EN2{j-1}(:,1) = abs(XN{j}.f(:,1)-XN{j-1}.f(:,1));
                                    EN2{j-1}(:,2) = abs(XN{j}.f(:,2)-XN{j-1}.f(:,2));
                                end
                            end
                            
                            hold on;
                            plot(msh.f.Xv,EN1{2}.f.f(:,2),'-r');
                            plot(msh.f.Xv,EN2{2}    (:,2),'-k');
                           
                            %  > Plot...
                            %Fig_4_1D.WrapUp_Fig_4_1D(true,false,[4,5],msh.f.Xv,i,E1,E2);
                            
                            
                        otherwise
                            return;   
                    end
            end
        end
    end
end