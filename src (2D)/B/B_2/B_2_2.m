classdef B_2_2
    methods (Static)
        %% > Wrap-up B_2_2.
        function [pde] = WrapUp_B_2_2(msh,pde,np_x,np_y,ng,wf,vx,vy,gx,gy,bnd_fc,bnd_ff,ft,st)
            %  > Select weight function coefficients.
            a   = 1;
            b   = 2;
            %  > Solve...
            pde = B_2_2.Assemble_A_B(msh,pde,np_x,np_y,ng,wf,a,b,vx,vy,gx,gy,bnd_fc,bnd_ff,ft,st);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [pde] = Assemble_A_B(msh,pde,np_x,np_y,ng,wf,a,b,vx,vy,gx,gy,bnd_fc,bnd_ff,ft,st)
            % >> X*Tf.
            %    Remark: Tf=[A,B,C,D,E,F,G,H,I,J,...], where: Cell dependent coefficients: A,B,C,...,G.
            %                                                 Face dependent coefficients: H,i,J,...,(...).
            %  > Tf.
            [Tf_C,Tf_D] = B_2_1.Assemble_Tf(msh,np_x,np_y,ng,wf,a,b,bnd_fc,bnd_ff);
                                                             
            % >> Face 'i'...
            for i = 1:msh.f.NF
                %  > Convection/diffusion contribution(s).
                V_Tf_C{i}(1,:) = vx.*Tf_C{i};
                V_Tf_C{i}(2,:) = vy.*Tf_C{i};
                G_Tf_D{i}(1,:) = gx.*Tf_D{i}(1,:);
                G_Tf_D{i}(2,:) = gy.*Tf_D{i}(2,:);
                %  > Cell/Face indices.
                Phi_fc{i}      = [msh.s.c{:,i}];
                if ~isempty(msh.s.xy_v_f{i})
                    Phi_ff{i}  = [msh.s.f{:,i}];
                end
            end
            for i = 1:msh.c.NC
                for j = 1:size(msh.c.f.Sf{i},2)
                    %  > Face 'j' of cell 'i'.
                    face     {i}(j) = msh.c.f.f {i}(j);
                    %  > Convection/diffusion contribution(s).
                    V_Tf_C_Sf{i}{j} = msh.c.f.Sf{i}(:,j)'*V_Tf_C{face{i}(j)};
                    G_Tf_D_Sf{i}{j} = msh.c.f.Sf{i}(:,j)'*G_Tf_D{face{i}(j)};
                end
            end
            
            % >> A.
            A = B_2_2.Assemble_A(msh.c.NC,face,Phi_fc,V_Tf_C_Sf,G_Tf_D_Sf);
            % >> B.
            B = B_2_2.Assemble_B(msh,pde.an.fn.func,st,ng,face,Phi_fc,Phi_ff,V_Tf_C_Sf,G_Tf_D_Sf,bnd_ff,pde.an.bnd);
            
            % >> Solve PDE...
            if strcmpi(ft,'Implicit')
                %  > Flux reconstruction: Implicit.
                pde = B_2_2.PDE_Implicit(msh,pde,A,B);
            elseif strcmpi(ft,'Explicit')
                %  > Flux reconstruction: Explicit.
                pde = B_2_2.PDE_Implicit(msh,pde,A,B);
            else
                return;
            end
        end
        % >> 1.3.1. -------------------------------------------------------
        function [A] = Assemble_A(NC,face,Phi_fc,V_Tf_C_Sf,G_Tf_D_Sf)
            %  > Initialize.
            A = zeros(NC,NC);
            
            % >> A (Cell dependent coefficients).
            %  > ... for cell #i: Phi_f(j_1)=A*Phi(X)+B*Phi(Y)+C*Phi(Z)+...
            %                     Phi_f(j_2)=D*Phi(M)+E*Phi(N)+F*Phi(O)+...
            %                     Phi_f(...)=...
            for i = 1:NC
                for j = 1:length(face{i})
                    %  > Index in the j-direction.
                    numb_c   {i}(j)     = length(Phi_fc{face{i}(j)});
                    kc_i     {i}(j)     = 1;
                    kc_f     {i}(j)     = numb_c{i}(j);
                    k                   = kc_i{i}(j):kc_f{i}(j);
                    ijk_c    {i}{j}     = Phi_fc{face{i}(j)}(k);
                    %  > Convection/diffusion contribution(s).
                    A(i,ijk_c{i}{j}(k)) = A(i,ijk_c{i}{j}(k))+V_Tf_C_Sf{i}{j}(k);
                    A(i,ijk_c{i}{j}(k)) = A(i,ijk_c{i}{j}(k))-G_Tf_D_Sf{i}{j}(k);
                end
            end
        end
        % >> 1.3.2. -------------------------------------------------------
        function [B] = Assemble_B(msh,func,st,ng,face,Phi_fc,Phi_ff,V_Tf_C_Sf,G_Tf_D_Sf,bnd_ff,bnd_fv)
            %  > Initialize.
            B = zeros(msh.c.NC,1);
            
            % >> B (Face dependent coefficients).
            %  > Compute source term.
            F_Vol = B_1_2.Compute_SourceTerm(msh,func,st,ng);
            for i = 1:msh.c.NC
                %  > Boundary contribution(s).
                for j = 1:length(face{i})
                    %  > Index in the j-direction.
                    numb_c{i}(j) = length(Phi_fc{face{i}(j)});
                    numb_f{i}(j) = length(Phi_ff{face{i}(j)});
                    kf_i  {i}(j) = numb_c{i}(j)+1;
                    kf_f  {i}(j) = numb_c{i}(j)+numb_f{i}(j);
                    k            = kf_i{i}(j):kf_f{i}(j);
                    ijk_f {i}{j} = Phi_ff{face{i}(j)}(k-numb_c{i}(j));
                    %  > Convection/diffusion contribution(s).
                    if all(isempty(ijk_f{i}{j}))
                        %  > Skip cells that don't have faces whose sentil does not include boundary faces.
                        break;
                    else
                        if ~isempty(ijk_f{i}{j})
                            for l = 1:length(ijk_f{i}{j})
                                %  > Indices.
                                bnd_j{i}{j}(l) = find(ijk_f{i}{j}(l) == bnd_ff);
                                %  > Boundary values.
                                bnd_v{i}{j}(l) = bnd_fv(bnd_j{i}{j}(l));
                            end
                        end
                    end
                    B(i,1) = B(i,1)-V_Tf_C_Sf{i}{j}(k)*bnd_v{i}{j}';
                    B(i,1) = B(i,1)+G_Tf_D_Sf{i}{j}(k)*bnd_v{i}{j}';
                end
                %  > Source term contribution(s).
                B(i,1) = B(i,1)+F_Vol(i);
            end
        end
        % >> 1.3.3. -------------------------------------------------------
        function [X] = SetUp_Solver(A,B,str)
            if strcmpi(str,'backslash')
                X       = A\B; % Equivalent to: V*((U'*B)./s).
            elseif strcmpi(str,'Tikhonov')
                lmbd    = 1E-06;
                [U,S,V] = svd_lapack(A);
                s       = diag(S);
                X       = V*(s.*(U'*B)./(s.^2+lmbd));
            else
                return;
            end
        end
        % >> 1.3.4. -------------------------------------------------------
        %  > Implicit flux reconstruction.
        function [pde] = PDE_Implicit(msh,pde,A,B)
            %  > Phi_c.
            pde.Phi = B_2_2.SetUp_Solver(A,B,'backslash');
            %  > Error norms.
            i       = 1:msh.c.NC;
            X(i)    = abs(pde.an.blk(i)-pde.Phi(i));
            pde.E   = B_2_2.Compute_ErrorNorms(msh,X);
            %  > Face polynomial.
            %    B_2_2.Compute_FacePolynomial(msh,Phi_fc,Phi_ff,pde.Phi,pde.bnd.f,bnd_ff,Tf_C,Tf_D);
        end
        % >> 1.3.5. -------------------------------------------------------
        %  > Explicit flux reconstruction.
        function [pde] = PDE_Explicit(msh,pde,A,B)
            %  > Error norms.
            X     = abs(A*pde.blk.f-B)';
            pde.E = B_2_2.Compute_ErrorNorms(msh,X);
            %  > Face polynomial.
            %    B_2_2.Compute_FacePolynomial(msh,Phi_fc,Phi_ff,pde.blk.f,pde.bnd.f,bnd_ff,Tf_C,Tf_D);
        end
        % >> 1.4. ---------------------------------------------------------
        function [] = Compute_FacePolynomial(msh,Phi_fc,Phi_ff,Phi_cv,Phi_fv,bnd_ff,Tf_C,Tf_D)
            for i = 1:msh.f.NF
                %  > Stencil values (a posteriori).
                nc (i)       = length(Phi_fc{i});
                nf (i)       = length(Phi_ff{i});
                j  {i}       = 1:nc(i);
                Phi{i}(j{i}) = Phi_cv(Phi_fc{i}(j{i}));
                k  {i}       = nc(i)+1:nc(i)+nf(i);
                if ~isempty(k{i})
                    for l = 1:length(k{i})
                        Phi{i}(k{i}) = Phi_fv(Phi_ff{i}(k{i}(l)-nc(i)) == bnd_ff);
                    end
                end
                %  > Phi_f,gradPhi_f.
                Phi_f(i)       = Tf_C{i}*Phi{i}';
                gradPhi_f(:,i) = Tf_D{i}*Phi{i}';
            end
        end
        
        %% > 2. -----------------------------------------------------------
        function [E]   = Compute_ErrorNorms(msh,X)
            i          = 1:msh.c.NC;
            E.EA(i)    = X(i);
            E_iX{1}(i) = X(i).*msh.c.h(i);
            E_iX{2}(i) = X(i).^2.*msh.c.h(i).^2;
            E.EN{1}    = sum(E_iX{1})./sum(msh.c.h);
            E.EN{2}    = sum(sqrt(E_iX{2}))./sum(sqrt(msh.c.h.^2));
            E.EN{3}    = max(X);
        end
    end
end