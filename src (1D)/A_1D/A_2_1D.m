classdef A_2_1D
    methods (Static)
        %% > Wrap-up A_2 (1D).
        function [inp,msh] = WrapUp_A_2_1D(inp)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            eg   = inp.msh.eg;
            
            switch eg
                case '1'
                    h               = inp.msh.h;
                    [inp,NX_c,Xd_x] = A_2_1D.SquareMesh_Uniform(inp,Xv_i,Xv_f,h);
                case '2'
                    NX_v            = inp.msh.Nv;
                    [inp,NX_c,Xd_x] = A_2_1D.SquareMesh_NonUniform_1(inp,Xv_i,Xv_f,NX_v);
                otherwise
                    return;
            end
            %  > Number of cells/faces.
            msh.c.NC = NX_c;
            msh.f.NF = NX_c+1;
            %  > Cell vertex coordinates.
            msh.f.Xv = Xd_x;
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [inp,NX_c,Xd_x] = SquareMesh_Uniform(inp,Xv_i,Xv_f,h)
            % >> Vertex coordinates.
            %  > Xd.
            NX_c = round(1./h.*(Xv_f-Xv_i));
            Xd_x = linspace(Xv_i,Xv_f,NX_c);

            % >> Boundary conditions.
            inp.pr.t.EB = 'Dirichlet'; %  > East  boundary type.
            inp.pr.t.WB = 'Dirichlet'; %  > West  boundary type.
        end
        % >> 1.2. ---------------------------------------------------------
        function [inp,NX_c,Xd_x] = SquareMesh_NonUniform_1(inp,Xv_i,Xv_f,NX_v)
            % >> Local variables.
            %  > Nf_Unf: Normalized computational domain coordinate (uniform distribution).
            %  > Pt    : Stretching location in "domain percentage".
            %  > B     : Stretching parameter.
            Nf_X = inp.msh.s_nu.Nf_X;
            Ks_X = inp.msh.s_nu.Ks_X;
            
            % >> X.
            Nf_Unf_X = linspace(0,1,NX_v);
            Pt_X     = (Nf_X-0)./(1-0);
            B_X      = 1./(2.*Ks_X).*log((1+(exp(Ks_X)-1).*(Pt_X))./(1+(exp(-Ks_X)-1).*(Pt_X)));
            %  > NF_X.
            i        = 1:NX_v;
            NF_X     = (Nf_X-0).*(1+sinh(Ks_X.*(Nf_Unf_X(i)-B_X))./sinh(Ks_X.*B_X))+0;

            % >> Vertex coordinates.
            %  > Xd.
            Xd_x = NF_X(1:NX_v).*(Xv_f-Xv_i);
            NX_c = NX_v-1;

            % >> Boundary conditions.
            inp.pr.t.EB = 'Dirichlet'; %  > East  boundary type.
            inp.pr.t.WB = 'Dirichlet'; %  > West  boundary type.
        end
    end
end