classdef A_2
    methods (Static)
        %% > Wrap-up A_2.
        function [struct,msh] = WrapUp_A_2(inp)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            pt   = inp.msh.pt;
            eg   = inp.msh.eg;
            dm   = inp.msh.dm;
            
            % >> Select polygon type...
            switch pt
                case 'v'
                    switch eg
                        case '1'
                            h      = inp.msh.h;
                            struct = A_2.fft_Distmesh2D(dm,Xv_i,Xv_f,Yv_i,Yv_f,h);
                        case '2'
                            NX_v   = inp.msh.Nv(1);
                            NY_v   = inp.msh.Nv(2);
                            xy_v   = A_2.TriangleMesh_NonUniform(Xv_i,Xv_f,Yv_i,Yv_f,NX_v,NY_v);
                            struct = delaunayTriangulation(xy_v(:,1),xy_v(:,2));
                        otherwise
                            return;
                    end
                case 's'
                    switch eg
                        case '1'
                            h                = inp.msh.h;
                            [NX_v,NY_v,xy_v] = A_2.SquareMesh_Uniform(Xv_i,Xv_f,Yv_i,Yv_f,h);
                        case '2'
                            NX_v = inp.msh.Nv(1);
                            NY_v = inp.msh.Nv(2);
                            xy_v = A_2.SquareMesh_NonUniform(inp,Xv_i,Xv_f,Yv_i,Yv_f,NX_v,NY_v);
                        otherwise
                            return;
                    end
                    %  > Number of cells/vertices.
                    [numb_C,numb_V] = ...
                        deal((NX_v-1).*(NY_v-1),NX_v.*NY_v);
                    %  > Connectivity list.
                    CList                        = reshape(1:numb_V,NY_v,NX_v);
                    struct.ConnectivityList(:,1) = reshape(CList(1:NY_v-1,1:NX_v-1),numb_C,1); % > SW.
                    struct.ConnectivityList(:,2) = reshape(CList(1:NY_v-1,2:NX_v  ),numb_C,1); % > SE.
                    struct.ConnectivityList(:,3) = reshape(CList(2:NY_v  ,2:NX_v  ),numb_C,1); % > NE.
                    struct.ConnectivityList(:,4) = reshape(CList(2:NY_v  ,1:NX_v-1),numb_C,1); % > NW.
                    %  > Points.
                    struct.Points(:,1) = xy_v(:,1);
                    struct.Points(:,2) = xy_v(:,2);
                otherwise
                    return;
            end
            %  > Number of cells.
            msh.c.NC = size(struct.ConnectivityList,1);
            %  > Cell vertex coordinates.
            for i = 1:msh.c.NC
                msh.c.xy_v{i}(:,1) = struct.Points(struct.ConnectivityList(i,:),1);
                msh.c.xy_v{i}(:,2) = struct.Points(struct.ConnectivityList(i,:),2);
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [struct] = fft_Distmesh2D(eg,Xv_i,Xv_f,Yv_i,Yv_f,h)
            %  > Bounding box.
            BBOX = [[Xv_i,Yv_i];[Xv_f,Yv_f]]; %  > [[Xmin,Ymin];[Xmax,Ymax]].
            switch eg
                case '1'
                    %% > Example 1: Uniform mesh on square.
                    %  > A,B,C,D.
                    A  = [Xv_i,Yv_i];
                    B  = [Xv_f,Yv_i];
                    C  = [Xv_f,Yv_f];
                    D  = [Xv_i,Yv_f];
                    %  > Function handles.
                    fd = @(p) drectangle(p,Xv_i,Xv_f,Yv_i,Yv_f);
                    fh = @(p) ones(size(p,1),1);
                    %  > Call 'distmesh2d'.
                    [struct.Points,struct.ConnectivityList] = ...
                        distmesh2d(fd,fh,h,[BBOX(1,:);BBOX(2,:)],[A;D;B;C]);
                case '2'
                    %% > Example 2: Circle/Ellipse.
                    %  > (rx,ry).
                    rx = 1./2.*(Xv_f-Xv_i);
                    ry = 1./2.*(Yv_f-Yv_i);
                    %  > (X0,Y0).
                    X0 = Xv_i+rx;
                    Y0 = Yv_i+ry;
                    %  > Function handles.
                    fd = @(p) (p(:,1)-X0).^2/rx.^2+(p(:,2)-Y0).^2/ry.^2-1;
                    %  > Call 'distmesh2d'.
                    [struct.Points,struct.ConnectivityList] = ...
                        distmesh2d(fd,@huniform,h,[BBOX(1,:);BBOX(2,:)],[]);
                otherwise
                    return;
            end
        end
        % >> 1.2. ---------------------------------------------------------
        function [xy_v] = TriangleMesh_NonUniform(Xv_i,Xv_f,Yv_i,Yv_f,NX_v,NY_v)
            %  > (Xd,Yd).
            Xd = sort((Xv_f-Xv_i).*rand(NY_v,NX_v-2)+Xv_i,2,'ascend');
            Xd = cat(2,ones(NY_v,1).*Xv_i,Xd,ones(NY_v,1).*Xv_f);
            Yd = sort((Yv_f-Yv_i).*rand(NY_v-2,NX_v)+Yv_i,1,'ascend');
            Yd = cat(1,ones(1,NX_v).*Yv_i,Yd,ones(1,NX_v).*Yv_f);
            %  > (Xv,Yv).
            xy_v(:,1) = reshape(Xd,[],1);
            xy_v(:,2) = reshape(Yd,[],1);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [NX,NY,xy_v] = SquareMesh_Uniform(Xv_i,Xv_f,Yv_i,Yv_f,h)
            %  > (Xd,Yd).
            [Xd_x,Yd_y] = deal(Xv_i:h:Xv_f,Yv_i:h:Yv_f);
            [Xd,Yd]     = meshgrid(Xd_x,Yd_y);
            %  > (NX,NY).
            NX = length(Xd_x);
            NY = length(Xd_x);
            %  > (Xv,Yv).
            xy_v(:,1) = reshape(Xd,[],1);
            xy_v(:,2) = reshape(Yd,[],1);
        end
        % >> 2.2. ---------------------------------------------------------
        function [xy_v] = SquareMesh_NonUniform(inp,Xv_i,Xv_f,Yv_i,Yv_f,NX_v,NY_v)
            % >> Local variables.
            %  > Nf_Unf: Normalized computational domain coordinate (uniform distribution).
            %  > Pt    : Stretching location in "domain percentage".
            %  > B     : Stretching parameter.
            Nf_X = inp.msh.s_nu.Nf_X;
            Nf_Y = inp.msh.s_nu.Nf_Y;
            Ks_X = inp.msh.s_nu.Ks_X;
            Ks_Y = inp.msh.s_nu.Ks_Y;
            
            % >> X.
            Nf_Unf_X = linspace(0,1,NX_v);
            Pt_X     = (Nf_X-0)./(1-0);
            B_X      = 1./(2.*Ks_X).*log((1+(exp(Ks_X)-1).*(Pt_X))./(1+(exp(-Ks_X)-1).*(Pt_X)));
            %  > NF_X.
            i        = 1:NX_v;
            NF_X     = (Nf_X-0).*(1+sinh(Ks_X.*(Nf_Unf_X(i)-B_X))./sinh(Ks_X.*B_X))+0;
            % >> Y.
            Nf_Unf_Y = linspace(0,1,NY_v);
            Pt_Y     = (Nf_Y-0)./(1-0);
            B_Y      = 1./(2.*Ks_Y).*log((1+(exp(Ks_Y)-1).*(Pt_Y))./(1+(exp(-Ks_Y)-1).*(Pt_Y)));
            %  > NF_Y.
            j        = 1:NY_v;
            NF_Y     = (Nf_Y-0).*(1+sinh(Ks_Y.*(Nf_Unf_Y(j)-B_Y))./sinh(Ks_Y.*B_Y))+0;
            
            % >> Domain vertices.
            %  > (Xd,Yd).
            [Xd_p,Yd_p] = ...
                deal(NF_X(1:NX_v).*(Xv_f-Xv_i),NF_Y(1:NY_v).*(Yv_f-Yv_i));
            [Xd,Yd] = meshgrid(Xd_p,Yd_p);
            %  > (Xv,Yv).
            xy_v(:,1) = reshape(Xd,[],1);
            xy_v(:,2) = reshape(Yd,[],1);
        end
    end
end