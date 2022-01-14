classdef SubClass_1_3
    methods (Static)
        %% > Wrap-up SubClass_1_3.
        % >> --------------------------------------------------------------
        % >> 1.   Set 'struct' structure: generate mesh vertex connectivity.
        %  > 1.1. 'v': Face polygon: triangle.
        %  > 1.2. 's': Face polygon: square.
        % >> 2.   Determine grid properties.
        %  > 2.1. Wrap up cell (cell components).
        %  > 2.2. Wrap up face (face components).
        %         Remark: 'WrapUp_Cell' is called first since some of the components of this field help setting up 'WrapUp_Face'.
        % >> --------------------------------------------------------------
        function [msh] = WrapUp_1_3(inp,msh)
            % >> 1.
            [struct,msh] = SubClass_1_3.Set_struct(inp,msh);
            % >> 2.
            %  > 2.1.
            msh = SubClass_1_3.WrapUp_Cell(struct,msh);
            %  > 2.2.
            msh = SubClass_1_3.WrapUp_Face(inp,msh,struct);
        end
        
        %% > 1. -----------------------------------------------------------
        function [struct,msh] = Set_struct(inp,msh)
            % >> Local variables.
            NX_v = inp.msh.Nv(1);
            NY_v = inp.msh.Nv(2);
            T1   = inp.msh.T_1.t;
                        
            if strcmpi(T1,'v')
                % >> Triangles.
                struct = delaunayTriangulation(msh.d.xy_v(:,1),msh.d.xy_v(:,2));
            elseif strcmpi(T1,'s')
                % >> Squares.
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
                struct.Points(:,1) = msh.d.xy_v(:,1);
                struct.Points(:,2) = msh.d.xy_v(:,2);
            end
            %  > Number of cells.
            msh.c.NC = size(struct.ConnectivityList,1);
            %  > Cell vertex coordinates.
            for i = 1:msh.c.NC
                msh.c.xy_v{i}(:,1) = struct.Points(struct.ConnectivityList(i,:),1);
                msh.c.xy_v{i}(:,2) = struct.Points(struct.ConnectivityList(i,:),2);
            end
        end      
               
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [msh] = WrapUp_Cell(struct,msh)
            for i = 1:msh.c.NC
                %  > Cell volume.
                msh.c.vol   (i) = polyarea(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2));
                %  > Cell centroid.
                msh.c.mean(1,i) = mean(msh.c.xy_v{i}(:,1));
                msh.c.mean(2,i) = mean(msh.c.xy_v{i}(:,2));
            end
            %  > Cell neighbouring cells.
            msh = SubClass_1_3_1.Set_CellNeighbours(struct,msh);
        end
        % >> 2.2. ---------------------------------------------------------
        function [msh] = WrapUp_Face(inp,msh,struct)
            % >> Local variables.
            NT    = inp.fr.nt;
            Order = inp.fr.np;
            
            %  > Face coordinates.
            msh = SubClass_1_3_1.Set_DomainFaces(struct,msh);
            %  > Cell reference length.
            msh = SubClass_1_3_1.Set_ReferenceLength(msh);
            
            for j = 1:size(msh.f.xy_v,2)
                %  > Face centroid.
                msh.f.mean(1,j) = mean(msh.f.xy_v{j}(:,1));
                msh.f.mean(2,j) = mean(msh.f.xy_v{j}(:,2));
                %  > Face length.
                msh.f.len   (j) = pdist2(msh.f.xy_v{j}(1,:),msh.f.xy_v{j}(2,:));
            end
            %  > Face normals.
            msh = SubClass_1_3_1.Set_FaceNormals(msh);
            % >> Stencil.
            %  > Neighbours.
            msh = SubClass_1_3_2.Set_FaceNeighbours(msh,NT,Order);
        end
    end
end