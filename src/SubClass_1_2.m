classdef SubClass_1_2
    methods (Static)
        %% > Wrap up SubClass_1_2.
        function [msh] = WrapUp_1_2(inp,msh)
            %  ------------------------------------------------------------
            % >> 1. Wrap up mesh (determine all necessary components).
            %  ------------------------------------------------------------
            
            %  > 1.
            msh = SubClass_1_2.WrapUp_msh(inp,msh);
        end
        
        %% > Tools.
        %% > 1.) ----------------------------------------------------------
        function [msh] = WrapUp_msh(inp,msh)
            %  ------------------------------------------------------------
            % >> 1.   Generate mesh vertices/vertex connectivity (i.e. the indices of the vertices that compose a given cell).
            %  > 1.1. '^' -> Face polygon: triangle -> Delaunay triangulation.
            %  > 1.2. 's' -> Face polygon: square   -> Cartesian.
            % >> 2.   Determine mesh properties.
            %  > 2.1. Wrap up cell (cell components).
            %  > 2.2. Wrap up face (face components). 
            %         Remark: 'WrapUp_Cell' is called first since some of the components of this field help setting up 'WrapUp_Face'.
            %  > 2.3. Finalize... 
            %  ------------------------------------------------------------
            % >> Local variables.
            NX_v  = inp.msh.Nv(1);
            NY_v  = inp.msh.Nv(2);
            NT    = inp.fr.nt;
            Order = inp.fr.n;
            T1    = inp.msh.T_1.t;
                        
            %  > 1.
            if strcmpi(T1,'^')
                struct = delaunayTriangulation(msh.d.XY_v(:,1),msh.d.XY_v(:,2));               
            elseif strcmpi(T1,'s')
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
                struct.Points(:,1) = msh.d.XY_v(:,1);
                struct.Points(:,2) = msh.d.XY_v(:,2);
            end
            %  ------------------------------------------------------------
            %  > 3.2. 
            %  > Number of cells.
            msh.c.NC = size(struct.ConnectivityList,1);
            %  > Cell vertices coordinates.
            for i = 1:msh.c.NC
                msh.c.XY_v{i}(:,1) = struct.Points(struct.ConnectivityList(i,:),1);
                msh.c.XY_v{i}(:,2) = struct.Points(struct.ConnectivityList(i,:),2);
            end
            msh = SubClass_1_2.WrapUp_Cell(struct,msh);
            msh = SubClass_1_2.WrapUp_Face(struct,msh,NT,Order);
            %  ------------------------------------------------------------
            %  > 3.3.
            msh = SubClass_1_2.WrapUp_Finalize(msh);
        end       
        % >> 3.2.1.) ------------------------------------------------------
        function [msh] = WrapUp_Cell(struct,msh)
            for i = 1:msh.c.NC
                %  > Volume.
                msh.c.Vol   (i) = polyarea(msh.c.XY_v{i}(:,1),msh.c.XY_v{i}(:,2));
                %  > Centroid.
                msh.c.mean(1,i) = mean(msh.c.XY_v{i}(:,1));
                msh.c.mean(2,i) = mean(msh.c.XY_v{i}(:,2));
            end
            %  > Cell neighbours.
            msh = SubClass_1_3.Set_CellNeighbours(struct,msh);
        end
        % >> 3.2.2.) ------------------------------------------------------
        function [msh] = WrapUp_Face(struct,msh,NT,Order)
            %  > Faces' coordinates.
            msh = SubClass_1_3.Set_DomainFaces(struct,msh);

            for j = 1:size(msh.f.XY_v,2)
                %  > Centroid.
                msh.f.mean(1,j) = mean(msh.f.XY_v{j}(:,1));
                msh.f.mean(2,j) = mean(msh.f.XY_v{j}(:,2));
                %  > Length.
                msh.f.len   (j) = pdist2(msh.f.XY_v{j}(1,:),msh.f.XY_v{j}(2,:)); 
            end
            %  > Normals.
            msh = SubClass_1_3.Set_FaceNormals(msh);
            % >> Stencil.
            %  > Stencil neighbours.
            msh = SubClass_1_3.Set_FaceNeighbours(msh,NT,Order);
            %  > Stencil limits.
            msh = SubClass_1_3.Set_Limits(msh);
        end        
        % >> 3.2.3.) ------------------------------------------------------
        function [msh] = WrapUp_Finalize(msh)
            %  > Reference length.
            msh = SubClass_1_3.Set_ReferenceLength(msh);
        end
    end
end