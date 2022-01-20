classdef A_3
    methods (Static)
        %% > Wrap-up A_3.
        % >> --------------------------------------------------------------
        % >> 1. Set 'struct' structure.
        % >> 2. Determine grid properties.
        % >> 3. Set up stencil.
        % >> --------------------------------------------------------------
        function [msh] = WrapUp_A_3(inp,msh)
            % >> Local variables.
            nt   = inp.fr.nt;
            p    = inp.fr.np;
            nl   = 1./2.*(p+1);  
            et_1 = inp.fr.et_1;
            et_2 = inp.fr.et_2;
            
            % >> 1.
            [struct,msh] = A_3.Set_struct(inp,msh);
            % >> 2.
            msh = A_3_1.WrapUp_A_3_1(struct,msh);
            % >> 3.
            msh = A_3_2.WrapUp_A_3_2(msh,nt,nl,p,et_1,et_2);
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
    end
end