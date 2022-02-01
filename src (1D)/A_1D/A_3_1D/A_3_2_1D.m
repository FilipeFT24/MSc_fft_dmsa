classdef A_3_2_1D
    methods (Static)
        %% > Wrap-up A_3_2 (1D).
        function [msh] = WrapUp_A_3_2_1D(msh,np)
            %  > Auxiliary variables.
            Xc = msh.c.Xc;
            NC = msh.c.NC;
            Xf = msh.f.Xv;
            NF = msh.f.NF;
            
            %  > Set up stencil for face 'i'...
            for i = 1:msh.f.NF
                %  > Cell/face indices.
                [msh.s.c{i},msh.s.f{i}] = A_3_2_1D.SetUp_Stencil_f(i,np,Xc,NC,Xf,NF);
                %  > Cell/face coordinates.
                [msh.s.x_v_c{i},msh.s.x_v_f{i}] = ...
                    A_3_2_1D.Compute_Coordinates_cf(Xc,Xf,msh.s.c{i},msh.s.f{i});
                msh.s.x_v_t {i} = ...
                    A_3_2_1D.Compute_Coordinates_tt(msh.s.x_v_c{i},msh.s.x_v_f{i});
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [sc,sf] = SetUp_Stencil_f(i,np,Xc,NC,Xf,NF)
            %  > Add faces...
            np_2 =  np./2;
            if i <= np_2
                sf = 1;
            elseif i >= NF-np_2+1
                sf = NF;
            else
                sf = [];
            end
            %  > Add cells...
            
            df     = Xf(i)-Xc;
            [~,ic] = min(abs(df));
            vc     = Xc(ic);
            
            if i <= np_2
                sc = 1:1:np-1;
            elseif i >= NF-np_2+1
                sc = NC-np+2:1:NC;
            else
                sc = i-np_2:1:i+np_2-1;
            end
        end
        % >> 1.2. ---------------------------------------------------------
        function [xc,xf] = Compute_Coordinates_cf(Xc,Xf,sc,sf)
            xc = Xc(sc);
            xf = Xf(sf);
        end
        %  > 1.3. ---------------------------------------------------------
        function [xt] = Compute_Coordinates_tt(st_c,st_f)
            xt = sort([st_c,st_f]);
        end
    end
end