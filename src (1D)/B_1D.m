classdef B_1D
    methods(Static)
        %% > Wrap-up B (1D).
        function [msh,pde] = WrapUp_B_1D(inp,msh)
            v         = inp.pr.v;
            g         = inp.pr.g;
            bnd_w     = inp.pr.w;
            bnd_e     = inp.pr.e;
            bnd       = [string(bnd_w),string(bnd_e)];
            p_adapt   = inp.fr.p_adapt;
            t1(1)     = inp.fr.type_1.v;
            t1(2)     = inp.fr.type_1.g;
            t2(1)     = inp.fr.type_2.v;
            t2(2)     = inp.fr.type_2.g;

            [msh,pde] = B_1D.SetUp(msh,v,g,"exp",bnd,p_adapt,t1,t2);
        end
        
        %% > Auxiliary functions.
        % >> 1. -----------------------------------------------------------
        function [pde,A,B,stl,s] = Initialize(msh,v,g,ft,t1,t2)
            A          = zeros(msh.c.NC);
            B          = zeros(msh.c.NC,1);
            pde        = B_1_1D.WrapUp_B_1_1D(msh,v,g,ft);
            B          = B+pde.FV;
            stl.p(1,:) = repelem(t2(1),msh.f.NF);
            stl.p(2,:) = repelem(t2(2),msh.f.NF);
            stl.s(1,:) = 1:msh.f.NF;
            stl.s(2,:) = 1:msh.f.NF;
            stl.t(1,:) = repelem(t1(1),msh.f.NF);
            stl.t(2,:) = repelem(t1(2),msh.f.NF);
            s.c        = cell (2,msh.f.NF);
            s.f        = cell (2,msh.f.NF);
            s.xt       = cell (2,msh.f.NF);
            s.Ls       = zeros(2,msh.f.NF);
            s.bnd      = cell (2,msh.f.NF);
            s.Tf       = cell (2,msh.f.NF);
            s.xf       = cell (2,msh.f.NF);
        end
        % >> 2. -----------------------------------------------------------
        function [msh,pde] = SetUp(msh,v,g,ft,bnd,p_adapt,tv,tg)
            % >> Initialize.
            [pde,A,B,stl,s] = B_1D.Initialize(msh,v,g,ft,tv,tg);
            % >> Set up problem...
            switch p_adapt
                case false
                    %  > ...solve PDE.
                    [s,xn] = B_2_1D.Update_stl(msh,s,stl.p,stl.s,stl.t,A,B,pde.sn,bnd,v,g);
                    [e,xn] = B_2_1D.Update_pde(msh,xn,pde.sn,s);
                case true
                    %  > ...solve PDE while...
                    i = 1;
                    while i < 5
                        [s,xn] = B_2_1D.Update_stl(msh,s,stl.p,stl.s,stl.t,A,B,pde.sn,bnd,v,g);
                        [e,xn] = B_2_1D.Update_pde(msh,xn,pde.sn,s);
                        [stl]  = B_2_1D.Select_f(stl,e);
                        i      = i+1;
                    end
                    s.stl = stl;
                otherwise
                    return;
            end
            %  > ...update structures.
            [msh,pde] = B_1D.Update(msh,pde,s,xn,e);
        end
        % >> 3. -----------------------------------------------------------
        function [msh,pde] = Update(msh,pde,s,xn,e)
            msh.s = s;
            msh   = Tools_1D.Sort_struct(msh);
            pde.x = xn; 
            pde.e = e; 
        end
    end
end