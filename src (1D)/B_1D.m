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
            allow_odd = inp.fr.allow_odd;
            n         = inp.fr.n;
            ee        = inp.fr.test_ee;
            t1(1)     = inp.fr.type_1.v;
            t1(2)     = inp.fr.type_1.g;
            t2(1)     = inp.fr.type_2.v;
            t2(2)     = inp.fr.type_2.g;
            
            [msh,pde] = B_1D.SetUp(msh,v,g,"sin",bnd,p_adapt,allow_odd,n,ee,t1,t2);
        end
        
        %% > Auxiliary functions.
        % >> 1. -----------------------------------------------------------
        function [pde,A,B,stl,s] = Initialize(msh,v,g,ft,t1,t2)
            A             = zeros(msh.c.NC);
            B             = zeros(msh.c.NC,1);
            pde           = B_1_1D.WrapUp_B_1_1D(msh,v,g,ft);
            B             = B+pde.f.st;
            stl.p{1}(:,1) = repelem(t2(1),msh.f.NF);
            stl.p{2}(:,1) = repelem(t2(2),msh.f.NF);
            stl.s{1}(:,1) = 1:msh.f.NF;
            stl.s{2}(:,1) = 1:msh.f.NF;
            stl.t{1}(:,1) = repelem(t1(1),msh.f.NF);
            stl.t{2}(:,1) = repelem(t1(2),msh.f.NF);
            s.c           = cell (2,msh.f.NF);
            s.f           = cell (2,msh.f.NF);
            s.xt          = cell (2,msh.f.NF);
            s.Ls          = zeros(2,msh.f.NF);
            s.bnd         = cell (2,msh.f.NF);
            s.Tf          = cell (2,msh.f.NF);
            s.xf          = cell (2,msh.f.NF);
        end
        % >> 2. -----------------------------------------------------------
        function [msh,pde] = SetUp(msh,v,g,ft,bnd,p_adapt,ao,n,ee,tv,tg)
            % >> Select...
            switch ee
                case false
                    % >> Initialize.
                    [pde,A,B,stl,s] = B_1D.Initialize(msh,v,g,ft,tv,tg);
                    % >> Set up problem...
                    switch p_adapt
                        case false
                            %  > ...solve PDE.
                            [stl,s,x,ea] = B_2_1_1D.Update_stl(msh,stl,s,A,B,pde.a,bnd,v,g);
                            [e,x]        = B_2_1_1D.Update_pde(msh,pde.a,s,x,ea,v,g);
                        case true
                            %  > ...solve PDE while...
                            i = 0;
                            while i < n
                                [stl,s,x,ea] = B_2_1_1D.Update_stl(msh,stl,s,A,B,pde.a,bnd,v,g);
                                [e,x]        = B_2_1_1D.Update_pde(msh,pde.a,s,x,ea,v,g);
                                if i+1 ~= n
                                    [stl] = B_2_2_1D.Select_fCD(ao,stl,e.c.c,e.t);
                                end
                                i = i+1;
                            end
                        otherwise
                            return;
                    end
                    %  > ...update structures.
                    [msh,pde] = B_1D.Update(msh,pde,stl,s,x,e);
                case true
                otherwise
                    return;
            end
        end
        % >> 3. -----------------------------------------------------------
        function [msh,pde] = Update(msh,pde,stl,s,xn,e)
            s.stl = stl;
            msh.s = s;
            msh   = Tools_1D.Sort_struct(msh);
            pde.x = xn;
            pde.e = e;
        end
    end
end