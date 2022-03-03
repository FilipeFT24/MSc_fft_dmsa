classdef B_1D
    methods(Static)
        %% > Wrap-up B (1D).
        function [msh,pde] = WrapUp_B_1D(inp,msh)
            %  > Initialize.
            [obj,pde,s,stl] = B_1D.Initialize(inp,msh);
            
            switch obj.ee
                case false
                    %  > 'Standard' and 'p-adaptative' runs w/ analytic values.
                    [msh,pde] = B_1D.Run_p(obj,msh,pde,s,stl);
                case true
                    %  > Check error estimators.
                    B_2_3_1D.SetUp_ee(obj,msh,pde,s);
                otherwise
                    return;
            end
        end
        
        %% > Auxiliary functions.
        % >> 0. -----------------------------------------------------------
        function [obj] = SetUp_obj(inp)
            %  > Auxiliary variables.
            bnd_w       = inp.pr.w;
            bnd_e       = inp.pr.e;
            t1(1)       = inp.fr.type_1.v;
            t1(2)       = inp.fr.type_1.g;
            t2(1)       = inp.fr.type_2.v;
            t2(2)       = inp.fr.type_2.g;
            
            %  > 'obj' structure.
            obj.v       = inp.pr.v;
            obj.g       = inp.pr.g;
            obj.bnd     = [string(bnd_w),string(bnd_e)];
            obj.ft      = inp.pr.ft;
            obj.p_adapt = inp.fr.p_adapt;
            obj.ao      = inp.fr.allow_odd;
            obj.n       = inp.fr.n;
            obj.gf      = inp.fr.gf;
            obj.t1      = t1;
            obj.t2      = t2;
            obj.ee      = inp.fr.test_ee;
        end
        % >> 1. -----------------------------------------------------------
        %  > Initialize problem (for all tests).
        function [obj,pde,s,stl] = Initialize(inp,msh)
            %  > obj.
            obj = B_1D.SetUp_obj(inp);
            %  > pde.
            pde = B_1_1D.WrapUp_B_1_1D(msh,obj.v,obj.g,obj.ft);
            %  > A/B.
            i = 1:2;
            for j = i
                A{j} = zeros(msh.c.NC);
                B{j} = zeros(msh.c.NC,1);
            end
            %  > s.
            s.A     = A;
            s.B     = B;
            s.c     = cell(2,msh.f.NF);
            s.f     = cell(2,msh.f.NF);
            s.bnd_i = cell(2,msh.f.NF);
            s.bnd_v = cell(2,msh.f.NF);
            s.xf    = cell(2,msh.f.NF);
            s.xt    = cell(2,msh.f.NF);
            %  > stl.
            switch obj.ee
                case false
                    stl   = A_2_1D.Initialize_stl(msh,obj.t1,obj.t2);
                case true
                    stl.p = cell(1,2);
                    stl.s = cell(1,2);
                    stl.t = cell(1,2);
                otherwise
                    return;
            end
        end
        % >> 2. -----------------------------------------------------------
        %  > 2.1 ----------------------------------------------------------
        %  > Set up 'p-standard' and 'p-adaptative' runs.
        function [msh,pde] = Run_p(obj,msh,pde,s,stl)
            switch obj.p_adapt
                case false
                    %  > 'p-standard' run.
                    [pde,s,stl] = B_2_1_1D.SetUp_p(obj,msh,pde,s,stl);
                case true
                    %  > 'p-adaptative' run.
                    [pde,s,stl] = B_2_1_1D.SetUp_p(obj,msh,pde,s,stl);
                    %  > While...
                    i = 0;
                    while 1
                        %  > Adapt...
                        if i ~= obj.n
                            [stl]       = B_2_2_1D.Adapt_p(obj,stl,pde.e);
                            [pde,s,stl] = B_2_1_1D.SetUp_p(obj,msh,pde,s,stl);
                        else
                            break;
                        end
                        i = i+1;
                    end
                otherwise
                    return;
            end
            %  > Update structures.
            [msh,pde] = B_1D.Update(msh,pde,s,stl);
            %  > Plot...
            Fig_1_1D.WrapUp_Fig_1_1_1D(1,msh,pde);
            Fig_1_1D.WrapUp_Fig_1_2_1D(0,msh,pde);
        end
        % >> 3. -----------------------------------------------------------
        %  > Update 'msh' and 'pde' structures.
        function [msh,pde] = Update(msh,pde,s,stl)
            s.stl = stl;
            msh.s = s;
            msh   = Tools_1D.Sort_struct(msh);
        end
    end
end