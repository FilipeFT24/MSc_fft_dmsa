classdef B_1D
    methods(Static)
        %% > Wrap-up B (1D).
        function [msh,pde] = WrapUp_B_1D(inp,msh)
            %  > Initialize.
            [pde,s,stl] = B_1D.Initialize(inp,msh);
            
            switch inp.ee.test
                case false
                    %  > 'Standard' and 'p-adaptative' runs w/ analytic values.
                    [msh,pde] = B_1D.Run_p(inp,msh,pde,s,stl);
                case true
                    %  > Check error estimators.
                    [msh,pde] = B_1D.SetUp_TX(inp,msh,pde,s);
                otherwise
                    return;
            end
        end
        
        %% > Auxiliary functions.
        % >> 1. -----------------------------------------------------------
        %  > Initialize problem (for all tests).
        function [pde,s,stl] = Initialize(inp,msh)
            %  > 'pde'.
            pde = B_1_1D.Update_pde(inp,msh);
            %  > A/B.
            i = 1:2;
            for j = i
                A{j} = zeros(msh.c.NC);
                B{j} = zeros(msh.c.NC,1);
            end
            %  > 's'.
            s.A     = A;
            s.B     = B;
            s.c     = cell(2,msh.f.NF);
            s.f     = cell(2,msh.f.NF);
            s.bnd_i = cell(2,msh.f.NF);
            s.bnd_v = cell(2,msh.f.NF);
            s.xf    = cell(2,msh.f.NF);
            s.xt    = cell(2,msh.f.NF);
            s.Inv   = cell(2,msh.f.NF);
            s.vg    = [inp.pv.v(1),-inp.pv.v(2)];
            %  > 'stl'.
            switch inp.ee.test
                case false
                    s.p   = A_2_1D.Compute_p(inp.ps.p,inp.ps.s);
                    stl   = A_2_1D.Initialize_stl(msh,inp.ps.p,inp.ps.s);
                case true
                    s.p   = A_2_1D.Compute_p(inp.ee.p,inp.ee.s);
                    stl.p = cell(1,2);
                    stl.s = cell(1,2);
                    stl.t = cell(1,2);
                otherwise
                    return;
            end
        end
        % >> 2. -----------------------------------------------------------
        %  > Set up 'p-standard' and 'p-adaptative' runs.
        function [msh,pde] = Run_p(inp,msh,pde,s,stl)
            switch inp.pa.adapt
                case false
                    %  > 'p-standard' run.
                    [msh,pde] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                    %  > Plot...
                    if inp.pl.all
                        Fig_1_1D.Plot(msh,pde);
                        Fig_2_1D.Plot(msh,pde,inp.pl.tt);
                    end
                case true
                    %  > 'p-adaptative' run.
                    [msh,pde] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                    %  > While...
                    i = 0;
                    s = msh.s;
                    while 1
                        %  > Adapt...
                        if i ~= inp.pa.n
                            stl       = B_2_2_1D.Adapt_p(inp,stl,pde.e);
                            [msh,pde] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                        else
                            break;
                        end
                        i = i+1;
                    end
                    %  > Plot...
                    if inp.pl.all
                        Fig_1_1D.Plot(msh,pde);
                    end
                otherwise
                    return;
            end
        end
        % >> 3. -----------------------------------------------------------
        %  > Set up other tests.
        function [msh,pde] = SetUp_TX(inp,msh,pde,s)
            if inp.ee.flag(1)
                pde = B_2_3_1D.T1(inp,msh,pde,s);
            end
            if inp.ee.flag(2)
                pde = B_2_3_1D.T2(inp,msh,pde,s);
            end
        end
    end
end