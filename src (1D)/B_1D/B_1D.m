classdef B_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        function [msh,pde] = Set_B(inp,msh)
            %  > Initialize...
            [pde,s,stl] = B_1D.Initialize(inp,msh);
            %  > Set up 'standard' and 'p-adaptative' runs.
            [msh,pde]   = B_1D.Run_p(inp,msh,pde,s,stl);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
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
                    stl   = A_2_1D.Initialize_stl(msh,inp.ps.p,inp.ps.t);
                case true
                    s.p   = A_2_1D.Compute_p(inp.ee.p,inp.ee.t);
                    stl.p = cell(1,2);
                    stl.s = cell(1,2);
                    stl.t = cell(1,2);
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Set up 'p-standard' and 'p-adaptative' runs.
        function [msh,pde] = Run_p(inp,msh,pde,s,stl)
            switch inp.pa.adapt
                case false
                    %  > 'p-standard' run.
                    [msh,pde] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                    %  > Plot...
                    if inp.pl.all
                        Fig_1_1D.Plot(msh,pde);
                        Fig_2_1D.Plot(inp,msh,pde);
                    end
                case true
                    %  > 'p-adaptative' run.
                    [msh,pde] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                    %  > While...
                    i = 0;
                    while 1
                        %  > Adapt...
                        if i ~= inp.pa.n
                            [pde,stl] = B_2_2_1D.Adapt_p(inp,pde,stl);
                            [msh,pde] = B_2_1_1D.SetUp_p(inp,msh,pde,msh.s,stl);
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
    end
end