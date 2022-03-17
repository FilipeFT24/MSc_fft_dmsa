classdef B_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        %  > Initialize problem and set up 'standard'/'p-adaptative' runs.
        function [msh,pde] = Set_B(inp,msh)
            [pde,s,stl] = B_1D.Initialize(inp,msh);
            [msh,pde]   = B_1D.Run_p(inp,msh,pde,s,stl);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Initialize problem.
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
            stl     = A_2_1D.Initialize_stl(msh,inp.ps.p,inp.ps.t);
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
                        Fig_2_1D.Plot(msh,pde,inp.pl.tt);
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