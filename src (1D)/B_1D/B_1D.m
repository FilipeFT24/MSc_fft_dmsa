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
        %    Initialize problem (for all tests).
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
            %  > 'stl'.
            switch inp.ee.test
                case false
                    stl   = A_2_1D.Initialize_stl(msh,inp.ps.p,inp.ps.s);
                case true
                    stl.p = cell(1,2);
                    stl.s = cell(1,2);
                    stl.t = cell(1,2);
                otherwise
                    return;
            end
        end
        % >> 2. -----------------------------------------------------------
        %    Set up 'p-standard' and 'p-adaptative' runs.
        function [msh,pde] = Run_p(inp,msh,pde,s,stl)
            switch inp.pa.adapt
                case false
                    %  > 'p-standard' run.
                    [pde,s] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                case true
                    %  > 'p-adaptative' run.
                    [pde,s] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                    %  > While...
                    i = 0;
                    while 1
                        %  > Adapt...
                        if i ~= inp.pa.n
                            [stl]   = B_2_2_1D.Adapt_p(inp,stl,pde.e);
                            [pde,s] = B_2_1_1D.SetUp_p(inp,msh,pde,s,stl);
                        else
                            break;
                        end
                        i = i+1;
                    end
                otherwise
                    return;
            end
            msh = Tools_1D.Set_msh(msh,stl,s);
            %  > Plot...
            Fig_1_1D.Plot(msh,pde);
        end
        % >> 3. -----------------------------------------------------------
        %    Set up error estimators run.
        function [msh,pde] = SetUp_TX(inp,msh,pde,s)
            % >> 1.
            %    Truncated terms' magnitude (w/ analytic derivatives).
            if inp.ee.flag(1)
                [msh,pde] = B_2_3_1D.T_1(inp,msh,pde,s,[5,5]); 
            end
            %  > 2.
            if inp.ee.flag(2)
                %  > Set schemes order/type (v/g).
                n        = 2;
                inp_lo.p = [1,1];     % > LO.
                inp_ho.p = [3,3];     % > HO.
                inp_lo.s = ["c","c"]; % > LO.
                inp_ho.s = ["c","c"]; % > HO.
                for i = 1:n
                    p_lo(i) = A_2_1D.Compute_p(inp_lo.p(i),inp_lo.s(i));
                    p_ho(i) = A_2_1D.Compute_p(inp_ho.p(i),inp_ho.s(i));
                end
                inp_lo.o  = p_lo;
                inp_lo.nt = p_ho-p_lo;
                %  > Compute runcated terms' magnitude (w/ higher-order solution).
                B_2_3_1D.EE_2(inp,inp_lo,inp_ho,msh,pde,s);
            end
            %  > 3.
            if inp.ee.flag(3)
                pde = B_2_3_1D.T_3(inp,msh,pde,s);
            end
        end
    end
end