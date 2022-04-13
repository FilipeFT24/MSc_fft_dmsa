classdef Fig_V1_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(0,exp);
            x.a = 1;  %  > Convection(1)/Diffusion(2).
            x.b = 6; %  > Face index(f).
            
%             
%             face   = msh.f.Xv(x.b);
%             xc_loc = msh.c.Xc(obj.s.c{x.b,x.a});
%             xc_fac = [0.09,0.11];
%             cell   = [1,2,5,6];
%             fac    = [3,4];
%             for i  = 1:6
%                 D(cell,i) = (xc_loc-face).^(i-1);
%             end
%             for i  = 1:6
%                 if i == 1
%                     D(fac ,i) = 0;
%                 else
%                     D(fac ,i) = (i-1).*(xc_fac-face).^(i-2);
%                 end
%             end
%             phi_s(cell,1) = obj.x{1}.nv.x.c(obj.s.c{x.b,x.a}); 
%             val = [0.0806386944967543;0.0768327862357099];
%             phi_s(fac ,1) = val;
%            
%             
%             
%             
%             
%             
%             
% %             Ds = D([1,3,4,6],[1,2,3,4]);
% %             phis_s = phi_s([1,3,4,6]);
% %             
%             %coeffs{1} = 0;%inv(Ds)*phis_s;
            coeffs = 0;%inv(D)*phi_s;
            
            
            
            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V1_2_1D.Plot_1_1(msh.c.Xc,msh.f.Xv,obj,x,fig,0,coeffs);
                subplot(1,2,2);
                Fig_V1_2_1D.Plot_1_2(msh.c.Xc,msh.f.Xv,obj,x,fig,0,coeffs);
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. ------------------------------------------------------- 
        function [func] = Tools_1(obj,a,b,coeffs)
            for i = 1:size(obj.x,2)
                %  > func_f(A).
                c         = length(obj.x{i}.cf.a{b,a});
                x         = Fig_V1_2_1D.Set_var(a,c);
                func.a{i} = obj.f.ps{a}(x.d)*obj.x{i}.cf.a{b,a}(x.e);
                %  > func_f(p).
                c         = length(obj.x{i}.cf.x{b,a});
                x         = Fig_V1_2_1D.Set_var(a,c);
                
                %if i == 2
                %    func.f{i} = obj.f.ps{a}(x.d)*coeffs(x.e);
                %else
                   %func.f{i} = obj.f.ps{a}(x.d)*coeffs{1}(x.e);
                   func.f{i} = obj.f.ps{a}(x.d)*obj.x{i}.cf.x{b,a}(x.e);
                %end
                
            end
        end
        %  > 2.1.2. -------------------------------------------------------
        function [x] = Set_var(a,c)
            switch a
                case 1
                    d = 1:c;
                    e = d;
                case 2
                    d = 1:c-1;
                    e = d+1;
                otherwise
                    return;
            end
            x.d = d;
            x.e = e;
        end
        % >> 2.2. ---------------------------------------------------------
        %  > \tau_f(x).
        function [V] = Tools_2(obj,a,b,f,y,coeffs)
            func  = Fig_V1_2_1D.Tools_1(obj,a,b,coeffs); 
            for i = 1:size(obj.x,2)-1
                %  > d\tau_f(A).
                x(1)       = func.a{i+1}-func.a{i};
                tau_f.a{i} = matlabFunction(abs(x(1)));
                %  > d\tau_f(p).
                x(2)       = func.f{i+1}-func.f{i};
                tau_f.p{i} = matlabFunction(abs(x(2)));
                %  > d\tau_f(d).
                x(3)       = func.a{i+1}-func.f{i+1};
                x(4)       = func.a{i}  -func.f{i};
                tau_f.d{i} = matlabFunction(abs(x(3)-x(4)));
            end
            V = Fig_V1_2_1D.Tools_4([tau_f.a,tau_f.p,tau_f.d],f,y);
        end
        % >> 2.3. ---------------------------------------------------------
        %  > dfunc_f(x).
        function [V] = Tools_3(obj,a,b,f,y,coeffs)
            func  = Fig_V1_2_1D.Tools_1(obj,a,b,coeffs);
            for i = 1:size(obj.x,2)
                %  > f-func_f(p).
                %d_f{i} = matlabFunction(abs(diff(diff(diff(diff(func.f{i+1}))))));
                %d_f{i+1} = matlabFunction(abs(diff(diff(diff(diff(func.a{i+1}))))));
                d_f{i} = matlabFunction(abs(obj.f.fh{a}-func.f{i}));
            end
            V = Fig_V1_2_1D.Tools_4(d_f,f,y);
        end
        % >> 2.4. ---------------------------------------------------------
        function [V] = Tools_4(v,f,y)
            syms x;
            for i = 1:size(v,2)
                switch nargin(v{i})
                    case 0
                        V{1,i} = @(x)(feval(v{i})).*x.^(0);
                        V{2,i} = repelem(feval(v{i}),length(y))';
                    case 1
                        V{1,i} = v{i}(x);
                        V{2,i} = v{i}(y);
                    case 2
                        V{1,i} = v{i}(f,x);
                        V{2,i} = v{i}(f,y);
                    otherwise
                        return;
                end
            end
        end
 
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        function [] = Plot_1_1(Xc,Xv,obj,x,fig,zoom,coeffs)
            %  > Auxiliary variables.
            fid = "V1_2_1D (1)";
            a   = x.a;
            b   = x.b;
            c   = Xv(b);
            X1  = sort([Xc(b-4:b+3);c]);
            
            %  > Select variables.
            M1(1,:)    = ["-",":","--"];
            M1(2,:)    = ["o","v","s"];
            L1         = Fig_V1_2_1D.Set_Legend_1(a);
            V1         = Fig_V1_2_1D.Tools_2(obj,a,b,c,X1,coeffs);
            %  > Plot variables.
            [L1,P1,Y1] = Fig_Tools_1D.Var_2(fig,M1,L1,X1,V1);
            f          = xline([c,c],'Color',fig.C(4,:),'Linewidth',fig.LW./10,'Linestyle',"-");
            
            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L1,P1,0,X1,Y1,[-1,1],1);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 3.1.2. -------------------------------------------------------
        function [L] = Set_Legend_1(a)
            S    = Fig_Tools_1D.Set_str_1(a);
            L{1} = join(["$|\bar{\tau}_{f^{\left(a\right)}}^{",S,"}\left(x\right)|$"]);
            L{2} = join(["$|\bar{\tau}_{f^{\left(p\right)}}^{",S,"}\left(x\right)|$"]);
            L{3} = join(["$|\bar{\tau}_{f^{\left(d\right)}}^{",S,"}\left(x\right)|$"]);
        end 
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        function [] = Plot_1_2(Xc,Xv,obj,x,fig,zoom,coeffs)
            %  > Auxiliary variables.
            fid = "V1_2_1D (2)";
            a   = x.a;
            b   = x.b;
            c   = Xv(b);
            X1  = sort([Xc(b-4:b+3);c]);
            
            %  > Select variables.
            M1(1,:)    = ["-",":"];
            M1(2,:)    = ["o","v"];
            L1         = Fig_V1_2_1D.Set_Legend_2(a);
            V1         = Fig_V1_2_1D.Tools_3(obj,a,b,c,X1,coeffs);
            %  > Plot variables.
            [L1,P1,Y1] = Fig_Tools_1D.Var_2(fig,M1,L1,X1,V1);
            f          = xline([c,c],'Color',fig.C(4,:),'Linewidth',fig.LW./10,'Linestyle',"-");
            
            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L1,P1,0,X1,Y1,[-5,1],1);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 3.2.2. -------------------------------------------------------
        function [L] = Set_Legend_2(a)
            S    = Fig_Tools_1D.Set_str_2(a);
            L{1} = join(["$|",S,"_{f^{\left(a\right)}}^{\left(p^{-}\right)}\left(x\right)-",S,"_{f^{\left(p\right)}}^{\left(p^{-}\right)}\left(x\right)|$"]);
            L{2} = join(["$|",S,"_{f^{\left(a\right)}}^{\left(p^{+}\right)}\left(x\right)-",S,"_{f^{\left(p\right)}}^{\left(p^{+}\right)}\left(x\right)|$"]);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Create function handle of \phi and \nabla\phi with 25 terms (polynomial shape).
        function [func] = pol_shape()
            syms x f;
            n       = 25;
            l_1     = 1:n;
            func{1} = (x-f).^(l_1-1);
            l_2     = 1:n-1;
            func{2} = l_2.*(x-f).^(l_2-1);
        end
    end    
end