%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format short; beep off; warning('off');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories();
%  > Run...
run_1 = true;
% >> ----------------------------------------------------------------------
if run_1
    [inp,msh] = A_1D.WrapUp_A_1D(0.0025);
    [msh,pde] = B_1D.WrapUp_B_1D(inp,msh);
end

% for i = 1:2
%     for j = 1:msh.f.NF
%         p(j,i) = A_2_1D.Compute_p(msh.s.stl.p{i}(j),msh.s.stl.t{i}(j));
%     end
% end
% figure(100);
% hold on;
% plot(msh.f.Xv,p(:,1),'-or');
% plot(msh.f.Xv,p(:,2),'-b');
% figure(101);
% hold on;
% plot(msh.f.Xv,pde.a.f(:,2),'-r'); %1a derivada
% plot(msh.f.Xv,pde.a.f(:,3),'-b'); %2a
% plot(msh.f.Xv,0.1.*pde.a.f(:,4),'-k'); %3a
% plot(msh.f.Xv,0.001.*pde.a.f(:,5),'-m'); %4a

% >> ----------------------------------------------------------------------

