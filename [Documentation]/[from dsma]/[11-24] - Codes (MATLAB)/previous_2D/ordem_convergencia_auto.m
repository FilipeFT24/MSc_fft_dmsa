clear all
clc
close all

global u_convec_x u_convec_y gamma_diff numero_condicao cell_side Lref;

save_results=true;
erro_fitting_poly = false;

max_iter = 10;

cellnum = round(linspace(10,80,max_iter));

u_convec_x = 1;
u_convec_y = 2;
gamma_diff = -0.01;

iter=0;
vector_save2=[];
vector_save4=[];
vector_save6=[];
vector_save8=[];


while iter<max_iter
iter = iter +1

cell_side=cellnum(iter);

[Lref, norma1_phi, erro_phi_max, erro_medio_poly, erro_max_poly]=solver2d_function('WLS_2');
vector_save2(iter,:)=[Lref norma1_phi erro_phi_max erro_medio_poly erro_max_poly];

[Lref, norma1_phi, erro_phi_max, erro_medio_poly, erro_max_poly]=solver2d_function('WLS_4');
vector_save4(iter,:)=[Lref norma1_phi erro_phi_max erro_medio_poly erro_max_poly];

[Lref, norma1_phi, erro_phi_max, erro_medio_poly, erro_max_poly]=solver2d_function('WLS_6');
vector_save6(iter,:)=[Lref norma1_phi erro_phi_max erro_medio_poly erro_max_poly];

[Lref, norma1_phi, erro_phi_max, erro_medio_poly, erro_max_poly]=solver2d_function('WLS_8');
vector_save8(iter,:)=[Lref norma1_phi erro_phi_max erro_medio_poly erro_max_poly];

pause(0.01)
end


if erro_fitting_poly

figure(1)
loglog(vector_save2(:,1), vector_save2(:,2), 's-','color','r','MarkerFaceColor', 'r')
hold on
loglog(vector_save2(:,1), vector_save2(:,4), 's-','color','r')
loglog(vector_save4(:,1), vector_save4(:,2), 'o-','color','g','MarkerFaceColor', 'g')
loglog(vector_save4(:,1), vector_save4(:,4), 'o-','color','g')
loglog(vector_save6(:,1), vector_save6(:,2), '^-','color','b','MarkerFaceColor', 'b')
loglog(vector_save6(:,1), vector_save6(:,4), '^-','color','b')
loglog(vector_save8(:,1), vector_save8(:,2), 'd-','color','k','MarkerFaceColor', 'k')
loglog(vector_save8(:,1), vector_save8(:,4), 'd-','color','k')

figure(2)
loglog(vector_save2(:,1), vector_save2(:,3), 's-','color','r','MarkerFaceColor', 'r')
hold on
loglog(vector_save2(:,1), vector_save2(:,5), 's-','color','r')
loglog(vector_save4(:,1), vector_save4(:,3), 'o-','color','g','MarkerFaceColor', 'g')
loglog(vector_save4(:,1), vector_save4(:,5), 'o-','color','g')
loglog(vector_save6(:,1), vector_save6(:,3), '^-','color','b','MarkerFaceColor', 'b')
loglog(vector_save6(:,1), vector_save6(:,5), '^-','color','b')
loglog(vector_save8(:,1), vector_save8(:,3), 'd-','color','k','MarkerFaceColor', 'k')
loglog(vector_save8(:,1), vector_save8(:,5), 'd-','color','k')

else

figure(1)
loglog(vector_save2(:,1), vector_save2(:,2), 's-','color','r','MarkerFaceColor', 'r')
hold on
loglog(vector_save4(:,1), vector_save4(:,2), 'o-','color','g','MarkerFaceColor', 'g')
loglog(vector_save6(:,1), vector_save6(:,2), '^-','color','b','MarkerFaceColor', 'b')
loglog(vector_save8(:,1), vector_save8(:,2), 'd-','color','k','MarkerFaceColor', 'k')

figure(2)
loglog(vector_save2(:,1), vector_save2(:,3), 's-','color','r','MarkerFaceColor', 'r')
hold on
loglog(vector_save4(:,1), vector_save4(:,3), 'o-','color','g','MarkerFaceColor', 'g')
loglog(vector_save6(:,1), vector_save6(:,3), '^-','color','b','MarkerFaceColor', 'b')
loglog(vector_save8(:,1), vector_save8(:,3), 'd-','color','k','MarkerFaceColor', 'k')



end


%% rectas das varias ordens
for i=1:2

if i==1
    filename=sprintf('sin_conv_%s_%s_dif_%s_norm1', strrep(num2str(u_convec_x),'.','dot'), strrep(num2str(u_convec_y),'.','dot'), strrep(num2str(gamma_diff),'.','dot'));
else
    filename=sprintf('sin_conv_%s_%s_dif_%s_normmax', strrep(num2str(u_convec_x),'.','dot'), strrep(num2str(u_convec_y),'.','dot'), strrep(num2str(gamma_diff),'.','dot'));
end
    
figure(i)
ponto_ini_recta = round(max_iter/2+1);

b2=vector_save2(ponto_ini_recta,i+1)/(vector_save2(ponto_ini_recta,1)^2);
fun2 = @(x) (x.^2)*b2;
 
b4=vector_save4(ponto_ini_recta,i+1)/(vector_save4(ponto_ini_recta,1)^4);
fun4 = @(x) b4*x.^4;

b6=vector_save6(ponto_ini_recta,i+1)/(vector_save6(ponto_ini_recta,1)^6);
fun6 = @(x) b6*x.^6;

b8=vector_save8(ponto_ini_recta,i+1)/(vector_save8(ponto_ini_recta,1)^8);
fun8 = @(x) b8*x.^8;


x = [0.1 vector_save2(ponto_ini_recta,1) 1/cellnum(end)];
y2 = fun2(x);
y4 = fun4(x);
y6 = fun6(x);
y8 = fun8(x);


loglog(x, y2,'--','color','r')
loglog(x, y4,'--','color','g')
loglog(x, y6,'--','color','b')
loglog(x, y8,'--','color','k')

if erro_fitting_poly 
legend({'\bf 2^{nd}','\bf Poly 2^{nd}' ,'\bf 4^{th}','\bf Poly 4^{th}','\bf 6^{th}', '\bf Poly 6^{th}','\bf 8^{th}', '\bf Poly 8^{th}'},'Location', 'East','FontSize',12)
else
legend({'\bf 2^{nd}','\bf 4^{th}', '\bf 6^{th}', '\bf 8^{th}'},'Location', 'East','FontSize',12)
end

xlabel('\bf h_{ref}','fontsize',12)
if i==1
ylabel('\bf ||e||_{1}','fontsize',12)
else
ylabel('\bf ||e||_{\infty}','fontsize',12)  
end
% set(gca,'fontsize',12)
set(gca, 'XDir','reverse')
%axis tight
set(gca,'box','off')
legend boxoff 

% yticks([10^-15 10^-12 10^-9 10^-6 10^-3 10^0 10^3 10^6])
% ylim([0.1*min(vector_save8(:,2)) 10*max(vector_save2(:,3))])
% xlim([10^-3 10^-1])

if save_results
    print(filename, '-dpng')
    print(filename, '-deps')
    savefig(filename)
end

end

beep