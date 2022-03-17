clear all
clc
close all

global u_convec gamma_diff numero_condicao cell_num fl0 Lref;

save_results=true;

max_iter = 15;

cellnum = round(linspace(10,400,max_iter));

u_convec = 1;
gamma_diff = -1;

iter=0;
vector_save2=[];
vector_save4=[];
vector_save6=[];
vector_save8=[];


while iter<max_iter
iter = iter +1

cell_num=cellnum(iter);

[numero_condicao,fl0,Peclet_number, Lref, norma1_phi, erro_phi_max]=solver1ddc_function('2ndOrder');
vector_save2(iter,:)=[Lref norma1_phi erro_phi_max];

[numero_condicao,fl0,Peclet_number, Lref, norma1_phi, erro_phi_max]=solver1ddc_function('4thOrder');
vector_save4(iter,:)=[Lref norma1_phi erro_phi_max];

[numero_condicao,fl0,Peclet_number, Lref, norma1_phi, erro_phi_max]=solver1ddc_function('6thOrder');
vector_save6(iter,:)=[Lref norma1_phi erro_phi_max];

[numero_condicao,fl0,Peclet_number, Lref, norma1_phi, erro_phi_max]=solver1ddc_function('8thOrder');
vector_save8(iter,:)=[Lref norma1_phi erro_phi_max];

end


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



%% rectas das varias ordens
for i=1:2

if i==1
    filename=sprintf('exp_conv_%s_dif_%s_norm1', strrep(num2str(u_convec),'.','dot'), strrep(num2str(gamma_diff),'.','dot'));
else
    filename=sprintf('exp_conv_%s_dif_%s_normmax', strrep(num2str(u_convec),'.','dot'), strrep(num2str(gamma_diff),'.','dot'));
end
    
figure(i)
ponto_ini_recta = 8;

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

legend({'\bf 2^{nd}','\bf 4^{th}', '\bf 6^{th}', '\bf 8^{th}'},'Location', 'East','FontSize',12)
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

yticks([10^-15 10^-12 10^-9 10^-6 10^-3 10^0 10^3 10^6])
ylim([0.1*min(vector_save8(:,2)) 10*max(vector_save2(:,3))])
xlim([10^-3 10^-1])

if save_results
    print(filename, '-dpng')
    print(filename, '-deps')
    savefig(filename)
end

end