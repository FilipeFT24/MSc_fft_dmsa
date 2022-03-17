clc
close all


figure(1)
loglog(vector_save2_1(:,1), vector_save2_1(:,2), 's-','color','r','MarkerFaceColor', 'r')
 hold on
loglog(vector_save4_1(:,1), vector_save4_1(:,2), 'o-','color','g','MarkerFaceColor', 'g')

loglog(vector_save6_1(:,1), vector_save6_1(:,2), '^-','color','b','MarkerFaceColor', 'b')
loglog(vector_save8_1(:,1), vector_save8_1(:,2), 'd-','color','k','MarkerFaceColor', 'k')

loglog(vector_save2_2(:,1), vector_save2_2(:,2), 's--','color','r','MarkerFaceColor', 'r')
loglog(vector_save4_2(:,1), vector_save4_2(:,2), 'o--','color','g','MarkerFaceColor', 'g')
loglog(vector_save6_2(:,1), vector_save6_2(:,2), '^--','color','b','MarkerFaceColor', 'b')
loglog(vector_save8_2(:,1), vector_save8_2(:,2), 'd--','color','k','MarkerFaceColor', 'k')

figure(2)
loglog(vector_save2_1(:,1), vector_save2_1(:,3), 's-','color','r','MarkerFaceColor', 'r')
hold on
loglog(vector_save4_1(:,1), vector_save4_1(:,3), 'o-','color','g','MarkerFaceColor', 'g')

loglog(vector_save6_1(:,1), vector_save6_1(:,3), '^-','color','b','MarkerFaceColor', 'b')
loglog(vector_save8_1(:,1), vector_save8_1(:,3), 'd-','color','k','MarkerFaceColor', 'k')


loglog(vector_save2_2(:,1), vector_save2_2(:,3), 's--','color','r','MarkerFaceColor', 'r')
loglog(vector_save4_2(:,1), vector_save4_2(:,3), 'o--','color','g','MarkerFaceColor', 'g')
loglog(vector_save6_2(:,1), vector_save6_2(:,3), '^--','color','b','MarkerFaceColor', 'b')
loglog(vector_save8_2(:,1), vector_save8_2(:,3), 'd--','color','k','MarkerFaceColor', 'k')




%% rectas das varias ordens
for i=1:2

if i==1
    filename=sprintf('norm1');
else
    filename=sprintf('norm_max');
end
    
 figure(i)
% ponto_ini_recta = round(max_iter/2+1);
% 
% b2=vector_save2(ponto_ini_recta,i+1)/(vector_save2(ponto_ini_recta,1)^2);
% fun2 = @(x) (x.^2)*b2;
%  
% b4=vector_save4(ponto_ini_recta,i+1)/(vector_save4(ponto_ini_recta,1)^4);
% fun4 = @(x) b4*x.^4;
% 
% b6=vector_save6(ponto_ini_recta,i+1)/(vector_save6(ponto_ini_recta,1)^6);
% fun6 = @(x) b6*x.^6;
% 
% b8=vector_save8(ponto_ini_recta,i+1)/(vector_save8(ponto_ini_recta,1)^8);
% fun8 = @(x) b8*x.^8;
% 
% 
% x = [0.1 vector_save2(ponto_ini_recta,1) 1/cellnum(end)];
% y2 = fun2(x);
% y4 = fun4(x);
% y6 = fun6(x);
% y8 = fun8(x);
% 
% 
% loglog(x, y2,'--','color','r')
% loglog(x, y4,'--','color','g')
% loglog(x, y6,'--','color','b')
% loglog(x, y8,'--','color','k')


legend({'\bf 2^{nd}','\bf 4^{th}', '\bf 6^{th}', '\bf 8^{th}'},'Location', 'East','FontSize',12)

% if erro_fitting_poly 
% legend({'\bf 2^{nd}','\bf Poly 2^{nd}' ,'\bf 4^{th}','\bf Poly 4^{th}','\bf 6^{th}', '\bf Poly 6^{th}','\bf 8^{th}', '\bf Poly 8^{th}'},'Location', 'East','FontSize',12)
% else
% legend({'\bf 2^{nd}','\bf 4^{th}', '\bf 6^{th}', '\bf 8^{th}'},'Location', 'East','FontSize',12)
% end

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

if true
    print(filename, '-dpng')
    print(filename, '-deps')
    savefig(filename)
end

end
