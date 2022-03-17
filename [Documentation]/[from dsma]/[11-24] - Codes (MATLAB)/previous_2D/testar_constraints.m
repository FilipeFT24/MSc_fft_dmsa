function []=testar_constraints(M,T_border,constrained_source,c,Q_border,soluti,order)

global face_num G faces face_bound 


%% Sem constraint


for face_v2 = 1:face_num

    if face_bound(face_v2,1)==1
        
      
        
        
xf=faces(face_v2,1);
yf=faces(face_v2,2);

 erro_totalunc(face_v2) = 0;
  erro_totalgradxunc(face_v2) =0;
   erro_totalgradyunc(face_v2) =0; 

 x= M{face_v2}*soluti{face_v2}'; % unconstrained



for qq=1:order/2
    
    xc=G{face_v2}(qq,1);
    yc=G{face_v2}(qq,2);  
    
    if order ==6
phi(face_v2)=PolyReconstruction6thOrder(x,xc,yc,xf,yf,'poly');
gradx(face_v2)=PolyReconstruction6thOrder(x,xc,yc,xf,yf,'xflux'); 
grady(face_v2)=PolyReconstruction6thOrder(x,xc,yc,xf,yf,'yflux'); 
    elseif order==8
phi(face_v2)=PolyReconstruction8thOrder(x,xc,yc,xf,yf,'poly');
gradx(face_v2)=PolyReconstruction8thOrder(x,xc,yc,xf,yf,'xflux'); 
grady(face_v2)=PolyReconstruction8thOrder(x,xc,yc,xf,yf,'yflux'); 
    end

analitico_phi(face_v2)=SolutionDiffusion('sin',xc,yc,'anal');
analitico_grad_x(face_v2) = SolutionDiffusion('sin',xc,yc,'xflux');
analitico_grad_y(face_v2)= SolutionDiffusion('sin',xc,yc,'yflux');


 erro_totalunc(face_v2) = erro_totalunc(face_v2)+2/order*(abs(phi(face_v2)-analitico_phi(face_v2)));
  erro_totalgradxunc(face_v2) =erro_totalgradxunc(face_v2)+ 2/order*(abs(gradx(face_v2)-analitico_grad_x(face_v2)));
   erro_totalgradyunc(face_v2) =erro_totalgradyunc(face_v2)+ 2/order*(abs(grady(face_v2)-analitico_grad_y(face_v2))); 
end
    end
end



%% Com constraint
% for face_v2 = 1:20
% 
%     if face_bound(face_v2,1)==1
%         
%         
%         
%         
% xf=faces(face_v2,1);
% yf=faces(face_v2,2);
% 
%  erro_totalcon(face_v2) = 0;
%   erro_totalgradxcon(face_v2) =0;
%    erro_totalgradycon(face_v2) =0; 
% 
% x= T_border{face_v2}*soluti{face_v2}'+constrained_source{face_v2}; % constrained
% 
%  
%  
% for qq=1:order/2
%     
%     xc=G{face_v2}(qq,1);
%     yc=G{face_v2}(qq,2);  
% 
%     if order ==6
% phi(face_v2)=PolyReconstruction6thOrder(x,xc,yc,xf,yf,'poly');
% gradx(face_v2)=PolyReconstruction6thOrder(x,xc,yc,xf,yf,'xflux'); 
% grady(face_v2)=PolyReconstruction6thOrder(x,xc,yc,xf,yf,'yflux'); 
%     elseif order==8
% phi(face_v2)=PolyReconstruction8thOrder(x,xc,yc,xf,yf,'poly');
% gradx(face_v2)=PolyReconstruction8thOrder(x,xc,yc,xf,yf,'xflux'); 
% grady(face_v2)=PolyReconstruction8thOrder(x,xc,yc,xf,yf,'yflux'); 
%     end
% 
% analitico_phi(face_v2)=SolutionDiffusion('sin',xc,yc,'anal');
% analitico_grad_x(face_v2) = SolutionDiffusion('sin',xc,yc,'xflux');
% analitico_grad_y(face_v2)= SolutionDiffusion('sin',xc,yc,'yflux');
% 
% 
%  erro_totalcon(face_v2) = erro_totalcon(face_v2)+2/order*(abs(phi(face_v2)-analitico_phi(face_v2)));
%   erro_totalgradxcon(face_v2) =erro_totalgradxcon(face_v2)+ 2/order*(abs(gradx(face_v2)-analitico_grad_x(face_v2)));
%    erro_totalgradycon(face_v2) =erro_totalgradycon(face_v2)+ 2/order*(abs(grady(face_v2)-analitico_grad_y(face_v2))); 
% end
%     end
% end




%% constrained
%    media_phicon=mean(erro_totalcon)
%    media_gradxcon=mean(erro_totalgradxcon)
%    media_gradycon=mean(erro_totalgradycon)
%    
%    max_phicon=max(erro_totalcon)
%    max_gradxcon=max(erro_totalgradxcon)
%    max_gradycon=max(erro_totalgradycon)
%    
%    soma_phicon=sum(erro_totalcon)
%    soma_gradxcon=sum(erro_totalgradxcon)
%    soma_gradycon=sum(erro_totalgradycon)
   
   
%% unconstrained   
    media_phiunc=mean(erro_totalunc)
   media_gradxunc=mean(erro_totalgradxunc)
   media_gradyunc=mean(erro_totalgradyunc)
   
   max_phiunc=max(erro_totalunc)
   max_gradxunc=max(erro_totalgradxunc)
   max_gradyunc=max(erro_totalgradyunc)
   
   soma_phiunc=sum(erro_totalunc)
   soma_gradxunc=sum(erro_totalgradxunc)
   soma_gradyunc=sum(erro_totalgradyunc)
   


pause(0.1)
end