function Ainv=matrix_inverter(a)

% A = [0.8147    0.0975    0.1576    0.1419    0.6557
%     0.9058    0.2785    0.9706    0.4218    0.0357
%     0.1270    0.5469    0.9572    0.9157    0.8491
%     0.9134    0.9575    0.4854    0.7922    0.9340
%     0.6324    0.9649    0.8003    0.9595    0.6787];

A=a;
tamanho=size(A);

Ainv=[];

I=eye(tamanho(1));
 LL=[];
 UU=[];
setup.type = 'ilutp';
setup.milu = 'col';
setup.droptol = 10^-1;
[LL,UU] = ilu(sparse(A),setup);

for i=1:tamanho(1)
b = I(:,i);


[aprox, flag,relres]=bicgstab(A,b,10^-13,50000, LL, UU);
i 
flag 
relres
 Ainv(:,i)=aprox;
end

%  [r,c] = size(a);
% 
% b = eye(r);
% 
% 
% for j = 1 : r
%     for i = j : r
%         if a(i,j) ~= 0
%             for k = 1 : r
%                 s = a(j,k); a(j,k) = a(i,k); a(i,k) = s;
%                 s = b(j,k); b(j,k) = b(i,k); b(i,k) = s;
%             end
%             t = 1/a(j,j);
%             for k = 1 : r
%                 a(j,k) = t * a(j,k);
%                 b(j,k) = t * b(j,k);
%             end
%             for L = 1 : r
%                 if L ~= j
%                     t = -a(L,j);
%                     for k = 1 : r
%                         a(L,k) = a(L,k) + t * a(j,k);
%                         b(L,k) = b(L,k) + t * b(j,k);
%                     end
%                 end
%             end  
%         else
%             
%             a
%             pause(1)
%             
%             
%         end
%          break
%     end
% 
%     if a(i,j) == 0
%         disp('Warning: Singular Matrix')
%         b = 'error';
%         return
%     end
% end
% 
% Ainv=b;
% 
% 
% % [L2,U,P] = lu(a);
% % 
% %  Ainv=inv(U)*inv(L2);
% % 
% % teste= Ainv*(P*a);
% 
% [U,S,V] = svd(a);
%    Ainv=V*inv(S)*U';



end
