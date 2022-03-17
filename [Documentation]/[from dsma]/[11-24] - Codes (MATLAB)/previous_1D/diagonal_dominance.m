function factor=diagonal_dominance(A)

tamanho = size(A);

a = sum(abs(A(round(tamanho(1)/2),:)));
b = abs(A(round(tamanho(1)/2),round(tamanho(1)/2)));
factor1 = 2*b-a;

a = sum(abs(A(1,:)));
b = abs(A(1,1));
factor2 = 2*b-a;

a = sum(abs(A(tamanho(1),:)));
b = abs(A(tamanho(1),tamanho(1)));
factor3 = 2*b-a;

factor = max(abs([factor1 factor2 factor3]));
factor = 1.1*(factor);

end