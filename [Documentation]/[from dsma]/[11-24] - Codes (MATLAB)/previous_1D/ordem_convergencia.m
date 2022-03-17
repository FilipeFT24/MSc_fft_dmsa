clc
clear all
close all



%% 2ª ordem
href_uni_exp = [3.125000E-02 7.812500E-03 3.906250E-03 9.765625E-04]; 
erro_uni_exp = [2.161906E-03 1.326343E-04 3.317221E-05 2.072693E-06]; %malha uniforme exponencial

href_nuni_exp = [3.125000E-02 7.812500E-03 3.906250E-03 9.765625E-04];
erro_nuni_exp = [3.607757E-03 1.428607E-04 3.518088E-05 2.401025E-06];%malha nao uniforme exponencial

href_uni_sin = [6.250000E-02 1.562500E-02 3.906250E-03 9.765625E-04 2.500000E-04];
erro_uni_sin = [1.905711E-02 1.171591E-03 7.310265E-05 4.568508E-06 2.994038E-07];%malha  uniforme sin

href_nuni_sin = [1.562500E-02 7.812500E-03 3.906250E-03 9.765625E-04];
erro_nuni_sin = [1.275669E-03 3.179957E-04 7.984329E-05 5.030109E-06];%malha nao uniforme sin
%% 4ªordem
href_uni_exp4 = [6.25*10^-2 3.125*10^-2 9.765625E-04 6.250000E-04]; 
erro_uni_exp4 = [5.541045*10^-4 3.619314*10^-5 3.659273E-11 8.163930E-12];

href_nuni_sin4 = [3.125000E-02 7.812500E-03 1.111111E-03];
erro_nuni_sin4 = [8.620397E-05 1.793936E-07 6.980842E-11];%malha nao uniforme sin

%% 6ªordem
href_uni_exp6 = [6.250000E-02 3.125000E-02 1.562500E-02 7.812500E-03 3.906250E-03 %{9.765625E-04}%
    ]; 
erro_uni_exp6 = [3.164634E-04 3.453163E-06 5.216922E-08 8.076907E-10 1.262032E-11 %7.181325E-12
    ];

href_nuni_sin6 = [6.250000E-02 3.125000E-02 1.562500E-02 7.812500E-03 3.906250E-03 %1.000000E-03
    ];
erro_nuni_sin6 = [9.985171E-05 6.631159E-06 4.930998E-08 2.465071E-10 3.816441E-12 %5.086554E-14
    ];%malha nao uniforme sin
%% 8ªordem
href_uni_exp8 = [6.250000E-02 3.125000E-02 1.562500E-02 7.812500E-03 3.906250E-03]; 
erro_uni_exp8 = [1.030251E-03 6.165797E-07 2.728913E-09 8.770699E-12 1.038236E-13];

href_nuni_sin8 = [6.250000E-02 3.125000E-02 1.562500E-02 7.812500E-03];
erro_nuni_sin8 = [1.161590E-04 5.025493E-07 1.009145E-09 1.720992E-12];%malha nao uniforme sin

%% caso Analitico 1 (0.1 -0.1) 6ªordem explicito
href_case1_6 = [6.250000E-02 3.125000E-02 2.222222E-02 1.562500E-02];
erro_case1_6 = [9.090512E-10 1.464473E-11 1.911227E-12 1.588700E-13];

%% caso Analitico 1 (0.1 -0.1) 6ªordem implicito+diagonal dom
href_case1_6imp = [6.250000E-02 3.125000E-02 1.562500E-02 7.812500E-03];
erro_case1_6imp = [1.038054E-09 5.279065E-10 2.047806E-09 8.177142E-09];

%% caso Analitico 2 (2.5 -0.1) 8ªordem implicito
href_case2_8 = [6.250000E-02 3.125000E-02 1.562500E-02 7.812500E-03 3.906250E-03 3.333333E-03];
erro_case2_8 = [3.620526E-04 6.697933E-06 6.639659E-08 4.378811E-10 2.325382E-12 6.736316E-13];
%% caso Analitico 2 (2.5 -0.1) 8ªordem implicito diagonal dom
href_case2_8imp = [6.250000E-02 3.125000E-02 1.562500E-02 7.812500E-03 3.906250E-03];
erro_case2_8imp = [3.620526E-04 6.697933E-06 6.647796E-08 1.183861E-09 2.949271E-09];
%% 2d 4ºordem
href_sin_4 = [0.1 0.05 0.0333333333 0.025 2.000000E-02 1.428571E-02 1.000000E-02 6.666667E-03];
erro_sin_4 = [2.2729*10^-2 1.5330*10^-3 3.2198*10^-4 1.0608*10^-4 4.490664E-05 1.215457E-05 3.013571E-06 6.110220E-07];
%% rectas das varias ordens
b2=erro_uni_exp(1)/(href_uni_exp(1)^2);
fun2 = @(x) (x.^2)*b2;
% 
b4=erro_sin_4(1)/(href_sin_4(1)^4);
fun4 = @(x) b4*x.^4;

% b6=erro_uni_exp6(1)/(href_uni_exp6(1)^6);
% fun6 = @(x) b6*x.^6;

b6=erro_case1_6(1)/(href_case1_6(1)^6);
fun6 = @(x) b6*x.^6;

% b8=erro_uni_exp8(2)/(href_uni_exp8(2)^8);
% fun8 = @(x) b8*x.^8;

b8=erro_case2_8imp(1)/(href_case2_8imp(1)^8);
fun8 = @(x) b8*x.^8;

x = [0.1  0.05 0.01  0.005 0.001 0.0005 0.0001];
y2 = fun2(x);
y4 = fun4(x);
y6 = fun6(x);
y8 = fun8(x);




%% Plots
figure()
loglog(href_uni_exp, erro_uni_exp,'o')
hold on
% %loglog(href_nuni_exp, erro_nuni_exp,'*')
% %loglog(href_uni_sin, erro_uni_sin,'o')
% loglog(href_nuni_sin, erro_nuni_sin,'*')
% 
% loglog(href_uni_exp4, erro_uni_exp4,'o')
% loglog(href_nuni_sin4, erro_nuni_sin4,'*')
% 
% loglog(href_uni_exp6, erro_uni_exp6,'o')
% loglog(href_nuni_sin6, erro_nuni_sin6,'*')
% 
% loglog(href_uni_exp8, erro_uni_exp8,'o')
% loglog(href_nuni_sin8, erro_nuni_sin8,'*')

% loglog(href_case1_6, erro_case1_6,'*')
% loglog(href_case2_8, erro_case2_8,'o')
% loglog(href_case1_6imp, erro_case1_6imp,'d')
loglog(href_sin_4, erro_sin_4,'x')

%loglog(x, x,'-')
loglog(x, y2,'-')
 loglog(x, y4,'-')
%loglog(x, y6,'-')
%loglog(x, y8,'-')

% legend('exp uniforme 2ªordem','sin nao-uniforme 2ªordem','exp uniforme 4ªordem','sin nao-uniforme 4ºordem','exp uniforme 6ªordem','sin nao-uniforme 6ºordem','exp uniforme 8ªordem','sin nao-uniforme 8ºordem','caso 1','caso 2', '2ºordem', '4ºordem','6ªordem','8ºordem','Location', 'Best')
legend('exp uniforme 2ªordem','sin 2D 4ºordem','2ªordem', '4ªordem','Location', 'Best')