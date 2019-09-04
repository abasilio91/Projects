function result = eval_objective(x,y)

% x is the project variables vector 
% y is the experimental data matrix

koh=x(1);
Eoh=x(2);
koc=x(3);
Eoc=x(4);
koL=x(5);
EoL=x(6);
fra1=x(7);
fra2=x(8);
fra3=x(9);
n1=x(10);
n2=x(11);
n3=x(12);
  
time = y(:,1);
mexp  = y(:,2);
dmexp =(-1)* y(:,3);
xexp  = y(:,4);
dydt  = (-1)*y(:,5);
Temperature = y(:,6);

options = odeset('RelTol',1e-5,'AbsTol',1e-5);

% solver for the solution of the ODE system
sol = ode23s(@DEQ,[0 time(length(time))],[0 0 0],options,koh,Eoh,koc,Eoc,koL,EoL,n1,n2,n3);

% associating the solution to each pseudocomponent curve
Y(:,1) = deval(sol,time,1);
Y(:,2) = deval(sol,time,2);
Y(:,3) = deval(sol,time,3);
Yt=fra1*Y(:,1)'+fra2*Y(:,2)'+fra3*Y(:,3)';

% Numerical differentiation of the solution
mH = 10.32687*ones(1,length(mexp))'-Y(:,1).*(10.32687-3.45199);
mC = 10.32687*ones(1,length(mexp))'-Y(:,2).*(10.32687-3.45199);
mL = 10.32687*ones(1,length(mexp))'-Y(:,3).*(10.32687-3.45199);
m = fra1*mH'+fra2*mC'+fra3*mL';

dyH=fra1*diff(Y(:,1))/5;
dyC=fra2*diff(Y(:,2))/5;
dyL=fra3*diff(Y(:,3))/5;
dy=dyH'+dyC'+dyL';

dH=(-1)*fra1*diff(mH)/5;
dC=(-1)*fra2*diff(mC)/5;
dL=(-1)*fra3*diff(mL)/5;
dm=dH'+dC'+dL';


erro = dy' - dydt(1:length(dy));
% erro = Yt' - xexp(1:length(Yt)); 
result = sum(erro.*erro); % objective function - minimization of the square error

function dY = DEQ(tt,yy,koh,Eoh,koc,Eoc,koL,EoL,n1,n2,n3)
% This is the convertion rate equation in terms of conversion (yy)
% considering the Arrhenius equation for the reaction constant as a
% function of the Temperature.
% dY = A *(1-Y)^n
% A = k*exp(-E/RT(t)). R stands for the universal constant for perfect
% gases

dY = zeros(3,1);
dY(1) = koh*exp(-Eoh/(8.314*(450.35+5/60*tt)))*(1-yy(1))^n1;
dY(2) = koc*exp(-Eoc/(8.314*(450.35+5/60*tt)))*(1-yy(2))^n2;
dY(3) = koL*exp(-EoL/(8.314*(450.35+5/60*tt)))*(1-yy(3))^n3;

