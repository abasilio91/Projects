function graphsDTG(X,y)

koh=X(1);
Eoh=X(2);
koc=X(3);
Eoc=X(4);
koL=X(5);
EoL=X(6);
fra1=X(7);
fra2=X(8);
fra3=X(9);
n1=X(10);
n2=X(11);
n3=X(12);

load data_carpel_5.txt

time = y(:,1);
mexp  = y(:,2);
dmexp =(-1)* y(:,3);
xexp  = y(:,4);
dydt  = (-1)*y(:,5);
Temperature = y(:,6);

% Confirming the numerical temperature ramping
% T=450.35+(5/60)*(time);
% figure(1)
% plot(time,Temperature,'b',time,T,'r');

options = odeset('RelTol',1e-5,'AbsTol',1e-5);

sol = ode23s(@DEQ,[0 time(length(time))],[0 0 0],options,koh,Eoh,koc,Eoc,koL,EoL,n1,n2,n3);

Y(:,1) = deval(sol,time,1);
Y(:,2) = deval(sol,time,2);
Y(:,3) = deval(sol,time,3);
Yt=fra1*Y(:,1)'+fra2*Y(:,2)'+fra3*Y(:,3)';

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

% FIGURES
figure(2)
plot(Temperature(1:length(dH)),dH,'g--',Temperature(1:length(dC)),dC,'k-.',Temperature(1:length(dyL)),dL,'r-',Temperature(1:length(dm)),dm,'b-',Temperature(1:10:length(dm)),dmexp(1:10:length(dm)),'k.')
ylabel('-dm/dt [mg/s]')
xlabel('Temperature [K]')
legend('Hemicellulose','Cellulose','Lignin', 'Total calculated', 'Experimental',1)

figure(3)
plot(Temperature(1:10:length(m)),mexp(1:10:length(m)),'k.',Temperature,m,'k-')
ylabel('Mass [mg]')
xlabel('Temperature [K]')
legend('Experimental','Calculated',1)
% 
% figure(4)
% plot(Temperature(1:length(dy)),dyH,'r-',Temperature(1:length(dy)),dyC,'b-',Temperature(1:length(dy)),dyL,'g-',Temperature(1:length(dy)),dy,'k-',Temperature(1:10:length(dy)),dydt(1:10:length(dy)),'k.')
% ylabel('-dy/dt')
% xlabel('Temperature [K]')
% legend('Hemicellulose','Cellulose','Lignin', 'Total calculated', 'Experimental',2)
% 
% figure(5)
% plot(Temperature(1:10:length(xexp)),xexp(1:10:length(xexp)),'r.',Temperature,Yt,'b-')
% ylabel('y')
% xlabel('Temperature [K]')
% legend('Experimental','Calculated',2)

fprintf(1,' ===============================================================================================================\n');
fra = fra1+fra2+fra3;
fprintf ('total fraction = %f\n', fra);

fprintf(1,' ===============================================================================================================\n');
nt=length(Yt);
resid=xexp-Yt';
SQR=sum(resid.*resid);
Syy=sum((xexp-mean(xexp)).*(xexp-mean(xexp)));
R2=(Syy-SQR)/Syy;
% fprintf(1,'R2 for X total= %f\n',R2);

% fprintf(1,' ===============================================================================================================\n');
mnt=length(m);
mresid=mexp-m';
mSQR=sum(mresid.*mresid);
mSyy=sum((mexp-mean(mexp)).*(mexp-mean(mexp)));
mR2=(mSyy-mSQR)/mSyy;
% fprintf(1,'R2 para m total= %f\n',mR2);
mF=((sum(m'.*m'))/2)/(mSQR/(mnt-2));
mFIT= (100*(mSQR/mnt)^0.5)/max(m);
fprintf(1,'FIT for TG total = %f\n',mFIT);

% fprintf(1,' ===============================================================================================================\n');
dmnt=length(dmexp);
dmresid=dmexp(1:length(dy))-dm';
dmSQR=sum(dmresid.*dmresid);
dmSyy=sum((dmexp-mean(dmexp)).*(dmexp-mean(dmexp)));
dmR2=(dmSyy-dmSQR)/dmSyy;
% fprintf(1,'R2 for dmdt total= %f\n',dmR2);
dmF=((sum(dm'.*dm'))/2)/(dmSQR/(dmnt-2));
dmFIT= (100*(dmSQR/dmnt)^0.5)/max(dm);
fprintf(1,'FIT for DTG total = %f\n',dmFIT);

% fprintf(1,' ===============================================================================================================\n');
dnt=length(dydt);
dresid=dydt(1:length(dy))-dy';
dSQR=sum(dresid.*dresid);
dvar=dSQR/(dnt-2);
dSyy=sum((dydt-mean(dydt)).*(dydt-mean(dydt)));
dR2=(dSyy-dSQR)/dSyy;
% fprintf(1,'R2 for dxdt total = %f\n',dR2);
dF=((sum(dy'.*dy'))/2)/(dSQR/(dnt-2));
dFIT= (100*(dSQR/dnt)^0.5)/max(dydt);
% fprintf(1,'FIT for dxdt total = %f\n',dFIT);

Hdnt=length(dydt(1:297));
Hdresid=dydt(1:297)-dy(1:297)';
HdSQR=sum(Hdresid.*Hdresid);
Hdvar=HdSQR/(Hdnt-2);
HdSyy=sum((dydt(1:297)-mean(dydt(1:297))).*(dydt(1:297)-mean(dydt(1:297))));
HdR2=(HdSyy-HdSQR)/HdSyy;
% fprintf(1,'R2 for dxdt Hemicelulose = %f\n',HdR2);
HdF=((sum(dy(1:297)'.*dy(1:297)'))/2)/(HdSQR/(Hdnt-2));
HdFIT= (100*(HdSQR/Hdnt)^0.5)/max(dydt(1:297));

Cdnt=length(dydt(298:415));
Cdresid=dydt(298:415)-dy(298:415)';
CdSQR=sum(Cdresid.*Cdresid);
Cdvar=CdSQR/(Cdnt-2);
CdSyy=sum((dydt(298:415)-mean(dydt(298:415))).*(dydt(298:415)-mean(dydt(298:415))));
CdR2=(CdSyy-CdSQR)/CdSyy;
% fprintf(1,'R2 for dxdt Celulose= %f\n',CdR2);
CdF=((sum(dy(298:415)'.*dy(298:415)'))/2)/(CdSQR/(Cdnt-2));
CdFIT= (100*(CdSQR/Cdnt)^0.5)/max(dydt(298:415));

Ldnt=length(dydt(416:957));
Ldresid=dydt(416:957)-dy(416:957)';
LdSQR=sum(Ldresid.*Ldresid);
Ldvar=LdSQR/(Ldnt-2);
LdSyy=sum((dydt(416:957)-mean(dydt(416:957))).*(dydt(416:957)-mean(dydt(416:957))));
LdR2=(LdSyy-LdSQR)/LdSyy;
% fprintf(1,'R2 for dxdt lignina = %f\n',LdR2);
LdF=((sum(dy(416:957)'.*dy(416:957)'))/2)/(LdSQR/(Ldnt-2));
LdFIT= (100*(LdSQR/Ldnt)^0.5)/max(dydt(416:957));

% ==========================================================================

 % OUT FILES
u = fopen('R2.dat', 'wt');
fprintf(u,'R2 x = %1.3f R2 dxdt = %1.3f ',R2,dR2);
fclose(u);

uH = fopen('parametersH.dat', 'wt');
fprintf(uH,'koh = %1.4e Eh = %1.4e N = %4.0f VAR = %1.4e',koh,Eoh,Hdnt,Hdvar);
fclose(uH);

uC = fopen('parametersC.dat', 'wt');
fprintf(uC,'koc = %1.4e Ec = %1.4e N = %4.0f VAR = %1.4e',koc,Eoc,Cdnt,Cdvar);
fclose(uC);

uL = fopen('parametersL.dat', 'wt');
fprintf(uL,'koL = %1.4e EL = %1.4e N = %4.0f VAR = %1.4e',koL,EoL,Ldnt,Ldvar);
fclose(uL);

H=[Temperature(1:297) xexp(1:297)];
save parallelsH.dat H -ASCII

C=[Temperature(298:415) xexp(298:415)];
save parallelsC.dat C -ASCII

L=[Temperature(416:957) xexp(416:957)];
save parallelsL.dat L -ASCII

Hs= fopen('Inferential_resultsH.dat','wt');
fprintf(Hs,'Variance = %1.4e R2 = %1.3f F = %6.3f FIT = %6.3f',Hdvar,HdR2,HdF,HdFIT);
fclose(Hs);

Cs= fopen('Inferential_resultsC.dat','wt');
fprintf(Cs,'Variance = %1.4e R2 = %1.3f F = %6.3f FIT = %6.3f',Cdvar,CdR2,CdF,CdFIT);
fclose(Cs);

Ls= fopen('Inferential_resultsL.dat','wt');
fprintf(Ls,'Variance = %1.4e R2 = %1.3f F = %6.3f FIT = %6.3f',Ldvar,LdR2,LdF,LdFIT);
fclose(Ls);

G=[Temperature(1:length(dm)) dH dC dL dm(1:length(dm))' dmexp(1:length(dm))];
xlswrite('dmorigin', G)

H=[Temperature(1:length(dm)) m(1:length(dm))' mexp(1:length(dm))];
xlswrite('morigin', H)

J=[Temperature(1:length(dy)) Yt(1:length(dy))' xexp(1:length(dy))];
xlswrite('xorigin', J)

L=[Temperature(1:length(dy)) dyH(1:length(dy)) dyC(1:length(dy)) dyL(1:length(dy)) dydt(1:length(dy)) dy(1:length(dy))'];
xlswrite('dyorigin', L)

function dY = DEQ(tt,yy,koh,Eoh,koc,Eoc,koL,EoL,n1,n2,n3)
dY = zeros(3,1);
dY(1) = koh*exp(-Eoh/(8.314*(450.35+5/60*tt)))*(1-yy(1))^n1;
dY(2) = koc*exp(-Eoc/(8.314*(450.35+5/60*tt)))*(1-yy(2))^n2;
dY(3) = koL*exp(-EoL/(8.314*(450.35+5/60*tt)))*(1-yy(3))^n3;
