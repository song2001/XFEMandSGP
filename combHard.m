function [stress, ep]=combHard(mp,deps,stressNP,epN)
% Inputs:
% mp = [lambda, mu, beta, H, Y0];
% stressN = [s11, s22, s33, t12, t23, t13];
% alphaN = [a11, a22, a33, a12, a23, a13];
%
stressN(1)=stressNP(1,1);
stressN(2)=stressNP(2,2);
stressN(3)=stressNP(3,3);
stressN(4)=stressNP(2,3);
stressN(5)=stressNP(1,3);
stressN(6)=stressNP(1,2);
stressN=stressN';
Iden = [1 1 1 0 0 0]'; 
two3 = 2/3; stwo3=sqrt(two3);                   %constants
mu=mp(1)/(2*(1+mp(2))); lam=mp(1)*mp(2)/((1+mp(2))*(1-2*mp(2)));
H=mp(3); Y0=mp(4);        %material properties
ftol = Y0*1E-6;                                 %tolerance for yield
D=[lam+2*mu lam      lam      0  0  0;
   lam      lam+2*mu lam      0  0  0;
   lam      lam      lam+2*mu 0  0  0;
   0        0        0        mu 0  0;
   0        0        0        0  mu 0;
   0        0        0        0  0  mu];
stresstr = stressN + D*deps;                    %trial stress
I1 = sum(stresstr(1:3));                        %trace(sigmatr)
str = stresstr - I1*Iden/3;                     %deviatoric stress
eta = str;                             %shifted stress
etat = sqrt(eta(1)^2 + eta(2)^2 + eta(3)^2 ...
          + 2*(eta(4)^2 + eta(5)^2 + eta(6)^2));%norm of eta
fyld = etat - stwo3*(Y0+H*epN);        %trial yield function
if fyld < ftol                                  %yield test
    stressP = stresstr; ep = epN;                %trial states are final 
    stress(1,1)=stressP(1);
    stress(2,2)=stressP(2);
    stress(3,3)=stressP(3);     
    stress(1,2)=stressP(6);
    stress(2,1)=stressP(6);
    stress(1,3)=stressP(5);
    stress(3,1)=stressP(5);
    stress(2,3)=stressP(4);
    stress(3,2)=stressP(4);
   
    return;
else
    gamma = fyld/(2*mu + two3*H);               %plastic consistency param
    ep = epN + gamma*stwo3;                     %updated eff. plastic strain
end
N = eta/etat;                                   %unit vector normal to f
stressP = stresstr - 2*mu*gamma*N;               %updated stress
    stress(1,1)=stressP(1);
    stress(2,2)=stressP(2);
    stress(3,3)=stressP(3);     
    stress(1,2)=stressP(6);
    stress(2,1)=stressP(6);
    stress(1,3)=stressP(5);
    stress(3,1)=stressP(5);
    stress(2,3)=stressP(4);
    stress(3,2)=stressP(4);  
    