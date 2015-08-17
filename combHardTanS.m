%
% Tangent stiffness for linear isotropic hardening model
%
function [Dtan]=combHardTanS(mp,D,deps,stressN,r)
% Inputs:
% mp = [E, mu, H, Y0];
% D = elastic stiffness matrix
% stressN = [s11, s22, t12, s33];
%
tol=1E-6;
mu=mp(1)/(2*(1+mp(2)));
K=mp(1)/(3*(1-2*mp(2)));
Sold=stressN;                               % Save stress at the beggining
St=stressN(4);
stressN(4)=[];
stresstr = stressN + D*deps;
stresstr(4)=St+D(2,1)*(deps(1)+deps(2));    % Obtain trial stress
Iden = [1 1 0 1]'; 
I1 = stresstr(1)+stresstr(2)+stresstr(4);       %trace(sigmatr)
eta = stresstr - I1*Iden/3;                     %deviatoric stress
etat=sqrt((3/2)*(eta(1)^2+eta(2)^2+eta(4)^2+2*(eta(3)^2))); %effective
N = eta/etat;                               % Trial flow direction
fyld = etat-mp(4)-r;                           % Yield condition
 if fyld<0
    Dtan = D;                                %elastic
 else
  r0=r;
  dp=0;
  res=1;    
   while abs(res)>tol    
   res=etat-3*mu*dp-r-mp(4);
   dp=dp+res/(3*mu+mp(3));
   r=r0+mp(3)*dp;
   end
 XR=(etat-3*mu*dp)/etat;
 Q=(1/(1+3*mu/mp(3))-XR)*3/2;
 Dtan(1,1)=2*mu*Q*N(1)*N(1)+(K-mu*XR*2/3)+2*mu*XR;
 Dtan(2,2)=2*mu*Q*N(2)*N(2)+(K-mu*XR*2/3)+2*mu*XR;
 Dtan(1,2)=2*mu*Q*N(1)*N(2)+(K-mu*XR*2/3);
 Dtan(2,1)=Dtan(1,2);
 Dtan(2,3)=2*mu*Q*N(2)*N(3);
 Dtan(3,2)=Dtan(2,3); 
 Dtan(1,3)=2*mu*Q*N(1)*N(3);
 Dtan(3,1)=Dtan(1,3);
 Dtan(3,3)=2*mu*Q*N(3)*N(3)+mu*XR;
 end    
end