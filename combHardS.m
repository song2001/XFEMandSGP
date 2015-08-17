%
% Linear combined isotropic/kinamtic hardening model
%
function [stress, dp, r]=combHardS(mp,D,deps,stressN,r)
% Inputs:
% mp = [lambda, mu, beta, H, Y0];
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
    stressT = stresstr; dp = 0; r=r;          %trial states are final
    stress(1,1)=stressT(1);
    stress(2,2)=stressT(2);
    stress(3,3)=stressT(4);
    stress(1,2)=stressT(3);
    stress(2,1)=stress(1,2);
    stress(1,3)=0.0;
    stress(2,3)=0.0;
    stress(3,1)=stress(1,3);
    stress(3,2)=stress(2,3);
 else
    r0=r;
    dp=0.0;
    res=1.0;    
      while abs(res)>tol    
       res=etat-3*mu*dp-r-mp(4);
       dp=dp+res/(3*mu+mp(3));
       r=r0+mp(3)*dp;
      end
    dpstrn=3/2*dp*eta/etat;                    % Plastic strain increments
    dpstrn(3)=2*dpstrn(3);                     % Engineering shear strains
    for i=1:3
    destrn(i)=deps(i)-dpstrn(i); 
    end                                        % Elastic strain increments
    destrn(4)=0-dpstrn(4); 
    Dnew=[D(1,1) D(1,2) D(1,3) D(1,2);
          D(2,1) D(2,2) D(2,3) D(1,2);
          D(3,1) D(3,2) D(3,3) 0;
          D(1,2) D(1,2) 0 D(2,2)];
    dstress=Dnew*destrn';
    stressT=Sold+dstress;
    stress(1,1)=stressT(1);
    stress(2,2)=stressT(2);
    stress(3,3)=stressT(4);
    stress(1,2)=stressT(3);
    stress(2,1)=stress(1,2);
    stress(1,3)=0.0;
    stress(2,3)=0.0;
    stress(3,1)=stress(1,3);
    stress(3,2)=stress(2,3);
 end
end