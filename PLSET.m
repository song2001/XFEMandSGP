function ETAN=PLSET(PROP)
%**********************************************************************
% Initialize elastic stiffness matrix
%
% ETAN  : Elastic stiffness matrix
%**********************************************************************
%%
  %
  LAM=PROP(1)*PROP(2)/((1+PROP(2))*(1-2*PROP(2)));
  MU=PROP(1)/(2*(1+PROP(2)));
  %
  %
  ETAN=[LAM+2*MU LAM       0;
        LAM      LAM+2*MU  0;
        0        0        MU];   
end