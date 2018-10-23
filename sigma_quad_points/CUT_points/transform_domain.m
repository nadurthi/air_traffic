function [XB,d]=transform_domain(Xb,bl,bu,Bl,Bu)

ns=size(Xb,1);
nx=size(Xb,2);

rb=0.5*(bl+bu);
%shifting to 0 mean
X0=Xb-repmat(rb,ns,1);

%scaling the points
db=bu-bl;
dB=Bu-Bl;

d=dB./db;

Xs=repmat(d,ns,1).*X0;

%shifting the mean

rB=0.5*(Bl+Bu);
% rmu=rB-rb;

XB=Xs+repmat(rB,ns,1);



end

