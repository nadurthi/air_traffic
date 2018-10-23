function [xu,Pu]=EKF(fn,F,h,H,xu,Pu,ym,Q,R,t0,dt,tf)
%this is discrete dynamics
global Js xs n q
Js=F;
xs=xu;
q=Q;

n=length(xu);
% forecast
% xf=f(xu);
[~,x]=ode45(fn,t0:dt:tf,xu');
xf=x(end,:)';
[~,P]=ode45(@cov_prop,t0:dt:tf,reshape(Pu,1,n*n));
Pf=reshape(P(end,:),n,n);
% Pf=F(xu)*Pu*F(xu)'+Q;

%update
if ym==-1
    xu=xf;
    Pu=Pf;
else
    K=Pf*H(xf)'*inv(H(xf)*Pf*H(xf)'+R);
    xu=xf+K*(ym-h(xf));
    Pu=(eye(n)-K*H(xf))*Pf;
end
end

function dp=cov_prop(t,p)
global Js xs n q
F=Js;
xu=xs;
Q=q;
P=reshape(p,n,n);
dpp=F(xu)*P+P*F(xu)'+Q;
dp=reshape(dpp,1,n*n);
dp=dp';
end