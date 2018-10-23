function [t,xk1]=ode45_discc(FF,t0,dt,tF,xk,Qsq)
xk=xk(:);
xk1=xk';
for i=t0:dt:tF-dt
    xk=FF(xk)+Qsq*randn(length(xk),1);
    xk1=vertcat(xk1,xk');

end
t=t0:dt:tF;
end