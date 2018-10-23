function dx=lorenz_cont(t,x)
sig=x(4);
p=x(5);

B=8/3;
dx(1,1)=sig*(x(2)-x(1));
dx(2,1)=p*x(1)-x(2)-x(1)*x(3);
dx(3,1)=x(1)*x(2)-B*x(3);
dx(4,1)=0;
dx(5,1)=0;

