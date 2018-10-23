function xk1=KIRB_UM_eg_dyn_disc(xk,T)
%  global T
% T=1;

xk1=[1,0,T,0;...
     0,1,0,T;
     0,0,1,0;
     0,0,0,1]*xk;
end