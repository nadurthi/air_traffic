function Jxk1=KIRB_CT_eg_dyn_jac_disc(xk,para)
 T=para(1);
omg=xk(end,1);

F=[cos(omg*T)*T/omg*xk(3)-sin(omg*T)/omg^2*xk(3)-sin(omg*T)*T/omg*xk(4)-(-1+cos(omg*T))/omg^2*xk(4);
   sin(omg*T)*T/omg*xk(3)-(1-cos(omg*T))/omg^2*xk(3)-cos(omg*T)*T/omg*xk(4)-(sin(omg*T))/omg^2*xk(4); 
   -sin(omg*T)*T*xk(3)-cos(omg*T)*T*xk(4);
  cos(omg*T)*T*xk(3)-sin(omg*T)*T*xk(4) ];

Jxk1=[1,0,sin(omg*T)/omg,-(1-cos(omg*T))/omg,F(1);...
     0,1,(1-cos(omg*T))/omg,sin(omg*T)/omg,F(2);
     0,0,cos(omg*T),-sin(omg*T),F(3);
     0,0,sin(omg*T),cos(omg*T),F(4);
     0,0,0,0,1];


end