function xk1=KIRB_CT_eg_dyn_disc(xk,T)
%  global T
% T=1;
xk=xk(:);
omg=xk(end,1);
try
xk1=[1,0,sin(omg*T)/omg,-(1-cos(omg*T))/omg,0;...
     0,1,(1-cos(omg*T))/omg,sin(omg*T)/omg,0;
     0,0,cos(omg*T),-sin(omg*T),0;
     0,0,sin(omg*T),cos(omg*T),0;
     0,0,0,0,1]*xk;
catch
    keyboard
end

% xk1=xk1';
end
%% discrete dynamic equations for MC
