function Jyk=KIRB_eg_meas_jac_disc2(xk1,para)

Jyk=[xk1(1)/sqrt(xk1(1)^2+xk1(2)^2),xk1(2)/sqrt(xk1(1)^2+xk1(2)^2),0,0,0;-(xk1(2)/(xk1(1)^2+xk1(2)^2)),xk1(1)/(xk1(1)^2+xk1(2)^2),0,0,0];
end