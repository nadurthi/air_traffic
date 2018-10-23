function xint=transform_points_uniform(xint)
lb_vent_rad = 65; ub_vent_rad = 150; a = 1/2*(lb_vent_rad+ub_vent_rad); b = 1/2*(ub_vent_rad-lb_vent_rad);

xint(:,1) = a + b*xint(:,1); % vent radius is uniform
%%

lb_vel = 45; ub_vel = 124;  a = 1/2*(lb_vel+ub_vel); b = 1/2*(ub_vel-lb_vel);% vent velocity is Glomap distributed

xint(:,2) = a + b*xint(:,2);
%%
lb1_gs = 1.5; lb2_gs = 3; ub1_gs = 2; ub2_gs = 5;

a = 1/2*(lb1_gs+ub2_gs); b = 1/2*(ub2_gs-lb1_gs);% grain size
xint(:,3) = a + b*xint(:,3);

%%
me_sigma = 1.9; sig_sigma = 0.6; 
lb_sig = me_sigma-sig_sigma; ub_sig = me_sigma+sig_sigma;

a = 1/2*(lb_sig+ub_sig); b = 1/2*(ub_sig-lb_sig);
xint(:,4) = a + b*xint(:,4);
end