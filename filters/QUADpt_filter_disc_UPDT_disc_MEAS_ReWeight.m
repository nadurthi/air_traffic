function [xu,Pu]=QUADpt_filter_disc_UPDT_disc_MEAS_ReWeight(model,mu_ut,P_ut,ym,para)

fx=model.fx;
hx=model.hx;
Q=model.Q;
R=model.R;
n=model.fn;
para_dt=model.para_dt;


      qd_pts=@conjugate_dir_gausspts_till_8moment;
 

    [x,w]=qd_pts(mu_ut,P_ut);
    
    [mu_ut_sf,P_ut_sf]=prop_mean_cov_points_discr_modified(x,w,fx,para_dt,model);

      P_ut_sf=P_ut_sf+Q;
      
%    Check if meas is avail at this time step
        if ym==-1234
        xu=mu_ut_sf;
        Pu=P_ut_sf;
        return
        end
        
         %obs forecast
         model.dynamics='discrete';
        
    [x,w]=qd_pts(mu_ut_sf,P_ut_sf);
    [mu_ut_of,P_ut_of]=prop_mean_cov_points_discr_modified(x,w,hx,para_dt,model);

    P_ut_of=P_ut_of+R;

    %cross cov
    [x,w]=qd_pts(mu_ut_sf,P_ut_sf);
    Pcc=0;
    for i=1:1:length(w)
        Pcc=Pcc+w(i)*(x(i,:)'-mu_ut_sf)*(hx(x(i,:)')-mu_ut_of)';
    end
    %kalman gain
    K=Pcc/P_ut_of;
    %update
    xu=mu_ut_sf+K*(ym-mu_ut_of);
    Pu=P_ut_sf-K*P_ut_of*K';

end