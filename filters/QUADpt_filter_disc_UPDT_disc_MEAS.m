function [xu,Pu]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_ut,P_ut,ym,method,para)

fx=model.fx;
 

Q=model.Q;

n=model.fn;
para_dt=model.para_dt;
global kappa
    kappa=para;
    switch lower(method)
   case {'ckf'}
      qd_pts=@cubature_KF_points;
   case 'ut'
      qd_pts=@(m,P)UT_sigmapoints(m,P,2);
   case 'cut4'
      qd_pts=@conjugate_dir_gausspts;
   case 'cut6'
      qd_pts=@conjugate_dir_gausspts_till_6moment_scheme2;
   case 'cut8'
      qd_pts=@conjugate_dir_gausspts_till_8moment;
   case 'gh'
      qd_pts=@(m,P)GH_pts(m,P,para);
   otherwise
      error('smthg is wrong: DONT ask me what')
    end
%     mu_ut
%     diag(P_ut)

    [x,w]=qd_pts(mu_ut,P_ut);
    
    [mu_ut_sf,P_ut_sf]=prop_mean_cov_points_discr_modified(x,w,fx,para_dt,model);

      P_ut_sf=P_ut_sf+Q;
      
%    Check if meas is avail at this time step
        if ym==-1234
        xu=mu_ut_sf;
        Pu=P_ut_sf;
        return
        end
     R=model.R;  
     hx=model.hx;
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
%     keyboard
    xu=mu_ut_sf+K*(ym-mu_ut_of);
    Pu=P_ut_sf-K*P_ut_of*K';


end