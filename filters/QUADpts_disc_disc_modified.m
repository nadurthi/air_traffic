function [xu,Pu]=QUADpts_disc_disc_modified(fn,hn,mu_ut,P_ut,ym,Q,R,method,para)
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
      qd_pts=@(m,P)GH_points(m,P,para);
   otherwise
      error('smthg is wrong: DONT ask me what')
    end
    n=length(mu_ut);
    [x,w]=qd_pts(mu_ut,P_ut);
    
    [mu_ut_sf,P_ut_sf]=prop_mean_cov_points_discr_modified(x,w,fn);

      P_ut_sf=P_ut_sf+Q;
      
   
        if ym==-1234
        xu=mu_ut_sf;
        Pu=P_ut_sf;
         return
        end
         %obs forecast
     [x,w]=qd_pts(mu_ut_sf,P_ut_sf);
    [mu_ut_of,P_ut_of]=prop_mean_cov_points_discr_modified(x,w,hn);
%     mu_ut_of=mu_ut_of';
%     mm=length(mu_ut_of);
%     P_ut_of=reshape(P_ut_of,mm,mm);
    
    P_ut_of=P_ut_of+R;
%     keyboard
    %cross cov
    [x,w]=qd_pts(mu_ut_sf,P_ut_sf);
    Pcc=0;
    for i=1:1:length(w)
        Pcc=Pcc+w(i)*(x(i,:)'-mu_ut_sf)*(hn(x(i,:)')-mu_ut_of)';
    end
    %kalman gain
    K=Pcc*inv(P_ut_of);
    %update
  
    xu=mu_ut_sf+K*(ym-mu_ut_of);
    Pu=P_ut_sf-K*P_ut_of*K';
% %     inter_mu_ut_sf(end,:)=XXu';
%     inter_P_ut_sf(end,:)=reshape(PPu,1,n*n);
%     xu=inter_mu_ut_sf;
%     Pu=inter_P_ut_sf;
%     end
end