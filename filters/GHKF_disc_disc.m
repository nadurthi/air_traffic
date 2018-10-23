function [xu,Pu]=GHKF_disc_disc(fn,hn,mu_ut,P_ut,ym,Q,R,NN)

    [x,w]=conjugate_dir_gausspts(mu_ut,P_ut);
%     [x,w]=GH_points(mu_ut,P_ut,NN);
    [mu_ut_sf,P_ut_sf]=prop_mean_cov_points_discr(x,w,fn);
    P_ut_sf=P_ut_sf+Q;
    
    %obs forecast
        if ym==-1
            
        xu=mu_ut_sf;
        Pu=P_ut_sf;
         
        else
              CKF=P_ut_sf;
              [x,w]=conjugate_dir_gausspts(mu_ut_sf,P_ut_sf);
%      [x,w]=GH_points(mu_ut_sf,P_ut_sf,NN);
    [mu_ut_of,P_ut_of]=prop_mean_cov_points_discr(x,w,hn);
    P_ut_of=P_ut_of+R;
    %cross cov
    [x,w]=conjugate_dir_gausspts(mu_ut_sf,P_ut_sf);
%     [x,w]=GH_points(mu_ut_sf,P_ut_sf,NN);
    Pcc=0;
    for i=1:1:length(w)
        Pcc=Pcc+w(i)*(x(i,:)'-mu_ut_sf)*(hn(x(i,:)')-mu_ut_of)';
    end
    %kalman gain
    K=Pcc*inv(P_ut_of);
    %update
  
    xu=mu_ut_sf+K*(ym-mu_ut_of);
    Pu=P_ut_sf-K*P_ut_of*K';
    
    end
end