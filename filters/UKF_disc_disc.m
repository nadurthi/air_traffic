function [xu,Pu]=UKF_disc_disc(fn,hn,mu_ut,P_ut,ym,Q,R)
%     size(mu_ut)
%     size(P_ut)
    [x,w]=UT_sigmapoints(mu_ut,P_ut,2);
    [mu_ut_sf,P_ut_sf]=prop_mean_cov_points_discr(x,w,fn);
    P_ut_sf=P_ut_sf+Q;
    
    %obs forecast
        if ym==-1
            
        xu=mu_ut_sf;
        Pu=P_ut_sf;
         
       else
     [x,w]=UT_sigmapoints(mu_ut_sf,P_ut_sf,2);
    [mu_ut_of,P_ut_of]=prop_mean_cov_points_discr(x,w,hn);
    P_ut_of=P_ut_of+R;
    %cross cov
    [x,w]=UT_sigmapoints(mu_ut_sf,P_ut_sf,2);
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
