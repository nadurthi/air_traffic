function OUTPUT=Allfilters(model,filter,pf,time)
% have to still make it general : like taking in the control input , model
% takes time value into propagation etc. Needds to be revised. have to take
% the dt into the model dyns i.e here T

% global T

% T=time.dt;

%% __-------------------------------------------------------
%TRUTH INITIALIZATION of the model and filters
%     x0tr=model.x0tr;
    
% FIlter INIitalization
x_0f=filter.x0_filt_start;
P_0f=filter.P0_filt_start;
%     P0tr_UM=diag([50^2,50^2,100,100]);
 
%% %%%%%%%%%%%%START FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  NNN=1; %% Simulate over NNN iterations


xNNN_mc=zeros(time.nSteps,model.fn);
YNNN_mc=zeros(time.nSteps,model.hn);
xNNN_ukf=zeros(time.nSteps,model.fn);
xNNN_ckf=zeros(time.nSteps,model.fn);
xNNN_cut4=zeros(time.nSteps,model.fn);
xNNN_cut6=zeros(time.nSteps,model.fn);
xNNN_cut8=zeros(time.nSteps,model.fn);
xNNN_gh=zeros(time.nSteps,model.fn);
xNNN_mupf=zeros(time.nSteps,model.fn);
xNNN_mopf=zeros(time.nSteps,model.fn);
xNNN_ekf=zeros(time.nSteps,model.fn);
xNNN_gmm=zeros(time.nSteps,model.fn);

PNNN_ukf=zeros(time.nSteps,model.fn^2);
PNNN_ckf=zeros(time.nSteps,model.fn^2);
PNNN_cut4=zeros(time.nSteps,model.fn^2);
PNNN_cut6=zeros(time.nSteps,model.fn^2);
PNNN_cut8=zeros(time.nSteps,model.fn^2);
PNNN_gh=zeros(time.nSteps,model.fn^2);
PNNN_pf=zeros(time.nSteps,model.fn^2);
PNNN_ekf=zeros(time.nSteps,model.fn^2);
PNNN_gmm=zeros(time.nSteps,model.fn^2);

histdata=cell(time.nSteps,2);

% for j=1:1:NNN
    
    %%-------------------------------------------------------------
    %%monte carlo truth generation using dynamics
    %%_-----------------------------------------------------------
%     [t,x_mc1]=ode45_discc(@KIRB_UM_eg_dyn_disc,time.t0,time.dt,125,x0tr,1e-200);
%     [t,x_mc2]=ode45_discc(@KIRB_CT_eg_dyn_disc,t(end),time.dt,t(end)+90,[x_mc1(end,:)';1*pi/180],1e-200);
%     [t,x_mc3]=ode45_discc(@KIRB_UM_eg_dyn_disc,t(end),time.dt,t(end)+125,x_mc2(end,1:end-1)',1e-200);
%     [t,x_mc4]=ode45_discc(@KIRB_CT_eg_dyn_disc,t(end),time.dt,t(end)+30,[x_mc3(end,:)';-3*pi/180],1e-200);
%     [t,x_mc5]=ode45_discc(@KIRB_UM_eg_dyn_disc,t(end),time.dt,t(end)+125,x_mc4(end,1:end-1)',1e-200);
%     x_mc=[[x_mc1(1:end-1,1:4),zeros(size(x_mc1,1)-1,1)];x_mc2(1:end-1,1:5);[x_mc3(1:end-1,1:4),zeros(size(x_mc3,1)-1,1)];x_mc4(1:end-1,1:5);[x_mc5(1:end,1:4),zeros(size(x_mc5,1),1)]];
%     
%     %generating the measurement
%     ym=zeros(size(x_mc,1),model.hn);
%     for i=1:1:time.nSteps
%         ym(i,:)=(model.hx(x_mc(i,:)')+model.sR*randn(model.hn,1))';
%     end
    xNNN_mc(:,:)=filter.truth;
    YNNN_mc(:,:)=filter.ymeas;
    ym=filter.ymeas;
%         figure(1)
% plot(x_mc(:,1),x_mc(:,2),ym(:,1),ym(:,2),'k*')
%% ------------------------------------------------------------------------

    
    % UKF filter
    mu_ut=x_0f;
    P_ut=P_0f;
    xNNN_ukf(1,:)=mu_ut';
    PNNN_ukf(1,:)=reshape(P_ut,1,model.fn^2);
    
    % CKF filter
    mu_ckf=x_0f;
    P_ckf=P_0f;
    xNNN_ckf(1,:)=mu_ckf';
    PNNN_ckf(1,:)=reshape(P_ckf,1,model.fn^2);
 
    % CUT4 filter
    mu_cut4=x_0f;
    P_cut4=P_0f;
    xNNN_cut4(1,:)=mu_cut4';
    PNNN_cut4(1,:)=reshape(P_cut4,1,model.fn^2);
    
    % CUT6 filter
    mu_cut6=x_0f;
    P_cut6=P_0f;
    xNNN_cut6(1,:)=mu_cut6';
    PNNN_cut6(1,:)=reshape(P_cut6,1,model.fn^2);
    
    % CUT8 filter
    mu_cut8=x_0f;
    P_cut8=P_0f;
    xNNN_cut8(1,:)=mu_cut8';
    PNNN_cut8(1,:)=reshape(P_cut8,1,model.fn^2);
    
    % GMM filter
    if strcmp(filter.GMMF ,'true')
        mu_gmm=filter.GMM.mean;
        P_gmm=filter.GMM.cov;
        xNNN_gmm(1,:)=mu_gmm';
        PNNN_gmm(1,:)=reshape(P_gmm,1,model.fn^2);
        GMM=filter.GMM;
    end
    % ghKF filter
    mu_gh=x_0f;
    P_gh=P_0f;
    xNNN_gh(1,:)=mu_gh';
    PNNN_gh(1,:)=reshape(P_gh,1,model.fn^2);
       
       
       % PF
       mu_pf=x_0f;
       mo_pf=x_0f;
       P_pf=P_0f;
       PNNN_pf(1,:)=reshape(P_pf,1,model.fn^2);
       xNNN_mupf(1,:)=mu_pf';
       xNNN_mopf(1,:)=mo_pf';
       X_pf=repmat(mu_pf,1,pf.no_particles)+sqrtm(P_pf)*randn(model.fn,pf.no_particles);
       w_pf = ones(1, pf.no_particles) / pf.no_particles;
       [PF.x{1},PF.P{1},PF.minB{1},PF.maxB{1},tmp_X_pf] = getPFdata(X_pf,w_pf);
       histdata{1,1}=tmp_X_pf;
       histdata{1,2}=w_pf;
        %EKF
    mu_ekf=x_0f;
    P_ekf=P_0f;
    xNNN_ekf(1,:)=mu_ekf';
    PNNN_ekf(1,:)=reshape(P_ekf,1,model.fn^2); 
        
        
  for k=2:1:time.nSteps
disp([num2str(k),' of ',num2str(time.nSteps)] )
        
        if rem(k,filter.freq)==0
            zm=ym(k,:)';
        else
            zm=-1234; % some non-measurement number
        end
        
        %run filters that are only switched on
          if strcmp(filter.UKF ,'true')
              %tic
      [mu_ut,P_ut]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_ut,P_ut,zm,'ut',filter.paras_ukf_kappa);
              %time_ukf=toc;
          end
          
          if strcmp(filter.CKF ,'true')
              %tic
      [mu_ckf,P_ckf]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_ckf,P_ckf,zm,'ckf',0);
          %time_ckf=toc;
          end
          
          if strcmp(filter.CUT4KF ,'true')
              %tic
      [mu_cut4,P_cut4]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_cut4,P_cut4,zm,'cut4',0);
          %time_cut4=toc;
          end
          
          if strcmp(filter.CUT6KF ,'true')
              %tic
      [mu_cut6,P_cut6]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_cut6,P_cut6,zm,'cut6',0);
          %time_cut6=toc;
          end
          
          if strcmp(filter.CUT8KF ,'true')
              %tic
      [mu_cut8,P_cut8]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_cut8,P_cut8,zm,'cut8',0);
          %time_cut8=toc;
          end
          
          if strcmp(filter.GHKF ,'true')
              %tic
      [mu_gh,P_gh]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_gh,P_gh,zm,'gh',filter.paras_gh_pts); 
          %time_gh=toc;
          end
            
        if strcmp(filter.GMMF ,'true')
            %tic
%             [mu_gh,P_gh]=QUADpt_filter_disc_UPDT_disc_MEAS(model,mu_gh,P_gh,zm,'gh',filter.paras_gh_pts); 
            GMM=prior_prop_GMM(GMM,model);
            GMM=post_GMM(GMM,model,zm);
            GMM.w
            GMM.mu
            [mu_gmm,P_gmm]=GMM_mean_cov(GMM);
            %time_gh=toc;
        end

          
          
          if strcmp(filter.PF ,'true')
              %tic
      [X_pf,w_pf]=BOOTSTRAP_Pfilter_disc_UPDT_disc_MEAS(model,pf,X_pf,w_pf,zm);
          %time_pf=toc;
          end
          
          if strcmp(filter.EKF ,'true')
              %tic
      [mu_ekf,P_ekf]=EKF_disc(model,mu_ekf,P_ekf,zm);
          %time_ekf=toc;
          end
      
      %storing data for this run
    xNNN_ukf(k,:)=mu_ut';
    PNNN_ukf(k,:)=reshape(P_ut,1,model.fn^2);
      
    xNNN_ckf(k,:)=mu_ckf';
    PNNN_ckf(k,:)=reshape(P_ckf,1,model.fn^2);
   
    xNNN_cut4(k,:)=mu_cut4';
    PNNN_cut4(k,:)=reshape(P_cut4,1,model.fn^2);
    
    xNNN_cut6(k,:)=mu_cut6';
    PNNN_cut6(k,:)=reshape(P_cut6,1,model.fn^2);
    
    xNNN_cut8(k,:)=mu_cut8';
    PNNN_cut8(k,:)=reshape(P_cut8,1,model.fn^2);

    xNNN_gh(k,:)=mu_gh';
    PNNN_gh(k,:)=reshape(P_gh,1,model.fn^2);
    
    if strcmp(filter.GMMF ,'true')
        xNNN_gmm(k,:)=mu_gmm';
        PNNN_gmm(k,:)=reshape(P_gmm,1,model.fn^2);
    end
    
    xNNN_ekf(k,:)=mu_ekf';
    PNNN_ekf(k,:)=reshape(P_ekf,1,model.fn^2);
    
      
    [PFx,P_pf,PFminB,PFmaxB,tmp_X_pf] = getPFdata(X_pf, w_pf);
    histdata{k,1}=tmp_X_pf;
    histdata{k,2}=w_pf;
    
    xNNN_mupf(k,:)=PFx';
    PNNN_pf(k,:)=reshape(P_pf,1,model.fn^2);
    
    MD=zeros(model.fn,1);
    for mo=1:1:model.fn
    rang=PFminB(mo):(PFmaxB(mo)-PFminB(mo))/(pf.no_bins-1):PFmaxB(mo);
    if numel(rang)==0
       MD(mo,1)=PFminB(mo);
        continue
    end
    nb=histc(tmp_X_pf(mo,:),rang);
    [lkll,ind]=max(nb);

    MD(mo,1)=rang(ind);
    end
    
    xNNN_mopf(k,:)=MD';
    
    
    if strcmpi(pf.plotHists,'Yes')
        disp('plotting histograms')
        figure(1)
       for ffn=1:1:model.fn
           subplot(ceil(model.fn/3),3,ffn)
          histNorm({tmp_X_pf(ffn,:)}) 
          title(num2str(ffn))
       end
       pause(0.2)
       saveas(gcf,strcat('PFhist_',pf.plotHiststagname,'k_',num2str(k)),'png')
        
    end
  end
 OUTPUT.mc=xNNN_mc;
 OUTPUT.ymeas=YNNN_mc;
  OUTPUT.mopf=xNNN_mopf;
  OUTPUT.mupf=xNNN_mupf;
  OUTPUT.Ppf=PNNN_pf;
   OUTPUT.mu_gh=xNNN_gh;
   OUTPUT.mu_ekf=xNNN_ekf;
   OUTPUT.mu_gmm=xNNN_gmm;
   
   OUTPUT.mu_cut8=xNNN_cut8;
  OUTPUT.mu_cut6=xNNN_cut6;
  OUTPUT.mu_cut4=xNNN_cut4;
  OUTPUT.mu_ckf=xNNN_ckf;
  OUTPUT.mu_ukf=xNNN_ukf;
    OUTPUT.P_gh=PNNN_gh;
    OUTPUT.P_ekf=PNNN_ekf;
    OUTPUT.P_gmm=PNNN_gmm;
    
  OUTPUT.P_cut8=PNNN_cut8;
  OUTPUT.P_cut6=PNNN_cut6;
  OUTPUT.P_cut4=PNNN_cut4;
  OUTPUT.P_ckf=PNNN_ckf;
  OUTPUT.P_ukf=PNNN_ukf;
  if strcmp(filter.save_pf_data ,'true')
  OUTPUT.PFdata=histdata;
  else
      OUTPUT.PFdata=[];
  end

end

% % Averaging over all the runs
% est_fin_ukf=zeros(time.nSteps,3);
% est_fin_ckf=zeros(time.nSteps,3);
% est_fin_cut4=zeros(time.nSteps,3);
% est_fin_cut6=zeros(time.nSteps,3);
% est_fin_cut8=zeros(time.nSteps,3);
% est_fin_gh=zeros(time.nSteps,3);
% est_fin_mupf=zeros(time.nSteps,3);
% est_fin_mopf=zeros(time.nSteps,3);
% 
% for r=1:1:time.nSteps
%         for c=1:1:3
%            if c==1
%            est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_ukf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            est_fin_ckf(r,c) =sqrt(mean((xNNN_ckf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_ckf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            est_fin_cut4(r,c)=sqrt(mean((xNNN_cut4(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_cut4(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            est_fin_cut6(r,c)=sqrt(mean((xNNN_cut6(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_cut6(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            est_fin_cut8(r,c)=sqrt(mean((xNNN_cut8(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_cut8(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            est_fin_gh(r,c)=  sqrt(mean((xNNN_gh(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_gh(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            est_fin_mupf(r,c)=  sqrt(mean((xNNN_mupf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_mupf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            est_fin_mopf(r,c)=  sqrt(mean((xNNN_mopf(r,c,:)-xNNN_mc(r,c,:)).^2+(xNNN_mopf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2));
%            end
%            if c==2
%            est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_ukf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            est_fin_ckf(r,c) =sqrt(mean((xNNN_ckf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_ckf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            est_fin_cut4(r,c)=sqrt(mean((xNNN_cut4(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut4(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            est_fin_cut6(r,c)=sqrt(mean((xNNN_cut6(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut6(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            est_fin_cut8(r,c)=sqrt(mean((xNNN_cut8(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_cut8(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            est_fin_gh(r,c)=  sqrt(mean((xNNN_gh(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_gh(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            est_fin_mupf(r,c)=  sqrt(mean((xNNN_mupf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_mupf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            est_fin_mopf(r,c)=  sqrt(mean((xNNN_mopf(r,c+1,:)-xNNN_mc(r,c+1,:)).^2+(xNNN_mopf(r,c+2,:)-xNNN_mc(r,c+2,:)).^2));
%            end 
%            if c==3
%            est_fin_ukf(r,c)= sqrt(mean((xNNN_ukf(r,5,:)-xNNN_mc(r,5,:)).^2));
%            est_fin_ckf(r,c) =sqrt(mean((xNNN_ckf(r,5,:)-xNNN_mc(r,5,:)).^2));
%            est_fin_cut4(r,c)=sqrt(mean((xNNN_cut4(r,5,:)-xNNN_mc(r,5,:)).^2));
%            est_fin_cut6(r,c)=sqrt(mean((xNNN_cut6(r,5,:)-xNNN_mc(r,5,:)).^2));
%            est_fin_cut8(r,c)=sqrt(mean((xNNN_cut8(r,5,:)-xNNN_mc(r,5,:)).^2));
%            est_fin_gh(r,c)=  sqrt(mean((xNNN_gh(r,5,:)-xNNN_mc(r,5,:)).^2));
%            est_fin_mupf(r,c)=  sqrt(mean((xNNN_mupf(r,5,:)-xNNN_mc(r,5,:)).^2));
%            est_fin_mopf(r,c)=  sqrt(mean((xNNN_mopf(r,5,:)-xNNN_mc(r,5,:)).^2));
%            end 
%         end
% end
% save(model.name,'model','time','pf','filter_paras','histdata','est_fin_ukf','est_fin_ckf','est_fin_cut4','est_fin_cut6','est_fin_cut8','est_fin_gh','est_fin_mupf','est_fin_mopf','xNNN_mc','YNNN_mc','xNNN_ukf','xNNN_ckf','xNNN_cut4','xNNN_cut6','xNNN_cut8','xNNN_gh','xNNN_mupf','xNNN_mopf','PNNN_ukf','PNNN_ckf','PNNN_cut4','PNNN_cut6','PNNN_cut8','PNNN_gh','PNNN_pf')
% 




% %% Plotting
% t=time.tspan;
% figure(1)
% plot(t,est_fin_mupf(:,1),'-.',t,est_fin_mopf(:,1),':',t,est_fin_ckf(:,1),t,est_fin_ukf(:,1),t,est_fin_cut4(:,1),t,est_fin_cut6(:,1),t,est_fin_cut8(:,1),t,est_fin_gh(:,1))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% figure(2)
% plot(t,est_fin_mupf(:,2),'-.',t,est_fin_mopf(:,2),':',t,est_fin_ckf(:,2),t,est_fin_ukf(:,2),t,est_fin_cut4(:,2),t,est_fin_cut6(:,2),t,est_fin_cut8(:,2),t,est_fin_gh(:,2))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% figure(3)
% plot(t,(180/pi)*est_fin_mupf(:,3),'-.',t,(180/pi)*est_fin_mopf(:,3),':',t,(180/pi)*est_fin_ckf(:,3),t,(180/pi)*est_fin_ukf(:,3),t,(180/pi)*est_fin_cut4(:,3),t,(180/pi)*est_fin_cut6(:,3),t,(180/pi)*est_fin_cut8(:,3),t,(180/pi)*est_fin_gh(:,3))
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh')
% figure(4)
% plot(xNNN_mupf(:,1,end),xNNN_mupf(:,2,end),'-.',xNNN_mopf(:,1,end),xNNN_mopf(:,2,end),':',xNNN_ckf(:,1,end),xNNN_ckf(:,2,end),xNNN_ukf(:,1,end),xNNN_ukf(:,2,end),xNNN_cut4(:,1,end),xNNN_cut4(:,2,end),xNNN_cut6(:,1,end),xNNN_cut6(:,2,end),xNNN_cut8(:,1,end),xNNN_cut8(:,2,end),xNNN_gh(:,1,end),xNNN_gh(:,2,end),xNNN_mc(:,1,end),xNNN_mc(:,2,end),'--',YNNN_mc(:,1,end),YNNN_mc(:,2,end),'k*','LineWidth',2)
% legend('PF-mu','PF-mo','ckf','ukf','cut4','cut6','cut8','gh','MC','sensor')



