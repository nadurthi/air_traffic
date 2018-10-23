

%% ------------------------------------------------------------------------
% time
%--------------------------------------------------------------------------
time.t0 = 0;
time.dt = 0.1;
time.tf = 10;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);
%_________________________________
%% ------------------------------------------------------------------------
% model
T=time.dt;

Qf=1*diag([1,1,1,1,1]);

Qt=0*diag([0,0,2.4064*1e-5,2.4064*1e-5,0]);

% sigr=100*1e-3;
% sigth=10*1e-3;

% R=diag([sigr^2,sigth^2]);
% R=diag([sigth^2]);



model.fn = 5;               % state space dimensionality
model.fx = @lorenz_cont;
model.fx_dt =time.dt;
model.dynamics='continuous';

R=diag([0.5,0.5,0.5]);
model.hn =3;               % measurement dimensionality
model.hx =@(x,dt)[x(1);x(2);x(3)];
model.hx_jac=@(x,dt)[eye(model.hn),zeros(3,2)];


model.gx=@(x,dt)1;
model.Q =Qf;
model.sQ =sqrtm(Qf);
model.R = R;
model.sR=sqrtm(R);
model.Qtruth=Qt;
model.Qt_sq=sqrtm(Qt);
model.dt=time.dt;

model.propagate=1;
model.para_dt=0;

% model.x0tr=[6500.4,349.14,-1.8093,-6.7967,0.6932]';
% model.P0tr=diag([1e-5,1e-5,1e-5,1e-5,0]);
% [1.50887,-1.531271,25.46091,12,25]';
model.x0tr=[1.50887,-1.531271,25.46091,10,28]';
model.P0tr=diag([4,4,4,2,4]);


%% ------------------------------------------------------------------------
% particle filter settings
%--------------------------------------------------------------------------
pf.no_particles =100;
pf.no_bins = 100;               % when computing the histogram

pf.resample = true;             % resampling after measurement update
pf.neff = 0.8;                  % treshold for the number of effective samples

filter.paras_ukf_kappa=1;
filter.paras_gh_pts=2;

filter.freq=1; %% This is actually the number of 'dt' steps after which a meas updt is done

% start the filter too at the true estimate
x0f=model.x0tr;
P0f=model.P0tr;

filter.x0_filt_start=x0f;
filter.P0_filt_start=P0f;

%% Generate measurements

% use the model's true initial state to generate measurements
x0trr = model.x0tr;%mvnrnd(model.x0tr', model.P0tr);
[t,x_mc]=ode45(model.fx,time.tspan,x0trr',model.Qt_sq);
ym=zeros(time.nSteps,model.hn);
for ii=1:1:time.nSteps
    ym(ii,:)=(model.hx(x_mc(ii,:)')+(sqrtm(R)*randn(model.hn,1)))';
end

% start the filter from random point
filter.x0_filt_start = mvnrnd(model.x0tr', model.P0tr)';


%     filter.x0_filt_start =[1.50887,-1.531271,25.46091,12,25]';
filter.ymeas=ym;


%% running
% FIlter INIitalization

mu_pf=filter.x0_filt_start;
P_pf=filter.P0_filt_start;

% initial condition for the filter
X_pf=repmat(mu_pf,1,pf.no_particles)+sqrtm(P_pf)*randn(model.fn,pf.no_particles);
w_pf = ones(1, pf.no_particles) / pf.no_particles;



xNNN_mupf=zeros(time.nSteps,model.fn);
PNNN_pf=zeros(time.nSteps,model.fn^2);


for k=2:1:time.nSteps
    k
    if rem(k,filter.freq)==0
        zm=filter.ymeas(k,:)';
    else
        zm=-1234;
    end
    
    [X_pf,w_pf]=BOOTSTRAP_Pfilter_disc_UPDT_disc_MEAS(model,pf,X_pf,w_pf,zm);
    
    [PFx,P_pf,PFminB,PFmaxB,tmp_X_pf] = getPFdata(X_pf, w_pf);
    
    xNNN_mupf(k,:)=PFx';
    PNNN_pf(k,:)=reshape(P_pf,1,model.fn^2);
    
end

plot(time.tspan,xNNN_mupf(:,1),time.tspan,xNNN_mupf(:,1)+sqrt(PNNN_pf(:,1)),time.tspan,xNNN_mupf(:,1)-sqrt(PNNN_pf(:,1)) )