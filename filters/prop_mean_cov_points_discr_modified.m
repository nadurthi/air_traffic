function [mk,Pk]=prop_mean_cov_points_discr_modified(x,w,F,para_dt,model)


%%%%%%%%%discrete version%%%%%%%%%%


% x has columns= dim
% and rows as number of samples
[N,n]=size(x);
%n is dim of system
%N no. of samples

% F gives out a column vector of dim m
% Y=zeros(size(x));
Pk=[];
mk=[];

Y=[];
if strcmp(model.dynamics,'discrete')
    
    for i=1:1:N
%         keyboard
        Y(i,:)=F(x(i,:)',model.fx_dt);
    end
else
    for i=1:1:N
        [tt,xx]=ode45(@model.fx,[0,model.dt],x(i,:)');
        Y(i,:)=xx(end,:);
    end
    
end
% [~,m]=size(Y);
% propagating the mean
% mk=zeros(m,1);
% for i=1:1:N
%     mk=mk+w(i)*Y(i,:)';
% end

[N,n]=size(Y);
W=repmat(w,1,n);
mk=sum(W.*Y,1)';


MU=repmat(mk',N,1);
X=Y-MU;
Pk=X'*(W.*X);


end