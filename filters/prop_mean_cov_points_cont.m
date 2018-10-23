function [mk,Pk]=prop_mean_cov_points_cont(x,w,F,ti,dt,tf)


%%%%%%%%%continuous version%%%%%%%%%%


% x has columns= dim
% and rows as number of samples
[N,n]=size(x);
%n is dim of system
%N no. of samples
options=odeset('RelTol',1e-15,'AbsTol',1e-15);
% F gives out a column vector of dim n 
for i=1:1:N
    [t,xt]=ode45(F,ti:dt:tf,x(i,:)');
    Y(i,:)=xt(end,:);
end

% propagating the mean
mk=zeros(n,1);
for i=1:1:N
    mk=mk+w(i)*Y(i,:)';
end

%propagating the cov
Pk=0;
for i=1:1:N
    Pk=Pk+w(i)*(Y(i,:)'-mk)*(Y(i,:)'-mk)';
end
