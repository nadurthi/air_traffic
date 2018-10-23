function [mk,Pk]=prop_mean_cov_points_discr(x,w,F,f)


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
for j=1:1:f
    Y=[];
for i=1:1:N
    Y(i,:)=F(x(i,:)');
end

% [~,m]=size(Y);
% propagating the mean
% mk=zeros(m,1);
% for i=1:1:N
%     mk=mk+w(i)*Y(i,:)';
% end
[N,n]=size(Y);
    W=repmat(w,1,n);
    inter_mk=sum(W.*Y,1)';
    mk=vertcat(mk,inter_mk');

MU=repmat(inter_mk',N,1);
X=Y-MU;
Pk=vertcat(Pk,reshape(X'*(W.*X),1,n*n));
x=Y;
[N,n]=size(x);
end
end