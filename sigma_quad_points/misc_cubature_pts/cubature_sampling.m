clear
clc
N=10;
X=[];
for i=1:1:N
    p=rand(1,2);
    % now make it fully symmetric
    X=vertcat(X,[p(1),p(2)]);
    X=vertcat(X,[-p(1),-p(2)]);
    X=vertcat(X,[-p(1),p(2)]);
    X=vertcat(X,[p(1),-p(2)]);
    X=vertcat(X,[p(2),p(1)]);
     X=vertcat(X,[-p(2),-p(1)]);
      X=vertcat(X,[-p(2),p(1)]);
       X=vertcat(X,[p(2),-p(1)]);
end
plot(X(:,1),X(:,2),'bo')
w=ones(length(X),1)/length(X);
grid
axis([-1,1,-1,1])
axis square
m1=cal_moments_wrt_pts(X,w,4,eye(2))
X=rand(size(X));
X=-1+2*X;
m2=cal_moments_wrt_pts(X,w,4,eye(2))
