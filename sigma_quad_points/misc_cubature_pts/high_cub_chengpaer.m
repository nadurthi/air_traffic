function [x,w]=high_cub_chengpaer(mu,P)
n=length(mu);
x=[];
w=[];
% n is dimension
sp=[];
sm=[];
e=eye(n);
for k=1:1:n
    for l=1:1:n
        if k<l
        sp=vertcat(sp,sqrt(0.5)*(e(k,:)+e(l,:)));
        sm=vertcat(sm,sqrt(0.5)*(e(k,:)-e(l,:)));
        end
    end
end
x=zeros(1,n);
w=2/(n+2);

w=vertcat(w,1/(n+2)^2*ones(n*(n-1)/2,1));
x=vertcat(x,sqrt(n+2)*sp);

w=vertcat(w,1/(n+2)^2*ones(n*(n-1)/2,1));
x=vertcat(x,-sqrt(n+2)*sp);

w=vertcat(w,1/(n+2)^2*ones(n*(n-1)/2,1));
x=vertcat(x,sqrt(n+2)*sm);

w=vertcat(w,1/(n+2)^2*ones(n*(n-1)/2,1));
x=vertcat(x,-sqrt(n+2)*sm);

w=vertcat(w,(4-n)/(2*(n+2)^2)*ones(n,1));
x=vertcat(x,sqrt(n+2)*e);


w=vertcat(w,(4-n)/(2*(n+2)^2)*ones(n,1));
x=vertcat(x,-sqrt(n+2)*e);

A=sqrtm(P);
mu=mu(:);
for i=1:1:length(w)
    x(i,:)=A*x(i,:)'+mu;
end
