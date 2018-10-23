function [x,w]=cubature_KF_points(mu,P)
n=length(P);
A=sqrtm(P);
x=zeros(2*n,n);
r=sqrt(2*n/2);
for i=1:1:n
    x(i,:)=r*A(:,i)';
    x(i+n,:)=-r*A(:,i)';
end
w=(1/(2*n))*ones(2*n,1);
for i=1:1:2*n
    x(i,:)=x(i,:)+mu';
end
end
