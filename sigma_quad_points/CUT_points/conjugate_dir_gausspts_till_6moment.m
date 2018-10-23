function [X,w]=conjugate_dir_gausspts_till_6moment(mu,P)
%dimension n of the system
n=length(mu);
% r1=1.5;
% w1=0.2;
% r2=2;
% w2=0.02;
% r3=3;
% w3=0.002;
% r=2.4142;
% the number of points in this scheme are 2n+2^n+1 points
% x0=[r1,w1,r2,w2,r3,w3,r];
r=2.4142;
if n==2
x0=[1.7,0.2,3,0.02,5,0.002];
 x=fsolve(@D2sys,x0);
f=D2sys(x)
r1=x(1);
w1=x(2);
r2=x(3);
w2=x(4);
r3=x(5);
w3=x(6);

end
if n==3
%     r=2.7320;
%     options=optimset('MaxFunEvals',5000000,'MaxIter',10000);
% x0=[3,0.3/6,4,0.3/8,5,0.3/24];%,3];
% %  x=fsolve(@D3sys,x0,options);
%  x=fmincon(@(x)(2*x(1)^8*x(2)+8*x(3)^8*x(4)+8*x(5)^8*x(6)*(2+r^8)-105)^2,x0,[],[],[],[],[0,0,0,0,0,0],[50,1,50,1,50,1],@D3sys,options);
% r1=x(1);
% w1=x(2);
% r2=x(3);
% w2=x(4);
% r3=x(5);
% w3=x(6);
% %  r=x(7);
r=2;
r1=1.6881059940166;
r2=1.545223430876;
r3=2.4158283260926;
w1=0.111809;
w2=0.0178901;
w3=0.0000750444;

end
if n==4
x0=[sqrt(3),0.2,3,0.02,5,0.002,2.4142];
x=fsolve(@D4sys,x0);
r1=x(1);
w1=x(2);
r2=x(3);
w2=x(4);
r3=x(5);
w3=x(6);
r=x(7);
end
[u,s,v]=svd(P);
A=chol(s)*u;
A=A';
A=sqrtm(P);
X=zeros(2*n+2^n+n*2^n+1,n);
index=GenerateIndex(n,n*ones(1,n));
[roww,coll]=size(index);
dr=[];
for i=1:1:roww
if length(find(index(i,:)>2))==0
    dr=vertcat(dr,index(i,:));
end
end
for i=1:1:n+1
    if i==1
        X(i,:)=zeros(1,n);
    else
              
    X(i,:)=r1*A(:,i-1);
    X(i+n,:)=-r1*A(:,i-1);
        
    end
end
mo=-1*ones(1,n);
for i=1:1:2^n
    rr=mo.^dr(i,:);
    sig=0;
    for j=1:1:n
        sig=sig+rr(j)*A(:,j);
    end
    X(2*n+1+i,:)=r2*sig;
end

% for j=0:1:n-1
% for i=1:1:2^n
%     X(2*n+1+2^n+i+j*2^n,:)=X(2*n+1+i,:).*[ones(1,j),r,ones(1,n-j-1)];
% end
% end
mo=-1*ones(1,n);
p=0;
    for i=1:1:n*2^n
        p=p+1;
        
        rr=mo.^dr(p,:);
        if rem(i,2^n)==0
            p=0;
        end
        k=ceil(i/2^n)-1;
        pp=[ones(1,k),r,ones(1,n-k-1)];
    sig=0;
    for j=1:1:n
        sig=sig+pp(j)*rr(j)*A(:,j);
    end
     X(2*n+1+2^n+i,:)=r3*sig;
    end
w0=1-2*n*w1-2^n*w2-n*2^n*w3
r1
w1
r2
w2
r3
w3
r
w=[w0,w1*ones(1,2*n),w2*ones(1,2^n),w3*ones(1,n*2^n)]';
for i=1:1:n
    X(:,i)=X(:,i)+mu(i);
end
