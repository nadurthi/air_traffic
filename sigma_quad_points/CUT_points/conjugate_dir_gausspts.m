function [X,w]=conjugate_dir_gausspts(mu,P)
%dimension n of the system

n=length(mu);
% the number of points in this scheme are 2n+2^n+1 points
%x=[r1,w1,r2,w2]
% [x,fval]=fsolve(@(x)[2*x(1)^2*x(2)+2^n*x(3)^2*x(4)-1,2*x(1)^4*x(2)+2^n*x(3)^4*x(4)-3,2^n*x(3)^4*x(4)-1,2^n*x(3)^6*x(4)-3],x0);
% fval
% n=2;
% x=fmincon(@(x)(2*x(1)^2+x(2)^2-3)^2,[2,4],[],[],[],[],[0.1,0.1],[20,20],@moment_4th_ND_eqns);
if n==2 || n==1
    if n==1
        w0=0.5811010092660772;
        w1=0.20498484723245053;
        w2=0.00446464813451093;
        r1=1.4861736616297834;
        r2=3.2530871022700643;
        
        A=sqrtm(P);
X=zeros(2*n+2^n+1,n);
X(1,:)=zeros(1,1);

X(2,:)=r1*A;
X(3,:)=-r1*A;
X(4,:)=r2*A;
X(5,:)=-r2*A;


w=[w0,w1*ones(1,2),w2*ones(1,2)]';

    else
%     w0=0.3;
%     a2=(1-sqrt(1-2*w0))/2;
%     a1=(1-a2)/2;
% r1=1/sqrt(a1);
% r2=1/sqrt(a2);
% w1=a1^2;
% w2=a2^2/2^2;
w0=0.41553535186548973;
w1=0.021681819434216532;
w2=0.12443434259941118;
r1=2.6060099476935847;
r2=1.190556300661233;
A=sqrtm(P);
X=zeros(2*n+2^n+1,n);
X(1,:)=zeros(1,2);

X(2,:)=r1*A(:,1);
X(3,:)=-r1*A(:,1);
X(4,:)=r1*A(:,2);
X(5,:)=-r1*A(:,2);

X(6,:)=r2*(A(:,1)+A(:,2));
X(7,:)=-r2*(A(:,1)+A(:,2));
X(8,:)=r2*(A(:,1)-A(:,2));
X(9,:)=-r2*(A(:,1)-A(:,2));
w=[w0,w1*ones(1,4),w2*ones(1,4)]';
    end
else
%  w0=2/(2+n);
%  a1=1/(2+n);
%  a2=n/(2+n);
% a1=2/(n+2);
% a2=(n-2)/(n+2);
% r1=1/sqrt(a1);
% r2=1/sqrt(a2);
% w1=a1^2;
% w2=a2^2/2^n;
r1=sqrt((n+2)/2);
r2=sqrt((n+2)/(n-2));
w1=4/(n+2)^2;
w2=(n-2)^2/(2^n*(n+2)^2);
% [u,s,v]=svd(P);
% A=chol(s)*u;
% A=A';
A=sqrtm(P);
X=zeros(2*n+2^n,n);

%%%%%% conjugate directions %%%%%%%%%%%%%%
% index=GenerateIndex(n,n*ones(1,n));
% [roww,coll]=size(index);
% dr=[];
% for i=1:1:roww
% if length(find(index(i,:)>2))==0
%     dr=vertcat(dr,index(i,:));
% end
% end
dr=prod_conjugate_dir(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         X(1,:)=zeros(1,n);
  
for i=1:1:n
 
              
    X(i,:)=r1*A(:,i);
    X(i+n,:)=-r1*A(:,i);
        

end
% mo=-1*ones(1,n);
for i=1:1:2^n
%     r=mo.^dr(i,:);
    sig=0;
    for j=1:1:n
        sig=sig+dr(i,j)*A(:,j);
    end
    X(2*n+i,:)=r2*sig;
end
% w0=1-2*n*w1-2^n*w2;

w=[w1*ones(1,2*n),w2*ones(1,2^n)]';
end
for i=1:1:n
    X(:,i)=X(:,i)+mu(i);
end

if isreal(X)==1 
1;
else
  error('imag points ') 
end

if isreal(w)==1 
1;
else
  error('imag weights ') 
end

if min(w)<0 
  error('neg weights ') 
end