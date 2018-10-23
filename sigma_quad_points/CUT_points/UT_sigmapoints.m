function [x,w]=UT_sigmapoints(mu,P,M)
mu=mu(:);

% x=mu(:)';
% w=1;
% return
% global kappa
kappa=1;
n=length(P);
% m=2 is 2n+1 points

if M==2 
    if n<4
    k=3-n;

    else
        
%         k=kappa;
k=1;
    end
        
    x(1,:)=mu';
%     try
    w(1)=k/(n+k);
%     catch
%         keyboard
%     end
if min(eig(P))<=0
    min(eig(P))
end

    A=sqrtm((n+k)*P);
%     A=chol((n+k)*P);
   for i=1:1:n
       x(i+1,:)=(mu+A(:,i))';
       x(i+n+1,:)=(mu-A(:,i))';
       w(i+1)=1/(2*(n+k));
       w(i+n+1)=1/(2*(n+k));
   end
   w=w';
end

% m=4 is 4n+1 points

if M==4 
    a=normpdf(2)/normpdf(1);
    b=normpdf(3)/normpdf(1);
    x(1,:)=mu';
    w(1)=1-n*(1+a+b)/(1+4*a+9*b);
    A=sqrtm(P);
   for i=1:1:n
       x(i+1,:)=(mu+A(:,i))';
       x(i+n+1,:)=(mu-A(:,i))';
       x(i+2*n+1,:)=(mu+2*A(:,i))';
       x(i+3*n+1,:)=(mu-2*A(:,i))';
       w(i+1)=1/(2+8*a+18*b);
       w(i+n+1)=1/(2+8*a+18*b);
       w(i+2*n+1)=a/(2+8*a+18*b);
       w(i+3*n+1)=a/(2+8*a+18*b);
   end
   w=w';
end

% m=6 is 6n+1 points

if M==6 
    a=normpdf(2)/normpdf(1);
    b=normpdf(3)/normpdf(1);
    x(1,:)=mu';
    w(1)=1-n*(1+a+b)/(1+4*a+9*b);
    A=sqrtm(P);
   for i=1:1:n
       x(i+1,:)=(mu+A(:,i))';
       x(i+n+1,:)=(mu-A(:,i))';
       x(i+2*n+1,:)=(mu+2*A(:,i))';
       x(i+3*n+1,:)=(mu-2*A(:,i))';
       x(i+4*n+1,:)=(mu+3*A(:,i))';
       x(i+5*n+1,:)=(mu-3*A(:,i))';
       w(i+1)=1/(2+8*a+18*b);
       w(i+n+1)=1/(2+8*a+18*b);
       w(i+2*n+1)=a/(2+8*a+18*b);
       w(i+3*n+1)=a/(2+8*a+18*b);
       w(i+4*n+1)=b/(2+8*a+18*b);
       w(i+5*n+1)=b/(2+8*a+18*b);
   end
   w=w';
end
