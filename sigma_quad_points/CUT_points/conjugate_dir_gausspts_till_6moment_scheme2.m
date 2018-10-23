function [X,w]=conjugate_dir_gausspts_till_6moment_scheme2(mu,P)
%dimension n of the system
n=length(mu);
if n==5
    r1=2.121320343559643;
    r2=1.1338934190276817;
    r3=2.9999999999999996;

w1=0.03292181069958846;
w2=0.014703360768175577;
w3=0.000685871056241427;

X1=r1*[eye(n);-eye(n)];
X2=r2*general_conj_axis(n,n);
X3=r3*general_conj_axis(n,2);
X0=zeros(1,n);
%*********************************************
X=[X0;X1;X2;X3];
w=[w1*ones(size(X1,1),1);w2*ones(size(X2,1),1);w3*ones(size(X3,1),1)];
w=[1-sum(w);w];

end    
% 
if n<7 && n~=5
    
    options=optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000);
[x,f]=fsolve(@(x)moment_6th_ND_eqns(x,n),[2,3,4]',options);
r1=x(1);
r2=x(2);
r3=x(3);
w=[(7-(n-1))/r1^6;1/(2^n*r2^6);1/(2*r3^6)];
w1=w(1);
w2=w(2);
w3=w(3);
X1=r1*[eye(n);-eye(n)];
X2=r2*general_conj_axis(n,n);
X3=r3*general_conj_axis(n,2);
X0=zeros(1,n);
%*********************************************
X=[X0;X1;X2;X3];
w=[w1*ones(size(X1,1),1);w2*ones(size(X2,1),1);w3*ones(size(X3,1),1)];
w=[1-sum(w);w];


end



if n==7
    

r1=2.5512003554818197;
r2=0.964263097900639;
r3 =2.3255766977088315;
w1=0.01269406283896717; 
w2=0.004859445930542121;
w3 =0.000395089978993786;

X1=r1*[eye(n);-eye(n)];
X2=r2*general_conj_axis(n,n);
X3=r3*general_conj_axis(n,3);
X0=zeros(1,n);
%*********************************************
X=[X0;X1;X2;X3];
w=[w1*ones(size(X1,1),1);w2*ones(size(X2,1),1);w3*ones(size(X3,1),1)];
w=[1-sum(w);w];

end


if n==8
    
r1=2.449489742783178;
r2=1;
r3=2.449489742783178;
w1=0.013888888888888888;
w2=0.00234375;
w3=0.0002314814814814815;


X1=r1*[eye(n);-eye(n)];
X2=r2*general_conj_axis(n,n);
X3=r3*general_conj_axis(n,3);
X0=zeros(1,n);
%*********************************************
X=[X0;X1;X2;X3];
w=[w1*ones(size(X1,1),1);w2*ones(size(X2,1),1);w3*ones(size(X3,1),1)];
w=[1-sum(w);w];

end


if n==9
    
w1 =0.015076391098114098;
w2 = 0.0011342717964254396; 
w3 =0.0001572731368706344;
r1 = 2.3439073215294153;
r2= 1.023262223053077;
r3 =2.5342864499001747;


X1=r1*[eye(n);-eye(n)];
X2=r2*general_conj_axis(n,n);
X3=r3*general_conj_axis(n,3);
X0=zeros(1,n);
%*********************************************
X=[X0;X1;X2;X3];
w=[w1*ones(size(X1,1),1);w2*ones(size(X2,1),1);w3*ones(size(X3,1),1)];
w=[1-sum(w);w];

end


%% Tranformation of the points
A=sqrtm(P);
for i=1:1:length(w)
    X(i,:)=(A*X(i,:)'+mu)';
end

%% discard if neagative weight
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
end