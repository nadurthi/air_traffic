function [X,w]=conjugate_dir_gausspts_till_8moment(mu,P)
%dimension n of the system
n=length(mu);

if n==6
   if exist('CUT86D.mat')==2
       load('CUT86D')
       X=XXcut8;
       w=WWcut8;
       A=sqrtm(P);
       for i=1:1:length(w)
       X(i,:)=A*X(i,:)'+mu;
       end
       return
   end
       
end

if n==2
A=sqrtm(P);
n1=A(:,1);
n2=A(:,2);
a1 = 0.23379853497231115;
a2 = 1.3867121461809555;
a4 = 0.288547926731565; 
a3 =0.771286446121831; 



w1 = 0.04382264267013926;
w2 = 0.1405096621714662;
w3 = 0.0009215768861610588;
w4=0.01240953967762697; 


r1=1/sqrt(a1);
r2=1/sqrt(a2);
r3=1/sqrt(a3);
r4=1/sqrt(a4);

h=3;
X(1,:)=[0,0];

X(2,:)=r1*n1;
X(3,:)=-r1*n1;
X(4,:)=r1*n2;
X(5,:)=-r1*n2;

X(6,:)=r2*(n1+n2);
X(7,:)=r2*(n1-n2);
X(8,:)=-r2*(n1+n2);
X(9,:)=-r2*(n1-n2);

X(10,:)=r4*(n1+n2);
X(11,:)=r4*(n1-n2);
X(12,:)=-r4*(n1+n2);
X(13,:)=-r4*(n1-n2);

X(14,:)=r3*(n1+h*n2);
X(15,:)=r3*(n1-h*n2);
X(16,:)=r3*(-n1+h*n2);
X(17,:)=r3*(-n1-h*n2);
X(18,:)=r3*(h*n1+n2);
X(19,:)=r3*(h*n1-n2);
X(20,:)=r3*(-h*n1+n2);
X(21,:)=r3*(-h*n1-n2);

w0=1-4*w1-4*w2-4*w4-8*w3;

w=[w0,w1*ones(1,4),w2*ones(1,4),w4*ones(1,4),w3*ones(1,8)]';
for i=1:1:n
    X(:,i)=X(:,i)+mu(i);
end
return
end


if n==3

a1 = 0.1966319276379789;
a2 = 1.9427321792767849;
a3 = 0.2944016021306629; 
a4 = 0.41171525390196356; 
a6 = 0.5866854673078312;


w1 = 0.024631993437193266;
w2 = 0.08151009408908164;
w3 = 0.009767235524166815;
w4=  0.00577248937435553; 
w6 = 0.000279472936899139;

r1=1/sqrt(a1);
r2=1/sqrt(a2);
r3=1/sqrt(a3);
r4=1/sqrt(a4);
r6=1/sqrt(a6);

h=2.74;
%%%%%%%%%%%%%% generating directions and corresponding points***********
A=sqrtm(P);
X=zeros(2*n+2*2^n+2*n*(n-1)+n*2^n+1,n);
%*******generating the CA direction***********
index=GenerateIndex(n,n*ones(1,n));
[roww,coll]=size(index);
dr=[];
for i=1:1:roww
if length(find(index(i,:)>2))==0
    dr=vertcat(dr,index(i,:));
end
end
dr
%***************  PA   ******************************
for i=1:1:n+1
    if i==1
        X(i,:)=zeros(1,n);
    else
              
    X(i,:)=r1*A(:,i-1);
    X(i+n,:)=-r1*A(:,i-1);
        
    end
end
%**************** CA - Space diagnols**************
mo=-1*ones(1,n);
for i=1:1:2^n
    rr=mo.^dr(i,:);
    sig=0;
    for j=1:1:n
        sig=sig+rr(j)*A(:,j);
    end
    X(2*n+1+i,:)=r2*sig;
    X(2*n+1+2^n+i,:)=r4*sig;
end
%*******generating the Plane Diagonal direction***********

index=GenerateIndex(n,n*ones(1,n));

dr=[];
for i=3:1:n
index(find(index==i))=0;
end

[roww,coll]=size(index);
for i=1:1:roww
if length(find(index(i,:)==0))==n-2
    dr=vertcat(dr,index(i,:));
end
end
[rowwdr,coll]=size(dr);
drr=dr(1,:);
for i=1:1:rowwdr
    [rdr,coll]=size(drr);
    dd=0;
    for j=1:1:rdr
        dd(j)=sum(abs(drr(j,:)-dr(i,:)));
    end
    
    if length(find(dd==0))==0
        drr=[drr;dr(i,:)];
    end
end
drr(find(drr==2))=-1;
%*********************************************
    for i=1:1:2*n*(n-1)
    sig=0;
    for j=1:1:n
        sig=sig+drr(i,j)*A(:,j);
    end
  
     X(2*n+1+2*2^n+i,:)=r3*sig;
    end
    

% *********   space multisector  ************
index=GenerateIndex(n,n*ones(1,n));
[roww,coll]=size(index);
dr=[];
for i=1:1:roww
if length(find(index(i,:)>2))==0
    dr=vertcat(dr,index(i,:));
end
end
mo=-1*ones(1,n);
p=0;
    for i=1:1:n*2^n
        p=p+1;
        
        rr=mo.^dr(p,:);
        if rem(i,2^n)==0
            p=0;
        end
        k=ceil(i/2^n)-1;
        pp=[ones(1,k),h,ones(1,n-k-1)];
    sig=0;
    for j=1:1:n
        sig=sig+pp(j)*rr(j)*A(:,j);
    end
     X(2*n+1+2*2^n+2*n*(n-1)+i,:)=r6*sig;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0=1-2*n*w1-2^n*w2-2^n*w4-2*n*(n-1)*w3-n*2^n*w6;

w=[w0,w1*ones(1,2*n),w2*ones(1,2^n),w4*ones(1,2^n),w3*ones(1,2*n*(n-1)),w6*ones(1,n*2^n)]';
for i=1:1:n
    X(:,i)=X(:,i)+mu(i);
end
return
end


if n==4

a1 = 0.20629093125597198;
a2 = 1.585407549822809;
a3 = 0.2851818320349825; 
a4 = 0.5660749628608427; 
a6 = 0.7889090077901053;


w1 = 0.01811008737283111;
w2 = 0.032063273384586845;
w3 = 0.006614353755080834;
w4=0.003489906522946932; 
w5=0.0006510416666666666; 
w6 = 0.00025218336987488566;

r1=1/sqrt(a1);
r2=1/sqrt(a2);
r3=1/sqrt(a3);
r4=1/sqrt(a4);
r5=2;
r6=1/sqrt(a6);

h=3;

end


if n==5

a1 = 0.1866956121576737;
a2 = 1.420294749459945;
a3 = 0.29836021128926843; 
a4 = 0.5123685659872401; 
a6 = 0.8065591548429262;


w1 = 0.010529034221546607;
w2 = 0.015144019639537572;
w3 = 0.0052828996967816825;
w4=0.0010671298950159158; 
w5=0.0006510416666666666; 
w6 = 0.00013776017592074394;

r1=1/sqrt(a1);
r2=1/sqrt(a2);
r3=1/sqrt(a3);
r4=1/sqrt(a4);
r5=2;
r6=1/sqrt(a6);

h=3;
end

if n==6
    
a1 = 0.16666666666666666;
a3 = 0.3333333333333333;
a2 = 1.251685733072443;
a4 = 0.42609204470533507;
a6 = 0.8333333333333334;

w1 = 0.006172839506172839; 
w2 = 0.006913443044833937;
w3 =0.004115226337448559; 
w4=0.0002183265828666806;
w5=0.0006510416666666666; 
w6 = 0.00007849171328446504;

r1=1/sqrt(a1);
r2=1/sqrt(a2);
r3=1/sqrt(a3);
r4=1/sqrt(a4);
r6=1/sqrt(a6);
r5=2;

h=3;
end
% [u,s,v]=svd(P);
% A=chol(s)*u;
% A=A';

%%%%%%%%%%%%%% generating directions and corresponding points***********
A=sqrtm(P);
X=zeros(2*n+2*2^n+2*n*(n-1)+4*n*(n-1)*(n-2)/3+n*2^n+1,n);
%*******generating the CA direction***********
index=GenerateIndex(n,n*ones(1,n));
[roww,coll]=size(index);
dr=[];
for i=1:1:roww
if length(find(index(i,:)>2))==0
    dr=vertcat(dr,index(i,:));
end
end
%***************  PA   ******************************
for i=1:1:n+1
    if i==1
        X(i,:)=zeros(1,n);
    else
              
    X(i,:)=r1*A(:,i-1);
    X(i+n,:)=-r1*A(:,i-1);
        
    end
end
%**************** CA - Space diagnols**************
mo=-1*ones(1,n);
for i=1:1:2^n
    rr=mo.^dr(i,:);
    sig=0;
    for j=1:1:n
        sig=sig+rr(j)*A(:,j);
    end
    X(2*n+1+i,:)=r2*sig;
    X(2*n+1+2^n+i,:)=r4*sig;
end
%*******generating the Plane Diagonal direction***********

index=GenerateIndex(n,n*ones(1,n));

dr=[];
for i=3:1:n
index(find(index==i))=0;
end

[roww,coll]=size(index);
for i=1:1:roww
if length(find(index(i,:)==0))==n-2
    dr=vertcat(dr,index(i,:));
end
end
[rowwdr,coll]=size(dr);
drr=dr(1,:);
for i=1:1:rowwdr
    [rdr,coll]=size(drr);
    dd=0;
    for j=1:1:rdr
        dd(j)=sum(abs(drr(j,:)-dr(i,:)));
    end
    
    if length(find(dd==0))==0
        drr=[drr;dr(i,:)];
    end
end
drr(find(drr==2))=-1;
%*********************************************
    for i=1:1:2*n*(n-1)
    sig=0;
    for j=1:1:n
        sig=sig+drr(i,j)*A(:,j);
    end
  
     X(2*n+1+2*2^n+i,:)=r3*sig;
    end
    
%*********   3-Subspace bisectors**********
index=GenerateIndex(n,n*ones(1,n));
for i=n:-1:4
     index(find(index==i))=0;
end
index(find(index==3))=1;

[roww,coll]=size(index);
dr=[];
for i=1:1:roww
if length(find(index(i,:)==0))==n-3
    dr=vertcat(dr,index(i,:));
end
end
[rowwdr,coll]=size(dr);
drr=dr(1,:);
for i=1:1:rowwdr
    [rdr,coll]=size(drr);
    dd=0;
    for j=1:1:rdr
        dd(j)=sum(abs(drr(j,:)-dr(i,:)));
    end
    
    if length(find(dd==0))==0
        drr=[drr;dr(i,:)];
    end
end
drr(find(drr==2))=-1;
    for i=1:1:4*n*(n-1)*(n-2)/3
    sig=0;
    for j=1:1:n
        sig=sig+drr(i,j)*A(:,j);
    end
  
     X(2*n+1+2*2^n+2*n*(n-1)+i,:)=r5*sig;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *********   space multisector  ************
index=GenerateIndex(n,n*ones(1,n));
[roww,coll]=size(index);
dr=[];
for i=1:1:roww
if length(find(index(i,:)>2))==0
    dr=vertcat(dr,index(i,:));
end
end
mo=-1*ones(1,n);
p=0;
    for i=1:1:n*2^n
        p=p+1;
        
        rr=mo.^dr(p,:);
        if rem(i,2^n)==0
            p=0;
        end
        k=ceil(i/2^n)-1;
        pp=[ones(1,k),h,ones(1,n-k-1)];
    sig=0;
    for j=1:1:n
        sig=sig+pp(j)*rr(j)*A(:,j);
    end
     X(2*n+1+2*2^n+2*n*(n-1)+4*n*(n-1)*(n-2)/3+i,:)=r6*sig;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0=1-2*n*w1-2^n*w2-2^n*w4-2*n*(n-1)*w3-4*n*(n-1)*(n-2)/3*w5-n*2^n*w6;

w=[w0,w1*ones(1,2*n),w2*ones(1,2^n),w4*ones(1,2^n),w3*ones(1,2*n*(n-1)),w5*ones(1,4*n*(n-1)*(n-2)/3),w6*ones(1,n*2^n)]';
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