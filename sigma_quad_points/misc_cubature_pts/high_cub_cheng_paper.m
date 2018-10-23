function [x,w]=high_cub_cheng_paper(mu,P,N)
n=length(mu);

m=(N-1)/2;

p=0:1:m*10^(n-1);
p=p(:);
M=zeros(length(p),n);
for i=1:1:length(p)
    s=p(i);
    for j=n-1:-1:0
        M(i,n-j)=floor(s/10^(j)) ;
        s=s-10^(j)*M(i,n-j);
    end
end
ind=find(sum(M,2)==m);

M=M(ind,:);
np=length(ind);

up=sqrt(M/m);
uf=@(pi)sqrt(pi/m);

% s=sym('s',[1,n]);
wp=zeros(np,1);
Ie=eye(n);
for pp=1:1:np
    S=[1,zeros(1,n)];
    for i=1:n
        for j=0:1:M(pp,i)-1
            c=(uf(M(pp,i))^2-uf(j)^2);
            P2=[1/c,2*Ie(i,:);-uf(j)^2/c,zeros(1,n)];
            S=multiply_polyND(S,P2);
        end
    end
    wp(pp)=spherical_weight_cheng(S);
    cp=length(find(M(pp,:)>0));
    wp(pp)=2^(-cp)*wp(pp);
    
    
end
xs=[];
ws=[];


V=general_conj_axis(n,n);
nv=size(V,1);
for i=1:1:np
    X=repmat(up(i,:),nv,1).*V;
    X=simplify_polyND([zeros(nv,1),X]);
    X=X(:,2:end);
    W=wp(i)*ones(size(X,1),1);
    
    xs=vertcat(xs,X);
    ws=vertcat(ws,W);
end
ws=ws/pi^(n/2);

[ xr, wr ] = gen_laguerre_compute ( m-1, (n-1)/2 );
xr=sqrt(xr);
if m==1  %3
    xr=sqrt(n/2);
    wr=1/(2)*gamma(n/2);
end
if m==2  % 5
    xr=[0,sqrt(n/2+1)];
    wr=[1/(n+2)*gamma(n/2),n/(2*(n+2))*gamma(n/2)];
end
if m==3 % 7
    xr=[0.5*sqrt(4+2*n+2*sqrt(4+2*n)),0.5*sqrt(4+2*n-2*sqrt(4+2*n))];
    wr=[(1/4-0.5/sqrt(4+2*n))*gamma(n/2),(1/4+0.5/sqrt(4+2*n))*gamma(n/2)];
end
if m==4 % 9
    a=0.5*gamma(n/2);
    b=0.5*gamma(n/2+1);
    c=0.5*gamma(n/2+2);
    d=0.5*gamma(n/2+3);
    e=0.5*gamma(n/2+4);
    rr2=(-c*d + b*e + sqrt(-3*c^2*d^2 + 4*b*d^3 + 4*c^3*e - 6*b*c*d*e +b^2*e^2))/(2*(-c^2 + b*d));
    rr3=(c*d - b*e + sqrt(-3 *c^2 *d^2 + 4 *b* d^3 + 4 *c^3 *e - 6 *b* c* d *e + b^2*e^2))/(2*c^2 - 2*b*d);
    B=[1,1,1;0,rr2,rr3;0,rr2^2,rr3^2];
    W=inv(B)*[a;b;c];
    w1=W(1);
    w2=W(2);
    w3=W(3);
    r2=sqrt(rr2);
    r3=sqrt(rr3);
    xr=[0,r2,r3];
    wr=[w1,w2,w3];
end

% wr=0.5*wr;
% wr=wr/sum(wr);
% x=zeros(length(ws)*length(wr),n);
% w=zeros(length(ws)*length(wr),1);
x=zeros(1,n);
w=0;

k=1;
cnt=0;

for j=1:1:length(wr)
    for i=1:1:length(ws)
        if xr(j)==0
            if cnt==0
                x(k,:)=xs(i,:)*xr(j);
                cnt=1;
            end
            w(k)=w(k)+ws(i);
            if i==length(ws)
                w(k)=w(k)*wr(j);
            k=k+1;
            end
        else
            x(k,:)=xs(i,:)*xr(j);
            w(k)=ws(i)*wr(j);
            k=k+1;
        end
    end
end
x=sqrt(2)*x;
w=w(:);
% 
% ind=find(sum(abs(x-repmat(x(k,:),size(x,1),1)),2)==0);
% if length(ind)>1
%     x(ind(2:end),:)=[];
%     w(ind(2:end))=[];
% end




end