function X=scaled_conj_axis(n,h)
g=general_conj_axis(n,n);
X=[];
% X=zeros(n*2^n,n);
% cnt=1;
for i=1:1:n
    p=g;
    p(:,i)=p(:,i)*h;
    X=vertcat(X,p);
% X(cnt:cnt+2^n,1:n)=h*g(:,i);
% cnt=cnt+2^n+1;
end
