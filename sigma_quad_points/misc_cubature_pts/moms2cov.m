function P=moms2cov(y,m2)
N=size(y,2);
n=size(y,1);
P=zeros(N);
I=eye(N);
for i=1:1:N
    for j=1:1:N
         z=repmat(sum(I([i,j],:),1),n,1);
            ind=find(sum(abs(z-y),2)==0);
        P(i,j)=m2(ind);
    end
end