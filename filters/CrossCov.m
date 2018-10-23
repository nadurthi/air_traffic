function Pcc=CrossCov(x,mx,z,mz,w)
mx=mx(:);
mz=mz(:);
[N,n]=size(x);
Pcc=0;
for i=1:1:N
    Pcc=Pcc+w(i)*(x(i,:)'-mx)*(z(i,:)'-mz)';
end


end