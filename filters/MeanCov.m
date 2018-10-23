function [mk,Pk]=MeanCov(x,w)
[N,n]=size(x);

W=repmat(w,1,n);
mk=sum(W.*x,1)';

MU=repmat(mk',N,1);
Z=x-MU;
Pk=Z'*(W.*Z);

end


