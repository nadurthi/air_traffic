function [muk1,Pk1]=prop_mean_cov(x,w,F)
%F takes in a column vector and gives out column vector
n=length(w);
%calculating the mean
muk1=0;
for i=1:1:n
    muk1=muk1+w(i)*F(x(i,:)');
end

%calculating the covariance
Pk1=0;
for i=1:1:n
    Pk1=Pk1+w(i)*(F(x(i,:)')-muk1)*(F(x(i,:)')-muk1)';
end
end
