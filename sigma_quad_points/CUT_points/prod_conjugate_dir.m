function dr=prod_conjugate_dir(n)

if n==1
    dr=[1;-1];
else
p=prod_conjugate_dir(n-1);
dr=zeros(2^(n),n);
dr(1:1:2^(n-1),1)=1;
dr(1:1:2^(n-1),2:1:n)=p;
dr(2^(n-1)+1:1:2^n,:)=-dr(1:1:2^(n-1),:);
end
end
    