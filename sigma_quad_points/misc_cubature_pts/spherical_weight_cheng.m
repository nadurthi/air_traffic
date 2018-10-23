function wp=spherical_weight_cheng(S)

[m,d]=size(S);
n=d-1;
wp=0;
for c=1:1:m
    C=S(c,1);
    K=S(c,2:end);
    
w=1;
for i=1:1:n
w=w*gamma((K(i)+1)/2);
end
w=2*w/gamma((sum(K)+n)/2);

wp=wp+C*w;
end