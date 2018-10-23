function X=unifpts_trans(X,bdd_low,bdd_up)

n=size(X,2);
mu=(bdd_low+bdd_up)/2;
h=-bdd_low+bdd_up;
for i=1:1:n
    X(:,i)=(h(i)/2)*X(:,i)+mu(i);
end