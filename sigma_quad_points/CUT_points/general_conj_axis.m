function D=general_conj_axis(n,m)
%n is the dimension of system
%m is the number of ones 1(s) I need in each vector
dr=prod_conjugate_dir(m);
C = nchoosek(1:n,m);
[rc,cc]=size(C);
[rdr,cdr]=size(dr);

D=[];
for i=1:1:nchoosek(n,m)
    A=zeros(2^m,n);
    for j=1:1:m
        A(:,C(i,j))=dr(:,j);
    end
    D=vertcat(D,A);
end