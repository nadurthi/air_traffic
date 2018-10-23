function D=third_conj_axis(n)
dr=prod_conjugate_dir(3);
C = nchoosek(1:n,3);
[rc,cc]=size(C);
[rdr,cdr]=size(dr);
D=zeros(4*n*(n-1)*(n-2)/3,n);
for i=0:rc:rdr*rc-rc
    for j=1:1:rc
   D(i+j,C(j,1))=dr(floor(i/rc)+1,1);
   D(i+j,C(j,2))=dr(floor(i/rc)+1,2);
   D(i+j,C(j,3))=dr(floor(i/rc)+1,3);
    end
end

    

