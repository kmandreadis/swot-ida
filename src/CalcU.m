function U = CalcU (D)

M=D.nR*D.nt;
N=D.nR*(D.nt-1);
 
 u=[zeros(1,D.nt-1)
    tril(ones(D.nt-1))];
 
U=zeros(M,N);          
          
 for i=1:D.nR,
     U(D.nt*(i-1)+[1:D.nt],(D.nt-1)*(i-1)+[1:D.nt-1])=u;
 end

 return