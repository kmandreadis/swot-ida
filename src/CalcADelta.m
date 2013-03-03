function Delta = CalcADelta(nr,nt)

Row=0;

N=nr*(nt-1);
M=nr*nt;
Delta=nan(N,M);

for r=1:nr,    
        
    for i=1:nt-1,
        Row=Row+1;
        
        Delta(Row,:)=zeros(M,1);

        t1=i; t2=i+1;
        Col1=t1+(r-1)*nt; 
        Col2=t2+(r-1)*nt;
                
        Delta(Row,Col1)=-1;
        Delta(Row,Col2)=1;
    end
end

return