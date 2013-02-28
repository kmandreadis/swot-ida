function B = CalcB(nr,nt)

Row=0;

N=nr*(nt-1);
M=nr*nt;
B=nan(N,M);

for r=1:nr,    
        
    for i=1:nt-1,
        Row=Row+1;
        
        B(Row,:)=zeros(M,1);

        t1=i; t2=i+1;
        Col1=t1+(r-1)*nt; 
        Col2=t2+(r-1)*nt;
                
        B(Row,Col1)=.5;
        B(Row,Col2)=.5;
    end
end

return