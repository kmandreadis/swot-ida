function Delta = CalcDelta(nr,nt,L)

Row=0;

N=nr*(nt-1);
M=nr*nt;
Delta=nan(N,M);

for r=1:nr,
    [a,b,c]=GetABC(r,nr,L);            
    
    for i=1:nt-1,
        t1=i; t2=i+1;
        
        Row=Row+1;
        
        Delta(Row,:)=zeros(M,1);
                
        if r==1,%upstream
            Col1=t1;
            Col2=t2;
            Col3=t1+nt;
            Col4=t2+nt;
            Delta(Row,Col1:Col2)=-c/2;
            Delta(Row,Col3:Col4)=b/2;
        elseif r==nr,
            Col1=nt+t1;
            Col2=nt+t2;
            Col3=2*nt+t1;
            Col4=2*nt+t2;
            Delta(Row,Col1:Col2)=-a/2;
            Delta(Row,Col3:Col4)=-c/2;            
        else            
            Col1=t1;
            Col2=t2;
            Col3=nt+t1;
            Col4=nt+t2;
            Col5=2*nt+t1;
            Col6=2*nt+t2;
            Delta(Row,Col1:Col2)=-a/2;
            Delta(Row,Col3:Col4)=-c/2;
            Delta(Row,Col5:Col6)=b/2;
        end
    end
end

return