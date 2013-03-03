function [a,b,c]=GetABC(r,nR,L)

if r==1,
    a=nan;
    b=2/(L(r)+L(r+1));
    c=2/(L(r)+L(r+1));
elseif r==nR,    
    a=2/(L(r)+L(r-1));
    b=nan;
    c=-2/(L(r)+L(r-1)); %this - to agree with notation in the paper
else
    a=1/(L(r)+L(r-1));
    b=1/(L(r)+L(r+1));
    c=(L(r-1)-L(r+1))/( (L(r)+L(r-1))*(L(r)+L(r+1)) ); 
end    

return