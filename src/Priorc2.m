function [y]=Priorc2(c2,c2min)

ymax=2./-c2min;

yv=-ymax./c2min.*c2+ymax;

yv(c2>0)=0;
yv(c2<c2min)=0;
y=prod(yv);

return