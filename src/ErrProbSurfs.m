function ErrProbSurfs(r,Chain,Truth,Obs,D,Prior,C,R,Estimate)

start=Chain.Nburn+1;

mA0=mean(Chain.thetaA0(:,start:end));
mn=mean(Chain.thetan(:,start:end));

A0=Chain.thetaA0(r,start:end);
n=Chain.thetan(r,start:end);

values=[A0; n;];
bins1=25;
bins2=26; %do this to make sure I'm mapping everything correctly
[n1, x1]  = hist(values(1,:),bins1);
[n2, x2]  = hist(values(2,:),bins2);
delta_x1 = x1(2)-x1(1);
delta_x2 = x2(2)-x2(1);
% 2. Initialize a 2-D matrix for the 2-D histogram
n2d      = zeros(length(x1), length(x2)); 
% 3. For each row, find the indices of the X_1 values which fall into that row.
%    Compute a histogram for the X_2 values, and put it in the 2-D histogram for that row.
for i = 1:length(x1),
   ind  = (values(1,:) > x1(i)-delta_x1/2) & (values(1,:) <= x1(i)+delta_x1/2);
   n2d(i,1:length(x2)) = hist(values(2,ind), x2);
end
pdf = n2d./(sum(sum(n2d)));

figure(9)
contour(x1,x2,pdf');
title('PDF from the chain')

for i=1:Chain.N,
    QerrChain(i)=mean(sqrt(mean((Truth.Q(:,:)-Chain.thetaQ(:,:,i)).^2)));
end

F=TriScatteredInterp(A0',n',QerrChain(start:end)');

x1s=min(x1):5:max(x1);
x2s=min(x2):.0005:max(x2);

[X1,X2]=meshgrid(x1s,x2s);

QerrHat=F(X1,X2);

figure(10)
contour(X1,X2,QerrHat)

set(gca,'FontSize',14)
hold on
plot(mean(Prior.meanA0),mean(Prior.meann),'rx',mean(Estimate.A0hat),mean(Estimate.nhat),'ro',...
    mean(Truth.A0),mean(Truth.n),'r+','LineWidth',2)
hold off

xlim( [  min([mean(Truth.A0) mean(Prior.meanA0) min(x1)])*.99 max([mean(Truth.A0) mean(Prior.meanA0) max(x1)])*1.01])
ylim([min([mean(Truth.n) mean(Prior.meann) min(x2)])*.99 max([mean(Truth.n) mean(Prior.meann) max(x2)])*1.01])
title(['Total Q Error, Reach=' num2str(r)])
colorbar

[Delta,DeltaA,B,~,~,~,qhatv] =  InitializeMetropolis (D,C,Prior,R);

f=nan(length(x1),length(x2));
ObjFuncR=nan(length(x1),length(x2));
Qerr=nan(length(x1),length(x2));


for i=1:length(x1),
    for j=1:length(x2),
        A0in(1:D.nR,1)=x1(i).*(Estimate.A0hat./mean(Estimate.A0hat));
        nin(1:D.nR,1)=x2(j).*(Estimate.nhat./mean(Estimate.nhat));
        [f(i,j),dQdxv,dAdtv]=CalcLklhd(Obs,A0in,nin,D,Prior,Delta,DeltaA,B,qhatv);

        A0v=reshape((A0in*ones(1,D.nt))',D.nR*D.nt,1);
        nv=reshape((nin*ones(1,D.nt))',D.nR*D.nt,1);        
        Qv=1./nv.*(A0v+Obs.dAv).^(5/3).*Obs.wv.^(-2/3).*sqrt(Obs.Sv);
        
        Errs=reshape(dQdxv+dAdtv,D.nt-1,D.nR)';
        ObjFuncR(i,j)=sqrt(mean(Errs(r,:).^2));
        Qs=reshape(Qv,D.nt,D.nR)';
        Qerr(i,j)=sqrt(mean( (Qs(r,:)-Truth.Q(r,:)).^2 ));
    end
end

figure(11)
contour(x1,x2,ObjFuncR')
set(gca,'FontSize',14)
hold on
plot(Prior.meanA0(r),Prior.meann(r),'rx',Estimate.A0hat(r),Estimate.nhat(r),'ro',...
    Truth.A0(r),Truth.n(r),'r+','LineWidth',2)

hold off

xlim( [  min([Truth.A0(r) Prior.meanA0(r) min(x1)])*.99 max([Truth.A0(r) Prior.meanA0(r) max(x1)])*1.01])
ylim([min([Truth.n(r) Prior.meann(r) min(x2)])*.99 max([Truth.n(r) Prior.meann(r) max(x2)])*1.01])
title(['Objective function, Reach=' num2str(r)])
colorbar

figure(12)
contour(x1,x2,Qerr')
set(gca,'FontSize',14)
hold on
plot(Prior.meanA0(r),Prior.meann(r),'rx',Estimate.A0hat(r),Estimate.nhat(r),'ro',...
    Truth.A0(r),Truth.n(r),'r+','LineWidth',2)

hold off

xlim( [  min([Truth.A0(r) Prior.meanA0(r) min(x1)])*.99 max([Truth.A0(r) Prior.meanA0(r) max(x1)])*1.01])
ylim([min([Truth.n(r) Prior.meann(r) min(x2)])*.99 max([Truth.n(r) Prior.meann(r) max(x2)])*1.01])
title(['Qerr for this reach, Reach=' num2str(r)])
colorbar

figure(13)
contour(x1,x2,f')
set(gca,'FontSize',14)
hold on
plot(mean(Prior.meanA0),mean(Prior.meann),'rx',mean(Estimate.A0hat),mean(Estimate.nhat),'ro',...
    mean(Truth.A0),mean(Truth.n),'r+','LineWidth',2)

hold off

xlim( [  min([mean(Truth.A0) mean(Prior.meanA0) min(x1)])*.99 max([mean(Truth.A0) mean(Prior.meanA0) max(x1)])*1.01])
ylim([min([mean(Truth.n) mean(Prior.meann) min(x2)])*.99 max([mean(Truth.n) mean(Prior.meann) max(x2)])*1.01])
title('Objective function all reaches')
colorbar


return