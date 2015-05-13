function [Q]=RunBjerklieAlg(D,Obs,PriorMAF,iSlopeVar)

ManningFunc=@(p,x) (1./x(:,1).* ( (x(:,2)-p).*x(:,3) ).^(5/3) .* x(:,4) .* x(:,5).^(1/2)    );
ManningAvgFunc=@(p,x) (mean( 1./x(:,1).* ( (x(:,2)-p).*x(:,3) ).^(5/3) .* x(:,4) .* x(:,5).^(1/2)  )  );
H0g=0;

for r=1:D.nR,
    Sa=mean(Obs.S(r,:));

    Wa=mean(Obs.w(r,:));
    Ha=mean(Obs.h(r,:));
    chiW=std(Obs.w(r,:))/Wa;
    chiH=std(Obs.h(r,:))/Ha;

    N=1.06*(chiW/chiH)^-1.11;
    b=1-1/(1+N);

    c1=0.85;
    x1=2.257+1.308*log10(chiH)+0.99*log10(chiW)+0.435*log10(Sa);
    na=0.22*Sa^0.18;

    n=c1.*( Obs.w(r,:).*Obs.h(r,:)./Wa./Ha ).^x1 .* na;


    %columsn of inputs
    %1 ~ roughness coefficient
    %2 ~ water elevation
    %3 ~ b (invariant)
    %4 ~ width
    %5 ~ slope (invariant)

    if iSlopeVar,
        X=[n' Obs.h(r,:)' b*ones(D.nt,1) Obs.w(r,:)' Obs.S(r,:)' ];
    else
        X=[n' Obs.h(r,:)' b*ones(D.nt,1) Obs.w(r,:)' Sa*ones(D.nt,1) ];
    end    
    
    H0=nlinfit(X,PriorMAF,ManningAvgFunc,H0g);

    Qr=ManningFunc(H0,X);
    
    if abs(mean(Qr)-PriorMAF)/PriorMAF > 1E-3,
        disp('Uh Oh! Didn''t hit target flow using nlinfit...')
    end
    
    Qall(:,r)=Qr;
end

Q=mean(Qall,2);

return