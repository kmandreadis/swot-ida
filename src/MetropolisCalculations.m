function [C] = MetropolisCalculations(Prior,D,Obs,jmp,C,R,DAll,AllObs)

[Delta,DeltaA,B,C,~,thetaun,thetauq,thetauQbar,R]=InitializeMetropolis (D,C,Prior,R);

%6.1) initial probability calculations
thetauA0=( Prior.meanQbase.*Prior.meann.*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
thetauA0(thetauA0<jmp.A0min)=jmp.A0min(thetauA0<jmp.A0min)*1.1; %sketchy
pu1=exp(-0.5.*(thetauA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetauA0-Prior.meanA0));
pu2=exp(-0.5.*(thetaun-Prior.meann)'*diag(Prior.stdn.^-2)*(thetaun-Prior.meann));
if C.Estimateq,
    pu3=exp(-0.5.*(thetauq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetauq-Prior.meanq));
end
v=(Prior.covQbar*Prior.meanQbar)^2;
[mu,sigma] = logninvstat(Prior.meanQbar,v);
pu4=lognpdf(thetauQbar,mu,sigma);

thetauQb=Prior.meanQbase;

[fu,dQdx,dAdt]=CalcLklhd(Obs,thetauA0,thetaun,D,Prior,Delta,DeltaA,B,thetauq);

if Prior.meanq==-1,
    Prior.meanq=dQdx+dAdt;
    C.thetaq(:,1)=Prior.meanq;
    thetauq=C.thetaq(:,1);
    pu3=exp(-0.5.*(thetauq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetauq-Prior.meanq));
end

%6.2) Validity check on A0 min: ensure no A=A0+dA values... this is now
%done in "ProcessPrior", and checks against entire Obs timeseries
% jmp.A0min=ceil(-min(Obs.dA,[],2));

%6.3) The loop
tic

jmp.stdQb=jmp.stdQbburn;

C.n_a1=0;
C.n_a2=0;
C.n_a3=0;

for i=1:C.N,
    if mod(i,C.Nburn/2)==0, 
        disp(['Iteration #' num2str(i) '/' num2str(C.N) '.']); 
    end    
    
    if i==C.Nburn,
        jmp.stdQb=jmp.stdQb;
    end
    
    thetavQb=thetauQb+jmp.stdQb.*R.z1(1,i);   
    thetavQb(thetavQb<jmp.Qbmin)=jmp.Qbmin(thetavQb<jmp.Qbmin); %could scalarize this line, but fine as is    
    thetavA0=( thetavQb.*thetaun.*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
    
    Qhatv = mean(mean(1./(thetaun*ones(1,DAll.nt)) .* ...
        ( (thetauA0-AllObs.A0Shift)*ones(1,DAll.nt)+AllObs.dA).^(5/3).*...
        AllObs.w.^(-2/3).*sqrt(AllObs.S) ));
        
    pv4=lognpdf(Qhatv,mu,sigma);
    
    pv1=exp(-0.5.*(thetavA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetavA0-Prior.meanA0));    
    if any(thetavA0<jmp.A0min),
        pv1=0;
    end
    fv=CalcLklhd(Obs,thetavA0,thetaun,D,Prior,Delta,DeltaA,B,thetauq);                         

%     MetRatio=exp(fv-fu)*pv1/pu1*pv4/pu4;
    MetRatio=exp(fv-fu)*pv4/pu4;
    if pv1==0,
        MetRatio=0;
    end
    if MetRatio>R.u1(i),
        C.n_a1=C.n_a1+1; %increment
        thetauQbar=thetavQb;pu4=pv4; fu=fv; %update u->v     
        thetauA0=thetavA0; pu1=pv1; %these are sort of "diagnostics"
    end    
    C.thetaQb(i)=thetauQbar;
    
    %n
    thetavn=thetaun+jmp.stdn.*R.z2(:,i);
    thetavn(thetavn<jmp.nmin)=jmp.nmin;
    pv2=exp(-0.5.*(thetavn-Prior.meann)'*diag(Prior.stdn.^-2)*(thetavn-Prior.meann));
    
    thetavA0=( thetauQbar.*thetavn.*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
    pv1=exp(-0.5.*(thetavA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetavA0-Prior.meanA0));        
    if any(thetavA0<jmp.A0min),
        pv1=0;
    end
    fv=CalcLklhd(Obs,thetavA0,thetavn,D,Prior,Delta,DeltaA,B,thetauq);    
    
%     MetRatio=exp(fv-fu)*pv2/pu2*pv1/pu1;
    MetRatio=exp(fv-fu)*pv2/pu2;
    if pv1==0,
        MetRatio=0;
    end
    if MetRatio>R.u2(i),
        C.n_a2=C.n_a2+1; %increment
        thetaun=thetavn; fu=fv; pu2=pv2; %update u->v  
        thetauA0=thetavA0; pu1=pv1; %these are sort of "diagnostics"
    end    
    C.thetan(:,i)=thetaun;    
    C.thetaA0(:,i)=thetauA0;
    
    if C.Estimateq,
        thetavq=thetauq+jmp.stdq.*R.z3(:,i);
        thetavq(thetavq<jmp.qmin)=jmp.qmin;
        pv3=exp(-0.5.*(thetavq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetavq-Prior.meanq));
        fv=CalcLklhd(Obs,thetauA0,thetaun,D,Prior,Delta,DeltaA,B,thetavq);    

        MetRatio=exp(fv-fu)*pv3/pu3;
        if MetRatio>R.u3(i),
            C.n_a3=C.n_a3+1; %increment
            thetauq=thetavq; fu=fv; pu3=pv3; %update u->v     
        end    
        C.thetaq(:,i)=thetauq;  
    end
    
    C.Like(i)=exp(fu);
    C.LogLike(i)=fu;
end
toc


disp(['Q base: Acceptance rate =' num2str(C.n_a1/C.N*100) ' pct.'])
disp(['n: Acceptance rate =' num2str(C.n_a2/C.N*100) ' pct.'])
if C.Estimateq,
    disp(['q: Acceptance rate =' num2str(C.n_a3/C.N*100) ' pct.'])
end

return