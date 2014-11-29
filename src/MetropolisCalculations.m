function [C] = MetropolisCalculations(Prior,D,Obs,jmp,C,R)

[Delta,DeltaA,B,C,~,thetaun,thetauq,thetauQb,R]=InitializeMetropolis (D,C,Prior,R);

%6.1) initial probability calculations
thetauA0=( Prior.meanQbase.*Prior.meann.*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
pu1=exp(-0.5.*(thetauA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetauA0-Prior.meanA0));
pu2=exp(-0.5.*(thetaun-Prior.meann)'*diag(Prior.stdn.^-2)*(thetaun-Prior.meann));
if C.Estimateq,
    pu3=exp(-0.5.*(thetauq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetauq-Prior.meanq));
end
v=(Prior.covQbase*Prior.meanQbase)^2;
[mu,sigma] = logninvstat(Prior.meanQbase,v);
pu4=prod(lognpdf(thetauQb,mu,sigma));


%temp... set up estimation of "c2" parameter, which corresponds to
%Bjerklie's "x1" parameter, to allow roughness to vary with stage... for
%now use a uniform distribution on 0,1
% thetauc2=Prior.meanc2;
thetauc2=zeros(D.nR,1);
% pu5=exp(-0.5.*(thetauc2-Prior.meanc2)'*diag(Prior.stdc2.^-2)*(thetauc2-Prior.meanc2));
pu5=Priorc2(thetauc2,Prior.meanc2);

C.thetac2=nan(D.nR,C.N);
C.thetac2(:,1)=thetauc2;

[fu,dQdx,dAdt]=CalcLklhd(Obs,thetauA0,thetaun,D,Prior,Delta,DeltaA,B,thetauq,thetauc2);

if Prior.meanq==-1,
    Prior.meanq=dQdx+dAdt;
    C.thetaq(:,1)=Prior.meanq;
    thetauq=C.thetaq(:,1);
    pu3=exp(-0.5.*(thetauq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetauq-Prior.meanq));
end

%6.2) Validity check on A0 min: ensure no A=A0+dA values
% jmp.A0min=max(jmp.A0min, ceil(-min(Obs.dA,[],2)) );
% While the above is now obsolete may need to add a similar check

%6.3) The loop
tic

jmp.stdQb=jmp.stdQbburn;

C.n_a1=0;
C.n_a2=0;
C.n_a3=0;
C.n_a5=0;

for i=1:C.N,
    if mod(i,C.Nburn/2)==0, 
        disp(['Iteration #' num2str(i) '/' num2str(C.N) '.']); 
    end    
    
    if i==C.Nburn,
        jmp.stdQb=jmp.stdQb;
    end
    
    %6.3.1: Qb
    thetavQb=thetauQb+jmp.stdQb.*R.z1(1,i);   
    thetavQb(thetavQb<jmp.Qbmin)=jmp.Qbmin(thetavQb<jmp.Qbmin); %could scalarize this line, but fine as is
    pv4=prod(lognpdf(thetavQb,mu,sigma));
    [~,nvar]=VariableRoughness(Obs,D,thetauc2,thetaun);
    thetavA0=( thetavQb.*nvar(:,1).*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
    pv1=exp(-0.5.*(thetavA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetavA0-Prior.meanA0));    
    
    fv=CalcLklhd(Obs,thetavA0,thetaun,D,Prior,Delta,DeltaA,B,thetauq,thetauc2);                         

    MetRatio=exp(fv-fu)*pv1/pu1*pv4/pu4;
    if MetRatio>R.u1(i),
        C.n_a1=C.n_a1+1; %increment
        thetauQb=thetavQb;pu4=pv4; fu=fv; %update u->v     
        thetauA0=thetavA0; pu1=pv1; %these are sort of "diagnostics"
    end    
    C.thetaQb(i)=thetauQb;
    
    %6.3.2: n
    thetavn=thetaun+jmp.stdn.*R.z2(:,i);
    thetavn(thetavn<jmp.nmin)=jmp.nmin;
    pv2=exp(-0.5.*(thetavn-Prior.meann)'*diag(Prior.stdn.^-2)*(thetavn-Prior.meann));
    
    [~,nvar]=VariableRoughness(Obs,D,thetauc2,thetavn);   
    thetavA0=( thetauQb.*nvar(:,1).*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
    pv1=exp(-0.5.*(thetavA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetavA0-Prior.meanA0));        
    
    fv=CalcLklhd(Obs,thetavA0,thetavn,D,Prior,Delta,DeltaA,B,thetauq,thetauc2);    
    
    MetRatio=exp(fv-fu)*pv2/pu2*pv1/pu1;
    if MetRatio>R.u2(i),
        C.n_a2=C.n_a2+1; %increment
        thetaun=thetavn; fu=fv; pu2=pv2; %update u->v  
        thetauA0=thetavA0; pu1=pv1; %these are sort of "diagnostics"
    end    
    C.thetan(:,i)=thetaun;    
    C.thetaA0(:,i)=thetauA0;
    
    %6.3.3: q
    if C.Estimateq,
        thetavq=thetauq+jmp.stdq.*R.z3(:,i);
        thetavq(thetavq<jmp.qmin)=jmp.qmin;
        pv3=exp(-0.5.*(thetavq-Prior.meanq)'*diag(Prior.stdq.^-2)*(thetavq-Prior.meanq));
        fv=CalcLklhd(Obs,thetauA0,thetaun,D,Prior,Delta,DeltaA,B,thetavq,thetauc2);    

        MetRatio=exp(fv-fu)*pv3/pu3;
        if MetRatio>R.u3(i),
            C.n_a3=C.n_a3+1; %increment
            thetauq=thetavq; fu=fv; pu3=pv3; %update u->v     
        end    
        C.thetaq(:,i)=thetauq;  
    end
    
    %6.3.4: c2
    thetavc2=thetauc2+jmp.stdc2.*R.z5(:,i);
%     pv5=exp(-0.5.*(thetavc2-Prior.meanc2)'*diag(Prior.stdc2.^-2)*(thetavc2-Prior.meanc2));
    pv5=Priorc2(thetavc2,Prior.meanc2);
    
    [~,nvar]=VariableRoughness(Obs,D,thetavc2,thetaun);   
    thetavA0=( thetavQb.*nvar(:,1).*Obs.w(:,1).^(2/3).*Obs.S(:,1).^(-.5) ).^(3/5);
    pv1=exp(-0.5.*(thetavA0-Prior.meanA0)'*diag(Prior.stdA0.^-2)*(thetavA0-Prior.meanA0));        

    
    fv=CalcLklhd(Obs,thetavA0,thetaun,D,Prior,Delta,DeltaA,B,thetauq,thetavc2);        
    MetRatio=exp(fv-fu)*pv5/pu5*pv1/pu1;
    if MetRatio>R.u5(i),
        C.n_a5=C.n_a5+1; %increment
        thetauc2=thetavc2; fu=fv; pu5=pv5; %update u->v  
    end    
    C.thetac2(:,i)=thetauc2;    
    
    
    C.Like(i)=exp(fu);
    C.LogLike(i)=fu;
end
toc


disp(['Q base: Acceptance rate =' num2str(C.n_a1/C.N*100) ' pct.'])
disp(['n: Acceptance rate =' num2str(C.n_a2/C.N*100) ' pct.'])
disp(['c2: Acceptance rate =' num2str(C.n_a5/C.N*100) ' pct.'])
if C.Estimateq,
    disp(['q: Acceptance rate =' num2str(C.n_a3/C.N*100) ' pct.'])
end

return