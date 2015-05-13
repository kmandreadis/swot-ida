% this script runs a single experiment, taking Run Directory as input
% April 2015

function RunExp(RunDir,ShowFigs)

% RunDir='.';

disp(['Running ' RunDir])

ObsFile=[RunDir '/SWOTobs.txt'];
[DAll,Obs] = ReadObs(ObsFile);

ParamFile=[RunDir '/params.txt'];
[Chain,Prior,jmp,R,Exp] = ReadParams(ParamFile,DAll);

TruthFile=[RunDir '/truth.txt'];
AllTruth=ReadTruth (TruthFile,DAll);

[D,Obs,AllObs,DAll,Truth]=SelObs(DAll,Obs,Exp,AllTruth);

[Obs] = CalcdA(D,Obs);
[AllObs] = CalcdA(DAll,AllObs);

[Prior,jmp,AllObs]=ProcessPrior(Prior,Obs,D,AllObs,jmp,DAll); 

[Obs,Prior] = GetCovMats(D,Obs,Prior);

%limit slopes to zero...
Obs.S(Obs.S<0)=0;
AllObs.S(AllObs.S<0)=0;

Chain=MetropolisCalculations(Prior,D,Obs,jmp,Chain,R,DAll,AllObs);

[Estimate,Chain]=CalculateEstimates (Chain,D,Obs,Prior,DAll,AllObs);

[Estimate] = FilterEstimate(Estimate,Chain,D,Obs);

Err=CalcErrorStats(AllTruth,Estimate,DAll);

Err=DispRMSEStats(Err,Truth,Prior,Estimate);

if ShowFigs,
    MakeFigs(D,Truth,Prior,Chain,Estimate,Err,AllTruth,DAll);
end

WriteSummary (R,Err,Estimate,RunDir);

% CheckBals(Truth,Obs,D,Prior,Chain,R,Estimate)

save([RunDir '/RunData.mat'])

return