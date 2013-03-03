%script to estimate A0, n, and (optionally) q
% M Durand, 2013. Contact: durand.8@osu.edu

clear all

addpath ./src %directory where algorithm scripts are located

RunDir='./AirSWOTSac'; %directory where inputs located, & outputs are saved

ObsFile=[RunDir '/SWOTobs.txt'];
[D,Obs] = ReadObs(ObsFile);

ParamFile=[RunDir '/params.txt'];
[Chain,Prior,jmp,R] = ReadParams(ParamFile,D);

TruthFile=[RunDir '/truth.txt'];
Truth=ReadTruth (TruthFile,D);

[Obs] = CalcdA(D,Obs);

[Obs,Prior] = GetCovMats(D,Obs,Prior);

Chain=MetropolisCalculations(Prior,D,Obs,jmp,Chain,R);

[Estimate,Chain]=CalculateEstimates (Chain,D,Obs,Prior);

Err=CalcErrorStats(Truth,Prior,Estimate);

MakeFigs(D,Truth,Prior,Chain,Estimate,Err);

WriteSummary (R,Err,Estimate,RunDir);

save([RunDir '/RunData.mat'])

rmpath ./src