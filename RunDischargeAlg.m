%script to estimate A0, n, and (optionally) q
% M Durand, 2013. Contact: durand.8@osu.edu

clear all

ObsFile='dat/SWOTobs.txt';
[D,Obs] = ReadObs(ObsFile);

ParamFile='dat/params.txt';
[Chain,Prior,jmp,R] = ReadParams(ParamFile,D);

TruthFile='dat/truth.txt';
Truth=ReadTruth (TruthFile,D);

[Obs] = CalcdA(D,Obs);

Chain=MetropolisCalculations(Prior,D,Obs,jmp,Chain,R);

[Estimate,Chain]=CalculateEstimates (Chain,D,Obs,Prior);

Err=CalcErrorStats(Truth,Prior,Estimate);

MakeFigs(D,Truth,Prior,Chain,Estimate,Err);

WriteSummary (R,Err,Estimate);

save('RunData.mat')
