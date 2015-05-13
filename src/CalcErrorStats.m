function [Err]=CalcErrorStats (AllTruth,E,DAll)

%based on Analysis/SWOT/Discharge/Pepsi Challenge/src/CalcErrorStats
Qt=mean(AllTruth.Q,1);
QhatAvg=mean(E.AllQ,1);

Stats.RMSE=sqrt(mean( (Qt-QhatAvg).^2 ) );
Stats.rRMSE=sqrt(mean( (  (Qt-QhatAvg)./Qt   ).^2 ) );

r=QhatAvg-Qt;
logr=log(QhatAvg)-log(Qt);

Stats.MSC=log(  sum((Qt-mean(Qt)).^2)/sum(r.^2) -2*2/DAll.nt  );
Stats.bias=mean(r);
Stats.stdresid=std(r);
Stats.meanLogRes=mean(logr);
Stats.stdLogRes=std(logr);
Stats.meanRelRes=mean(r./Qt);
Stats.stdRelRes=std(r./Qt);

Stats.Qbart=mean(Qt);

%Copy to error struct
Err.Stats=Stats;

return