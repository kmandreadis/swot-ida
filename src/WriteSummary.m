function WriteSummary (R,Err,Estimate,RunDir)

fid=fopen([RunDir '/summary.csv'],'w');
fprintf(fid,'Seed	, 	Q  relative uncertainty	,Q relative Error ,  sigma n /n	, Error in n ,	sigma A0/A0	, Error in A0	, varQ/Q_n	, varQ/Q_A0	, varQ/Q_dA	, varQ/Q_w	, varQ/Q_S	, varQ/Q	, stdQ/Q \n');
fprintf(fid,'%f,', R.Seed);
fprintf(fid,'%f,', mean(mean(Estimate.QstdPost./Estimate.QhatPost)) );
fprintf(fid,'%f,', mean(Err.QRelErrPost));
fprintf(fid,'%f,', mean(Estimate.stdnPost./Estimate.nhat) );
fprintf(fid,'%f,', Err.RelErrN );
fprintf(fid,'%f,', mean(Estimate.stdA0Post./Estimate.A0hat) );
fprintf(fid,'%f,', Err.RelErrA0 );

for i=1:5,
    fprintf(fid,'%f,',mean(Estimate.QerrVarSum(:,i)));
end

fprintf(fid,'%f,',sum(mean(Estimate.QerrVarSum)) );
fprintf(fid,'%f,',sqrt(sum(mean(Estimate.QerrVarSum))) );

fclose(fid);

return