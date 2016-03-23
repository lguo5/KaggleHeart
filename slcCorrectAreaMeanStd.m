function [areasOut,mskAdj]=slcCorrectAreaMeanStd(areasIn,nstd)

areasOut=areasIn;
areasMean=mean(areasIn);
areasStd=std(areasIn);
mskAdj = areasIn>areasMean+nstd*areasStd | areasIn<areasMean-nstd*areasStd;
areasOut(mskAdj)=areasMean;