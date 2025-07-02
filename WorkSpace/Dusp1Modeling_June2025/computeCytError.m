function CytError = computeCytError(x,indsPars,SSAModel_100,extendedMod)
SSAModel_100.parameters(indsPars,2) = num2cell(x);
ssaSoln_100 = SSAModel_100.solve;
CytError = compareDistPlots(ssaSoln_100,extendedMod,[1:12],6,2);
CytError = CytError + compareDistPlots(ssaSoln_100,extendedMod,[1:12],5,1);
