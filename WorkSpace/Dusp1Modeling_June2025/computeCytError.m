function [TotError,ssaSoln] = computeCytError(x,indsPars,Model,saveBest,saveName,useGR,GRModel)
arguments
    x
    indsPars
    Model
    saveBest = false
    saveName = ''
    useGR = false
    GRModel = []
end
Model.parameters(indsPars,2) = num2cell(x);
ssaSoln = Model.solve;
[CytError,KSCyt] = compareDistPlots(ssaSoln,Model,6,2);
[NucError,KSNuc] = compareDistPlots(ssaSoln,Model,5,1);
TotError = CytError + NucError;
if useGR
    [GRErrorCyt,KSGRCyt] = compareDistPlots(ssaSoln,GRModel,3,1,[0:40],[10,11],5,7);
    [GRErrorNuc,KSGRNuc] = compareDistPlots(ssaSoln,GRModel,4,2,[0:40],[10,11],5,7);
    TotError = CytError + NucError + GRErrorCyt + GRErrorNuc;
end

if saveBest
    load(saveName,'minerr')
    if TotError<minerr
        minerr = TotError
        bestx = x;
        save(saveName,'minerr','bestx')
        
        f=figure(600); clf; set(f,'Name','Nuclear RNA')
        makeCytDistPlots(ssaSoln,Model,600,5,1,[0:5:300],true,'cdf');
        for i=1:length(KSNuc)
            subplot(3,4,i);
            title(['Err=',num2str(KSNuc(i))])   
        end
        f=figure(601); clf; set(f,'Name','Cytoplasmic RNA')
        makeCytDistPlots(ssaSoln,Model,601,6,2,[0:5:300],true,'cdf');
        for i=1:length(KSCyt)
            subplot(3,4,i);
            title(['Err=',num2str(KSCyt(i))])   
        end
        f=figure(602); clf; set(f,'Name','Nuc vs. Cyto RNA')
        makeNucCytScatterPlots(ssaSoln,Model,602,[5,6],[1,2],true);
        drawnow

        if useGR
            f=figure(603);clf; set(f,'Name','Cyto GR')
            makeCytDistPlots(ssaSoln,GRModel,601,3,1,[0:40],true,'cdf',[10,11],5,7);
            f=figure(603);clf; set(f,'Name','Nuc GR')
            makeCytDistPlots(ssaSoln,GRModel,601,4,2,[0:40],true,'cdf',[10,11],5,7);
        end
    end
end

