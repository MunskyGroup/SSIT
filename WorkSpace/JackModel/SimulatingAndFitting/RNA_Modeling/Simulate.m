function output = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot)
    a = 53*ones(Ncells,1);
    b = 39*ones(Ncells,1);
    u = 2*pi*rand(Ncells,1);
    posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data

    fileName = sprintf('kon_%.3g_koff_%.3g_w_%.3g_kex_%.3g_kr_%.3g_D_%.3g_gam_%.3g_N_%d.csv', ...
                            kon, koff, w, kex, kr, D, gam, NCells);
    fileName = strrep(fileName, '.', 'p');
    
    stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,NCells,fileName)
    
    [intensTS,binCounts,abEllipse,areaEllipse,ecentricity,...
        binTS,binCountsDim,binCountsBright,hisCountsBright,...
        hisCountsDim,XcRot,iTS,XCVcRot,XYfit] = ...
        processSpotPositionData(fileName,[],a,b,0);
    
    output.intensTS = intensTS;
    output.binCounts = binCounts;
    output.abEllipse = abEllipse;
    output.areaEllipse = areaEllipse;
    output.ecentricity = ecentricity;
    output.binTS = binTS;
    output.binCountsDim = binCountsDim;
    output.binCountsBright = binCountsBright;
    output.hisCountsBright = hisCountsBright;
    output.hisCountsDim = hisCountsDim;
    output.XcRot = XcRot;
    output.iTS = iTS;
    output.XCVcRot = XCVcRot;
    output.XYfit = XYfit;
end
