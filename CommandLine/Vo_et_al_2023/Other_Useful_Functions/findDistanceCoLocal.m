%% findDistanceCoLocal
clear all
timeArray = [0];

psf_z=350;               % Theoretical size of the PSF emitted by a [rna] spot in the z plan, in nanometers
psf_yx=160;              % Theoretical size of the PSF emitted by a [rna] spot in the yx plan, in nanometers
voxel_size_z=500;        % Microscope conversion px to nanometers in the z axis.
voxel_size_yx=160;       % Microscope conversion px to nanometers in the xy axis.
% distScale = ([voxel_size_yx/psf_yx, voxel_size_yx/psf_yx, voxel_size_z/psf_z]);
distScale = [1 1 1];

dThreshVec = linspace(0,5,21);

for iD = 1:length(dThreshVec)
    dThresh = dThreshVec(iD);
    time = 0;%timeArray(iTime+1);
    STR0 = ['Huy_intensity_data_correct/complete__MS2_CY5_time_0_int_550_400.csv'];
    X0 = importdata(STR0);
    CellData = findColocalzeData(X0,X0,distScale,dThresh,time,false);

    fColocalMCP(iD) = mean([CellData{:,6}]./[CellData{:,2}],'omitnan');
    fColocalFISH(iD) = mean([CellData{:,6}]./[CellData{:,3}],'omitnan');
    fColocalTotal(iD) = mean([CellData{:,6}]./[CellData{:,9}],'omitnan');
    fDoubleMCP(iD) =  mean([CellData{:,7}]./([CellData{:,7}]+[CellData{:,2}]),'omitnan');
    fDoubleFISH(iD) = mean([CellData{:,8}]./([CellData{:,7}]+[CellData{:,3}]),'omitnan');

end
figure(1);
plot(dThreshVec,fColocalTotal,'o');


function CellData = findColocalzeData(X0,X1,distScale,dThresh,time,makeFigures)

iType = find(strcmp(X0.colheaders,'spot_type'));
iXpos = find(strcmp(X0.colheaders,'x'));
iYpos = find(strcmp(X0.colheaders,'y'));
iZpos = find(strcmp(X0.colheaders,'z'));
iCellNum = find(strcmp(X0.colheaders,'cell_id'));
iSpotNum = find(strcmp(X0.colheaders,'spot_id'));
iClusterSize = find(strcmp(X0.colheaders,'cluster_size'));
iFragmented = strcmp(X0.colheaders,'is_cell_fragmented');
iNucIntCh0 = find(strcmp(X0.colheaders,'nuc_int_ch_0'));
iNucIntCh1 = find(strcmp(X0.colheaders,'nuc_int_ch_1'));
iNucIntCh2 = find(strcmp(X0.colheaders,'nuc_int_ch_2'));
iNucIntCh3 = find(strcmp(X0.colheaders,'nuc_int_ch_3'));

jSpots0 = X0.data(:,iType)==0&X0.data(:,iFragmented)~=1;
jSpots1 = X1.data(:,iType)==1&X1.data(:,iFragmented)~=1;

allSpots0 = X0.data(jSpots0,:);
allSpots1 = X1.data(jSpots1,:);

jNoSpots = X0.data(:,iType)==-1;
noSpots = X0.data(jNoSpots,:);

nCells = max(max(X0.data(:,iCellNum)),max(X1.data(:,iCellNum)))+1;

keepCells = zeros(1,nCells,'logical');

CellData = table;
for iCell = 0:nCells-1
    jThisCell0 = (allSpots0(:,iCellNum)==iCell);
    thisCell0 =  allSpots0(jThisCell0,[iXpos,iYpos,iZpos,iClusterSize,iNucIntCh0,iNucIntCh1,iNucIntCh2,iNucIntCh3]);
    nSpotsThisCell0 = sum(jThisCell0);

    jThisCell1 = (allSpots1(:,iCellNum)==iCell);
    thisCell1 =  allSpots1(jThisCell1,[iXpos,iYpos,iZpos,iClusterSize,iNucIntCh0,iNucIntCh1,iNucIntCh2,iNucIntCh3]);
    nSpotsThisCell1 = sum(jThisCell1);

    % Cell Nucleus intensity
    try
        nucIntens0 = round(thisCell0(1,5));
        nucIntens1 = round(thisCell0(1,6));
        nucIntens2 = round(thisCell0(1,7));
        nucIntens3 = round(thisCell0(1,8));
        keepCells(iCell+1)=1;
    catch
        try
            nucIntens0 = round(thisCell1(1,5));
            nucIntens1 = round(thisCell1(1,6));
            nucIntens2 = round(thisCell1(1,7));
            nucIntens3 = round(thisCell1(1,8));
            keepCells(iCell+1)=1;
        catch
            jThisCellNone = noSpots(:,iCellNum)==iCell;
            thisCellNone =  noSpots(jThisCellNone,[iXpos,iYpos,iZpos,iClusterSize,iNucIntCh0,iNucIntCh1,iNucIntCh2,iNucIntCh3]);
            try 
                nucIntens0 = round(thisCellNone(1,5));
                nucIntens1 = round(thisCellNone(1,6));
                nucIntens2 = round(thisCellNone(1,7));
                nucIntens3 = round(thisCellNone(1,8));
                keepCells(iCell+1)=1;
            catch
                nucIntens0 =NaN;
                nucIntens1 =NaN;
                nucIntens2 =NaN;
                nucIntens3 =NaN;
            end
        end
    end

    if nSpotsThisCell1==1&&thisCell1(1,1)==-1
        nSpotsThisCell1=0;
        thisCell1=[];
    end
    if nSpotsThisCell0==1&&thisCell0(1,1)==-1
        nSpotsThisCell0=0;
        thisCell0=[];
    end

    % Colocalization detection
    D = squareform(pdist([thisCell0(:,1:3);thisCell1(:,1:3)].*distScale))<=dThresh;

    D01 = D(1:nSpotsThisCell0,nSpotsThisCell0+1:end);
    coLocalSpots1 = sum(D01,1)~=0;
    nSpots1Only = nSpotsThisCell1-sum(coLocalSpots1);

    coLocalSpots0 = sum(D01,2)~=0;
    nSpots0Only = nSpotsThisCell0-sum(coLocalSpots0);

    nSpotsBoth = min(sum(coLocalSpots0),sum(coLocalSpots1));

    % Double spots detection
    D00 = D(1:nSpotsThisCell0,1:nSpotsThisCell0);
    nDoubleSpots0 = sum(sum(D00))-nSpotsThisCell0;
    D11 = D(nSpotsThisCell0+1:end,nSpotsThisCell0+1:end);
    nDoubleSpots1 = sum(sum(D11))-nSpotsThisCell1;

    nSpotsTotal = nSpots1Only + nSpots0Only + nSpotsBoth;

    % TS detection and quantification
    jClusters0 = thisCell0(:,4)>=5;
    jClusters1 = thisCell1(:,4)>=5;
    numberTS0 = sum(jClusters0);
    numberTS1 = sum(jClusters1);

    tsIntensities0 = sort(thisCell0(jClusters0,4),'descend')';
    tsIntensities1 = sort(thisCell1(jClusters1,4),'descend')';

    if isempty(tsIntensities0)
        tsIntensities0=zeros(1,4);
    elseif length(tsIntensities0)<4
        tsIntensities0(end+1:4) = 0;
    elseif length(tsIntensities0)>4
        tsIntensities0 =tsIntensities0(1:4);
    end
    if isempty(tsIntensities1)
        tsIntensities1=zeros(1,4);
    elseif length(tsIntensities1)<4
        tsIntensities1(end+1:4) = 0;
    elseif length(tsIntensities1)>4
        tsIntensities1 =tsIntensities1(1:4);
    end
    
%     if ~isempty(thisCell0)
%     keepCells(iCell+1) = true;
    CellData(iCell+1,:) = {iCell,...
        nSpotsThisCell0,nSpotsThisCell1,...
        nSpots0Only,nSpots1Only,...
        nSpotsBoth,...
        nDoubleSpots0,nDoubleSpots1,...
        nSpotsTotal,...
        numberTS0,numberTS1,...
        tsIntensities0(1),tsIntensities0(2),tsIntensities0(3),tsIntensities0(4),...
        tsIntensities1(1),tsIntensities1(2),tsIntensities1(3),tsIntensities1(4),...
        time,...
        nucIntens0,nucIntens1,nucIntens2,nucIntens3};

    %             % Plot cell images
    if makeFigures
        iCell
        hold off
        scatter3(thisCell0(~coLocalSpots0,1),thisCell0(~coLocalSpots0,2),thisCell0(~coLocalSpots0,3)); hold on
        scatter3(thisCell1(~coLocalSpots1,1),thisCell1(~coLocalSpots1,2),thisCell1(~coLocalSpots1,3));
        scatter3(thisCell0(coLocalSpots0,1),thisCell0(coLocalSpots0,2),thisCell0(coLocalSpots0,3),'ko'); hold on
        scatter3(thisCell1(coLocalSpots1,1),thisCell1(coLocalSpots1,2),thisCell1(coLocalSpots1,3),'kx'); hold on

        plot3(thisCell0(jClusters0,1),thisCell0(jClusters0,2),thisCell0(jClusters0,3),'co','MarkerSize',20,'MarkerFaceColor','c'); hold on
        plot3(thisCell1(jClusters1,1),thisCell1(jClusters1,2),thisCell1(jClusters1,3),'mo','MarkerSize',15,'MarkerFaceColor','m'); hold on
        pause
    end
end

CellData = CellData(keepCells,:);
end


