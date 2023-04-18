%% findCoLocalSpotsNucIntens
% This functions loads the processed smiFISH data and determines how many
% spots are in each color channel and how many are collocalized.  It also
% makes a table of the spot numbers and nuclear fluorescnece intensities
% for each cell. 
clear all
timeArray = [0,18,300];

psf_z=350;               % Theoretical size of the PSF emitted by a [rna] spot in the z plan, in nanometers
psf_yx=160;              % Theoretical size of the PSF emitted by a [rna] spot in the yx plan, in nanometers
voxel_size_z=500;        % Microscope conversion px to nanometers in the z axis.
voxel_size_yx=160;       % Microscope conversion px to nanometers in the xy axis.

% the code allows for different scaling of the x-, y-, and z- dimensions
% for the distance calculations, but here we will use pixels.
% distScale = ([voxel_size_yx/psf_yx, voxel_size_yx/psf_yx, voxel_size_z/psf_z]);
distScale = [1 1 1];

% makeFigures = false;
% tArray = [400:50:550];

% Distance threshold for collocalization detection.  Larger numbers are
% more permissive and will results in mor collocalized spots.  We choose 2
% here to allow for the spots to be separated by one pixel in each
% direction.  Allowable distances would therefore be: x<[1,1,1] and
% x<([2,0,0],[0,2,0],[0,0,2])
dThresh = 2;
allCellData = table;

for iTime = 0:2
    time = timeArray(iTime+1);
    switch iTime
        case 0
%             STR0 = ['Huy_intensity_data_correct/complete__MS2_CY5_time_0_int_550_400.csv'];
            STR0 = ['Huy_intensity_data_correct/dataframe_MS2-CY5_Cyto543_560_woStim.csv'];
        case 1
%             STR0 = ['Huy_intensity_data_correct/complete__MS2_CY5_time_0_int_550_400.csv'];
            STR0 = ['Huy_intensity_data_correct/dataframe_MS2-CY5_Cyto543_560_18minTPL_5uM.csv'];
        case 2
%             STR0 = ['Huy_intensity_data_correct/complete__MS2_CY5_time_300_int_550_400.csv'];
%             STR0 = ['Huy_intensity_data_correct/complete__MS2_CY5_time_300_int_550_400_removeDuplicateCell.csv'];
            STR0 = ['Huy_intensity_data_correct/dataframe_MS2-CY5_Cyto543_560_5hTPL_5uM.csv'];
    end
    X0 = importdata(STR0);
    CellData = findColocalzeData(X0,X0,distScale,dThresh,time,false);

    allCellData = [allCellData;CellData];

    total(iTime+1) = sum([CellData{:,9}]);
    both(iTime+1) = sum([CellData{:,6}]);

end
allCellData.Properties.VariableNames = {'CellNum',...
    'number_spots_type_0','number_spots_type_1',...
    'number_spots_type_0_only','number_spots_type_1_only',...
    'number_spots_type_0_1_both',...
    'nDoubleSpots0','nDoubleSpots1',...
    'total_spots',...
    'numberTS0','numberTS1',...
    'tsIntensities0_1','tsIntensities0_2','tsIntensities0_3','tsIntensities0_4',...
    'tsIntensities1_1','tsIntensities1_2','tsIntensities1_3','tsIntensities1_4',...
    'time',...
    'nucIntens0','nucIntens1','nucIntens2','nucIntens3'};

FN = ['Huy_intensity_data_correct/NucAndSpotClassification_dTS_',num2str(dThresh,3),'.csv'];
% writetable(allCellData, FN)


for tt = [0,18,300]
    disp(['At time t=',num2str(tt)])
    J = allCellData.time==tt;
    nCells = sum(J)
    nSpotsTotal = sum(allCellData.total_spots(J))
    nMCPOnly = sum(allCellData.number_spots_type_0_only(J))
    nFISHOnly = sum(allCellData.number_spots_type_1_only(J))
    nBoth = sum(allCellData.number_spots_type_0_1_both(J))
    fMCPOnly = sum(allCellData.number_spots_type_0_only(J))/sum(allCellData.total_spots(J))
    fFISHOnly = sum(allCellData.number_spots_type_1_only(J))/sum(allCellData.total_spots(J))
    fBoth = sum(allCellData.number_spots_type_0_1_both(J))/sum(allCellData.total_spots(J))
end


% [a,iT0] = max(sum(both,3)./sum(total,3));
% [a,iT1] = max(a);
% T0 = tArray(iT0(iT1))
% T1 = tArray(iT1)
%
% %% Make Figures of cells
% T0 = 550;
% T1 = 400;
%
% iTime = 0;
% time = timeArray(iTime+1);
% STR0 = ['Huy_paper_data/MS2_CY5_100X__nuc_70__cyto_0__psfz_350__psfyx_160/ts_',num2str(T0),'_',num2str(T0),'_time_',num2str(iTime),'.csv'];
% X0 = importdata(STR0);
% STR1 = ['Huy_paper_data/MS2_CY5_100X__nuc_70__cyto_0__psfz_350__psfyx_160/ts_',num2str(T1),'_',num2str(T1),'_time_',num2str(iTime),'.csv'];
% X1 = importdata(STR1);
%
% findColocalzeData(X0,X1,distScale,dThresh,time,true);
%

function CellData = findColocalzeData(X0,X1,distScale,dThresh,time,makeFigures)

iType = find(strcmp(X0.colheaders,'spot_type'));
iXpos = find(strcmp(X0.colheaders,'x'));
iYpos = find(strcmp(X0.colheaders,'y'));
iZpos = find(strcmp(X0.colheaders,'z'));
iCellNum = find(strcmp(X0.colheaders,'cell_id'));
iImagNum = find(strcmp(X0.colheaders,'image_id'));
iSpotNum = find(strcmp(X0.colheaders,'spot_id'));
iClusterSize = find(strcmp(X0.colheaders,'cluster_size'));
iFragmented = strcmp(X0.colheaders,'is_cell_fragmented');
iNucIntCh0 = find(strcmp(X0.colheaders,'nuc_int_ch_0'));
iNucIntCh1 = find(strcmp(X0.colheaders,'nuc_int_ch_1'));
iNucIntCh2 = find(strcmp(X0.colheaders,'nuc_int_ch_2'));
iNucIntCh3 = find(strcmp(X0.colheaders,'nuc_int_ch_3'));

nImagID = max(X0.data(:,iImagNum));
cellNums = [];
maxCellNum = max(unique([X0.data(:,iCellNum);X1.data(:,iCellNum)]));
for i=0:nImagID
    for j = 1:length(cellNums)
        K0 = find((X0.data(:,iImagNum)==i&X0.data(:,iCellNum)==cellNums(j)));
        K1 = find((X1.data(:,iImagNum)==i&X1.data(:,iCellNum)==cellNums(j)));
        if ~isempty(K0)||~isempty(K1)
            maxCellNum=maxCellNum+1;
            X0.data(K0,iCellNum)=maxCellNum;
            X1.data(K1,iCellNum)=maxCellNum;
        end
    end
    cellNums = [cellNums;unique(X0.data(X0.data(:,iImagNum)==i,iCellNum))];
end

jSpots0 = X0.data(:,iType)==0&X0.data(:,iFragmented)~=1;
jSpots1 = X1.data(:,iType)==1&X1.data(:,iFragmented)~=1;

allSpots0 = X0.data(jSpots0,:);
allSpots1 = X1.data(jSpots1,:);

jNoSpots = X0.data(:,iType)==-1;
noSpots = X0.data(jNoSpots,:);

nCellsVector = unique([X0.data(X0.data(:,iFragmented)~=1,iCellNum);X1.data(X1.data(:,iFragmented)~=1,iCellNum)]);
nCells = length(nCellsVector);
% nCells = max(max(X0.data(:,iCellNum)),max(X1.data(:,iCellNum)))+1;

keepCells = zeros(1,nCells,'logical');

CellData = table;
for iCellL = 1:nCells
    iCell = nCellsVector(iCellL);
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
        keepCells(iCellL)=1;
    catch
        try
            nucIntens0 = round(thisCell1(1,5));
            nucIntens1 = round(thisCell1(1,6));
            nucIntens2 = round(thisCell1(1,7));
            nucIntens3 = round(thisCell1(1,8));
            keepCells(iCellL)=1;
        catch
            jThisCellNone = noSpots(:,iCellNum)==iCell;
            thisCellNone =  noSpots(jThisCellNone,[iXpos,iYpos,iZpos,iClusterSize,iNucIntCh0,iNucIntCh1,iNucIntCh2,iNucIntCh3]);
            try 
                nucIntens0 = round(thisCellNone(1,5));
                nucIntens1 = round(thisCellNone(1,6));
                nucIntens2 = round(thisCellNone(1,7));
                nucIntens3 = round(thisCellNone(1,8));
                keepCells(iCellL)=1;
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
    D = squareform(pdist([thisCell0(:,1:3);thisCell1(:,1:3)].*distScale))<dThresh;

    D01 = D(1:nSpotsThisCell0,nSpotsThisCell0+1:end);
    coLocalSpots1 = sum(D01,1)~=0;
    nSpots1Only = nSpotsThisCell1-sum(coLocalSpots1);

    coLocalSpots0 = sum(D01,2)~=0;
    nSpots0Only = nSpotsThisCell0-sum(coLocalSpots0);

    nSpotsBoth = min(sum(coLocalSpots0),sum(coLocalSpots1));

    % Double spots detection
    D00 = D(1:nSpotsThisCell0,1:nSpotsThisCell0);
    DoubleSpots0 = sum(D00,1)~=1;
    nDoubleSpots0 = sum(DoubleSpots0)/2;

    D11 = D(nSpotsThisCell0+1:end,nSpotsThisCell0+1:end);
    DoubleSpots1 = sum(D11,1)~=1;
    nDoubleSpots1 = sum(DoubleSpots1)/2;

%     nSpotsBoth = nSpotsBoth - max(sum(DoubleSpots1.*coLocalSpots1),sum(DoubleSpots0.*coLocalSpots0'));
%     nSpots1Only = nSpots1Only - sum(DoubleSpots1.*(~coLocalSpots1));
%     nSpots0Only = nSpots0Only - sum(DoubleSpots0.*(~~coLocalSpots0'));

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
    CellData(iCellL,:) = {iCell,...
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


