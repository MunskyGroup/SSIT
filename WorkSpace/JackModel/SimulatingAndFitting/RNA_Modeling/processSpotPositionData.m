function [intensTS,binCounts,abEllipse,areaEllipse,ecentricity,binTS,binCountsDim,binCountsBright,hisCountsBright,hisCountsDim,XcRot,iTS,XCVcRot,XYfit] = ...
    processSpotPositionData(dataFile,maskFile,a,b,makePlots,nBins,minTS)
arguments
    dataFile
    maskFile = [];
    a =[];
    b =[];
    makePlots = false;
    nBins = 11;
    minTS =4;
end

if ~isempty(maskFile)
    %% Load Mask Data
    mask = load(maskFile);
elseif isempty(a)||isempty(b)
    error('need to provide masks or (a,b)')
end

%% Load Spot Data
% dataFile = '../../../data/dataframes_final/dataframe_MS2-CY5_Cyto543_560_woStim.csv';
% maskFile = '../../../data/mask_polygon/polygons_wo.mat';
X0 = importdata(dataFile);
iType = find(strcmp(X0.colheaders,'spot_type'));
iXpos = find(strcmp(X0.colheaders,'x'));
iYpos = find(strcmp(X0.colheaders,'y'));
iZpos = find(strcmp(X0.colheaders,'z'));
iCellNum = find(strcmp(X0.colheaders,'cell_id'));
iImagNum = find(strcmp(X0.colheaders,'image_id'));
iSpotNum = find(strcmp(X0.colheaders,'spot_id'));
iSpotIntens = find(strcmp(X0.colheaders,'spot_int_ch_3'));
iClusterSize = find(strcmp(X0.colheaders,'cluster_size'));
iFragmented = strcmp(X0.colheaders,'is_cell_fragmented');
if sum(iFragmented)>0
    jSpots0 = X0.data(:,iType)==0&X0.data(:,iFragmented)~=1;
else
    jSpots0 = X0.data(:,iType)==0;
end
allSpots0 = X0.data(jSpots0,:);
nCells = max(allSpots0(:,iCellNum));

%% Initialize data fields for processing results
intensTS = zeros(1,nCells+1);
numberSpots = zeros(1,nCells+1);
binTS = zeros(1,nCells+1);
areaEllipse = zeros(1,nCells+1);
abEllipse = zeros(2,nCells+1);
binCounts = zeros(nBins,nCells+1);
histBins = 20;
hisCountsBright = zeros(nBins,histBins);
hisCountsDim = zeros(nBins,histBins);

for iCell= 0:nCells  %25

    %% Look at specific cell
    jThisCell0 = (allSpots0(:,iCellNum)==iCell);
    thisCell0 =  allSpots0(jThisCell0,[iXpos,iYpos,iZpos,iClusterSize,iSpotIntens]);
    nSpotsThisCell0 = sum(jThisCell0);
    X = thisCell0(:,[2,1]);

    %% Get index of TS
    if ~isempty(X)
        [intensTS(iCell+1),iTS] = max(thisCell0(:,4));
        numberSpots(iCell+1) = size(thisCell0,1);


        %% Find Convex Hull using Mask Polygon
        if ~isempty(maskFile)
            maskXY = squeeze(mask.polygons_wo(iCell+1,:,:));
            XConvHull = convhull(maskXY);
            JConvHull = unique(XConvHull(:));
            XCV = maskXY(JConvHull,:);
        else
            u=-pi:0.1:pi;
            if length(a)==1
                XCV(:,2)=a*cos(u);
                XCV(:,1)=b*sin(u);
            else
                XCV(:,2)=a(iCell+1)*cos(u);
                XCV(:,1)=b(iCell+1)*sin(u);
            end
        end

        %% Center the mask at Zero
        cen = mean(XCV);
        Xc = X-cen;
        XCVc = XCV-cen;

        %% Rotate so major axis of mask is in first dimension.
        [~,~,V] = svd(XCVc);
        v1 = V(:,1)';
        v1 = v1/norm(v1);
        v23 = null(v1);
        PHI = [v1',v23];
        XcRot = Xc*PHI;
        XCVcRot = XCVc*PHI;

        %% Flip image so TS is on left side
        if (intensTS(iCell+1)>=minTS)&&(XcRot(iTS,1)>0)
            XcRot(:,1)=-XcRot(:,1);
            XCVcRot(:,1)=-XCVcRot(:,1);
        end


        if makePlots
            figure(1); clf;
            plot(XcRot(:,1),XcRot(:,2),'o'); hold on
            plot(XCVcRot(:,1),XCVcRot(:,2),'k--','lineWidth',1.5);
            if intensTS(iCell+1)>=minTS
                plot(XcRot(iTS,1),XcRot(iTS,2),'rx','MarkerSize',12,'lineWidth',2);  %scatter
            end
        end

        %% Fit convex hull with an ellipse
        beta = (XCVcRot.^2)\ones(length(XCVcRot),1);
        u=-pi:0.1:pi;
        XYfit=[1/sqrt(beta(1))*cos(u);1/sqrt(beta(2))*sin(u)]';

        if makePlots
            plot(XYfit(:,1),XYfit(:,2),'-','lineWidth',1.5)
        end
        %         h = legend('RNA Spots','Cell Mask','TS','Ellipse Approx.');
        %         set(h,'FontSize',20);
        %         title(['RNA Spatial Representation for Cell 18'],'FontSize',20)

        %% Identify bright and dim spots
        % 5 = the index of the intensity of spots.
        Jdim = thisCell0(:,5)<=median(thisCell0(:,5));
        Jbright = thisCell0(:,5)>median(thisCell0(:,5));

        %% Divide into bins along major axis
        binEdges = linspace(min(XYfit(:,1)),max(XYfit(:,1)),nBins+1);
        binEdges(1) = -inf; binEdges(end) = inf;

        for i=1:nBins
            binCounts(i,iCell+1) = sum((XcRot(:,1)>=binEdges(i)&XcRot(:,1)<=binEdges(i+1)));
            binCountsDim(i,iCell+1) = (sum((XcRot(Jdim,1)>=binEdges(i)&XcRot(Jdim,1)<=binEdges(i+1))));
            binCountsBright(i,iCell+1) = (sum((XcRot(Jbright,1)>=binEdges(i)&XcRot(Jbright,1)<=binEdges(i+1))));
        end

        if makePlots
            figure(3);
            stairs(binCounts(:,iCell+1)/sum(binCounts(:,iCell+1))); hold on;
        end

        %% Find Bin Position of TS for Experimental Data
        % Find XTS and XTSc
        if intensTS(iCell+1)>=minTS
            XTScRot = XcRot(iTS,:);

            % Find bin for TS
            for i=1:nBins
                if (XTScRot(1,1)>=binEdges(i)&&XTScRot(1,1)<=binEdges(i+1))
                    binTS(iCell+1) = i;
                    break
                end
            end
        else
            binTS(iCell+1)=NaN;
        end

        %% Record some features of the ellispe
        abEllipse(:,iCell+1) = 1./sqrt(beta);
        a = max(abEllipse(:,iCell+1));
        b = min(abEllipse(:,iCell+1));
        areaEllipse(iCell+1) = pi*sqrt(a*b);
        ecentricity(iCell+1) = sqrt(1-b^2/a^2);

    else
        abEllipse(:,iCell+1) = [NaN;NaN];
        areaEllipse(iCell+1) = NaN;
        ecentricity(iCell+1) = NaN;
        binCounts(1:nBins,iCell+1)=NaN;
        binTS(iCell+1)=NaN;
    end


    % Find the distance from the mRNAs to the TS.
    % To be counted the cell must have a TS and one or more mRNA.
    if (intensTS(iCell+1)>=minTS)&&(size(X,1)>=2)
        Dists = squareform(pdist(X));
        JTS = iTS(1,1);
        Dists2TS = Dists(JTS,:);
        hisCountsBright(iCell+1,:) = histcounts(Dists2TS(Jbright),linspace(0,100,histBins+1)); %bright spots only
        hisCountsDim(iCell+1,:) = histcounts(Dists2TS(Jdim),linspace(0,100,histBins+1)); %dim spots only

        if makePlots
                figure(4)
                subplot(1,2,1)
                stairs(linspace(0,100,histBins+1),[0,hisCountsBright(iCell+1,:)]) %bright spots only
                title('Bright Spots')
                xlabel('Distance away from TS')
                ylabel('Number of RNA')
                xlim([0 100])
                
                subplot(1,2,2)
                stairs(linspace(0,100,histBins+1),[0,hisCountsDim(iCell+1,:)]) %bright spots only
                title('Dim Spots')
                xlabel('Distance away from TS')
                ylabel('Number of RNA')
                xlim([0 100])

                figure(5); hold off;
                plot(linspace(0,100,histBins+1),...
                    cumsum([0,hisCountsBright(iCell+1,:)])/sum(hisCountsBright(iCell+1,:)),'LineWidth',3); hold on;
                plot(linspace(0,100,histBins+1),...
                    cumsum([0,hisCountsDim(iCell+1,:)])/sum(hisCountsDim(iCell+1,:)),'LineWidth',3); hold on;
                legend('Bright Spots','Dim Spots')
        end

    else
        hisCountsBright(iCell+1,:) = NaN; %bright spots only
        hisCountsDim(iCell+1,:) = NaN; %bright spots only
    end

    

end



