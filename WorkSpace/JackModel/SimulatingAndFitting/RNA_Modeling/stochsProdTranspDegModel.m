function stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

nSpots = 0;
for iCell = Ncells:-1:1 
    %                           kon,koff,w,kex,kr,D,gam
    %     [spotsX,spotsY] = simRNAposn(
%     kon = 1e-4;
%     koff = 2.6e-4;
%     w = 0.1;
%     kex = 10;
%     kr = 20;
%     %     [spotsX,spotsY] = simRNAposn(1e-4,2.6e-4,0.1,10,20,[0.001,0.08,0.08],[0.04;0.001;0.001]);
% 
%     D = [0.01,0.8,0.8];
%     gam =[0.04;0.001;0.001];
%     a = 55;
%     b = 40;
% 
%     u = 2*pi*rand;
%     posnTS=[a*cos(u),b*sin(u)].*rand(1,2);
%     makePlot = 0;
    
    [spotsX,spotsY] = simRNAposn(kon,koff,w,kex,kr,D,gam,posnTS(iCell,:),makePlot,a(iCell),b(iCell));
    
    for i=1:size(spotsX,2)
        spotsX1 = spotsX(:,i);
        spotsY1 = spotsY(:,i);
        spotsX1 = spotsX1(~isnan(spotsX1));
        spotsY1 = spotsY1(~isnan(spotsY1));
        if i>1
            cell_id = (iCell-1)*ones(size(spotsX1));
            spot_ch = 0*ones(size(spotsX1));
            spot_ty = 0*ones(size(spotsX1));
            spot_int_ch_0 = (size(spotsX,2)+1-i)*ones(size(spotsX1));
            clus_size = 0*ones(size(spotsX1));
        else % TS
            clus_size = length(spotsX1);
            spotsX1 = mean(spotsX1,1);
            spotsY1 = mean(spotsY1,1);
            cell_id = (iCell-1);
            spot_ch = 0;
            spot_ty = 0;
            spot_int_ch_0 = clus_size;
        end

        spot_id = [nSpots+1:nSpots+length(spotsX1)]';
        nSpots=nSpots+length(spotsX1);

        %export data as a csv file
        Tnew{iCell,i} = [spotsX1,spotsY1,zeros(size(spotsX1)),cell_id,spot_int_ch_0,spot_id,spot_ty,clus_size];
    end
end  

TabData = zeros(nSpots,8);
nSpots = 0;
for iCell = 1:Ncells
    for i=1:size(spotsX,2)
        m = size(Tnew{iCell,i},1);
        TabData(nSpots+1:nSpots+m,:)=Tnew{iCell,i};
        nSpots = nSpots+m;
    end
end

T = table(TabData(:,1),TabData(:,2),TabData(:,3),TabData(:,4),TabData(:,5),TabData(:,6),TabData(:,7),TabData(:,8),...
    'VariableNames', {'x','y','z','cell_id','spot_int_ch_3','spot_id','spot_type','cluster_size'});

writetable(T,fileName)
