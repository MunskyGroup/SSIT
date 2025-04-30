%% Function for Simulated RNA
function [spotsX, spotsY] = simRNAposn(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b)
arguments
    kon
    koff
    w
    kex
    kr
    D
    gam
    posnTS = [-1,0.2];    % X/Y positions of Transcription sites.
    makePlot = false
    a = 3;
    b = 2;
end

%% Simulation parameters
t = 0;
tfin = 1000;
plotTime = linspace(0,tfin,40);
iPlot = 1;
colors = {'r','b','c','m'};

%% Cell ellipse parameters
u=-pi:0.1:pi;
xEllipse=a*cos(u);
yEllipse=b*sin(u);
if makePlot
    plot(xEllipse,yEllipse,'+')
end

%% Transcription site positions and states
nTS = 1;

sTS = [-1 1 0 0;...
    1 -1 -1 1;...
    0 0  1 -1];
M = size(sTS,2);
N = size(sTS,1);
STS = sTS;
for i=2:nTS
    STS = blkdiag(STS,sTS);
end
funPropsTS = @(x)[kon*x(1);koff*x(2);w*x(2);kex*x(3)];

%% Randomize initial condition based on CME solution.
if rand<kon/(kon+koff)
    xTS = [0;1;0];
else
    xTS = [1;0;0];
end

%% Transcription events;
funPropsTr = @(x)kr*x(3);

%% mRNA Events
nDegSteps = length(gam);
% sRNA = -eye(nDegSteps);
% sRNA(2:end,1:end-1) = sRNA(2:end,1:end-1)+eye(nDegSteps-1);

%%
spotsX = nan*ones(2000,nDegSteps);
spotsY = nan*ones(2000,nDegSteps);

nSpots = zeros(length(gam),1);

recalcTS = 1;
recalcRNA = 1;

while t<tfin
    %% Propensities
    if recalcTS
        for i=1:nTS
            % Gene Transition events
            propsTS((i-1)*M+1:i*M,1) =  funPropsTS(xTS((i-1)*N+1:i*N));
            % Gene Transition events
            propsKR(i,1) =  funPropsTr(xTS((i-1)*N+1:i*N));
        end
        recalcTS=0;
    end

    % Degradation events
    if recalcRNA
        propsRNA = gam.*nSpots;
        recalcRNA=0;
    end

    % Choose next time
    propsAll = [propsTS;propsKR;propsRNA];
    w0 = sum(propsAll);
    delt = -log(rand)/w0;

    if t+delt<tfin
        t = t + delt;
        % Choose next reaction
        j = 1;
        cdf = propsAll(j);
        r = rand*w0;
        while cdf<r
            j=j+1;
            cdf = cdf+propsAll(j);
        end

        if j<=M*nTS
            % TS reactions
            xTS = xTS + STS(:,j);
            recalcTS=1;

        elseif j<=(M*nTS+nTS)
            % Production events
            spotsX(nSpots(1)+1,1)=posnTS(j-M*nTS,1);
            spotsY(nSpots(1)+1,1)=posnTS(j-M*nTS,2);
            nSpots(1) = nSpots(1)+1;
            recalcRNA=1;
        else
            % Degradation events
            ideg = j - (M*nTS+nTS);
            removeSpot = randi(nSpots(ideg));
            keepSpots = [1:removeSpot-1,removeSpot+1:nSpots(ideg)];
            % add one of the next spot type
            if ideg<nDegSteps
                spotsX(nSpots(ideg+1)+1,ideg+1)=[spotsX(removeSpot,ideg)];
                spotsY(nSpots(ideg+1)+1,ideg+1)=[spotsY(removeSpot,ideg)];
                nSpots(ideg+1) = nSpots(ideg+1)+1;
            end
            % remove one of the chosen spot type
            spotsX(1:nSpots(ideg),ideg)=[spotsX(keepSpots,ideg);NaN];
            spotsY(1:nSpots(ideg),ideg)=[spotsY(keepSpots,ideg);NaN];
            nSpots(ideg) = nSpots(ideg)-1;
            recalcRNA=1;
        end
    else
        delt = tfin-t;
        break
    end

    if max(~isnan(spotsX(1,:)))>0
        for i=1:nDegSteps
            spotsX1 = spotsX(:,i);
            spotsY1 = spotsY(:,i);
            spotsX1(~isnan(spotsX1)) = spotsX1(~isnan(spotsX1))+D(i)*sqrt(delt)*randn(size(spotsX1(~isnan(spotsX1))));
            spotsY1(~isnan(spotsY1)) = spotsY1(~isnan(spotsY1))+D(i)*sqrt(delt)*randn(size(spotsY1(~isnan(spotsY1))));
            spotsX(:,i) = spotsX1;
            spotsY(:,i) = spotsY1;
        end
    end

    [spotsX,spotsY] = reflectSpots(spotsX,spotsY,a,b);

    if makePlot
        while t>plotTime(iPlot)
            clf;
            if makePlot
                plot(xEllipse,yEllipse,'+'); hold on
                plot(posnTS(:,1),posnTS(:,2),'cd','MarkerEdgeColor','c','MarkerSize',15)
            end
            for i=1:nDegSteps
                spotsX1 = spotsX(:,i);
                spotsY1 = spotsY(:,i);
                if makePlot
                    plot(spotsX1(~isnan(spotsX1)),spotsY1(~isnan(spotsY1)),'o','MarkerFaceColor',colors{i});

                end
            end
            drawnow
            iPlot = iPlot+1;
            title(['t =',num2str(plotTime(iPlot))])
        end
    else
        while t>plotTime(iPlot)
            for i=1:nDegSteps
                spotsX1 = spotsX(:,i);
                spotsY1 = spotsY(:,i);
            end
            iPlot = iPlot+1;
        end
    end
end
% % Save final data
% for i=1:nDegSteps
%     spotsX1 = spotsX(:,i);
%     spotsY1 = spotsY(:,i);
%     spotsX1 = spotsX1(~isnan(spotsX1));
%     spotsY1 = spotsY1(~isnan(spotsY1));
%     if i==1
%         % spot number, x, y, spot intensity, cluster size
%         v = [mean(spotsX1),mean(spotsY1),length(spotsX1),length(spotsX1)];
%     else
%         v = [v;spotsX1,spotsY1,(nDegSteps+1-i)*ones(size(spotsX1)),zeros(size(spotsX1))];
%     end
% end
% v = [[0:size(v,1)-1]',v];

    function [X,Y] = reflectSpots(X,Y,a,b)
        J = find(~isnan(X)&(X.^2/a^2+Y.^2/b^2>1));

        if ~isempty(J)
            theta = acos((X(J)/a)./(sqrt(X(J).^2/a^2+Y(J).^2/b^2)));
            theta = theta.*sign(Y(J));
            newr = sqrt(2-(X(J).^2/a^2+Y(J).^2/b^2));
            X(J) = real(a*newr.*cos(theta));
            Y(J) = real(b*newr.*sin(theta));
        end
        %     xEllipse=a*cos(u);
        %     yEllipse=b*sin(u);
    end


end