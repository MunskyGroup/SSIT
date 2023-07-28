AfspFull = ssit.FspMatrix(propensities, stateSpace, constraintCount, varNames);
jacFull = AfspFull.createSingleMatrix(outputTimes(1));
jacFull = jacFull(1:end-constraintCount,1:end-constraintCount);

solVecFull = zeros(stateCount, 1);
solVecFull(1:size(initStates,2)) = initProbs;

% n = min(length(jacFull),100);
% 
% % [PHI,D] = eig(full(jacFull-diag(sum(jacFull))));
% % [PHI,D] = eig(full(jacFull));
% [PHI,D] = eigs((jacFull),n,0);
% % [PHI,D] = eigs((jacFull-diag(sum(jacFull))),47,0,'Tolerance',1e-16);
% [rD,I] = sort(real(diag(D)),'descend');
% PHI = PHI(:,I(1:n));
% % PHI = orth([solVecFull,jacFull*solVecFull,jacFull*(jacFull*solVecFull),PHI]);
% 
% nKr = 1;
% V = zeros(size(solVecFull,1),nKr);
% V(:,1) = solVecFull;
% for i=2:nKr
%     V(:,i) = jacFull*V(:,i-1);
% end
% PHI = orth([V,PHI]);

PHI = modRedTransformMatrices.phi;

% PHI = orth([PHI]);
jacEig = PHI'*jacFull*PHI;
solVecEig = PHI'*solVecFull;


% Compare later time
figure(2); clf
ta = [0,0.1,1,3,10];
for i=1:length(ta)
    t = ta(i);
    Pfull = expm(jacFull*t)*solVecFull;
    Pred = real(PHI*expm(jacEig*t)*solVecEig);
    Pred = Pred/sum(Pred);
    Pred = max(0,Pred);
    subplot(length(ta),2,2*i-1); scatter(Pfull,Pred)
    subplot(length(ta),2,2*i); plot(Pfull); hold on;  plot(Pred,'--');
end



%%
modRedTransformMatrices
whos PH*
figure(3)
plot(abs(PHI(:,1)));
hold on
plot(abs(modRedTransformMatrices.phi(:,1)),'--');

figure(4)
plot(abs(PHI(:,5)));
hold on
plot(abs(modRedTransformMatrices.phi(:,5)),'--');


