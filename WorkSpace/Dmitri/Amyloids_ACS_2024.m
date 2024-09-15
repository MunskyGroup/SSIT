clear all
addpath(genpath('../../src'));

Model = SSIT;
Model.species = {};
Model.initialCondition = [];
Model.propensityFunctions = {};
Model.propensityFunctions = {};
Model.stoichiometry = [];
Model.parameters = {};
Model.tSpan = linspace(0,2,21);

%% Define limits of nucleus size, reaction order of secondary nucleation:
MIN_NUCLEUS_SIZE = 3;
assert(MIN_NUCLEUS_SIZE > 1, ...
    'Minimum nucleus size must be more than a monomer!')
MAX_NUCLEUS_SIZE = 10;
assert(MAX_NUCLEUS_SIZE >= MIN_NUCLEUS_SIZE, ...
    'Maximum nucleus size must be no smaller than minimum size!')
nc = MIN_NUCLEUS_SIZE;
n2 = 2;
assert(n2 > 0, 'Secondary nucleation must proceed via a nucleus!')

%% Define initial monomer population:
INITIAL_MONOMER_POPULATION = 100;
assert(INITIAL_MONOMER_POPULATION >= MIN_NUCLEUS_SIZE, ...
    'There must be sufficiently many monomers to form one nucleus')

%% Define whether we are in "aggregate mass" mode:
% Experimentally, this is what is typically measured; namely, the total
% aggregate mass concentration as a function of time. The aggregate mass is
% defined as M(t) = sum(i*xi) for all i in [nc, MAX_NUCLEUS_SIZE].

aggregate_mass_mode = true;

%% Add species: xk denotes a k-mer, from 1 up to the maximum considered.

for i = 1 : MAX_NUCLEUS_SIZE
    if i == 1
        Model = Model.addSpecies(['x', num2str(i)], ...
            INITIAL_MONOMER_POPULATION);
    else
        Model = Model.addSpecies(['x', num2str(i)], 0);
    end
end

if aggregate_mass_mode
    Model = Model.addSpecies('M', 0);
end

%% Add parameters

Model.stoichiometry = [];

Model = Model.addParameter({ ...
    'kn', 0.1; ... % Primary nucleation
    'k2', 0.2; ... % Surface-catalyzed secondary nucleation
    'kplus', 0.3; ... % Elongation (at end of linear chain)
    'kminus',0.1 ... % Fragmentation (at any interior point of linear chain)
    });

%% Add reactions

stoichTempl = zeros(size(Model.species, 1), 1);

% Primary nucleation reaction
stoichPN = stoichTempl;
stoichPN(1) = -1 * nc; % n_c monomers are consumed to make the nucleus
stoichPN(nc) = 1; % One oligomer of size n_c is produced
if aggregate_mass_mode
    % The total aggregate mass increases by nc (one aggregate of mass nc):
    stoichPN(end) = nc;
end
Model = Model.addReaction(['kn*x1^', num2str(nc)], stoichPN);

% Fragmentation reactions
% Any nucleus of length [nc, MAX_NUCLEUS_SIZE] can be fragmented at any
% interior point in the chain. This will split the existing oligomer into
% two oligomers whose combined length will equal that of the original 
% chain.

for l = nc : MAX_NUCLEUS_SIZE
    for i = 1 : floor(l/2)
        stoichF = stoichTempl;
        stoichF(i) = 1; % One oligomer of length i is produced
        stoichF(l-i) = 1; % One oligomer of length l-i is produced
        stoichF(l) = -1; % One oligomer of length l is consumed
        if aggregate_mass_mode
            % The total aggregate mass changes based on the oligomers
            % produced and consumed and whether they are sufficiently
            % large to be nuclei. Since l >= nc, we know that the mass will
            % decrease by that much, but it may also increase based upon
            % the sizes of the product nuclei:

            delta_aggregate_mass = -l;
            if i >= nc
                delta_aggregate_mass = delta_aggregate_mass + i;
            end
            if (l-i) >= nc
                delta_aggregate_mass = delta_aggregate_mass + (l-i);
            end
            stoichF(end) = delta_aggregate_mass;
        end

        % The propensity will depend on the number of locations in the
        % chain where an equivalent fragmentation could occur. For an
        % exact bifurcation, there is only a single suitable point, e.g.
        % only one cut will split a 6-mer into two 3-mers. For inexact
        % bifurcations, e.g. a 7-mer into one 2- and one 5-mer, there are
        % two cuts:

        if i < (l/2)
            propensity = ['2*kminus*x', num2str(l)];
        else
            propensity = ['kminus*x', num2str(l)];
        end

        Model = Model.addReaction(propensity, stoichF);
    end % for i = 1 : floor(l/2) ...
end % for l = nc : MAX_NUCLEUS_SIZE ...

% Secondary nucleation reactions
% Any nucleus of length [nc, MAX_NUCLEUS_SIZE - n2] can grow via the 
% attachment of another nucleus at any interior point in the chain.
% As in the case of fragmentation reactions, the number of possible
% attachment sites depends on the length of the existing chain. 

for l = nc : (MAX_NUCLEUS_SIZE - n2)
    stoichSN = stoichTempl;      
    stoichSN(l) = -1; % One oligomer of length l is consumed
    stoichSN(1) = -n2; % n2 monomers are consumed
    stoichSN(l+n2) = 1; % One oligomer of length l+n2 is produced
    if aggregate_mass_mode
        % The total aggregate mass increases by n2 (one aggregate of mass
        % l is replaced by one of mass l+n2):
        stoichSN(end) = n2;
    end

    % The propensity will depend on the number of locations in the
    % chain where an equivalent attachment could occur. We assume that
    % chains are linear, before and after; therefore, there is
    % symmetry throughout, and the number of equivalent attachment 
    % sites equals the number of interior monomers, e.g. a 7-mer will
    % have five interior monomers.

    propensity = ['k2*(x1^', num2str(n2), ')*x', num2str(l), ...
        '*', num2str(l-2)];

    Model = Model.addReaction(propensity, stoichSN);
end % for l = nc : (MAX_NUCLEUS_SIZE - n2) ...

% Elongation reactions
% Any nucleus of length [nc, MAX_NUCLEUS_SIZE) can grow via the 
% attachment of a monomer at either endpoint of the chain. Regardless of
% the length of the existing chain, there are exactly two possible points
% of attachment (since the existing nucleus is at minimum a dimer, it has
% two distinct ends).

for l = nc : (MAX_NUCLEUS_SIZE - 1)
    stoichE = stoichTempl;      
    stoichE(l) = -1; % One oligomer of length l is consumed
    stoichE(1) = -1; % One monomer is consumed
    stoichE(l+1) = 1; % One oligomer of length l+1 is produced
    if aggregate_mass_mode
        % The total aggregate mass increases by 1 (one aggregate of mass
        % l is replaced by one of mass l+1):
        stoichE(end) = 1;
    end

    propensity = ['2*kplus*x1*x', num2str(l)];

    Model = Model.addReaction(propensity, stoichE);
end % for l = nc : (MAX_NUCLEUS_SIZE - 1) ...

Model.summarizeModel

%% Multi-model fitting:

replicates = 1:4;
initial_monomer_concentration = [5, 4, 3.5, 3, 2.5, 2, 1.75, 1.5, 1.35, 1.1];
initial_monomer_population = 1000 * initial_monomer_concentration;
ModelList = cell(1, length(initial_monomer_population));
ParameterList = cell(1, length(initial_monomer_population));

%% Create models and load data:

% We will create a series of models that will all share the same
% parameters. The difference will be that they will different initial
% monomoer populations.

for counter = 1 : length(initial_monomer_population)
    % Clone the template model:

    newModel = Model;
    
    % Update the initial monomer population:

    newModel.initialCondition(1) = ...
        initial_monomer_population(counter);

    % For now, we will only use data from the first replicate:

    column_to_use = ...
        ['Ab42_rep1_', num2str(initial_monomer_concentration(counter)), 'uM'];
    newModel = newModel.loadData('AmyloidData.csv', ...
        {'M', column_to_use});

    % Add the model to the list of all models, and add a corresponding
    % entry of all parameters to that list:

    ModelList{counter} = newModel;
    ParameterList{counter} = 1:length(newModel.parameters);
end
%%
newModel.solve
%% Set fitting options:
fitAlgorithm = 'fminsearch';
fitOptions = optimset('Display', 'final', 'MaxIter', 500);

%% Combine all models:

combinedModelDependent = SSITMultiModel(ModelList, ParameterList);
combinedModelDependent = combinedModelDependent.initializeStateSpaces;
allParsDependent = ([Model.parameters{:,2}]);
allParsDependent = combinedModelDependent.maximizeLikelihood(...
    allParsDependent, fitOptions, fitAlgorithm);
combinedModelDependent = combinedModelDependent.updateModels(allParsDependent);