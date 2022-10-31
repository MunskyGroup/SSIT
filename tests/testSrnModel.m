stoichiometry = [1 0; -1 0;0 1;0 -1]';
propensity_strs = {
    'kr/(1 + kf*x2)',...
    'gr*x1',...
    'kp*x1',...
    'gp*x2'
    };
parameters = containers.Map({'kr', 'kf', 'gr', 'kp', 'gp'}, {1.0E-4, 11.0, 1.0E-3, 1.0E-2, 1.0E-1});

testSpace = ssit.FiniteStateSet([0, 0]', stoichiometry);
fConstraints = @(x) [
-x(1,:);
-x(2,:);
x(1,:);
x(2,:)
];
testSpace = testSpace.expand(fConstraints, [0,0,10,10]');
model = ssit.SrnModel(stoichiometry, propensity_strs, {'kr', 'kf', 'gr', 'kp', 'gp'}, []);
%% Test propensity creation
propensities = model.createPropensities(parameters);

kr_ = parameters("kr");
kf_ = parameters("kf");
gr_ = parameters("gr");
kp_ = parameters("kp");
gp_ = parameters("gp");
true_propensities = cell(4,1);
true_propensities{1} = @(x) kr_./(1 + kf_.*x(2,:));
true_propensities{2} = @(x) gr_.*x(1,:);
true_propensities{3} = @(x) kp_.*x(1,:);
true_propensities{4} = @(x) gp_.*x(2,:);

for ir = 1:4
    propensityEvaluations = propensities{ir}.evaluate(0, testSpace.states);
    propensityExactValues = true_propensities{ir}(testSpace.states);
    assert(prod(propensityExactValues == propensityEvaluations));
end
%% Test symbolic differentiation
dprops = model.findPropensityDerivativesSymbolic(parameters);
assert(prod(size(dprops) == [4, 5]));
trueSparsityPattern = [
    1 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1
];
for i = 1:4
    for j = 1:5
        if trueSparsityPattern(i,j) == 1
            assert(isa(dprops{i,j}, 'ssit.Propensity'));
        else 
            assert(isempty(dprops{i,j}));
        end
    end
end