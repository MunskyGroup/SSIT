% In this test, we will use a 6x6 state space of the toggle switch to check
% that the FiniteStateSet object is working properly.

stoichMatrix = [   1  0;
                 1  0;
                -1  0;
                 0  1;
                 0  1;
                 0 -1;]';

t_final = 100;
x0 = [0 0]';

constraintFunctions = @(x) [-x(1,:); -x(2,:); x(1,:); x(2,:)];
constraintBounds = [ 0, 0, 2, 2]';
b_constraints_next = [ 0, 0, 5, 5]';

%% Test state exploration from the initial state
my_set = ssit.FiniteStateSet(x0, stoichMatrix);
my_set = my_set.expand(constraintFunctions, constraintBounds);

% Test for the total number of states added
assert(size(my_set.states, 2) == 9);
% Test for the inclusion of states
pass = true;
for i = 0:2
    for j = 0:2
        x = [i; j];
        if (~ismember(x', my_set.states', 'rows'))
            pass = false;
            break;
        end
    end
end
assert(pass);

%% Test state exploration from an existing set of states
my_set = ssit.FiniteStateSet(x0, stoichMatrix);
my_set = my_set.expand(constraintFunctions, constraintBounds);
my_set = my_set.expand(constraintFunctions, b_constraints_next);

% Test for the total number of states added
assert(size(my_set.states, 2) == 36);
% Test for the inclusion of states
pass = true;
for i = 0:5
    for j = 0:5
        x = [i; j];
        if (~ismember(x', my_set.states', 'rows'))
            pass = false;
            break;
        end
    end
end
assert(pass);
%%
disp('Finite State Set tests passed!');