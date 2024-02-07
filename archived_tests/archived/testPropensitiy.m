w1str = 'x1*x1'; 
w2str = '(1+sin(t))*x2';
w3str = '1/(1 + cos(t)*x2)';

propensity1 = ssit.Propensity.createFromSym(str2sym(w1str), [1 0]', 1);
propensity2 = ssit.Propensity.createFromSym(str2sym(w2str), [0 -1]', 2);
propensity3 = ssit.Propensity.createFromSym(str2sym(w3str), [0 1]', 3);

% Test that the propensities have correct types
assert(propensity1.isTimeDependent == false);
assert(propensity2.isTimeDependent == true);
assert(propensity3.isTimeDependent == true);
assert(propensity2.isFactorizable == true);
assert(propensity3.isFactorizable == false);

% Test the evaluation with both t and x
w1 = str2func(['@(t,x1,x2) ' w1str]);
w2 = str2func(['@(t,x1,x2) ' w2str]);
w3 = str2func(['@(t,x1,x2) ' w3str]);
pass = true;
for t = 0:0.1:1
for x1 = 0:20
    for x2 = 0:20
        err = abs(w1(t,x1,x2) - propensity1.eval(t, [x1 x2]')) + ...
            abs(w2(t,x1,x2) - propensity2.eval(t, [x1 x2]')) + ...
            abs(w3(t,x1,x2) - propensity3.eval(t, [x1 x2]'));
        if (err > eps)
            pass = false;
            break;
        end
    end
end
end
assert(pass);
%% Test the evaluation with only x
w1 = str2func(['@(x1,x2) ' w1str]);
w2 = str2func(['@(x1,x2) ' 'x2']);
pass = true;
for x1 = 0:20
    for x2 = 0:20
        err = abs(w1(x1,x2) - propensity1.evalStateFactor([x1 x2]')) + ...
            abs(w2(x1,x2) - propensity2.evalStateFactor([x1 x2]'));
        if (err > eps)
            pass = false;
            break;
        end
    end
end
assert(pass);
%%
fprintf('Propensity tests passed!\n');
