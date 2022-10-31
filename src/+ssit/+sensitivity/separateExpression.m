function [ft, fx, isFactorizable] = separateExpression(expr)
    ft = [];
    fx = [];
    isFactorizable = [];
    % Check input
    left_bracket_loc = strfind(expr, '{');
    right_bracket_loc = strfind(expr, '}');

    if (length(left_bracket_loc) > length(right_bracket_loc))
        fprintf('Syntax error, missing ''}'' in one of the propensities.\n');
        return;
    elseif (length(left_bracket_loc) < length(right_bracket_loc))
        fprintf('Syntax error, missing ''{'' in one of the propensities.\n');
        return;
    end

    if (length(left_bracket_loc) > 2)
        fprintf('We currently only support separating propensity into two factor.\n');
        return;
    end

    if (isempty(left_bracket_loc))
        % Is it a function of only x?
        if (~contains(expr, 't'))
            isFactorizable = true;
            ft = '1';
            fx = expr;
            return;
        elseif (~contains(expr, 'x'))
            isFactorizable = true;
            ft = expr;
            fx = '1';
            return;
        else
            isFactorizable = false;
            ft = [];
            fx = [];
            return;
        end
    end
    %
    isFactorizable = true;
    f1 = expr(left_bracket_loc(1)+1:right_bracket_loc(1)-1);
    f2 = expr(left_bracket_loc(2)+1:right_bracket_loc(2)-1);

    if (contains(f1, 't') || contains(f1, 'I'))
        ft = f1;
        fx = f2;
    else
        ft = f2;
        fx = f1;
    end
    %     ft = strip(ft, 'left', '{');
    %     ft = strip(ft, 'right', '}');
    %     fx = strip(fx, 'left', '{');
    %     fx = strip(fx, 'right', '}');
end