function summarizeModelReactions(Species,Reactions,Inputs,ax)
%% SSIT.summarizeModel - Prints a summary of an SSIT model:
%% Species, Reactions (with Stoichiometric updates),
%% Model Parameters
%
% Input:  SSIT model
%
% Output:  Summary text to screen
%
% Example:  Model.summarizeModel
arguments
    Species
    Reactions
    Inputs
    ax
end

nSpecies = length(Species);
nRxn = size(Reactions,1);
S = zeros(nSpecies,nRxn);                      % Assigns a zeros ReactionsTabOutputs.stoichMatrixiometry matrix for the three species to be filled in with the following for-loops
for iRxn=1:nRxn                                       % For-loop over all rows in the Reactions Table
    %% Build ReactionsTabOutputs.stoichMatrixiometry matrix
    Rcts = Reactions{iRxn,2};
    Prods = Reactions{iRxn,3};
    for iSp = 1:nSpecies
        % Reactants
        J = strfind(Rcts,Species{iSp,1});
        J1 = strfind(Rcts,'(');
        J2 = strfind(Rcts,')');
        for j=1:length(J)
            S(iSp,iRxn)=S(iSp,iRxn)-str2num(Rcts(J1(j)+1:J2(j)-1));      % Assign reaction effect on stoichiometry
        end

        % Products
        J = strfind(Prods,Species{iSp,1});
        J1 = strfind(Prods,'(');
        J2 = strfind(Prods,')');
        for j=1:length(J)                                % For-loop over the number of products found by J
            S(iSp,iRxn)=S(iSp,iRxn)+str2num(Prods(J1(j)+1:J2(j)-1)); % Assign product effect on stoichiometry
        end
    end
end

% Assume `ax` is a valid axes handle
cla(ax);  % clear previous content
axis(ax, 'off');  % hide axes
hold(ax, 'on');

% Start compiling text lines
textLines = {};
textLines{end+1} = '\textbf{Reactions:}';

for iR = 1:nRxn
    % Reactants
    s=S(:,iR);
    jReactant = find(s<0);
    jProd = find(s>0);
    if isempty(jReactant)
        reactTxt = 'NULL';
    else
        if s(jReactant(1))==-1
            reactTxt = Species{jReactant(1)};
        else
            reactTxt = [num2str(-s(jReactant(1))),Species{jReactant(1)}];
        end
        for i = 2:length(jReactant)
            if s(jReactant(i))==-1
                reactTxt = [reactTxt,' + ',Species{jReactant(i)}];
            else
                reactTxt = [reactTxt,' + ',num2str(-s(jReactant(i))),Species{jReactant(i)}];
            end
        end
    end
    if isempty(jProd)
        prodTxt = 'NULL';
    else
        if s((1))==1
            prodTxt = Species{jProd(1)};
        else
            prodTxt = [num2str(s(jProd(1))),Species{jProd(1)}];
        end
        for i = 2:length(jProd)
            if s(jProd(i))==1
                prodTxt = [prodTxt,' + ',Species{jProd(i)}];
            else
                prodTxt = [prodTxt,' + ',num2str(s(jProd(i))),Species{jProd(i)}];
            end
        end
    end

    syms(symvar(Reactions{iR,4}));
    propstext = latex(simplify(eval(Reactions{iR,4})));
    
    textLines{end+1} = sprintf('s_{%d}: \\quad %s \\rightarrow %s \\quad | \\quad w_{%d}(x): \\quad %s', ...
    iR, reactTxt, prodTxt, iR, propstext);
    

end

% Input Signals
if ~isempty(Inputs)
    textLines{end+1} = '-----------';
    textLines{end+1} = '\textbf{Input Signals:}';
    nI = size(Inputs, 1);
    for i = 1:nI
        syms(symvar(Inputs{i,1}));
        syms(symvar(Inputs{i,2}));
        InputsTxt1 = latex(simplify(eval(Inputs{i,1})));
        InputsTxt2 = latex(simplify(eval(Inputs{i,2})));
        textLines{end+1} = sprintf('%s(t) = %s', InputsTxt1, InputsTxt2);
    end
end

% Parameters
% textLines{end+1} = ' ';
% textLines{end+1} = '\textbf{Model Parameters:}';
% paramNames = obj.parameters(:,1);
% for i = 1:size(paramNames,1)
%     pname = paramNames{i};
%     pval = obj.parameters{i,2};
%     textLines{end+1} = sprintf('%s = %.4g', pname, pval);
% end

% Display all lines as LaTeX text on the axes
nLines = length(textLines);
y = 1; dy = 2;  % line spacing
for i = 1:nLines
    text(ax, 0, y, ['$$', textLines{i}, '$$'], ...
        'Interpreter', 'latex', 'FontSize', 14, ...
        'VerticalAlignment', 'top', 'Units', 'normalized');
    y = y - dy * 0.05;
end


end

