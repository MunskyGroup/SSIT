function handleOut = matlabFunctionSSIT(...
    symbolicExpression, VarNames, fn, isSparse, logVarsReps)
% This function writes m files for SSIT propensity functions.
arguments
    symbolicExpression sym
    VarNames cell
    fn (1, 1) string {mustBeNonzeroLengthText}
    isSparse (1, 1) logical = false;
    logVarsReps cell = {};
end

[path, name, ext] = fileparts(fn);
% Ensure that the filename ends in .m if it doesn't already.
if ~endsWith(ext, ".m")
    ext = ext + ".m";
end
fn = path + filesep + name + ext;
[~, name, ~] = fileparts(fn);

fid = fopen(fn, 'w');

txt = "function result = " + name + "(in1,in2,in3,in4)\r\n";
fprintf(fid, txt);

fprintf(fid, ...
    "%%This was automatically generated on " + ...
    string(datetime("now")) + ".\r\n");
fprintf(fid,'%%\r\n');

fprintf(fid,'arguments\r\n');
fprintf(fid,'  in1\r\n');
fprintf(fid,'  in2\r\n');
fprintf(fid,'  in3 = [];\r\n');
fprintf(fid,'  in4 = [];\r\n');

fprintf(fid,'end\r\n');
fprintf(fid,'\r\n');

% TODO - Diagnose bug that occurs if treating symbolic expression and
% logVarsReps as strings. The following error message is obtained when
% % trying to create the ToggleSwitch model:
% Operands to the short-circuit AND (&&) and OR (||) operators must be convertible to logical scalars. Use
% the ANY or ALL functions to reduce array operands to logical scalars. For elementwise operations, use the
% AND (&) and OR (|) operators instead.
charSymb = char(symbolicExpression);
for j = 1:length(VarNames)
    vName = VarNames{j};
    for i = 1:length(vName)
        charVn = char(vName(i));
        if contains(charSymb,charVn)||...
                (~isempty(logVarsReps)&&max(contains(logVarsReps(:,1),charVn)))
            if isSparse
                sparsityString = "sparse";
            else
                sparsityString = "";
            end
            fprintf(fid, string(charVn) + "=" + sparsityString + ...
                "(in" + j + "(" + i + ",:));\r\n");            
        end
    end
    fprintf(fid,'\r\n');
end


% TODO - Need to figure better way to create sparse matrix outputs.
% if isSparse
%     [I,J] = find(symbolicExpression~=0);
%     % txt = ['result = sparse(',num2str(I),'],[',num2str(J),']',char((symbolicExpression(symbolicExpression~=0))),';\r\n'];  
%     txt = ['result = sparse([',num2str(I'),'],[',num2str(J'),'],',char(symbolicExpression(sub2ind(size(symbolicExpression),I,J))'),');\r\n'];  
% else

    txt = "result = " + charSymb + ";\r\n";
% end

if ~isempty(logVarsReps)
    % Replace first string in function name
    for i = 1:size(logVarsReps,1)
        txt = strrep(txt, [logVarsReps{i,2}], ...
            "(" + logVarsReps{i,1} + ")");
    end

end

txt = stripOmitnanClause(txt);

% Ensure that all multiplication, division, and exponentiation is
% element-wise

multDivideSym = {"*", "/", "^"};
for j = 1:length(multDivideSym)
    txt = strrep(txt, multDivideSym{j}, "." + multDivideSym{j});
end
% The above might have introduced some double-dots that need to be removed.
txt = strrep(txt, "..", "."); 

fprintf(fid,txt);
fclose(fid);

handleOut = str2func(name);

end

function txtOut = stripOmitnanClause(txtIn)
% Remove optional max/min-style omitnan clause:
%   ..., 'omitnan', <any expression>)
% where <any expression> may itself contain nested parentheses.

txtOut = txtIn;
searchStart = 1;

while true
    if searchStart > numel(txtOut)
        break;
    end

    % Match omitnan with optional quote wrappers, e.g.:
    % 'omitnan', "omitnan", or "'omitnan'".
    [matchStarts, matchEnds] = regexpi(txtOut(searchStart:end), ...
        '(?<![A-Za-z0-9_])["'']*omitnan["'']*(?![A-Za-z0-9_])');
    if isempty(matchStarts)
        break;
    end
    omitLoc = searchStart + matchStarts(1) - 1;
    omitEnd = searchStart + matchEnds(1) - 1;

    depthAtOmit = parenDepthAt(txtOut,omitLoc);

    % Find the comma that starts the omitnan argument chunk.
    startRemove = omitLoc;
    while startRemove > 1
        c = txtOut(startRemove-1);
        if c == ','
            startRemove = startRemove - 1;
            break;
        elseif c == '(' || c == ';' || c == char(10) || c == char(13)
            % No valid argument prefix found; skip this occurrence.
            startRemove = [];
            break;
        else
            startRemove = startRemove - 1;
        end
    end

    if isempty(startRemove) || startRemove <= 0
        searchStart = omitEnd + 1;
        continue;
    end

    % Scan forward to the matching closing parenthesis of the current call.
    [closingParenLoc, foundParen] = findMatchingCallClose(txtOut,omitLoc,depthAtOmit);
    if ~foundParen
        searchStart = omitEnd + 1;
        continue;
    end

    % Remove from the comma before omitnan through the char before ')'.
    txtOut = [txtOut(1:startRemove-1), txtOut(closingParenLoc:end)];
    searchStart = max(startRemove,1);
end

end

function depth = parenDepthAt(txt,idx)
% Parenthesis depth at position idx (ignoring quoted strings).
depth = 0;
inSingle = false;
inDouble = false;

limit = idx - 1;
i = 1;
while i <= limit
    c = txt(i);

    if inSingle
        if c == ''''
            if i < numel(txt) && txt(i+1) == ''''
                i = i + 1;
            else
                inSingle = false;
            end
        end
        i = i + 1;
        continue;
    end

    if inDouble
        if c == '"'
            inDouble = false;
        end
        i = i + 1;
        continue;
    end

    if c == ''''
        inSingle = true;
    elseif c == '"'
        inDouble = true;
    elseif c == '('
        depth = depth + 1;
    elseif c == ')'
        depth = max(depth - 1,0);
    end

    i = i + 1;
end

end

function [closeLoc, found] = findMatchingCallClose(txt,startIdx,baseDepth)
% Find ')' that closes the function call containing startIdx.
closeLoc = numel(txt);
found = false;

depth = baseDepth;
inSingle = false;
inDouble = false;
i = startIdx;

while i <= numel(txt)
    c = txt(i);

    if inSingle
        if c == ''''
            if i < numel(txt) && txt(i+1) == ''''
                i = i + 1;
            else
                inSingle = false;
            end
        end
        i = i + 1;
        continue;
    end

    if inDouble
        if c == '"'
            inDouble = false;
        end
        i = i + 1;
        continue;
    end

    if c == ''''
        inSingle = true;
    elseif c == '"'
        inDouble = true;
    elseif c == '('
        depth = depth + 1;
    elseif c == ')'
        if depth == baseDepth
            closeLoc = i;
            found = true;
            return;
        end
        depth = max(depth - 1,0);
    end

    i = i + 1;
end

end