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

txt = strrep(txt,', ''omitnan'', false','');
txt = strrep(txt,', \"''omitnan''\", false','');

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