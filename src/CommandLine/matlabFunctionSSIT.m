function handleOut = matlabFunctionSSIT(symbolicExpression,VarNames,fn,isSparse,logVarsReps)
% This function writes m files for SSIT propensity functions.
arguments
    symbolicExpression
    VarNames
    fn
    isSparse = false;
    logVarsReps = {};
end

if ~strcmp(fn(end-1:end),'.m')
    fn(end+1:end+2) = '.m';
end

jFilName = max([strfind(fn,'/'),strfind(fn,'\')]);

fid = fopen(fn, 'w');

txt = ['function result = ',fn(jFilName+1:end-2),'(in1,in2,in3,in4)\r\n'];
fprintf(fid, txt);

fprintf(fid,['%%This is an automatically generated on ',datestr(now),'.\r\n']);
fprintf(fid,'%%\r\n');

fprintf(fid,'arguments\r\n');
fprintf(fid,'  in1\r\n');
fprintf(fid,'  in2\r\n');
fprintf(fid,'  in3 = [];\r\n');
fprintf(fid,'  in4 = [];\r\n');

fprintf(fid,'end\r\n');
fprintf(fid,'\r\n');

charSymb = char(symbolicExpression);
for j = 1:length(VarNames)
    vName = VarNames{j};
    for i = 1:length(vName)
        cVn = char(vName(i));
        if contains(charSymb,cVn)||...
                (~isempty(logVarsReps)&&max(contains(logVarsReps(:,1),cVn)))
            if isSparse
                fprintf(fid,[cVn,'=sparse(in',num2str(j),'(',num2str(i),',:));\r\n']);
            else
                fprintf(fid,[cVn,'=in',num2str(j),'(',num2str(i),',:);\r\n']);
            end
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
    txt = ['result = ',char(symbolicExpression),';\r\n'];
% end

if ~isempty(logVarsReps)
    % Replace first string in function name
    for i = 1:size(logVarsReps,1)
        txt = strrep(txt,[logVarsReps{i,2}],['(',logVarsReps{i,1},')']);
    end

end

txt = strrep(txt,', ''omitnan'', false','');
txt = strrep(txt,', \"''omitnan''\", false','');

multDivideSym = {'*','/','^'};
for j = 1:3
    txt = strrep(txt,multDivideSym{j},['.',multDivideSym{j}]);
end
txt = strrep(txt,'..','.');


fprintf(fid,txt);
fclose(fid);

handleOut = str2func(fn(jFilName+1:end-2));



