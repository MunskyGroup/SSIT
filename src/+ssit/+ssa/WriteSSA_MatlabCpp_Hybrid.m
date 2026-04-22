function WriteGPUSSA(k,w,S,tprint,fun_name)
Nspec = size(S,1);
Nrxn = size(S,2);
Nt = length(tprint);
funName = char(fun_name);

gpuIn = cell(1,Nspec);
ssaIn = cell(1,Nspec);
ssaOut = cell(1,Nspec*Nt);
idx = 1;
for i = 1:Nspec
    gpuIn{i} = ['x',num2str(i),'_0_GPU'];
    ssaIn{i} = ['x',num2str(i)];
    for j = 1:Nt
        ssaOut{idx} = ['x',num2str(i),'_',num2str(j)];
        idx = idx + 1;
    end
end

fileID = fopen([funName,'.m'],'w');
if fileID < 0
    error('WriteGPUSSA:fileOpenFailed','Unable to create %s.m',funName);
end
cleanupObj = onCleanup(@() fclose(fileID));

fprintf(fileID,'function [X]=%s(x0,N_run,parametersIn,useGPU)\n',funName);
fprintf(fileID,'%% This is an automatically generated MATLAB SSA Program.\n');
fprintf(fileID,'arguments\n');
fprintf(fileID,'   x0\n');
fprintf(fileID,'   N_run\n');
fprintf(fileID,'   parametersIn\n');
fprintf(fileID,'   useGPU=''CPU''\n');
fprintf(fileID,'end\n\n');

fprintf(fileID,'Nspec = %d;\n',Nspec);
fprintf(fileID,'Nt = %d;\n',Nt);
fprintf(fileID,'X = zeros(Nspec,Nt,N_run);\n\n');

fprintf(fileID,'if strcmp(useGPU,''GPU'')\n');
fprintf(fileID,'   g = parallel.gpu.RandStream(''Philox4x32-10'',''Seed'',0);\n');
fprintf(fileID,'   parallel.gpu.RandStream.setGlobalStream(g);\n');
for i = 1:Nspec
    fprintf(fileID,'   x%d_0_GPU = x0(%d)*gpuArray.ones(1,N_run);\n',i,i);
end
fprintf(fileID,'   %s_SSA_GPU = @(%s)%s_SSA(%s,parametersIn);\n',...
    funName,strjoin(gpuIn,','),funName,strjoin(gpuIn,','));
fprintf(fileID,'   [%s] = arrayfun(@%s_SSA_GPU,%s);\n',...
    strjoin(ssaOut,','),funName,strjoin(gpuIn,','));
for i = 1:Nspec
    for j = 1:Nt
        fprintf(fileID,'   X(%d,%d,:) = gather(x%d_%d);\n',i,j,i,j);
    end
end

fprintf(fileID,'elseif strcmp(useGPU,''Parallel'')\n');
fprintf(fileID,'   ssaCpuFun = str2func([mfilename,''_SSA_mex'']);\n');
fprintf(fileID,'   x0Col = x0(:);\n');
fprintf(fileID,'   p = gcp(''nocreate'');\n');
fprintf(fileID,'   if isempty(p)\n');
fprintf(fileID,'      nGroups = 1;\n');
fprintf(fileID,'   else\n');
fprintf(fileID,'      nGroups = p.NumWorkers;\n');
fprintf(fileID,'   end\n');
fprintf(fileID,'   nGroups = max(1,min(nGroups,N_run));\n');
fprintf(fileID,'   chunkSize = ceil(N_run/nGroups);\n');
fprintf(fileID,'   groupStarts = 1:chunkSize:N_run;\n');
fprintf(fileID,'   nGroups = numel(groupStarts);\n');
fprintf(fileID,'   groupResults = cell(1,nGroups);\n');
fprintf(fileID,'   parfor ig = 1:nGroups\n');
fprintf(fileID,'      iStart = groupStarts(ig);\n');
fprintf(fileID,'      iEnd = min(N_run,iStart+chunkSize-1);\n');
fprintf(fileID,'      nLocal = iEnd - iStart + 1;\n');
fprintf(fileID,'      groupResults{ig} = ssaCpuFun(x0Col,nLocal,parametersIn);\n');
fprintf(fileID,'   end\n');
fprintf(fileID,'   for ig = 1:nGroups\n');
fprintf(fileID,'      iStart = groupStarts(ig);\n');
fprintf(fileID,'      iEnd = min(N_run,iStart+chunkSize-1);\n');
fprintf(fileID,'      X(:,:,iStart:iEnd) = groupResults{ig};\n');
fprintf(fileID,'   end\n');

fprintf(fileID,'elseif strcmp(useGPU,''Series'')\n');
fprintf(fileID,'   ssaCpuFun = str2func([mfilename,''_SSA_mex'']);\n');
fprintf(fileID,'   x0Col = x0(:);\n');
fprintf(fileID,'   X = ssaCpuFun(x0Col,N_run,parametersIn);\n');
fprintf(fileID,'else\n');
fprintf(fileID,'   error(''Unknown useGPU mode: %%s'',useGPU);\n');
fprintf(fileID,'end\n');
fprintf(fileID,'end\n\n');

fprintf(fileID,'function [%s] = %s_SSA(%s,parametersIn)\n',...
    strjoin(ssaOut,','),funName,strjoin(ssaIn,','));
fprintf(fileID,'%% First we define the parameters.\n');
for i = 1:length(k)
    fprintf(fileID,'k%d=parametersIn(%d);\n',i,i);
end
fprintf(fileID,'\n');

fprintf(fileID,'%% Initialize the time and state.\n');
fprintf(fileID,'t=%.17g;\n',tprint(1));
for i = 1:Nspec
    fprintf(fileID,'x%dnew=x%d;\n',i,i);
end
fprintf(fileID,'\n');

for it = 1:Nt
    fprintf(fileID,'tstop = %.17g;\n',tprint(it));
    fprintf(fileID,'while t<tstop\n');
    for j = 1:Nspec
        fprintf(fileID,'  x%d=x%dnew;\n',j,j);
    end

    for ir = 1:Nrxn
        fprintf(fileID,'  w%d=%s;\n',ir,w{ir});
    end
    fprintf(fileID,'  w0=0');
    for ir = 1:Nrxn
        fprintf(fileID,'+w%d',ir);
    end
    fprintf(fileID,';\n');
    fprintf(fileID,'  t = t-1/w0*log(rand);\n');
    fprintf(fileID,'  r2w0=rand*w0;\n');

    for ir = 1:Nrxn
        if ir == 1
            fprintf(fileID,'  if r2w0<w1');
        else
            fprintf(fileID,'  elseif r2w0<w1');
        end
        for jr = 2:ir
            fprintf(fileID,'+w%d',jr);
        end
        fprintf(fileID,'\n');
        for js = 1:Nspec
            if S(js,ir) ~= 0
                fprintf(fileID,'    x%dnew=x%d+(%g);\n',js,js,S(js,ir));
            end
        end
    end
    fprintf(fileID,'  end\n');
    fprintf(fileID,'end\n');
    for js = 1:Nspec
        fprintf(fileID,'x%d_%d=x%d;\n',js,it,js);
    end
    fprintf(fileID,'\n');
end
fprintf(fileID,'end\n');

clear cleanupObj;

% Compile CPU SSA kernel used by Parallel/Series branches.
ssit.ssa.WriteCppSSA(k,w,S,tprint,[funName,'_SSA_mex']);