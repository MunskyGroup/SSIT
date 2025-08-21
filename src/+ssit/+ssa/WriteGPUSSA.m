function WriteGPUSSA(k,w,S,tprint,fun_name)
Nspec = size(S,1); % Number of species.
Nrxn = size(S,2);  % Number of reactions.
Nt = length(tprint); % Length of time at which to print results.

fileID = fopen([fun_name,'.m'],'w');
% Name of m-file to be written.

txt = ['function [X]=',fun_name,'(x0,N_run,useGPU)\r\n'];
fprintf(fileID,txt);
fprintf(fileID,'%%This is an automatically generated MATLAB SSA Program.\r\n');
fprintf(fileID,'%%The tools used to generate this file were covered at the.\r\n');
fprintf(fileID,'%%2016 q-bio Summer School at the Colorado State University.\r\n');
% Write Header for m-file.

fprintf(fileID,'arguments\r\n');
fprintf(fileID,'   x0\r\n');
fprintf(fileID,'   N_run\r\n');
fprintf(fileID,'   useGPU=CPU\r\n');
fprintf(fileID,'end\r\n\r\n');

fprintf(fileID,['Nspec = ',num2str(Nspec),'; %% Number of species.\r\n']);
fprintf(fileID,['Nt = ',num2str(Nt),'; %% Number of time points.\r\n']);

% fprintf(fileID,'n_run = N_run/N_split; %% Number of runs per CPU.\r\n');
fprintf(fileID,'X = zeros(Nspec,Nt,N_run); %% Initialize matrix of results. \r\n');

fprintf(fileID,'if strcmp(useGPU,''GPU'')\r\n');
fprintf(fileID,'   g = parallel.gpu.RandStream(''Philox4x32-10'',''Seed'',0); %% Set seed for RNG on GPU.\r\n');
fprintf(fileID,'   parallel.gpu.RandStream.setGlobalStream(g); %% Apply RNG seed to GPU.\r\n');

for i=1:Nspec
    fprintf(fileID,['   x',num2str(i),'_0_GPU = x0(',num2str(i),')*gpuArray.ones(1,N_run); %% Specific Initial Conditions.\r\n']);
end

% fprintf(fileID,'for i = 1:N_split  %% Run loop over different GPU cards.\r\n');
for i=1:Nspec
    for j=1:Nt
        if i==1&&j==1
            txt = '   [x1_1';
        else
            txt = [txt,',x',num2str(i),'_',num2str(j)];
        end
    end
    if i==1
        txt2 = 'x1_0_GPU';
    else
        txt2 = [txt2,',x',num2str(i),'_0_GPU'];
    end
end
txt = [txt,'] = arrayfun(@',fun_name,'_SSA,',txt2,');\r\n'];
fprintf(fileID,txt);

for i=1:Nspec
    for j=1:Nt
        txt = ['   X(',num2str(i),',',num2str(j),',:) =gather(x',num2str(i),'_',num2str(j),');\r\n'];
        fprintf(fileID,txt);
    end
end

fprintf(fileID,'elseif strcmp(useGPU,''Parallel'')\r\n');
for i=1:Nspec
    fprintf(fileID,['   x',num2str(i),'_0 = x0(',num2str(i),'); %% Specific Initial Conditions.\r\n']);
end

fprintf(fileID,'  parfor i = 1:N_run\r\n');

for i=1:Nspec
    for j=1:Nt
        if i==1&&j==1
            txt = '    [x1_1';
        else
            txt = [txt,',x',num2str(i),'_',num2str(j)];
        end
    end
    if i==1
        txt2 = 'x1_0';
    else
        txt2 = [txt2,',x',num2str(i),'_0'];
    end
end
txt3 = [txt,'] = ',fun_name,'_SSA(',txt2,');\r\n'];
fprintf(fileID,txt3);
txt4 = ['    X(:,:,i) = reshape(',txt(3:end),'],[Nt,Nspec])'';\r\n'];
fprintf(fileID,txt4);
fprintf(fileID,'  end\r\n');

fprintf(fileID,'elseif strcmp(useGPU,''Series'')\r\n');
fprintf(fileID,'  for i = 1:N_run\r\n');
for i=1:Nspec
    fprintf(fileID,['   x',num2str(i),'_0 = x0(',num2str(i),'); %% Specific Initial Conditions.\r\n']);
end
fprintf(fileID,txt3);
fprintf(fileID,txt4);
fprintf(fileID,'  end\r\n');

fprintf(fileID,'end\r\n');

% fprintf(fileID,'end\r\n');
fprintf(fileID,'\r\n\r\n');

%%

for i=1:Nspec
    for j=1:Nt
        if i==1&&j==1
            txt = 'function [x1_1';
        else
            txt = [txt,',x',num2str(i),'_',num2str(j)];
        end
    end
    if i==1
        txt2 = 'x1';
    else
        txt2 = [txt2,',x',num2str(i)];
    end
end
txt = [txt,'] = ',fun_name,'_SSA(',txt2,')\r\n'];
fprintf(fileID,txt);


fprintf(fileID,['%%First we define the parameters.\r\n']);
for i=1:length(k)
    txt = ['k',num2str(i),'=',num2str(k(i)),';\r\n'];
    fprintf(fileID,txt);
end
fprintf(fileID,'\r\n');

fprintf(fileID,['%%Initialize the time.\r\n']);
fprintf(fileID,['t=0;\r\n']);

fprintf(fileID,['%%Start the SSA.\r\n']);

for it = 1:Nt
    fprintf(fileID,['tstop = ',num2str(tprint(it)),';   %%Next time to print results.\r\n']);
    fprintf(fileID,['while t<tstop   %%Next time to print results.\r\n']);
    
    txt2 = '  w0=0';
    for i=1:Nrxn
        txt = strcat({'  w'},{num2str(i)},{'='},w{i},{';\r\n'});
        fprintf(fileID,txt{1});
        txt2 = [txt2,'+w',num2str(i)];
    end
    txt2 = [txt2,';\r\n'];
    fprintf(fileID,txt2);
    fprintf(fileID,'  t = t-1/w0*log(rand);\r\n');
    fprintf(fileID,'  if t<=tstop\r\n');
    fprintf(fileID,'    r2w0=rand*w0;\r\n');
    
    for i=1:Nrxn
        if i==1
            txt = '    if r2w0<w1';
        else
            txt = '    elseif r2w0<w1';
        end
        for j=2:i
            txt=[txt,'+w',num2str(j)];
        end
        txt = [txt,'\r\n'];
        fprintf(fileID,txt);
        
        for j=1:Nspec
            if S(j,i)~=0
                txt = ['      x',num2str(j),'=x',num2str(j),'+(',num2str(S(j,i)),');\r\n'];
                fprintf(fileID,txt);
            end
        end
    end
    fprintf(fileID,'    end\r\n');
    fprintf(fileID,'  end\r\n');
    fprintf(fileID,'end\r\n');
    
    for i=1:Nspec
        txt = ['x',num2str(i),'_',num2str(it),'=x',num2str(i),';\r\n'];
        fprintf(fileID,txt);
    end
    fprintf(fileID,'\r\n');
end
