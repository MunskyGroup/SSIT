function [smpl,accept,value,bestfound,bestObjToNow] = metropolisHastingsSample(start,nsamples,varargin)
% MHSAMPLE Generate Markov chain using Metropolis-Hasting algorithm 
%   SMPL = MHSAMPLE(START,NSAMPLES,'pdf',PDF,'proppdf',PROPPDF,'proprnd',PROPRND)
%   draws NSAMPLES random samples from a target stationary distribution PDF
%   using the Metropolis-Hasting algorithm. START is a row vector
%   containing the start value of the Markov Chain. NSAMPLES is an integer
%   specifying the number of samples to be generated. PDF, PROPPDF, and
%   PROPRND are function handles created using @. PROPPDF defines the
%   proposal distribution density and PROPRND defines the random number
%   generator for the proposal distribution. PDF and PROPRND take one
%   argument as an input and this argument has the same type and size as
%   START. PROPPDF takes two arguments as inputs and both arguments have
%   the same type and size as START. SMPL is a column vector or matrix
%   containing the samples. If log density function is preferred, 'pdf' and
%   'proppdf' can be replaced with 'logpdf' and 'logproppdf'. The density
%   functions used in Metropolis-Hasting algorithm are not necessarily 
%   normalized.  If proppdf or logproppdf q(x,y) satisfies q(x,y) = q(y,x),
%   i.e. the proposal distribution is symmetric, MHSAMPLE implements
%   Random Walk Metropolis-Hasting sampling.  If proppdf or logproppdf
%   q(x,y) satisfies q(x,y) = q(x), i.e. the proposal distribution is
%   independent of current values, MHSAMPLE implements Independent
%   Metropolis-Hasting sampling.
%
%   The proposal distribution q(x,y) gives the probability density for
%   choosing x as the next point when y is the current point.  It is
%   sometimes written as q(x|y) in the literature.
%  
%   SMPL = MHSAMPLE(...,'symmetric',SYM) draws NSAMPLES random samples from
%   a target stationary distribution PDF using the Metropolis-Hasting
%   algorithm. SYM is a logical value, which indicates whether the proposal
%   distribution is symmetric. The default value is false, which
%   corresponds to the asymmetric proposal distribution. If the SYM is
%   true, i.e. the proposal distribution is symmetric, PROPPDF and
%   LOGPROPDF are optional. 
% 
%   SMPL = MHSAMPLE(...,'burnin',K) generate a Markov chain with values
%   between the starting point and the K-th point omitted, but keep points
%   after that. K is a non-negative integer. The default value of K is 0. 
%
%   SMPL = MHSAMPLE(...,'thin',M) generate a markov chain with M-1 out
%   of M values omitted in the generated sequence. M is a positive integer.
%   The default value is 1.
%
%   SMPL = MHSAMPLE(...,'nchain',N) generate N Markov chains using the
%   Metropolis-Hasting algorithm. N is a positive integer. The default
%   value of N is 1. SMPL is a matrix containing the samples. The last
%   dimension contains the indices for individual chains.
%
%   [SMPL,ACCEPT] = MHSAMPLE(...) also returns ACCEPT as the acceptance
%   rate of the proposed distribution. ACCEPT is a scalar if a single chain
%   is generated and is a vector if multiple chains are generated.
%
%  Example:  
%    Estimate the second order moment of a Gamma distribution using
%    Independent Metropolis-Hasting sampling: 
%        alpha = 2.43;
%        beta = 1;
%        pdf = @(x) gampdf(x,alpha,beta);     %target distribution
%        proppdf = @(x,y) gampdf ...          % proposal pdf
%                   (x,floor(alpha),floor(alpha)/alpha);
%        proprnd = @(x) sum ...               % proposal random sampler
%                   (exprnd(floor(alpha)/alpha,floor(alpha),1));
%        nsamples = 5000;
%        smpl = mhsample(1,nsamples,'pdf',pdf,'proprnd', ...
%                   proprnd,'proppdf',proppdf);
%        xxhat = cumsum(smpl.^2)./(1:nsamples)';
%        plot(1:nsamples,xxhat)
%
%    Generate random samples from N(0,1) using Random Walk
%    Metropolis-Hasting sampling: 
%        delta = .5;
%        pdf = @(x) normpdf(x);                     %target distribution
%        proppdf = @(x,y) unifpdf(y-x,-delta,delta);% proposal pdf
%        proprnd = @(x) x + rand*2*delta - delta;   % proposal random sampler
%        nsamples = 15000;
%        x = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd,'symmetric',1);
%        histfit(x,50)
%
%  See also: SLICESAMPLE, RAND.
 
% Copyright 2005-2011 The MathWorks, Inc.

% if nargin<2  % start and nsamples are required.
%     error(message('stats:mhsample:LessArg'));
% end
 

% parse the information in the name/value pairs 
pnames = {'pdf' ,'logpdf', 'proppdf','logproppdf','proprnd', ...
    'burnin','thin','symmetric','nchains','progress','saveFileName'};
dflts =  {[] [] [],[],[] 0,1,false,1,false,'tmpMHResults'};
[pdf,logpdf,proppdf,logproppdf, proprnd,burnin,thin,sym,nchain,progress,saveFileName] = ...
       internal.stats.parseArgs(pnames, dflts, varargin{:});

if exist(saveFileName,'file')||exist([saveFileName,'.mat'],'file')
    load(saveFileName,'value','smpl','bestfound'); 
    k = find(value~=0,1,'last');
    start = smpl(k,:);
    burnin=0;
    iStart = k*thin+1;
    disp(['Found an interrupted MH run -- resuming from broken chain at step ',num2str(iStart)])
else
    iStart = 1-burnin;
end

% pdf or logpdf for target distribution has to be provided.
if isempty(pdf)&& isempty(logpdf)
    error(message('stats:mhsample:BadTarget'));
end
 
% pdf or logpdf for proposal distribution has to be provided.
if ~sym
    if (isempty(proppdf)&& isempty(logproppdf))
        error(message('stats:mhsample:NoProppdf'));
    end
end
if isempty(proprnd)
    error(message('stats:mhsample:NoProprnd'));
end 
 
% error checks for the functional handles
if ~isempty(pdf)
   checkFunErrs('pdf',pdf,start);
end
if ~isempty(logpdf)
   checkFunErrs('logpdf',logpdf,start);
end
if ~isempty(proppdf)
   checkFunErrs('proppdf',proppdf,start);
end
if ~isempty(logproppdf)
   checkFunErrs('logproppdf',logproppdf,start);
end
if ~isempty(proprnd)
   checkFunErrs('proprnd',proprnd,start);
end
 
%error checks for burnin and thin
if (burnin<0) || burnin~=round(burnin)
    error(message('stats:mhsample:BadBurnin'));
end
if (thin<=0)|| thin~=round(thin)
    error(message('stats:mhsample:BadThin'));
end
 
% log density is preferred for numerical stability
if ~isempty(pdf) && isempty(logpdf)
    logpdf = @(x) mylog(pdf(x));
end
if ~isempty(proppdf) && isempty(logproppdf)
    logproppdf = @(x,y) mylog(proppdf(x,y));
end     
if ~sym 
    if any (logpdf(start)==-Inf | logproppdf(start,start) == -Inf)
        error(message('stats:mhsample:BadStart'))
    end
else 
    if any (logpdf(start)==-Inf) 
            error(message('stats:mhsample:BadStartSym'))
    end
end
outclass = superiorfloat(start); % single or double

% Put the replicates dimension second.
distnDims = size(start,2);
if ~exist('smpl','var')
    smpl = zeros([nsamples,nchain,distnDims],outclass);
end

if ~exist('value','var')
    value = zeros([nsamples,nchain,1],outclass);
end
 
x0 = start;  %x0  is the place holder for the current value
accept =zeros(nchain,1,outclass);    
% Metropolis-Hasting Algorithm.
U = log(rand(nchain,nsamples*thin+burnin));

logpdf_x0 = logpdf(x0);   %Brian -- Do not recompute OBJ.

if ~exist('bestfound','var')
    bestfound = x0;
    bestObjToNow = logpdf_x0;
else
    bestObjToNow = logpdf(bestfound);
end  

if progress
    D = waitbar(0,'Metropolis Hastings Progress');
end

for i = iStart:nsamples*thin
    y = proprnd(x0); % sample from proposal dist'n
    if ~sym
        q1 = logproppdf(x0,y);
        q2 = logproppdf(y,x0);
        % this is a generic formula.
        logpdf_y = logpdf(y);  %Brian -- Do not recompute OBJ.
        rho = (q1+logpdf_y)-(q2+logpdf_x0);%Brian -- Do not recompute OBJ.
    else
        %save the evaluation time for symmetric proposal dist'n
        logpdf_y = logpdf(y);         %Brian -- Do not recompute OBJ.
        if isfinite(logpdf_y)
            rho = logpdf_y-logpdf_x0; %Brian -- Do not recompute OBJ.
        else
            rho = -inf; %Brian -- reject pars that lead to NaN objective funs.
        end
    end
    % Accept or reject the proposal.
    Ui = U(:,i+burnin);
    acc = Ui<=min(rho,0);
    if acc&&~isfinite(logpdf_y)
        error('not finite but accepted!!')
    end
    
    if acc&&logpdf_y>bestObjToNow 
        bestObjToNow = logpdf_y;
        bestfound = y;
    end

    x0(acc,:) = y(acc,:); % preserves x's shape.
    accept = accept+(acc);
    
    if isfinite(logpdf_y) %Brian -- fix bug that puts NaN values in chain.
        logpdf_x0 = logpdf_x0.*(1-acc) + logpdf_y.*acc;  %Brian -- Do not recompute OBJ.
    end
    if i>0 && mod(i,thin)==0 % burnin and thin
        smpl(i/thin,:,:) = x0;
        value(i/thin,:,1) = logpdf_x0;
    end
    
    if progress&&mod(i,20)==0
        waitbar((i+burnin)/(burnin+nsamples*thin),D,...
            ['Running MetHast, Acceptance = ',num2str(accept/(i+burnin))]);
        figure(123);clf;plot(value(value~=0));
    elseif mod(i,ceil(nsamples*thin/20))==0
        disp(['Completed ',num2str(i),' out of ',num2str(nsamples*thin),' samples. Acceptance = ',num2str(accept/(i+burnin))])
    end

    if mod(i,100)==0
        save(saveFileName,'value','smpl','bestfound')
        % if progress
        %     disp(['n=',num2str(i),'; acc=',num2str(accept/(i+burnin)),'. ',saveFileName])
        % end
    end
end
if progress
    close(D);
end
 
% Accept rate can be used to optimize the choice of scale parameters in
% random walk MH sampler. See for example Roberts, Gelman and Gilks (1997).
accept = accept/(nsamples*thin+burnin); 
 
% Move the replicates dimension to the end to make samples easier to
% manipulate.
smpl = permute(smpl,[1 3 2]);
 
%-------------------------------------------------
function  y = mylog(x)
% my log function is to define to avoid the warnings. 
y = -Inf(size(x));
y(x>0) = log(x(x>0));
 
%----------------------------------------------------
function checkFunErrs(type,fun,param)
%CHECKFUNERRS Check for errors in evaluation of user-supplied function
if isempty(fun), return; end
 
try
    switch type
    case 'logproppdf'
        out=fun(param,param);
    case 'proppdf'
        out=fun(param,param);
    case 'logpdf'
        out=fun(param);
    case 'pdf'
        out=fun(param);
    case 'proprnd'
        out=fun(param);
    end
catch ME
    m = message('stats:mhsample:FunEvalError',type,func2str(fun));
    throwAsCaller(addCause(MException(m.Identifier,'%s',getString(m)),ME));
end
switch type
    case 'logproppdf'
        if any(isnan(out))||any(isinf(out) & out>0 )
            error(message('stats:mhsample:NonfiniteLogproppdf'));
        end
    case 'proppdf'
        if any(~isfinite(out))
            error(message('stats:mhsample:NonfiniteProppdf'));
        end     
        if any(out<0)
            error(message('stats:mhsample:NegativeProppdf'));
        end 
    case 'logpdf'
        if any(isnan(out))||any(isinf(out) & out>0 )
            error(message('stats:mhsample:NonfiniteLogpdf'));
        end
    case 'pdf'
        if any(~isfinite(out))
            error(message('stats:mhsample:NonfinitePdf'));
        end    
        if any(out<0)
            error(message('stats:mhsample:NegativePdf'));
        end 
    case 'proprnd'
        if any(~isfinite(out))
            error(message('stats:mhsample:NonfiniteProprnd'));
        end     
end


