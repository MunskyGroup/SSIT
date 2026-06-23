%  [w, err, hump] = mexpv( t, A, v, tol, m )
%  MEXPV computes an approximation of w = exp(t*A)*v using Krylov
%  subspace projection techniques. This is a customised version for
%  Markov Chains. This means that a check is done within this code to
%  ensure that the resulting vector w is a probability vector, i.e.,
%  w must have all its components in [0,1], with sum equal to 1.
%  This check is done at some expense and the user may try EXPV
%  which is cheaper since it ignores probability constraints.
%
%  IMPORTANT: The check assumes that the transition rate matrix Q
%             satisfies Qe = 0, where e = (1,...,1)'. Don't use MEXPV
%             if this condition does not hold. Use EXPV instead.
%             MEXPV/EXPV require A = Q', i.e., the TRANSPOSE of Q.
%             Failure to remember this leads to wrong results.
%
%
%  MEXPV does not compute the matrix exponential in isolation but
%  instead, it computes directly the action of the exponential operator
%  on the operand vector. This way of doing so allows for addressing
%  large sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = mexpv( t, A, v )
%  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
%
%  [w, err] = mexpv( t, A, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = mexpv( t, A, v, tol )
%  overrides default tolerance.
%
%  [w, err] = mexpv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace.
%
%  [w, err, hump] = expv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace,
%  and renders an approximation of the `hump'.
%
%  The hump is defined as:
%          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
%  It is used as a measure of the conditioning of the matrix exponential
%  problem. The matrix exponential is well-conditioned if hump = 1,
%  whereas it is poorly-conditioned if hump >> 1.  However the solution
%  can still be relatively fairly accurate even when the hump is large
%  (the hump is an upper bound), especially when the hump and
%  ||w(t)||/||v|| are of the same order of magnitude (further details in
%  reference below). Markov chains are usually well-conditioned problems.
%
%  Example:
%  --------
%    % generate a transition rate matrix
%    n = 100;
%    A = rand(n);
%    for j = 1:n
%	 sumj = 0;
%        for i = 1:n
%            if rand < 0.5, A(i,j) = 0; end;
%            sumj = sumj + A(i,j);
%        end;
%	 A(j,j) = A(j,j)-sumj;
%    end;
%    v = eye(n,1);
%    A = sparse(A); % invaluable for a large and sparse matrix.
%
%    tic
%    [w,err] = expv(1,A,v);
%    toc
%
%    disp('w(1:10) ='); disp(w(1:10));
%    disp('err =');     disp(err);
%
%    tic
%    w_matlab = expm(full(A))*v;
%    toc
%
%    disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
%    gap = norm(w-w_matlab)/norm(w_matlab);
%    disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);
%
%  In the above example, n could have been set to a larger value,
%  but the computation of w_matlab will be too long (feel free to
%  discard this computation).
%
%  See also EXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials.
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

function [wout, err, hump, Time_array_out, P_array, P_lost, tryagain, te, ye] = ...
    mexpv_modified_2( t, Ain, vin, tol, m, N_prt, Time_array,fspTol,SINKS,...
    tNow,fspErrorCondition,resetSparsity,fixedEvents,maxTimeStep,clipStates)

[n] = size(Ain,1);
if nargin == 3
    tol = 1.0e-7;
    m = min(n,30);
    %% Changes by Brian Munsky
    N_prt=1;
    Time_array=[0 t];
    fspTol=1;
    SINKS = [];
    %%
end
if nargin == 4
    m = min(n,30);
    %% Changes by Brian Munsky
    N_prt=1;
    Time_array=[0 t];
    fspTol=1;
    SINKS = [];
    %%
end
%% Changes by Brian Munsky
if nargin == 5
    N_prt=1;
    Time_array=[0 t];
    fspTol=1;
    SINKS = [];
end
if nargin ==6
    Time_array=linspace(0,t,N_prt);
    fspTol=1;
    SINKS = [];
end
if nargin>=7
    N_prt = length(Time_array);
end
if nargin ==7
    fspTol=1;
    SINKS = [];
end
if nargin ==8
    tol=min(tol,fspTol/10);
    SINKS = [];
end
if nargin ==9
    tol=min(tol,fspTol/10);
end
if nargin<12
    resetSparsity = 0;
end
if nargin<13
    fixedEvents = {};
end
%%
maxTimeStep = 10;
% clipStates = false;

% R = 5e-2;
% Jclip = 1:length(vin);
wout = vin;
% if isempty(SINKS)
%     iCore = 1:length(vin);
% else
%     iCore = 1:(SINKS(1)-1);
% end

%%
anorm = norm(Ain,'inf');
mxrej = 10;  btol  = 1.0e-7;
gamma = 0.9; delta = 1.2;
mb    = m; t_out   = abs(t);
s_error = 0;
rndoff= anorm*eps;

k1 = 2; xm = 1/m; normv = norm(vin); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;

if t==0
    sgn = 1;
else
    sgn = sign(t); 
end

istep = 0;
n_prev_clip = 0;

w = wout; ye = w'; te = tNow;
hump = normv;

%% Changes by Brian Munsky
Time_array_out=[];
i_prt=1;
P_array = zeros(length(Time_array),length(vin));
%%
if tNow==Time_array(1)
    P_array(1,:)=wout';
    Time_array_out=tNow;
    i_prt=2;
end
%%

if ~isempty(fixedEvents)
    indNextFixedTime = find(fixedEvents.times>tNow,1,"first");
    nextFixedTime = fixedEvents.times(indNextFixedTime);
else
    nextFixedTime = inf;
end

% if ~clipStates
Acsr = build_csr(Ain);
% A = Ain;
[nr,nc] = size(Ain);
assert(nr==nc,'A must be square');
Jclip = 1:length(vin);
% end

while tNow < t_out && i_prt<=length(Time_array)
    istep = istep + 1;

    % if clipStates
    %     Asub = Ain(iCore, iCore);
    %     d = full(diag(Asub));
    %     inSet = false(length(iCore), 1);
    %     iSeed = find(wout(iCore) > 1e-10);
    %     inSet(iSeed) = true;
    %     frontier = iSeed(:);
    %     while ~isempty(frontier)
    %         newSet = false(length(iCore), 1);
    %         for fi = frontier'
    %             if d(fi) ~= 0
    %                 [jj, ~, vv] = find(Asub(:, fi));
    %                 nbrs = jj(wout(iCore(fi)) * abs(vv / d(fi)) > R & ~inSet(jj));
    %                 newSet(nbrs) = true;
    %             end
    %         end
    %         newSet(inSet) = false;
    %         frontier = find(newSet);
    %         inSet(frontier) = true;
    %     end
    %     Jcand = iCore(inSet);
    %     if ~isempty(SINKS)
    %         Jcand = unique([Jcand(:); SINKS(:)])';
    %     else
    %         Jcand = unique(Jcand(:))';
    %     end
    % 
    %     if length(Jcand) > 1000
    %         Jclip = Jcand;
    %         A = Ain(Jclip,Jclip);
    %         A = A - diag(sum(A));
    %         w = wout(Jclip);
    %     else
    %         Jclip = 1:length(vin);
    %         A = Ain;
    %         w = wout;
    %     end
    %     [nr,nc] = size(A);
    %     assert(nr==nc,'A must be square');
    %     n = nr;
    %     Acsr = build_csr(A);
    % end

    beta = norm(w);

    if n ~= n_prev_clip
        k1 = 2;
        mb = m;
        n_prev_clip = n;
    end

    % step to next print time, next step size, or next fixed event time
    t_step = min([Time_array(i_prt)-tNow,t_new,nextFixedTime-tNow,maxTimeStep]);
    k1_in = k1;
    mb_in = mb;
    t_step_in = t_step;

    % tic
    if nr<500
        orthDepth = m;%min(m,15);
    else
        orthDepth = min(m,30);
    end

    try
        [H,V,k1,mb,t_step] = ssit.fsp_ode_solvers.mexFunctionExpokit(n,m,w,beta,Acsr,btol,...
            Time_array(i_prt),tNow,k1_in,mb_in,t_step_in,orthDepth);
    catch
        [H,V,k1,mb,t_step] = ssit.fsp_ode_solvers.ExpensiveTask(n,m,w,beta,A,btol,...
            Time_array(i_prt),tNow,k1_in,mb_in,t_step_in,orthDepth);
    end
   
    if k1 ~= 0
        H(m+2,m+1) = 1;
        avnorm = norm(Ain*V(:,m+1));
    else
        avnorm = 0;
    end
    ireject = 0;

    % Define defaults for variables (needed for C conversion).
    F = zeros(size(H));
    err_loc = 0;

    while ireject <= mxrej
        mx = mb + k1;
        F = expm(sgn*t_step*H(1:mx,1:mx));
        if k1 == 0
            err_loc = btol;
            break;
        else
            phi1 = abs( beta*F(m+1,1) );
            phi2 = abs( beta*F(m+2,1) * avnorm );
            if phi1 > 10*phi2
                err_loc = phi2;
                xm = 1/m;
            elseif phi1 > phi2
                err_loc = (phi1*phi2)/(phi1-phi2);
                xm = 1/m;
            else
                err_loc = phi1;
                xm = 1/(m-1);
            end
        end
        if err_loc <= delta * t_step*tol
            break;
        else
            t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
            s = 10^(floor(log10(t_step))-1);
            t_step = ceil(t_step/s) * s;
            if ireject == mxrej
                %% Changes by Brian Munsky
                tryagain=1;
                P_lost=w(SINKS);
                %             time_of_excess = tNow/t;
                w=[];
                err=[];
                hump=[];
                Time_array_out=[];
                P_array=[];
                return
                %           error('The requested tolerance is too high.');
            end
            ireject = ireject + 1;
        end
    end
    mx = mb + max( 0,k1-1 );
    w = V(:,1:mx)*(beta*F(1:mx,1));
    beta = norm( w );
    hump = max(hump,beta);
    
    if resetSparsity
        WJ = w(1:SINKS(1)-1);
        WJ(WJ<1e-10) = 0;
        w(1:SINKS(1)-1) = WJ;
    end
    % neg = find(w < 0);
    % ineg = length(neg);
    % w(neg) = 0;
    % 
    wnorm = norm(w,1);
    % if ineg > 0
    %     w = (1/wnorm)*w;
    % end
    roundoff = abs(1.0d0-wnorm)/n;
    
    tNow = tNow + t_step;


    wout = zeros(size(vin));
    wout(Jclip) = w;

    if tNow == nextFixedTime
        % disp('fixed time reached')
        % implement fixed time effect
        wout = fixedEvents.matrices{fixedEvents.matrixInds(indNextFixedTime)}*wout;
        indNextFixedTime = find(fixedEvents.times>tNow,1,"first");
        nextFixedTime = fixedEvents.times(indNextFixedTime);
    end

    %% Changes by Brian Munsky
    tolCalc = fspTol*(tNow-fspErrorCondition.tInit)/(fspErrorCondition.tFinal-fspErrorCondition.tInit);
    if sum(wout(SINKS)) > tolCalc
        P_lost=wout(SINKS);
        err=[];
        hump=[];
        ye = wout';
        te = tNow-t_step;
        P_array = P_array(1:i_prt-1,:);       
        tryagain=0;
        return
    else
        ye = wout';
        te = tNow;
    end
    while i_prt<=length(Time_array)&&tNow>=Time_array(i_prt)
        Time_array_out = [Time_array_out;Time_array(i_prt)];
        P_array(i_prt,:) = wout';
        i_prt=i_prt+1;
    end
    %%
    
    t_new = gamma * t_step * (t_step*tol/err_loc)^xm;
    s = 10^(floor(log10(t_new))-1);
    t_new = ceil(t_new/s) * s;
    
    err_loc = max(err_loc,roundoff);
    err_loc = max(err_loc,rndoff);
    s_error = s_error + err_loc;
end
while i_prt<=length(Time_array)&&tNow>=Time_array(i_prt)
    Time_array_out = [Time_array_out;Time_array(i_prt)];
    P_array(i_prt,:) = wout';
    i_prt=i_prt+1;
end
err = s_error;
hump = hump / normv;
te = max(Time_array); ye = wout';
P_lost = wout(SINKS);

%% Changes by Brian Munsky
tryagain=0;
end

function [ir0, jc0, pr] = find_csc_like(S)
    % Returns 0-based CSC arrays similar to C sparse storage.
    [i,j,v] = find(S);                    % column-major order
    n = size(S,2);
    counts = accumarray(j,1,[n,1]);
    jc0 = [0; cumsum(counts)];            % 0-based column pointer
    ir0 = i - 1;                          % 0-based row indices
    pr = v;
end

function Acsr = build_csr(A)
    [row_idx, col_ptr, val] = find_csc_like(A);
    n = size(A,1);
    nnzA = numel(val);

    row_counts = accumarray(row_idx+1, 1, [n,1]);
    row_ptr = zeros(n+1,1);
    row_ptr(2:end) = cumsum(row_counts);

    col_ind = zeros(nnzA,1);
    val_csr = zeros(nnzA,1);

    next = row_ptr(1:end-1);
    for c = 1:n
        k0 = col_ptr(c) + 1;
        ak1 = col_ptr(c+1);
        for k = k0:ak1
            r = row_idx(k) + 1;
            dst = next(r) + 1;
            col_ind(dst) = c - 1;
            val_csr(dst) = val(k);
            next(r) = next(r) + 1;
        end
    end

    Acsr = struct('row_ptr', row_ptr, 'col_ind', col_ind, 'val', val_csr);
end