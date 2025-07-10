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

function [w, err, hump, Time_array_out, P_array, P_lost, tryagain, te, ye] = mexpv_modified_2( t, A, v, tol, m, N_prt, Time_array,fspTol,SINKS,tNow,fspErrorCondition)

[n] = size(A,1);
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
%%

anorm = norm(A,'inf');
mxrej = 10;  btol  = 1.0e-7;
gamma = 0.9; delta = 1.2;
mb    = m; t_out   = abs(t);
s_error = 0;
rndoff= anorm*eps;

k1 = 2; xm = 1/m; normv = norm(v); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;

if t==0
    sgn = 1;
else
    sgn = sign(t); 
end

istep = 0;

w = v; ye = w'; te = tNow;
hump = normv;

%% Changes by Brian Munsky
Time_array_out=[];
i_prt=1;
P_array = zeros(length(Time_array),length(w));
%%
if tNow==Time_array(1)
    P_array(1,:)=w';
    Time_array_out=tNow;
    i_prt=2;
end

while tNow < t_out && i_prt<=length(Time_array)
    istep = istep + 1;
    t_step = min(Time_array(i_prt)-tNow,t_new);
    V = zeros(n,m+1);
    H = zeros(m+2,m+2);
    
    V(:,1) = (1/beta)*w;
    for j = 1:m
        p = A*V(:,j);
        for i = 1:j
            H(i,j) = p'*V(:,i);    %8.832
            p = p-H(i,j)*V(:,i);
        end
        s = norm(p);
        if s < btol
            k1 = 0;
            mb = j;
            t_step = Time_array(i_prt)-tNow;
            break;
        end
        H(j+1,j) = s;
        V(:,j+1) = (1/s)*p;
    end
    if k1 ~= 0
        H(m+2,m+1) = 1;
        avnorm = norm(A*V(:,m+1));
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
    
    neg = find(w < 0);
    ineg = length(neg);
    w(neg) = 0;
       
    wnorm = norm(w,1);
    if ineg > 0
        w = (1/wnorm)*w;
    end
    roundoff = abs(1.0d0-wnorm)/n;
    
    tNow = tNow + t_step;
    %% Changes by Brian Munsky
    tolCalc = fspTol*(tNow-fspErrorCondition.tInit)/(fspErrorCondition.tFinal-fspErrorCondition.tInit);
    if sum(w(SINKS)) > tolCalc
        P_lost=w(SINKS);
        err=[];
        hump=[];
        ye = w';
        te = tNow-t_step;
        P_array=[P_array(1:i_prt-1,:)];
        tryagain=0;
        return
    else
        ye = w';
        te = tNow;
    end
    while i_prt<=length(Time_array)&&tNow>=Time_array(i_prt)
        Time_array_out = [Time_array_out;Time_array(i_prt)];
        P_array(i_prt,:) = w;
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
    P_array(i_prt,:) = w;
    i_prt=i_prt+1;
end
err = s_error;
hump = hump / normv;
te = max(Time_array); ye = w;
P_lost = w(SINKS);

%% Changes by Brian Munsky
tryagain=0;
