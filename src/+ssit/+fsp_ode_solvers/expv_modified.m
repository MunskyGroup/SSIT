%  [w, err, hump] = expv( t, A, v, tol, m )
%  EXPV computes an approximation of w = exp(t*A)*v for a
%  general matrix A using Krylov subspace  projection techniques.
%  It does not compute the matrix exponential in isolation but instead,
%  it computes directly the action of the exponential operator on the 
%  operand vector. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = expv( t, A, v )
%  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
%
%  [w, err] = expv( t, A, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = expv( t, A, v, tol )
%  overrides default tolerance.
%
%  [w, err, hump] = expv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace,
%  and renders an approximation of the `hump'.
%
%  The hump is defined as:
%          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
%  It is used as a measure of the conditioning of the matrix exponential
%  problem. The matrix exponential is well-conditioned if hump = 1,
%  whereas it is poorly-conditioned if hump >> 1. However the solution
%  can still be relatively fairly accurate even when the hump is large
%  (the hump is an upper bound), especially when the hump and 
%  ||w(t)||/||v|| are of the same order of magnitude (further details in 
%  reference below).
%

%  In the above example, n could have been set to a larger value,
%  but the computation of w_matlab will be too long (feel free to
%  discard this computation).
%
%  See also MEXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials. 
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

function [w, err, hump, Time_array_out, P_array, P_lost, tryagain, te, ye] = expv_modified( t, A, v, tol, m, N_prt, Time_array,fspTol,SINKS,tNow,fspErrorCondition)
arguments
    t
    A
    v
    tol = 1e-7
    m = 30
    N_prt = []
    Time_array = []
    fspTol = 1e-4
    SINKS = []
    tNow = 0
    fspErrorCondition = struct('tInit',0);
end
n = size(A,1);
% if nargin == 3
%   tol = 1.0e-7;
%   m = min(n,30);
% end
% if nargin == 4
%   m = min(n,30);
%   %% Changes by Brian Munsky
%   N_prt=1;
%   Time_array=[0 t];
%   fspTol=1;
%   SINKS = [];
% end
% %% Changes by Brian Munsky
% if nargin == 5
%     N_prt=1;
%     Time_array=[0 t];
%     fspTol=1;
%     SINKS = [];
% end
% if nargin ==6
%     Time_array=linspace(0,t,N_prt);
%     fspTol=1;
%     SINKS = [];
% end
% if nargin>=7
%     N_prt = length(Time_array);
% end
% if nargin ==7
%     fspTol=1;
%     SINKS = [];
% end
% if nargin ==8
%     tol=min(tol,fspTol/10);
%     SINKS = [];
% end
% if nargin ==9
%     tol=min(tol,fspTol/10);
% end

tryagain=0;
anorm = norm(A,'inf'); 
mxrej = 10;  btol  = 1.0e-7; 
gamma = 0.9; delta = 1.2; 
mb    = m; 
t_out = abs(t);
s_error = 0;
rndoff= anorm*eps;

k1 = 2; xm = 1/m; normv = norm(v); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s; 
sgn = sign(t); nstep = 0;

w = v; ye = w'; te = tNow;
hump = normv;

%% Changes by Brian Munsky
Time_array_out=[];
i_prt=1;
P_array = zeros(length(Time_array),length(w));
if tNow==Time_array(1)
    P_array(1,:)=w';
    Time_array_out=tNow;
    i_prt=2;
end
%%
while tNow < t_out && i_prt<=length(Time_array)
  nstep = nstep + 1;
  t_step = min(Time_array(i_prt)-tNow,t_new);
  V = zeros(n,m+1); 
  H = zeros(m+2,m+2);

  V(:,1) = (1/beta)*w;
  for j = 1:m
     p = A*V(:,j);
     for i = 1:j
        H(i,j) = dot(V(:,i),p);
        p = p-H(i,j)*V(:,i);
     end
     s = norm(p); 
     if s < btol
        k1 = 0;
        mb = j;
        t_step = t_out-tNow;
        break;
     end
     H(j+1,j) = s;
     V(:,j+1) = (1/s)*p;
  end
  if k1 ~= 0
     H(m+2,m+1) = 1;
     avnorm = norm(A*V(:,m+1)); 
  end
  ireject = 0;
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
%             Time_array_out=[];
%             P_array=P_array;
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

  tNow = tNow + t_step;
  %% Changes by Brian Munsky
  if max(w(SINKS))*length(SINKS)>fspTol*(tNow-fspErrorCondition.tInit)/(max(Time_array)-fspErrorCondition.tInit)
      P_lost=w(SINKS);
      err=[];
      hump=[];
      ye = w';
      te = tNow;
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

  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
end
%% Changes by Brian Munsky
while i_prt<=length(Time_array)&&tNow>=Time_array(i_prt)
    Time_array_out = [Time_array_out;Time_array(i_prt)];
    P_array(i_prt,:) = w;
    i_prt=i_prt+1;
end
err = s_error;
hump = hump / normv;
te = max(Time_array); ye = w;
P_lost = w(SINKS);
tryagain=0;
