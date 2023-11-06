% Two species Poisson Process
clear all

k10 = 5;
k11 = 2;
om1 = 1;
k20 = 7;
k21 = 3;
om2 = 1;

k1 = @(t)k10+k11*sin(om1*t);
k2 = @(t)k20+k21*sin(om2*t);
g1 = 1;
g2 = 0.5;

N = [50,80];  % Projection size
A0 = zeros(prod(N+1),prod(N+1));
A1 = zeros(prod(N+1),prod(N+1));
A2 = zeros(prod(N+1),prod(N+1));
C1 = zeros(N(1)+1,prod(N+1));
C2 = zeros(N(2)+1,prod(N+1));

for i1 = 0:N(1)
    for i2 = 0:N(2)
        k = i1*(N(2)+1)+i2+1;
        if i1>0
            A0(k,k) = A0(k,k)-g1*i1;
            A0(k-(N(2)+1),k) = g1*i1;
        end
        if i2>0
            A0(k,k) = A0(k,k)-g2*i2;
            A0(k-1,k) = g2*i2;
        end
        if i1<N(1)
            A1(k,k) = A1(k,k)-1;
            A1(k+(N(2)+1),k) = 1;
        end
        if i2<N(2)
            A2(k,k) = A2(k,k)-1;
            A2(k+1,k) = 1;
        end
        C1(i1+1,k) = 1;
        C2(i2+1,k) = 1;
    end
end

P0 = zeros(prod(N+1),1); 
P0(1) = 1;

t = 2.4;

A0 = sparse(A0);
A1 = sparse(A1);
A2 = sparse(A2);

A = @(t,x)A0+k1(t)*A1+k2(t)*A2;

Apat = A(0,rand(size(P0)))~=0;

ode_opts = odeset('Jacobian', A, 'Vectorized','on','JPattern',Apat,...
    'relTol',1e-8, 'absTol', 1e-10,'NonNegative',true);
rhs = @(t,x)A(t,x)*x;
tic
[tExport, yout] =  ode15s(rhs, [0,t/2,t], P0);
timeODE23s = toc;

Pt = yout(end,:)';

P1t = C1*Pt;
P2t = C2*Pt;

%% Check results

lam1 = k10/g1*(1-exp(-g1*t)) + ...
    (k11*om1*exp(-g1*t))/(g1^2 + om1^2) - (k11*(om1*cos(om1*t) - g1*sin(om1*t)))/(g1^2 + om1^2);
y1 = poisspdf([0:N(1)]',lam1);
plot([0:N(1)],P1t,[0:N(1)],y1,'r--')
err1 = sum(abs(y1-P1t))
% 
hold on
lam2 = k20/g2*(1-exp(-g2*t)) + ...
    (k21*om2*exp(-g2*t))/(g2^2 + om2^2) - (k21*(om2*cos(om2*t) - g2*sin(om2*t)))/(g2^2 + om2^2);
y2 = poisspdf([0:N(2)]',lam2);
plot([0:N(2)],P2t,[0:N(2)],y2,'r--')
err2 = sum(abs(y2-P2t))





