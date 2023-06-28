% Two species Poisson Process
clear all

k1 = 10;
k2 = 13;
g1 = 1;
g2 = 0.5;

N = [50,65];  % Projection size
A = zeros(prod(N+1),prod(N+1));
C1 = zeros(N(1)+1,prod(N+1));
C2 = zeros(N(2)+1,prod(N+1));

for i1 = 0:N(1)
    for i2 = 0:N(2)
        k = i1*(N(2)+1)+i2+1;
        if i1>0
            A(k,k) = A(k,k)-g1*i1;
            A(k-(N(2)+1),k) = g1*i1;
        end
        if i2>0
            A(k,k) = A(k,k)-g2*i2;
            A(k-1,k) = g2*i2;
        end
        if i1<N(1)
            A(k,k) = A(k,k)-k1;
            A(k+(N(2)+1),k) = k1;
        end
        if i2<N(2)
            A(k,k) = A(k,k)-k2;
            A(k+1,k) = k2;
        end
        C1(i1+1,k) = 1;
        C2(i2+1,k) = 1;
    end
end

P0 = zeros(prod(N+1),1); 
P0(1) = 1;

t = 2.4;
Pt = expm(A*t)*P0;

P1t = C1*Pt;
P2t = C2*Pt;
%% Check results
lam1 = k1/g1*(1-exp(-g1*t));
y1 = poisspdf([0:N(1)]',lam1);
plot([0:N(1)],P1t,[0:N(1)],y1,'r--')
err1 = sum(abs(y1-P1t))

hold on
lam2 = k2/g2*(1-exp(-g2*t));
y2 = poisspdf([0:N(2)]',lam2);
plot([0:N(2)],P2t,[0:N(2)],y2,'r--')
err2 = sum(abs(y2-P2t))





