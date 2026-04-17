function [H,V,k1,mb,t_step] = ExpensiveTask(n,m,w,beta,A,btol, ...
    Time_array_i_prt,tNow,k1,mb,t_step)

arguments
    n  % positive integer
    m  % positive integer
    w  % long dense non-negative vector
    beta % postive scalar
    A  % very large and very sparse double matrix 
    btol  % positive scalar
    Time_array_i_prt % positive scalar
    tNow % scalar
    k1
    mb
    t_step
end

V = zeros(n,m+1);
H = zeros(m+2,m+2);

V(:,1) = (1/beta)*w;
for j = 1:m
    p = A*V(:,j);
    p(abs(p)<1e-10) = 0;
    for i = 1:j
        H(i,j) = p'*V(:,i);    %8.832
        p = p-H(i,j)*V(:,i);
    end
    s = norm(p);
    if s < btol
        k1 = 0;
        mb = j;
        t_step = Time_array_i_prt-tNow;
        break;
    end
    H(j+1,j) = s;
    V(:,j+1) = (1/s)*p;
end
