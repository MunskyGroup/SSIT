clear
close all

%% Parameter Ranges
n_params = 5;
kon = 1e-3;
koff = 7.5e-5;
w = 0.0025;
kex = 750;
kr = 7.5e4;
D = [0.01,5,4];
gam =[0.035;0.0025;0.001];
makePlot = 0;
Ncells = 500;

kon_range = linspace(kon*0.7, kon*1.30, n_params);
koff_range = linspace(koff*0.7, koff*1.30, n_params);
w_range = linspace(w*0.7, w*1.30, n_params);
kex_range = linspace(kex*0.7, kex*1.30, n_params);
kr_range = linspace(kr*0.7, kr*1.30, n_params);

D_range = [linspace(D(1)*0.7, D(1)*1.30, n_params), linspace(D(2)*0.7, D(2)*1.30, n_params), linspace(D(3)*0.7, D(3)*1.30, n_params)];
gam_range = [linspace(gam(1)*0.7, gam(1)*1.30, n_params), linspace(gam(2)*0.7, gam(2)*1.30, n_params), linspace(gam(3)*0.7, gam(3)*1.30, n_params)];






%% Generate Simulated data
for i = 1:n_params
    for j = 1:n_params
        for k = 1:n_params
            for l = 1:n_params
                for m = 1:n_params
                    kon = kon_range(i);
                    koff = koff_range(j);
                    w = w_range(k);
                    kex = kex_range(l);
                    kr = kr_range(m);
                    D = [0.01,5,4];
                    gam =[0.035;0.0025;0.001];
                    a = 53*ones(Ncells,1);
                    b = 39*ones(Ncells,1);
                    u = 2*pi*rand(Ncells,1);
                    posnTS=[a.*cos(u),b.*sin(u)].*rand(Ncells,2); %position of TS for sim data

                    fileName = ['simulation_kon-', num2str(kon), '_koff-', num2str(koff), '_w-', num2str(w), '_kr-', num2str(kr), '_D1-', num2str(D(1)), '_D2-', num2str(D(2)), '_D3-', num2str(D(3)), '_gam1-', num2str(gam(1)), '_gam2-' num2str(gam(2)), '_gam3-', num2str(gam(3)), '.csv'];
                    
                    stochsProdTranspDegModel(kon,koff,w,kex,kr,D,gam,posnTS,makePlot,a,b,Ncells,fileName)

                end
            end
        end
    end
end















