function [] = exportBothFiguresSsa(app,trajs,hists)
arguments
    app
    trajs = true
    hists = true
end
%% Extract data to plot
samples = app.StochasticSimulationTabOutputs.samples;
%% This creates figures from the SSA tab into a Matlab figure
if isempty(samples)
    msgbox('Please generate the trajectories first.');
end
Nsim = app.SsaNumSimField.Value;
T_array = eval(app.PrintTimesEditField.Value);
cols = {'b','r','c','m','g','b--','r--','c--','m--','g--'};

%% Trajectory plot
if trajs
    figure()
    hold('off');
    legs = {};
    for isim = 1:Nsim
        for jSpec = 1:size(samples,1)
            if any(contains(app.SpeciestoShowListBox.Value,app.SSITModel.species{jSpec}))
                plot(T_array,samples(jSpec,:,isim),cols{jSpec}); hold('on');
                if isim==1
                    legs{end+1} = app.SSITModel.species{jSpec};
                end
            end
        end

        set(gca,'fontsize',20);
        title('Trajectories');
        xlabel('Time');
        ylabel('Response');
        legend(legs)
    end
    hold('off')
end

if hists
    % Histogram plot
    if Nsim>1
        figure()
        Nbin = ceil(sqrt(Nsim));                %%% Also need to check the bins here
        [~,j] = min(abs(T_array-app.SsaTimeSlider.Value));
        hold('off');
        legs = {};
        for jSpec = 1:size(samples,1)
            if any(contains(app.SpeciestoShowListBox.Value,app.SSITModel.species{jSpec}))
                if max(samples(1,j,:))<=Nbin
                    x = [0:max(samples(jSpec,j,:))];
                else
                    x = linspace(0,max(samples(jSpec,j,:)),Nbin);
                end

                [h] = hist(squeeze(samples(jSpec,j,:)),x);
                stairs(x,h/Nsim,cols{jSpec},'Linewidth',2);  hold('on');
                legs{end+1} = app.SSITModel.species{jSpec};

            end
        end
        set(gca,'fontsize',20);
        title(['histogram at t = ',num2str(T_array(j))]);
        xlabel('Molecule count');
        ylabel('Probability');
        legend(legs)
    end
end