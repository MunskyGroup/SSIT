function [] = exportHistogramSsa(app)
% Only creates the histogram as a matlab figure
samples = app.StochasticSimulationTabOutputs.samples;
if isempty(samples)
    SsaRunButtonPushed(app, event);
end
Nsim = app.SsaNumSimField.Value;
T_array = eval(app.PrintTimesEditField.Value);

% Histogram plot
if Nsim>1
    figure()
    Nbin = ceil(sqrt(Nsim));                        %%% Fix bin sizes here too!!!
    [~,j] = min(abs(T_array-app.SsaTimeSlider.Value));
    hold('off');
    if app.Ssax1CheckBox.Value
        if max(samples(1,j,:))<=Nbin
            x = [0:max(samples(1,j,:))];
        else
            x = linspace(0,max(samples(1,j,:)),Nbin);
        end
        [h] = hist(squeeze(samples(1,j,:)),x);
        stairs(x,h/Nsim,'b','Linewidth',2);  hold('on');
    end
    if app.Ssax2CheckBox.Value
        if max(samples(2,j,:))<=Nbin
            x = [0:max(samples(2,j,:))];
        else
            x = linspace(0,max(samples(2,j,:)),Nbin);
        end
        [h] = hist(squeeze(samples(2,j,:)),x);
        stairs(x,h/Nsim,'r','Linewidth',2);   hold('on');
    end
    if app.Ssax3CheckBox.Value
        if max(samples(3,j,:))<=Nbin
            x = [0:max(samples(3,j,:))];
        else
            x = linspace(0,max(samples(3,j,:)),Nbin);
        end
        [h] = hist(squeeze(samples(3,j,:)),x);
        stairs(x,h/Nsim,'g','Linewidth',2);   hold('on');
    end
    set(gca,'fontsize',20);
    title(['histogram at t = ',num2str(T_array(j))]);
    xlabel('Molecule count');
    ylabel('Probability');
    legend(char(app.NameTable.Data(1,2)),char(app.NameTable.Data(2,2)),char(app.NameTable.Data(3,2)))
end