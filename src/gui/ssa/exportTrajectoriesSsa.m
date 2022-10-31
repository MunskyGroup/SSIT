function [] = exportTrajectoriesSsa(app)
 % Only creates the Trajectoy plot as a matlab figure
 samples = app.StochasticSimulationTabOutputs.samples;
            if isempty(samples)
                SsaRunButtonPushed(app, event);
            end
            Nsim = app.SsaNumSimField.Value;
            T_array = eval(app.PrintTimesEditField.Value);
            
            % Trajectory plot
            figure()
            hold('off');
            for isim = 1:Nsim
                if app.Ssax1CheckBox.Value
                    plot(T_array,samples(1,:,isim),'b'); hold('on');
                    
                end
                if app.Ssax2CheckBox.Value
                    plot(T_array,samples(2,:,isim),'r'); hold('on');
                    
                end
                if app.Ssax3CheckBox.Value
                    plot(T_array,samples(3,:,isim),'g'); hold('on');
                    
                end
                set(gca,'fontsize',20);
                title('Trajectories');
                xlabel('Time');
                ylabel('Response');
                legend(char(app.NameTable.Data(1,2)),char(app.NameTable.Data(2,2)),char(app.NameTable.Data(3,2)));
            end
            hold('off')