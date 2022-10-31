function [] = plotMomentsInGui(app)

tspan = eval(app.MomentsPrintTimesField.Value);  % Pulls the time 
solutions = app.MomentTabOutputs.moments;

try
    %% plot average response
    legends = {};
    LG = {};
    hold(app.MomentsAxes,'off');
    Plts_to_make = [app.Momentsx1CheckBox.Value,app.Momentsx2CheckBox.Value,app.Momentsx3CheckBox.Value];
    cols = ['b','r','g'];
    
    % plot mean values
    n=3;
    var=[solutions(:,n+1)-solutions(:,1).^2, solutions(:,n+n+1)-solutions(:,2).^2, solutions(:,n+n+(n-1)+1)-solutions(:,3).^2];
    stdUpper=solutions(:,1:3)+sqrt(abs(var));
    stdLower=solutions(:,1:3)-sqrt(abs(var));
    for iplt=1:3
        if Plts_to_make(iplt)
            plot(app.MomentsAxes,tspan,solutions(:,iplt),cols(iplt));hold(app.MomentsAxes,'on');
            legends{end+1} = [char(app.NameTable.Data(iplt,2))];
            if app.MomentsShowVarianceCheckBox.Value == 1
                if all(solutions(:,iplt)==0)
                else
                    Lower=stdLower(:,iplt)';
                    Upper=stdUpper(:,iplt)';
                    mask = Upper > Lower;
                    fx=[tspan(mask), fliplr(tspan(mask))];
                    fy=[Lower(mask), fliplr(Upper(mask))];
                    std=fill(app.MomentsAxes,fx,fy,cols(iplt)); hold(app.MomentsAxes,'on');
                    alpha(std,0.3);
                    std.EdgeColor = 'none';
                    legends{end+1} = [char(app.NameTable.Data(iplt,2)),' \pm std'];
                end
            end
        end
    end
    legend(app.MomentsAxes,legends,'Location','best')
    
    %% plot histogram
    hold(app.MomentsHistAxes,'off');
    maxValue=max(max(solutions(:,1:3)));
    HistMax=maxValue+maxValue/2;
    x=linspace(0,HistMax,101);
    [~, timeIndex] = min(abs(tspan-app.MomentsTimeSlider.Value));
    for iplt=1:3
        if Plts_to_make(iplt)
            y=normpdf(x,solutions(timeIndex,iplt),sqrt(abs(var(timeIndex,iplt))));
            plot(app.MomentsHistAxes,x,y,cols(iplt)); hold(app.MomentsHistAxes,'on');
            legends{end+1} = [char(app.NameTable.Data(iplt,2))];
        end
        title(app.MomentsHistAxes,['Histogram at t = ',num2str(tspan(timeIndex))]);
    end
    legend(app.MomentsHistAxes,legends,'Location','best')
catch
    
end