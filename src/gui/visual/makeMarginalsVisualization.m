function makeMarginalsVisualization(app)
% This function creates graphs of the different marginals at the desired
% time point in the FSP 3 Species GUI. This will use the time selected from
% the time slider as the time point from each of the different combinations
% based on the species selected.

T_array = eval(app.FspPrintTimesField.Value);            
[~,j] = min(abs(T_array-app.FspTimeSlider.Value));
Z = squeeze(app.FspTabOutputs.solutions(j,:,:,:));


    figure()
    hold('off')
    if app.FspMarginalX1CheckBox.Value == 1
        stairs([0:size(Z,1)-1],squeeze(sum(sum(Z,2),3)),'b');  hold('on');
        title(sprintf('Marginals at time: t = %1.2f',T_array(j)));
        xlabel('Species Count')
        ylabel('Probability')
        
    end
    if app.FspMarginalX2CheckBox.Value == 1
        stairs([0:size(Z,2)-1],squeeze(sum(sum(Z,1),3)),'r');  hold('on');
        title(sprintf('Marginals at time: t = %1.2f',T_array(j)));
        xlabel('Species Count')
        ylabel('Probability')
    end
    if app.FspMarginalX3CheckBox.Value == 1
        stairs([0:size(Z,3)-1],squeeze(sum(sum(Z,1),2)),'g');  hold('on');
        title(sprintf('Marginals at time: t = %1.2f',T_array(j)));
        xlabel('Species Count')
        ylabel('Probability')
    end
    hold('off');



end


