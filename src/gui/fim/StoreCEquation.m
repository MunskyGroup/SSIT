function StoreCEquation(app)
P = permute(app.FspTabOutputs.solutions,[2 3 4 1]);
if ~app.x1CheckBox_10.Value
    C1 = ones(1,size(P,1));
else
    switch app.ImageProcessingErrorsDropDown.Value
        case 'Error Free' 
            max_x1 = size(P,1);
            C1 = eye(max_x1);
        case 'Binomial'
            max_x1 = size(P,1);
            for i = 0:max_x1-1
                C1(:,i+1) = binopdf(0:max_x1-1,i,app.ProbEditField.Value)';
            end
        case 'Poisson'  
            max_x1 = size(P,1);
            for i = 0:max_x1-1
                C1(:,i+1) = poisspdf(0:max_x1-1,app.MeanEditField.Value)';
            end
        case 'Gaussian'
            max_x1 = size(P,1);
            for i = 0:max_x1-1
                C1(:,i+1) = normpdf(0:max_x1-1,app.MeanEditField,Value,sqrt(app.VarianceEditField.Value))';
            end
    end
    if app.PerfectNoBinningCheckBox.Value
        B = 1;
    elseif app.UniformCheckBox.Value
        Nbins = min(app.ofBinsEditField.Value,size(P,1));
        B = zeros(Nbins,size(P,1));
        Edges = [linspace(0,size(P,1),Nbins+1)];
        BINS = discretize(0:size(P,1)-1,Edges);
        for i = 1:Nbins
            B(i,BINS==i) = 1;
        end
    elseif app.SpecifiedCheckBox.Value
        % todo
    end
    C1 = B.*C1;
end
   
if ~app.x2CheckBox_10.Value
    C2 = ones(1,size(P,2));
else
    switch app.ImageProcessingErrorsDropDown_4.Value
        case 'Error Free' 
            max_x2 = size(P,2);
            C2 = eye(max_x2);
        case 'Binomial'
            max_x2 = size(P,2);
            for i = 0:max_x2-1
                C2(:,i+1) = binopdf(0:max_x2-1,i,app.ProbEditField.Value)';
            end
        case 'Poisson'  
            max_x2 = size_Q23(P,2);
            for i = 0:max_x2-1
                C2(:,i+1) = poisspdf(0:max_x2-1,app.MeanEditField_2.Value)';
            end
        case 'Gaussian'
            max_x2 = size_Q23(P,2);
            for i = 0:max_x2-1
                C2(:,i+1) = normpdf(0:max_x2-1,app.MeanEditField_2,Value,sqrt(app.VarianceEditField_2.Value))';
            end
    end
    if app.PerfectNoBinningCheckBox.Value
        B = 1;
    elseif app.UniformCheckBox.Value
        Nbins = app.ofBinsEditField.Value;
        B = zeros(Nbins,size(P,2));
        Edges = [0,linspace(0,size(P,2),Nbins)];
        BINS = discretize(0:size(P,2)-1,Edges);
        for i = 1:Nbins
            B(i,BINS==i) = 1;
        end
    elseif app.SpecifiedCheckBox.Value
        % todo
    end
    C2 = B.*C2;
end
if ~app.x3CheckBox_10.Value
    C3 = ones(1,size(P,3));
else
    switch app.ImageProcessingErrorsDropDown_5.Value
        case 'Error Free' 
            max_x3 = size(P,3);
            C3 = eye(max_x3);
        case 'Binomial'
            max_x3 = size(P,3);
            for i = 0:max_x3-1
                C3(:,i+1) = binopdf(0:max_x3-1,i,app.ProbEditField_3.Value)';
            end
        case 'Poisson'  
            max_x3 = size(P,3);
            for i = 0:max_x3-1
                C3(:,i+1) = poisspdf(0:max_x3-1,app.MeanEditField_3.Value)';
            end
        case 'Gaussian'
            max_x3 = size(P,3);
            for i = 0:max_x3-1
                C3(:,i+1) = normpdf(0:max_x3-1,app.MeanEditField_3,Value,sqrt(app.VarianceEditField_3.Value))';
            end
    end
    if app.PerfectNoBinningCheckBox.Value
        B = 1;
    elseif app.UniformCheckBox.Value
        Nbins = app.ofBinsEditField.Value;
        B = zeros(Nbins,size(P,3));
        Edges = [linspace(0,size(P,3),Nbins)];
        BINS = discretize(0:size(P,3)-1,Edges);
        for i = 1:Nbins
            B(i,BINS==i) = 1;
        end
    elseif app.SpecifiedCheckBox.Value
        % todo
    end
    C3 = B.*C3;
end
 