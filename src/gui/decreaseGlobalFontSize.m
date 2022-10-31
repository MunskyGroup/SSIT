function [] = decreaseGlobalFontSize(app)
% This function increases the font size of all fields within the app.

nms = fields(app);
for i = 1:length(nms)
    try
        eval(['app.',nms{i},'.FontSize = app.',nms{i},'.FontSize - 1;']);
    catch
    end
end
