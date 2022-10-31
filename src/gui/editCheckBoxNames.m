function editCheckBoxNames(app)
% This function changes the names of all the species check boxes to match
% the table in the Reactions Tab.
newName = app.NameTable.Data(:,2);
fieldNames = fieldnames(app);
x1CheckBoxNames = fieldNames(contains(fieldNames,'X1CheckBox','IgnoreCase',true));
x2CheckBoxNames = fieldNames(contains(fieldNames,'X2CheckBox','IgnoreCase',true));
x3CheckBoxNames = fieldNames(contains(fieldNames,'X3CheckBox','IgnoreCase',true));
x1x2CheckBoxNames = fieldNames(contains(fieldNames,'X1X2CheckBox','IgnoreCase',true));
x1x3CheckBoxNames = fieldNames(contains(fieldNames,'X1X3CheckBox','IgnoreCase',true));
x2x3CheckBoxNames = fieldNames(contains(fieldNames,'X2X3CheckBox','IgnoreCase',true));
for i = 1:numel(x1CheckBoxNames)
    app.(x1CheckBoxNames{i}).Text = newName(1);
end
for i = 1:numel(x2CheckBoxNames)
    app.(x2CheckBoxNames{i}).Text = newName(2);
end
for i = 1:numel(x3CheckBoxNames)
    app.(x3CheckBoxNames{i}).Text = newName(3);
end 
for i = 1:numel(x1x3CheckBoxNames)
    app.(x1x3CheckBoxNames{i}).Text = strcat(newName(1),',',newName(3));
end 
for i = 1:numel(x2x3CheckBoxNames)
    app.(x2x3CheckBoxNames{i}).Text = strcat(newName(2),',',newName(3));
end 
for i = 1:numel(x1x2CheckBoxNames)
    app.(x1x2CheckBoxNames{i}).Text = strcat(newName(1),',',newName(2));
end 