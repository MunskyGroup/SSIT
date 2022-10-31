function conditionOnCallback(hObject,eventData,fieldId)

global conditionOnArray
str = get(hObject,'String');
conditionOnArray = string(conditionOnArray);
if str
    conditionOnArray(fieldId) = str;
else
    conditionOnArray(fieldId) = "0";
end