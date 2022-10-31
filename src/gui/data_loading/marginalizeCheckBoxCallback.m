function marginalizeCheckBoxCallback(hObject,eventData,checkboxId)

global marginalMatrix
rowInd = checkboxId(1);
colInd = checkboxId(2);
value = get(hObject,'Value');

if value
    marginalMatrix(rowInd,colInd) = 1;
else
    marginalMatrix(rowInd,colInd) = 0;
end