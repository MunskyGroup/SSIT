compile_app = true;

if exist('SSITGUI', 'file')
    answer = questdlg('Local SSITGUI.mlapp found. Do you want to overwrite? This will erase all unsaved changes you made to the GUI!!!');
    if strcmp(answer,'No') || strcmp(answer,'Cancel') || isempty(answer)
        compile_app = false;
    end
end

if compile_app
    fprintf('Packing app source files into a single app file...');
    copyfile('SSITGUI.mlapp',strrep(['SSITGUI_bak',date,'.mlapp'],'-','_'))
    zip('SSITGUI.mlapp', '*', 'appsrc');
    movefile('SSITGUI.mlapp.zip', 'SSITGUI.mlapp');
    fprintf('done.\n');
end