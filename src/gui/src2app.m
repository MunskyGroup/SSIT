compile_app = true;

if exist('src/gui/SSITGUI.mlapp', 'file')
    answer = questdlg('Local SSITGUI.mlapp found. Do you want to overwrite? This will erase all unsaved changes you made to the GUI!!!');
    if strcmp(answer,'No') || strcmp(answer,'Cancel') || isempty(answer)
        compile_app = false;
    end
end

if compile_app
    fprintf('Packing app source files into a single app file...');
    copyfile('src/gui/SSITGUI.mlapp',strrep(['src/gui/SSITGUI_bak',date,'.mlapp'],'-','_'))
    zip('src/gui/SSITGUI.mlapp', '*', 'src/gui/appsrc');
    movefile('src/gui/SSITGUI.mlapp.zip', 'src/gui/SSITGUI.mlapp');
    fprintf('done.\n');
end