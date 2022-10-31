compile_app = true;

if exist('SSIT_GUI', 'file')
    answer = questdlg('Local SSIT_GUI.mlapp found. Do you want to overwrite? This will erase all unsaved changes you made to the GUI!!!');
    if strcmp(answer,'No') || strcmp(answer,'Cancel') || isempty(answer)
        compile_app = false;
    end
end

if compile_app
    fprintf('Packing app source files into a single app file...');
    zip('SSIT_GUI.mlapp', '*', 'app_src');
    movefile('SSIT_GUI.mlapp.zip', 'SSIT_GUI.mlapp');
    fprintf('done.\n');
end