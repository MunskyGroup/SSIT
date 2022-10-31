function bgFitModel2DataMerged(app,SSITapp,fit_data_name,fit_results_file_name,Merged)
% Since the background fits will be started in new instance, they will
% need to know what paths to use for the codes.
addpath(genpath('src'));
ssit.parest.fitMergedModel2Data(app,SSITapp,fit_data_name,fit_results_file_name,Merged)