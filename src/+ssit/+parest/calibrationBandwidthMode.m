function mode = calibrationBandwidthMode(calibration)
% ssit.parest.calibrationBandwidthMode - Bandwidth mode recorded in a
% sequential-ABC calibration struct ('lcv' for calibrations made before
% the mode was recorded; 'fixed' for those with fixBandwidth = true).
if isfield(calibration, 'bandwidthMode')
    mode = calibration.bandwidthMode;
elseif isfield(calibration, 'fixBandwidth') && calibration.fixBandwidth
    mode = 'fixed';
else
    mode = 'lcv';
end
end
