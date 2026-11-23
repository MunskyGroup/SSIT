function tf = calibrationOutlierMatches(calibration, outlierFraction, outlierSupport)
% ssit.parest.calibrationOutlierMatches - True if an SSA likelihood
% calibration (SSIT.runSequentialABC) was made with the given
% outlier-robust likelihood settings (kdeOutlierFraction and
% kdeOutlierSupport). Calibrations without these fields were made
% without the outlier background (fraction 0).
frac = 0; support = [];
if isfield(calibration, 'kdeOutlierFraction'), frac = calibration.kdeOutlierFraction; end
if isfield(calibration, 'kdeOutlierSupport'), support = calibration.kdeOutlierSupport; end
tf = frac == outlierFraction && (outlierFraction == 0 || isequal(support(:)', outlierSupport(:)'));
end
