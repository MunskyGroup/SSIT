% build_expokit.m
% Compiles mexFunctionExpokit + expokitC into a MEX binary.
%
% Usage (from the MATLAB command window, in the folder containing both .c files):
%   build_expokit        % auto-detects platform
%   build_expokit debug  % adds -g flag for debugging symbols

function build_expokit(varargin)

    debug_flag = any(strcmpi(varargin, 'debug'));

    src_files = {'src/+ssit/+fsp_ode_solvers/mexFunctionExpokit.c', 'src/+ssit/+fsp_ode_solvers/expokitC.c'};

    % ----------------------------------------------------------------
    % Optimisation flags shared across platforms
    %   -O3          : full optimisation
    %   -ffast-math  : allows reassociation etc. (safe for this code)
    %   -march=native: tune for the build machine's CPU (remove for
    %                  portable binaries)
    % ----------------------------------------------------------------
    common_cflags = 'CFLAGS=$CFLAGS -O3 -ffast-math -march=native';

    if debug_flag
        common_cflags = 'CFLAGS=$CFLAGS -O0 -g';
    end

    % ----------------------------------------------------------------
    % Platform-specific BLAS linkage
    % ----------------------------------------------------------------
    if ismac
        % Apple Accelerate framework – always present on macOS.
        % -framework must go into LDFLAGS, not as a bare mex argument.
        blas_link = {'LDFLAGS=$LDFLAGS -framework Accelerate'};
        blas_inc  = {};

    elseif isunix
        % Linux: prefer MKL if available (e.g. Intel oneAPI), else OpenBLAS
        mkl_root = getenv('MKLROOT');
        if ~isempty(mkl_root)
            blas_inc  = {['-I' mkl_root '/include'], ...
                         'CFLAGS=$CFLAGS -DMKL_BLAS'};
            blas_link = {['-L' mkl_root '/lib/intel64'], ...
                         '-lmkl_rt', '-lpthread', '-lm', '-ldl'};
        else
            % OpenBLAS fallback
            blas_inc  = {};
            blas_link = {'-lopenblas', '-lpthread', '-lm'};
        end

    elseif ispc
        % Windows: link against MATLAB's own BLAS (lapack.lib ships with MATLAB)
        matlab_root = matlabroot;
        blas_inc    = {};
        blas_link   = {['-L' fullfile(matlab_root, 'extern', 'lib', computer('arch'))], ...
                       '-lmwblas'};
    else
        error('build_expokit:UnknownPlatform', 'Unsupported platform.');
    end

    % ----------------------------------------------------------------
    % Assemble and run the mex command
    % ----------------------------------------------------------------
    % cmd = [src_files, {common_cflags}, blas_inc, blas_link];
    cmd = [{'-R2018a'}, src_files, {common_cflags}, blas_inc, blas_link];
    fprintf('Building MEX: mex');
    for k = 1:numel(cmd), fprintf(' %s', cmd{k}); end
    fprintf('\n');

    mex(cmd{:});
    fprintf('Build successful.\n');
end