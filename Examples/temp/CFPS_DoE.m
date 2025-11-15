
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell-free protein synthesis (CFPS) mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries:
% clear
% close all
addpath(genpath('../src'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create an SSIT instance and call it 'CFPS_DoE':
CFPS_DoE = SSIT;    

%% Simple CFPS mechanistic model with DoE-relevant parameters

% Set species names for CFPS_DoE:
%   mRNA   : transcripts
%   Pimm   : immature (non-fluorescent) protein
%   P      : mature fluorescent protein (what your plate reader sees)
%   E      : lumped "energy" tokens (ATP/GTP equivalents)
CFPS_DoE.species = {'mRNA'; 'Pimm'; 'P'; 'E'};

% Set initial condition:
%   Start with no transcripts or protein, finite energy pool E0.
CFPS_DoE.initialCondition = [0; 0; 0; 1e6];   % [mRNA; Pimm; P; E]

% Set stoichiometry of reactions:
% Reactions:
%   1: ∅ + E    → mRNA + (E - 1)          (transcription, energy-consuming)
%   2: mRNA     → ∅                       (mRNA degradation)
%   3: mRNA + E → mRNA + Pimm + (E - 1)   (translation, energy-consuming)
%   4: Pimm     → P                       (protein maturation)
%   5: P        → ∅                       (protein degradation / loss)
%   6: ∅        → E                       (energy regeneration from 3-PGA, CFE, etc.)
CFPS_DoE.stoichiometry = [ ...
    1, -1,  0,  0,  0,  0;  ... % mRNA
    0,  0,  1, -1,  0,  0;  ... % Pimm
    0,  0,  0,  1, -1,  0;  ... % P
   -1,  0, -1,  0,  0,  1];     % E
%                 Reactions: 1, 2, 3, 4, 5, 6

% Optional: define "effective" transcription and translation rates as
% inputExpressions so you can encode DoE-identified interactions
% (e.g., Mg–PEG, 3-PGA–folinic, 3-PGA–spermidine) at the level of rate
% constants. These are time-independent here, but can depend on parameters.
CFPS_DoE.inputExpressions = { ...
    'k_tx_eff', ...
        ['k_tx0 * (1 + alpha_Mg_tx*(Mg - Mg_ref)) * ' ...
         '(1 + alpha_PEG_tx*(PEG - PEG_ref))']; ...
    'k_tl_eff', ...
        ['k_tl0 * (1 + alpha_3P_fol*(ThreePGA*Folinic - ThreePGA_ref*Folinic_ref)) * ' ...
         '(1 + alpha_3P_Sp*(ThreePGA*Spermidine - ThreePGA_ref*Spermidine_ref))'] ...
    };

% Set propensity functions (mass-action with effective rates):
%   a1 = k_tx_eff * E
%   a2 = gamma_m * mRNA
%   a3 = k_tl_eff * mRNA * E
%   a4 = k_mat * Pimm
%   a5 = gamma_p * P
%   a6 = k_reg * ThreePGA * CFE      (energy regeneration driven by 3-PGA & lysate)
CFPS_DoE.propensityFunctions = { ...
    'k_tx_eff * E';          ... % 1) transcription (energy-consuming)
    'gamma_m * mRNA';        ... % 2) mRNA degradation
    'k_tl_eff * mRNA * E';   ... % 3) translation (energy-consuming)
    'k_mat * Pimm';          ... % 4) protein maturation
    'gamma_p * P';           ... % 5) protein degradation / loss
    'k_reg * ThreePGA * CFE' ... % 6) energy regeneration from 3-PGA + lysate
    };

% Add parameters:
%   Core kinetic parameters:
%     k_tx0, k_tl0           : baseline transcription / translation rates
%     gamma_m                : mRNA degradation rate
%     k_mat                  : maturation rate
%     gamma_p                : protein loss rate
%     k_reg                  : baseline energy regeneration rate
%
%   DoE factors (buffer components, lysate concentration):
%     Mg, PEG, ThreePGA, Folinic, Spermidine, CFE, NAD, cAMP
%
%   Reference levels and interaction strengths (so you can encode linear
%   or weakly nonlinear effects consistent with your DSD/RSM fits):
%     Mg_ref, PEG_ref, ThreePGA_ref, Folinic_ref, Spermidine_ref
%     alpha_Mg_tx, alpha_PEG_tx, alpha_3P_fol, alpha_3P_Sp
CFPS_DoE.parameters = ({ ...
    % Core kinetic rates (example values; adjust/fit as needed)
    'k_tx0',          1e-4;  ... % baseline transcription rate per energy token
    'k_tl0',          1e-6;  ... % baseline translation rate per mRNA·energy
    'gamma_m',        1e-3;  ... % mRNA degradation rate
    'k_mat',          5e-3;  ... % immature → mature protein rate
    'gamma_p',        1e-4;  ... % mature protein loss rate
    'k_reg',          1e-3;  ... % baseline energy regeneration rate

    % DoE factor levels (set these from your actual experimental conditions)
    'Mg',             12;    ... % mM Mg-glutamate (example)
    'PEG',            3;     ... % % (w/v) PEG-8000 (example)
    'ThreePGA',       30;    ... % mM 3-PGA (example)
    'Folinic',        0.1;   ... % mM folinic acid (example)
    'Spermidine',     1;     ... % mM spermidine (example)
    'CFE',            0.4;   ... % µL CFE / µL reaction (effective lysate conc.)
    'NAD',            0.5;   ... % mM NAD (example; not explicitly used yet)
    'cAMP',           0.5;   ... % mM cAMP (example; not explicitly used yet)

    % Reference levels for interaction terms (e.g., center points of your RSM)
    'Mg_ref',         12;    ... % reference Mg (e.g., central DoE level)
    'PEG_ref',        3;     ... % reference PEG
    'ThreePGA_ref',   30;    ... % reference 3-PGA
    'Folinic_ref',    0.1;   ... % reference folinic acid
    'Spermidine_ref', 1;     ... % reference spermidine

    % Interaction strengths (dimensionless; to be fit from data / RSM)
    'alpha_Mg_tx',    0.0;   ... % Mg effect on transcription
    'alpha_PEG_tx',   0.0;   ... % PEG effect on transcription
    'alpha_3P_fol',   0.0;   ... % 3-PGA × folinic effect on translation
    'alpha_3P_Sp',    0.0    ... % 3-PGA × spermidine effect on lag/translation
    });

% Print a summary of CFPS model:
CFPS_DoE.summarizeModel
