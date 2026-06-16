%% example_SI_SBML
% Example script to demonstrate the loading of a model from SBML

Model = SSIT('Empty');
Model = Model.createModelFromSBML('../SBML_test_cases/00010/00010-sbml-l1v2.xml',true);
Model = Model.formPropensitiesGeneral('SBMEModel');
Model.Solutions = Model.solve;
Model.plotFSP(speciesNames=Model.species, plotType='meansAndDevs')