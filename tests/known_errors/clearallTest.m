classdef clearallTest < matlab.unittest.TestCase

    methods (Test)
        function clear(testCase1)
            for i=1:6
                
                clear; clc; close all
                addpath(genpath('../src'));
                
                %% From STEP1 of example_DUSP1_Regulation.m
                Model1 = SSIT;  
                Model1.species = {'offGene';'onGene';'rna'}; 
                Model1.initialCondition = [2;0;0];           
                Model1.propensityFunctions = {'kon*IGR*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
                Model1.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'}; 
                Model1.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; 
                Model1.parameters = ({'koff',0.014;'kon',0.002;'kr',1;'gr',0.004;...
                                      'a1',20;'r1',0.04;'r2',0.1}); 
                Model1.fspOptions.initApproxSS = true; 
                Model1.summarizeModel                  
                Model1 = Model1.formPropensitiesGeneral('ToyDUSP1Model');
            end
        end
        function clearall(testCase2)
            for i=1:6
                
                clear all; clc; close all
                addpath(genpath('../src'));
                
                %% From STEP1 of example_DUSP1_Regulation.m
                Model1 = SSIT;  
                Model1.species = {'offGene';'onGene';'rna'}; 
                Model1.initialCondition = [2;0;0];           
                Model1.propensityFunctions = {'kon*IGR*offGene';'koff*onGene';'kr*onGene';'gr*rna'};         
                Model1.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))*(t>0)'}; 
                Model1.stoichiometry = [-1,1,0,0;1,-1,0,0;0,0,1,-1]; 
                Model1.parameters = ({'koff',0.014;'kon',0.002;'kr',1;'gr',0.004;...
                                      'a1',20;'r1',0.04;'r2',0.1}); 
                Model1.fspOptions.initApproxSS = true; 
                Model1.summarizeModel                  
                Model1 = Model1.formPropensitiesGeneral('ToyDUSP1Model');
            end
        end
    end
end