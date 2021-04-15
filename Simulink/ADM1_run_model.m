% -------------------------------------------------------------------------
%
% Run different simplifications of a mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% Example: Anaerobic co-digestion of maize silage and cattle manure
%
% Implementation in Matlab/Simulink 
% R2019b (The MathWorks, Inc.)
%
% Level 2 S-Function
%
% Version 1.0
%
% https://github.com/soerenweinrich/ADM1
%
% Copyright (c) 2021 Sören Weinrich
% E-Mail: soeren.weinrich@dbfz.de
%
% Additional information (citation example):
%
% Weinrich, S.; Nelles, M. (2021).
% Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
% Model development and stoichiometric analysis. Bioresource Technology. 
% In press. https://doi.org/10.1016/j.biortech.2021.125124.
%
% -------------------------------------------------------------------------

clear all
clc

%% Initialise model
% Load standard model parameters
load('Model_data\ADM1_parameters.mat')
% Load experimental data
load('Model_data\ADM1_input_data.mat')
% Select model type: ADM1, ADM1_R1, ADM1_R2, ADM1_R3 or ADM1_R4
model = 'ADM1';
% Set time range in days
time_range = [0 100];
% Open Simulink model
open_system('ADM1') 

%% Run selcted model
% Set model type and parameters
h = getSimulinkBlockHandle('ADM1/Model');
set_param(h,'FunctionName',[model,'_mass']);
set_param(h,'Parameters',['initial.',model,'.Variables'',parameters.',model,'.Variables,system.Variables,3,0']);
% Set input composition
h = getSimulinkBlockHandle('ADM1/Maize silage');
set_param(h,'VariableName',['input.',model,'.Maize.Variables']);
h = getSimulinkBlockHandle('ADM1/Cattle manure');
set_param(h,'VariableName',['input.',model,'.Manure.Variables']);
h = getSimulinkBlockHandle('ADM1/Water');
set_param(h,'VariableName',['input.',model,'.Water.Variables']);
% Run simulation
out = sim('ADM1');

%% Set model output
num_t = size(out.tout,1);
output=output.(model);
output{1:num_t,:} = [out.tout out.yout(:,2:end)];
 
%% Plot model output
plot(output{:,1},output{:,2:end});
plotbrowser('on');
% Set legend
l = legend(output.Properties.VariableNames(2:end));
set(l,'Interpreter','none','visible','off');
% Set title
t = title(['Simulation results of the ',model]);
set(t,'Interpreter','none');
% Set axis labels
xlabel('Time [d]');
ylabel('Model output'); 

%% Clear variables 
clearvars h num_t t l time_range