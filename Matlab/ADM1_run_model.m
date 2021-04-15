% -------------------------------------------------------------------------
%
% Run different simplifications of a mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% Example: Anaerobic co-digestion of maize silage and cattle manure
%
% Implementation in Matlab 
% R2019b (The MathWorks, Inc.)
%
% Matlab ODE function
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

%% Run selected model 
switch model
    case 'ADM1'
        % Solve ODE of mass-based ADM1
        ode = ode15s(@(t,x) ADM1_mass(t,x,system.Variables,input.ADM1.Variables,parameters.ADM1.Variables),time_range,initial.ADM1.Variables);
        % Calculate model ouput
        num_t = size(ode.x,2);
        for i = 1:num_t
            [t(i),y(i,:)] = ADM1_mass_output(ode.x(1,i),ode.y(:,i),system.Variables,parameters.ADM1.Variables);
        end      
    case 'ADM1_R1'
        % Solve ODE of mass-based ADM1-R1
        ode = ode15s(@(t,x) ADM1_R1_mass(t,x,system.Variables,input.ADM1_R1.Variables,parameters.ADM1_R1.Variables),time_range,initial.ADM1_R1.Variables);
        % Calculate model ouput
        num_t = size(ode.x,2);
        for i = 1:num_t
            [t(i),y(i,:)] = ADM1_R1_mass_output(ode.x(1,i),ode.y(:,i),system.Variables,parameters.ADM1_R1.Variables);
        end
    case 'ADM1_R2'
        % Solve ODE of mass-based ADM1-R2
        ode = ode15s(@(t,x) ADM1_R2_mass(t,x,system.Variables,input.ADM1_R2.Variables,parameters.ADM1_R2.Variables),time_range,initial.ADM1_R2.Variables);
        % Calculate model ouput
        num_t = size(ode.x,2);
        for i = 1:num_t
            [t(i),y(i,:)] = ADM1_R2_mass_output(ode.x(1,i),ode.y(:,i),system.Variables,parameters.ADM1_R2.Variables);
        end      
    case 'ADM1_R3'
        % Solve ODE of mass-based ADM1-R3
        ode = ode15s(@(t,x) ADM1_R3_mass(t,x,system.Variables,input.ADM1_R3.Variables,parameters.ADM1_R3.Variables),time_range,initial.ADM1_R3.Variables);
        % Calculate model ouput
        num_t = size(ode.x,2);
        for i = 1:num_t
            [t(i),y(i,:)] = ADM1_R3_mass_output(ode.x(1,i),ode.y(:,i),system.Variables,parameters.ADM1_R3.Variables);
        end
    case 'ADM1_R4'
        % Solve ODE of mass-based ADM1-R4
        ode = ode15s(@(t,x) ADM1_R4_mass(t,x,system.Variables,input.ADM1_R4.Variables,parameters.ADM1_R4.Variables),time_range,initial.ADM1_R4.Variables);
        % Calculate model ouput
        num_t = size(ode.x,2);
        for i = 1:num_t
            [t(i),y(i,:)] = ADM1_R4_mass_output(ode.x(1,i),ode.y(:,i),system.Variables,parameters.ADM1_R4.Variables);
        end
    otherwise
        % Stop script if specific model name does not match available model structures
        warning(['Model name incorrect. The selected model name "',model,'" does not match available model structures.'])
        return
end

%% Set model output
output=output.(model);
output{1:num_t,:} = [t' y];

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
clearvars i num_t t y l time_range