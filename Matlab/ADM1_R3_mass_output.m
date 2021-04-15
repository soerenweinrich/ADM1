function [t,y] = ADM1_R3_mass_output(t,x,system,parameter) 

% -------------------------------------------------------------------------
%
% Simplified (mass-based)
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R3
%
% Implementation in Matlab 
% R2019b (The MathWorks, Inc.)
%
% Matlab ODE function 
%
% Version 1.1
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

% System parameters  
 
V_liq = system(1); 
V_gas = system(2); 
p_atm = system(3); 
 
% Model parameters 
 
K_H_ch4 =  parameter(1); 
K_H_co2 =  parameter(2); 
K_I_IN =  parameter(3); 
K_I_nh3 =  parameter(4); 
K_a_IN =  parameter(5); 
K_a_ac =  parameter(6); 
K_a_co2 =  parameter(7); 
K_ac =  parameter(8); 
K_w =  parameter(9); 
R =  parameter(10); 
T =  parameter(11); 
k_AB_IN =  parameter(12); 
k_AB_ac =  parameter(13); 
k_AB_co2 =  parameter(14); 
k_La =  parameter(15); 
k_ch =  parameter(16); 
k_dec =  parameter(17); 
k_li =  parameter(18); 
k_m_ac =  parameter(19); 
k_p =  parameter(20); 
k_pr =  parameter(21); 
pK_l_ac =  parameter(22); 
pK_u_ac =  parameter(23); 
p_h2o =  parameter(24); 
 
% Define algebraic equations 
 
S_nh4_i = x(4) - x(15); 
S_co2 = x(3) - x(14); 
phi = x(11) + S_nh4_i/17 - x(14)/44 - x(13)/60  - x(12); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_ch4 = x(16)*R*T/16; 
p_co2 = x(17)*R*T/44; 
p_gas = p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define output (states) 
 
for i = 1:17
    y(i) = x(i);  
end 
 
% Define output (algebraic components) 
 
y(18) = S_nh4_i; 
y(19) = S_co2; 
y(20) = phi; 
y(21) = S_H; 
y(22) = pH; 
y(23) = p_ch4; 
y(24) = p_co2; 
y(25) = p_gas; 
y(26) = q_gas; 