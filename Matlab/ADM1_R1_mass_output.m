function [t,y] = ADM1_R1_mass_output(t,x,system,parameter) 

% -------------------------------------------------------------------------
%
% Simplified (mass-based)
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R1
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
K_a_bu =  parameter(7); 
K_a_co2 =  parameter(8); 
K_a_pro =  parameter(9); 
K_a_va =  parameter(10); 
K_aa =  parameter(11); 
K_ac =  parameter(12); 
K_bu =  parameter(13); 
K_fa =  parameter(14); 
K_pro =  parameter(15); 
K_su =  parameter(16); 
K_va =  parameter(17); 
K_w =  parameter(18); 
R =  parameter(19); 
T =  parameter(20); 
k_AB_IN =  parameter(21); 
k_AB_ac =  parameter(22); 
k_AB_bu =  parameter(23); 
k_AB_co2 =  parameter(24); 
k_AB_pro =  parameter(25); 
k_AB_va =  parameter(26); 
k_La =  parameter(27); 
k_ch =  parameter(28); 
k_dec =  parameter(29); 
k_li =  parameter(30); 
k_m_aa =  parameter(31); 
k_m_ac =  parameter(32); 
k_m_bu =  parameter(33); 
k_m_fa =  parameter(34); 
k_m_pro =  parameter(35); 
k_m_su =  parameter(36); 
k_m_va =  parameter(37); 
k_p =  parameter(38); 
k_pr =  parameter(39); 
pK_l_aa =  parameter(40); 
pK_l_ac =  parameter(41); 
pK_u_aa =  parameter(42); 
pK_u_ac =  parameter(43); 
p_h2o =  parameter(44); 
 
% Define algebraic equations 
 
S_nh4_i = x(10) - x(29); 
S_co2 = x(9) - x(28); 
phi = x(22) + S_nh4_i/17 - x(28)/44 - x(27)/60 - x(26)/74 - x(25)/88 - x(24)/102 - x(23); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_ch4 = x(30)*R*T/16; 
p_co2 = x(31)*R*T/44; 
p_gas = p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define output (states) 
 
for i = 1:31
    y(i) = x(i);  
end 
 
% Define output (algebraic components)
 
y(32) = S_nh4_i; 
y(33) = S_co2; 
y(34) = phi; 
y(35) = S_H; 
y(36) = pH; 
y(37) = p_ch4; 
y(38) = p_co2; 
y(39) = p_gas; 
y(40) = q_gas; 