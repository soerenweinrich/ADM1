function [t,y] = ADM1_mass_output(t,x,system,parameter) 

% -------------------------------------------------------------------------
%
% Original (mass-based)
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1
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
K_H_h2 =  parameter(3); 
K_I_IN =  parameter(4); 
K_I_c4 =  parameter(5); 
K_I_fa =  parameter(6); 
K_I_nh3 =  parameter(7); 
K_I_pro =  parameter(8); 
K_a_IN =  parameter(9); 
K_a_ac =  parameter(10); 
K_a_bu =  parameter(11); 
K_a_co2 =  parameter(12); 
K_a_pro =  parameter(13); 
K_a_va =  parameter(14); 
K_aa =  parameter(15); 
K_ac =  parameter(16); 
K_bu =  parameter(17); 
K_fa =  parameter(18); 
K_h2 =  parameter(19); 
K_pro =  parameter(20); 
K_su =  parameter(21); 
K_va =  parameter(22); 
K_w =  parameter(23); 
R =  parameter(24); 
T =  parameter(25); 
k_AB_IN =  parameter(26); 
k_AB_ac =  parameter(27); 
k_AB_bu =  parameter(28); 
k_AB_co2 =  parameter(29); 
k_AB_pro =  parameter(30); 
k_AB_va =  parameter(31); 
k_La =  parameter(32); 
k_ch =  parameter(33); 
k_dec =  parameter(34); 
k_li =  parameter(35); 
k_m_aa =  parameter(36); 
k_m_ac =  parameter(37); 
k_m_bu =  parameter(38); 
k_m_fa =  parameter(39); 
k_m_h2 =  parameter(40); 
k_m_pro =  parameter(41); 
k_m_su =  parameter(42); 
k_m_va =  parameter(43); 
k_p =  parameter(44); 
k_pr =  parameter(45); 
pK_l_aa =  parameter(46); 
pK_l_ac =  parameter(47); 
pK_l_h2 =  parameter(48); 
pK_u_aa =  parameter(49); 
pK_u_ac =  parameter(50); 
pK_u_h2 =  parameter(51); 
p_h2o =  parameter(52); 
 
% Define algebraic equations 
 
S_nh4_i = x(11) - x(31); 
S_co2 = x(10) - x(30); 
phi = x(24) + S_nh4_i/17 - x(30)/44 - x(29)/60 - x(28)/74 - x(27)/88 - x(26)/102 - x(25); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_h2 = x(32)*R*T/2; 
p_ch4 = x(33)*R*T/16; 
p_co2 = x(34)*R*T/44; 
p_gas = p_h2 + p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define output (states) 
 
for i = 1:34 
    y(i) = x(i);  
end 
 
% Define output (algebraic components)
 
y(35) = S_nh4_i; 
y(36) = S_co2; 
y(37) = phi; 
y(38) = S_H; 
y(39) = pH; 
y(40) = p_h2; 
y(41) = p_ch4; 
y(42) = p_co2; 
y(43) = p_gas; 
y(44) = q_gas; 