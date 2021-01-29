function [t,y,rate,inhibition] = ADM1_R2_mass_output(t,x,s,input,parameter) 

% -------------------------------------------------------------------------
%
% Simplified and mass-based
% Anaerobic Digestion Model No. 1 (ADM1)
%
% ADM1-R2
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
% Citation example:
%
% Weinrich, S.; Nelles, M. (2021).
% Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
% Model development and stoichiometric analysis.
% Submitted to Bioresource Technology.
%
% -------------------------------------------------------------------------

% System parameters  
 
V_liq = s(1); 
V_gas = s(2); 
p_atm = s(4); 

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
K_ac =  parameter(11); 
K_bu =  parameter(12); 
K_pro =  parameter(13); 
K_va =  parameter(14); 
K_w =  parameter(15); 
R =  parameter(16); 
T =  parameter(17); 
k_AB_IN =  parameter(18); 
k_AB_ac =  parameter(19); 
k_AB_bu =  parameter(20); 
k_AB_co2 =  parameter(21); 
k_AB_pro =  parameter(22); 
k_AB_va =  parameter(23); 
k_La =  parameter(24); 
k_ch =  parameter(25); 
k_dec =  parameter(26); 
k_li =  parameter(27); 
k_m_ac =  parameter(28); 
k_m_bu =  parameter(29); 
k_m_pro =  parameter(30); 
k_m_va =  parameter(31); 
k_p =  parameter(32); 
k_pr =  parameter(33); 
pK_l_aa =  parameter(34); 
pK_l_ac =  parameter(35); 
pK_u_aa =  parameter(36); 
pK_u_ac =  parameter(37); 
p_h2o =  parameter(38); 
 
% Define algebraic equations 
 
S_nh4_i = x(7) - x(24); 
S_co2 = x(6) - x(23); 
phi = x(17) + S_nh4_i/17 - x(23)/44 - x(22)/60 - x(21)/74 - x(20)/88 - x(19)/102 - x(18); 
S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w); 
pH = -log10(S_H); 
p_ch4 = x(25)*R*T/16; 
p_co2 = x(26)*R*T/44; 
p_gas = p_ch4 + p_co2 + p_h2o ; 
q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm; 
 
% Define output (states) 
 
for i = 1:26 
    y(i) = x(i);  
end 
 
% Define output (algebraic components)
 
y(27) = S_nh4_i; 
y(28) = S_co2; 
y(29) = phi; 
y(30) = S_H; 
y(31) = pH; 
y(32) = p_ch4; 
y(33) = p_co2; 
y(34) = p_gas; 
y(35) = q_gas; 
 
% Define inhibtion functions 
 
inhibition(1) = x(7) / (x(7) + K_I_IN); 
inhibition(2) = 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2) / (S_H^(3/(pK_u_aa - pK_l_aa)) + 10^(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)); 
inhibition(3) = 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2) / (S_H^(3/(pK_u_ac - pK_l_ac)) + 10^(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)); 
inhibition(4) = K_I_nh3 / (K_I_nh3 + x(24)); 
 
% Define rate equations 
 
rate(1) = k_ch * x(9); 
rate(2) = k_pr * x(10); 
rate(3) = k_li * x(11); 
rate(4) = k_m_va * x(1) / (K_va + x(1))*x(13)*x(1)/(x(2) + x(1) + 1e-8) * inhibition(1) * inhibition(2); 
rate(5) = k_m_bu * x(2) / (K_bu + x(2))*x(14)*x(2)/(x(1) + x(2) + 1e-8) * inhibition(1) * inhibition(2); 
rate(6) = k_m_pro * x(3) / (K_pro + x(3))*x(15) * inhibition(1) * inhibition(2); 
rate(7) = k_m_ac * x(4) / (K_ac + x(4))*x(16) * inhibition(1) * inhibition(3) * inhibition(4); 
rate(8) = k_dec * x(12); 
rate(9) = k_dec * x(13); 
rate(10) = k_dec * x(14); 
rate(11) = k_dec * x(15); 
rate(12) = k_dec * x(16); 
rate(13) = k_AB_va*(x(19)  * (K_a_va + S_H) - K_a_va*x(1)); 
rate(14) = k_AB_bu*(x(20) * (K_a_bu + S_H) - K_a_bu * x(2)); 
rate(15) = k_AB_pro*(x(21) * (K_a_pro + S_H) - K_a_pro * x(3)); 
rate(16) = k_AB_ac*(x(22) * (K_a_ac + S_H) - K_a_ac * x(4)); 
rate(17) = k_AB_co2*(x(23) * (K_a_co2 + S_H) - K_a_co2 * x(6)); 
rate(18) = k_AB_IN*(x(24) * (K_a_IN + S_H) - K_a_IN * x(7)); 
rate(19) = k_La*(x(5) -16*(K_H_ch4*p_ch4)); 
rate(20) = k_La*(S_co2 - 44*(K_H_co2*p_co2)); 
