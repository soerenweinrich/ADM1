// Simplified and mass-based
// Anaerobic Digestion Model No. 1 (ADM1)
//
// ADM1-R2
//
// Implementation in Matlab/Simulink 
// R2019b (The MathWorks, Inc.)
//
// Level 2 S-Function
//
// Version 1.0
//
// https://github.com/soerenweinrich/ADM1
//
// Copyright (c) 2021 Sören Weinrich
// E-Mail: soeren.weinrich@dbfz.de
//
// Implementation based on the detailed C-Mex template
// provided by The MathWorks, Inc.
//
// Citation example:
//
// Weinrich, S.; Nelles, M. (2021).
// Systematic simplification of the Anaerobic Digestion Model No. 1 (ADM1) -
// Model development and stoichiometric analysis.
// Submitted to Bioresource Technology.

#define S_FUNCTION_NAME ADM1_R2_mass
#define S_FUNCTION_LEVEL 2 

#include "simstruc.h"
#include <math.h>

#define INITIAL(S)          ssGetSFcnParam(S,0)
#define PARAMETER(S)        ssGetSFcnParam(S,1)
#define SYSTEM(S)           ssGetSFcnParam(S,2)
#define NINPUTS(S)          ssGetSFcnParam(S,3)
#define PARAMETER_INPUT(S)  ssGetSFcnParam(S,4)
#define NOUTPUTS 1     

static void mdlInitializeSizes(SimStruct *S)
{
   int  i;  
   int_T input = mxGetPr(NINPUTS(S))[0];  

   if (!ssSetNumInputPorts(S, input)) return;  

   for (i = 0; i < input; i++) {  
      ssSetInputPortWidth(S, i,27);  
      ssSetInputPortDirectFeedThrough(S, i, 1);  
      ssSetInputPortRequiredContiguous(S, i, 1);  
      ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetInputPortOffsetTime(S, i, 0.0);  
      ssSetInputPortOverWritable(S, i, 0);  
   }

   if (!ssSetNumOutputPorts(S, NOUTPUTS)) return;  

   for (i = 0; i < NOUTPUTS; i++) {  
      ssSetOutputPortWidth(S, i,36);  
      ssSetOutputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetOutputPortOffsetTime(S, i, 0.0);  
   }

   ssSetNumContStates(    S, 26); 
   ssSetNumDiscStates(    S, 0); 
   ssSetNumSampleTimes(   S, 1); 
   ssSetNumSFcnParams(    S, 5); 
   ssSetNumRWork(         S, 0); 
   ssSetNumIWork(         S, 0); 
   ssSetNumPWork(         S, 0); 
   ssSetNumModes(         S, 0); 
   ssSetNumNonsampledZCs( S, 0);  
   
   ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE); 

   ssSetOptions(S, 0);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
   ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
   ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS 
#if defined(MDL_INITIALIZE_CONDITIONS) 

static void mdlInitializeConditions(SimStruct *S)
{
   real_T *x0 = ssGetContStates(S);
   int i;

   for (i = 0; i < 26; i++) {
      x0[i] = mxGetPr(INITIAL(S))[i];
   }
}

#endif 
#undef MDL_START  

static void mdlOutputs(SimStruct *S, int_T tid)
{
   real_T *y = ssGetOutputPortRealSignal(S,0);
   real_T *x = ssGetContStates(S);

   // System parameters

   double V_liq = mxGetPr(SYSTEM(S))[0];
   double V_gas = mxGetPr(SYSTEM(S))[1];
   double p_atm = mxGetPr(SYSTEM(S))[2];
   double q_in = 0;

   // Model parameters

   double K_H_ch4 =  mxGetPr(PARAMETER(S))[0];
   double K_H_co2 =  mxGetPr(PARAMETER(S))[1];
   double K_I_IN =  mxGetPr(PARAMETER(S))[2];
   double K_I_nh3 =  mxGetPr(PARAMETER(S))[3];
   double K_a_IN =  mxGetPr(PARAMETER(S))[4];
   double K_a_ac =  mxGetPr(PARAMETER(S))[5];
   double K_a_bu =  mxGetPr(PARAMETER(S))[6];
   double K_a_co2 =  mxGetPr(PARAMETER(S))[7];
   double K_a_pro =  mxGetPr(PARAMETER(S))[8];
   double K_a_va =  mxGetPr(PARAMETER(S))[9];
   double K_ac =  mxGetPr(PARAMETER(S))[10];
   double K_bu =  mxGetPr(PARAMETER(S))[11];
   double K_pro =  mxGetPr(PARAMETER(S))[12];
   double K_va =  mxGetPr(PARAMETER(S))[13];
   double K_w =  mxGetPr(PARAMETER(S))[14];
   double R =  mxGetPr(PARAMETER(S))[15];
   double T =  mxGetPr(PARAMETER(S))[16];
   double k_AB_IN =  mxGetPr(PARAMETER(S))[17];
   double k_AB_ac =  mxGetPr(PARAMETER(S))[18];
   double k_AB_bu =  mxGetPr(PARAMETER(S))[19];
   double k_AB_co2 =  mxGetPr(PARAMETER(S))[20];
   double k_AB_pro =  mxGetPr(PARAMETER(S))[21];
   double k_AB_va =  mxGetPr(PARAMETER(S))[22];
   double k_La =  mxGetPr(PARAMETER(S))[23];
   double k_ch =  mxGetPr(PARAMETER(S))[24];
   double k_dec =  mxGetPr(PARAMETER(S))[25];
   double k_li =  mxGetPr(PARAMETER(S))[26];
   double k_m_ac =  mxGetPr(PARAMETER(S))[27];
   double k_m_bu =  mxGetPr(PARAMETER(S))[28];
   double k_m_pro =  mxGetPr(PARAMETER(S))[29];
   double k_m_va =  mxGetPr(PARAMETER(S))[30];
   double k_p =  mxGetPr(PARAMETER(S))[31];
   double k_pr =  mxGetPr(PARAMETER(S))[32];
   double pK_l_aa =  mxGetPr(PARAMETER(S))[33];
   double pK_l_ac =  mxGetPr(PARAMETER(S))[34];
   double pK_u_aa =  mxGetPr(PARAMETER(S))[35];
   double pK_u_ac =  mxGetPr(PARAMETER(S))[36];
   double p_h2o =  mxGetPr(PARAMETER(S))[37];

   // Algebraic components

   double S_nh4_i;
   double S_co2;
   double phi;
   double S_H;
   double pH;
   double p_ch4;
   double p_co2;
   double p_gas;
   double q_gas;

   int i, n; 

   // Define input 

   n = ssGetNumInputPorts(S);  

   for (i = 0; i < n; i++) {
      const real_T *u = ssGetInputPortRealSignal(S,i);
      q_in = q_in + u[0]; 
   }

   // Define output (states) 

   y[0] = q_in; 

   for (i = 0; i < 26; i++){
      y[i+1] = x[i];
   }

   // Define algebraic equations

   S_nh4_i = x[6] - x[23];
   S_co2 = x[5] - x[22];
   phi = x[16] + S_nh4_i/17 - x[22]/44 - x[21]/60 - x[20]/74 - x[19]/88 - x[18]/102 - x[17];
   S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w);
   pH = -log10(S_H);
   p_ch4 = x[24]*R*T/16;
   p_co2 = x[25]*R*T/44;
   p_gas = p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define output (algebraic components)

   y[27] = S_nh4_i;
   y[28] = S_co2;
   y[29] = phi;
   y[30] = S_H;
   y[31] = pH;
   y[32] = p_ch4;
   y[33] = p_co2;
   y[34] = p_gas;
   y[35] = q_gas;
}

#define MDL_DERIVATIVES 
#undef MDL_UPDATE 
#if defined(MDL_DERIVATIVES) 

static void mdlDerivatives(SimStruct *S)
{
   real_T *dx = ssGetdX(S);
   real_T *x = ssGetContStates(S);

   // System parameters

   double V_liq = mxGetPr(SYSTEM(S))[0];
   double V_gas = mxGetPr(SYSTEM(S))[1];
   double p_atm = mxGetPr(SYSTEM(S))[2];
   double q_in;

   // Model parameters

   double K_H_ch4 =  mxGetPr(PARAMETER(S))[0];
   double K_H_co2 =  mxGetPr(PARAMETER(S))[1];
   double K_I_IN =  mxGetPr(PARAMETER(S))[2];
   double K_I_nh3 =  mxGetPr(PARAMETER(S))[3];
   double K_a_IN =  mxGetPr(PARAMETER(S))[4];
   double K_a_ac =  mxGetPr(PARAMETER(S))[5];
   double K_a_bu =  mxGetPr(PARAMETER(S))[6];
   double K_a_co2 =  mxGetPr(PARAMETER(S))[7];
   double K_a_pro =  mxGetPr(PARAMETER(S))[8];
   double K_a_va =  mxGetPr(PARAMETER(S))[9];
   double K_ac =  mxGetPr(PARAMETER(S))[10];
   double K_bu =  mxGetPr(PARAMETER(S))[11];
   double K_pro =  mxGetPr(PARAMETER(S))[12];
   double K_va =  mxGetPr(PARAMETER(S))[13];
   double K_w =  mxGetPr(PARAMETER(S))[14];
   double R =  mxGetPr(PARAMETER(S))[15];
   double T =  mxGetPr(PARAMETER(S))[16];
   double k_AB_IN =  mxGetPr(PARAMETER(S))[17];
   double k_AB_ac =  mxGetPr(PARAMETER(S))[18];
   double k_AB_bu =  mxGetPr(PARAMETER(S))[19];
   double k_AB_co2 =  mxGetPr(PARAMETER(S))[20];
   double k_AB_pro =  mxGetPr(PARAMETER(S))[21];
   double k_AB_va =  mxGetPr(PARAMETER(S))[22];
   double k_La =  mxGetPr(PARAMETER(S))[23];
   double k_ch =  mxGetPr(PARAMETER(S))[24];
   double k_dec =  mxGetPr(PARAMETER(S))[25];
   double k_li =  mxGetPr(PARAMETER(S))[26];
   double k_m_ac =  mxGetPr(PARAMETER(S))[27];
   double k_m_bu =  mxGetPr(PARAMETER(S))[28];
   double k_m_pro =  mxGetPr(PARAMETER(S))[29];
   double k_m_va =  mxGetPr(PARAMETER(S))[30];
   double k_p =  mxGetPr(PARAMETER(S))[31];
   double k_pr =  mxGetPr(PARAMETER(S))[32];
   double pK_l_aa =  mxGetPr(PARAMETER(S))[33];
   double pK_l_ac =  mxGetPr(PARAMETER(S))[34];
   double pK_u_aa =  mxGetPr(PARAMETER(S))[35];
   double pK_u_ac =  mxGetPr(PARAMETER(S))[36];
   double p_h2o =  mxGetPr(PARAMETER(S))[37];

   // Algebraic components

   double S_nh4_i;
   double S_co2;
   double phi;
   double S_H;
   double pH;
   double p_ch4;
   double p_co2;
   double p_gas;
   double q_gas;
   double inhibition[4];
   double rate[20];
   double process[26];
   double in_out[26];

   int i, j, n; 

   // Define algebraic equations

   S_nh4_i = x[6] - x[23];
   S_co2 = x[5] - x[22];
   phi = x[16] + S_nh4_i/17 - x[22]/44 - x[21]/60 - x[20]/74 - x[19]/88 - x[18]/102 - x[17];
   S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w);
   pH = -log10(S_H);
   p_ch4 = x[24]*R*T/16;
   p_co2 = x[25]*R*T/44;
   p_gas = p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define inhibtion functions

   inhibition[0] = x[6] / (x[6] + K_I_IN);
   inhibition[1] = pow(10,(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)) / (pow(S_H,(3/(pK_u_aa - pK_l_aa))) + pow(10,(-(3/(pK_u_aa - pK_l_aa))*(pK_l_aa+pK_u_aa)/2)));
   inhibition[2] = pow(10,(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)) / (pow(S_H,(3/(pK_u_ac - pK_l_ac))) + pow(10,(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)));
   inhibition[3] = K_I_nh3 / (K_I_nh3 + x[23]);

   // Define rate equations

   rate[0] = k_ch * x[8];
   rate[1] = k_pr * x[9];
   rate[2] = k_li * x[10];
   rate[3] = k_m_va * x[0] / (K_va + x[0])*x[12]*x[0]/(x[1] + x[0] + 1e-8) * inhibition[0] * inhibition[1];
   rate[4] = k_m_bu * x[1] / (K_bu + x[1])*x[13]*x[1]/(x[0] + x[1] + 1e-8) * inhibition[0] * inhibition[1];
   rate[5] = k_m_pro * x[2] / (K_pro + x[2])*x[14] * inhibition[0] * inhibition[1];
   rate[6] = k_m_ac * x[3] / (K_ac + x[3])*x[15] * inhibition[0] * inhibition[2] * inhibition[3];
   rate[7] = k_dec * x[11];
   rate[8] = k_dec * x[12];
   rate[9] = k_dec * x[13];
   rate[10] = k_dec * x[14];
   rate[11] = k_dec * x[15];
   rate[12] = k_AB_va*(x[18]  * (K_a_va + S_H) - K_a_va*x[0]);
   rate[13] = k_AB_bu*(x[19] * (K_a_bu + S_H) - K_a_bu * x[1]);
   rate[14] = k_AB_pro*(x[20] * (K_a_pro + S_H) - K_a_pro * x[2]);
   rate[15] = k_AB_ac*(x[21] * (K_a_ac + S_H) - K_a_ac * x[3]);
   rate[16] = k_AB_co2*(x[22] * (K_a_co2 + S_H) - K_a_co2 * x[5]);
   rate[17] = k_AB_IN*(x[23] * (K_a_IN + S_H) - K_a_IN * x[6]);
   rate[18] = k_La*(x[4] -16*(K_H_ch4*p_ch4));
   rate[19] = k_La*(S_co2 - 44*(K_H_co2*p_co2));

   // Define process equations

   process[0] =  0.15883 * rate[1] -10.1452 * rate[3];
   process[1] =  0.076292 * rate[0] +  0.20135 * rate[1] +  0.0092572 * rate[2] -10.9274 * rate[4];
   process[2] =  0.19032 * rate[0] +  0.046509 * rate[1] +  0.023094 * rate[2] +  6.9368 * rate[3] -14.4449 * rate[5];
   process[3] =  0.41 * rate[0] +  0.52784 * rate[1] +  1.7353 * rate[2] +  5.6494 * rate[3] +  14.0023 * rate[4] +  11.2133 * rate[5] -26.5447 * rate[6];
   process[4] =  0.047712 * rate[0] +  0.019882 * rate[1] +  0.18719 * rate[2] +  0.68644 * rate[3] +  0.87904 * rate[4] +  2.1242 * rate[5] +  6.7367 * rate[6] - rate[18];
   process[5] =  0.22553 * rate[0] +  0.18347 * rate[1] -0.64703 * rate[2] -2.6138 * rate[3] -3.0468 * rate[4] +  1.5366 * rate[5] +  18.4808 * rate[6] - rate[19];
   process[6] =  -0.013897 * rate[0] +  0.18131 * rate[1] -0.024038 * rate[2] -0.15056 * rate[3] -0.15056 * rate[4] -0.15056 * rate[5] -0.15056 * rate[6];
   process[7] =  -0.028264 * rate[0] -0.40923 * rate[1] -0.44342 * rate[2] -1.363 * rate[3] -1.7566 * rate[4] -1.2786 * rate[5] +  0.4778 * rate[6];
   process[8] =  - rate[0] +  0.18 * rate[7] +  0.18 * rate[8] +  0.18 * rate[9] +  0.18 * rate[10] +  0.18 * rate[11];
   process[9] =  - rate[1] +  0.77 * rate[7] +  0.77 * rate[8] +  0.77 * rate[9] +  0.77 * rate[10] +  0.77 * rate[11];
   process[10] =  - rate[2] +  0.05 * rate[7] +  0.05 * rate[8] +  0.05 * rate[9] +  0.05 * rate[10] +  0.05 * rate[11];
   process[11] =  0.092305 * rate[0] +  0.090036 * rate[1] +  0.15966 * rate[2] - rate[7];
   process[12] = rate[3] - rate[8];
   process[13] = rate[4] - rate[9];
   process[14] = rate[5] - rate[10];
   process[15] = rate[6] - rate[11];
   process[16] = 0;
   process[17] = 0;
   process[18] =  - rate[12];
   process[19] =  - rate[13];
   process[20] =  - rate[14];
   process[21] =  - rate[15];
   process[22] =  - rate[16];
   process[23] =  - rate[17];
   process[24] =  (V_liq/V_gas) * rate[18];
   process[25] =  (V_liq/V_gas) * rate[19];

   // Define input equations

   n = ssGetNumInputPorts(S); 

   for (i = 0; i < 26; i++) {
      in_out[i] = 0;  
   } 

   for (j = 0; j < n; j++) {
      const real_T *u = ssGetInputPortRealSignal(S,j); 
      q_in = u[0]; 

      in_out[0] = in_out[0] + q_in*(u[1] - x[0])/V_liq;
      in_out[1] = in_out[1] + q_in*(u[2] - x[1])/V_liq;
      in_out[2] = in_out[2] + q_in*(u[3] - x[2])/V_liq;
      in_out[3] = in_out[3] + q_in*(u[4] - x[3])/V_liq;
      in_out[4] = in_out[4] + q_in*(u[5] - x[4])/V_liq;
      in_out[5] = in_out[5] + q_in*(u[6] - x[5])/V_liq;
      in_out[6] = in_out[6] + q_in*(u[7] - x[6])/V_liq;
      in_out[7] = in_out[7] + q_in*(u[8] - x[7])/V_liq;
      in_out[8] = in_out[8] + q_in*(u[9] - x[8])/V_liq;
      in_out[9] = in_out[9] + q_in*(u[10] - x[9])/V_liq;
      in_out[10] = in_out[10] + q_in*(u[11] - x[10])/V_liq;
      in_out[11] = in_out[11] + q_in*(u[12] - x[11])/V_liq;
      in_out[12] = in_out[12] + q_in*(u[13] - x[12])/V_liq;
      in_out[13] = in_out[13] + q_in*(u[14] - x[13])/V_liq;
      in_out[14] = in_out[14] + q_in*(u[15] - x[14])/V_liq;
      in_out[15] = in_out[15] + q_in*(u[16] - x[15])/V_liq;
      in_out[16] = in_out[16] + q_in*(u[17] - x[16])/V_liq;
      in_out[17] = in_out[17] + q_in*(u[18] - x[17])/V_liq;
      in_out[18] = in_out[18] ;
      in_out[19] = in_out[19] ;
      in_out[20] = in_out[20] ;
      in_out[21] = in_out[21] ;
      in_out[22] = in_out[22] ;
      in_out[23] = in_out[23] ;
      in_out[24] = in_out[24] -x[24] * q_gas / V_gas;
      in_out[25] = in_out[25] -x[25] * q_gas / V_gas;
   }

   // Define differential equations

   dx[0] = in_out[0] + process[0];
   dx[1] = in_out[1] + process[1];
   dx[2] = in_out[2] + process[2];
   dx[3] = in_out[3] + process[3];
   dx[4] = in_out[4] + process[4];
   dx[5] = in_out[5] + process[5];
   dx[6] = in_out[6] + process[6];
   dx[7] = in_out[7] + process[7];
   dx[8] = in_out[8] + process[8];
   dx[9] = in_out[9] + process[9];
   dx[10] = in_out[10] + process[10];
   dx[11] = in_out[11] + process[11];
   dx[12] = in_out[12] + process[12];
   dx[13] = in_out[13] + process[13];
   dx[14] = in_out[14] + process[14];
   dx[15] = in_out[15] + process[15];
   dx[16] = in_out[16] + process[16];
   dx[17] = in_out[17] + process[17];
   dx[18] = in_out[18] + process[18];
   dx[19] = in_out[19] + process[19];
   dx[20] = in_out[20] + process[20];
   dx[21] = in_out[21] + process[21];
   dx[22] = in_out[22] + process[22];
   dx[23] = in_out[23] + process[23];
   dx[24] =  -x[24] * q_gas / V_gas + process[24];
   dx[25] =  -x[25] * q_gas / V_gas + process[25];
}

#endif 

static void mdlTerminate(SimStruct *S)
{
}

#ifdef	MATLAB_MEX_FILE    
#include "simulink.c"      
#else
#include "cg_sfun.h"      
#endif