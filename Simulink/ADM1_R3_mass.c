// Simplified and mass-based
// Anaerobic Digestion Model No. 1 (ADM1)
//
// ADM1-R3
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

#define S_FUNCTION_NAME ADM1_R3_mass
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
      ssSetInputPortWidth(S, i,18);  
      ssSetInputPortDirectFeedThrough(S, i, 1);  
      ssSetInputPortRequiredContiguous(S, i, 1);  
      ssSetInputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetInputPortOffsetTime(S, i, 0.0);  
      ssSetInputPortOverWritable(S, i, 0);  
   }

   if (!ssSetNumOutputPorts(S, NOUTPUTS)) return;  

   for (i = 0; i < NOUTPUTS; i++) {  
      ssSetOutputPortWidth(S, i,27);  
      ssSetOutputPortSampleTime(S, i, INHERITED_SAMPLE_TIME);  
      ssSetOutputPortOffsetTime(S, i, 0.0);  
   }

   ssSetNumContStates(    S, 17); 
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

   for (i = 0; i < 17; i++) {
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
   double K_a_co2 =  mxGetPr(PARAMETER(S))[6];
   double K_ac =  mxGetPr(PARAMETER(S))[7];
   double K_w =  mxGetPr(PARAMETER(S))[8];
   double R =  mxGetPr(PARAMETER(S))[9];
   double T =  mxGetPr(PARAMETER(S))[10];
   double k_AB_IN =  mxGetPr(PARAMETER(S))[11];
   double k_AB_ac =  mxGetPr(PARAMETER(S))[12];
   double k_AB_co2 =  mxGetPr(PARAMETER(S))[13];
   double k_La =  mxGetPr(PARAMETER(S))[14];
   double k_ch =  mxGetPr(PARAMETER(S))[15];
   double k_dec =  mxGetPr(PARAMETER(S))[16];
   double k_li =  mxGetPr(PARAMETER(S))[17];
   double k_m_ac =  mxGetPr(PARAMETER(S))[18];
   double k_p =  mxGetPr(PARAMETER(S))[19];
   double k_pr =  mxGetPr(PARAMETER(S))[20];
   double pK_l_ac =  mxGetPr(PARAMETER(S))[21];
   double pK_u_ac =  mxGetPr(PARAMETER(S))[22];
   double p_h2o =  mxGetPr(PARAMETER(S))[23];

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

   for (i = 0; i < 17; i++){
      y[i+1] = x[i];
   }

   // Define algebraic equations

   S_nh4_i = x[3] - x[14];
   S_co2 = x[2] - x[13];
   phi = x[10] + S_nh4_i/17 - x[13]/44 - x[12]/60  - x[11];
   S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w);
   pH = -log10(S_H);
   p_ch4 = x[15]*R*T/16;
   p_co2 = x[16]*R*T/44;
   p_gas = p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define output (algebraic components)

   y[18] = S_nh4_i;
   y[19] = S_co2;
   y[20] = phi;
   y[21] = S_H;
   y[22] = pH;
   y[23] = p_ch4;
   y[24] = p_co2;
   y[25] = p_gas;
   y[26] = q_gas;
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
   double K_a_co2 =  mxGetPr(PARAMETER(S))[6];
   double K_ac =  mxGetPr(PARAMETER(S))[7];
   double K_w =  mxGetPr(PARAMETER(S))[8];
   double R =  mxGetPr(PARAMETER(S))[9];
   double T =  mxGetPr(PARAMETER(S))[10];
   double k_AB_IN =  mxGetPr(PARAMETER(S))[11];
   double k_AB_ac =  mxGetPr(PARAMETER(S))[12];
   double k_AB_co2 =  mxGetPr(PARAMETER(S))[13];
   double k_La =  mxGetPr(PARAMETER(S))[14];
   double k_ch =  mxGetPr(PARAMETER(S))[15];
   double k_dec =  mxGetPr(PARAMETER(S))[16];
   double k_li =  mxGetPr(PARAMETER(S))[17];
   double k_m_ac =  mxGetPr(PARAMETER(S))[18];
   double k_p =  mxGetPr(PARAMETER(S))[19];
   double k_pr =  mxGetPr(PARAMETER(S))[20];
   double pK_l_ac =  mxGetPr(PARAMETER(S))[21];
   double pK_u_ac =  mxGetPr(PARAMETER(S))[22];
   double p_h2o =  mxGetPr(PARAMETER(S))[23];

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
   double inhibition[3];
   double rate[11];
   double process[17];
   double in_out[17];

   int i, j, n; 

   // Define algebraic equations

   S_nh4_i = x[3] - x[14];
   S_co2 = x[2] - x[13];
   phi = x[10] + S_nh4_i/17 - x[13]/44 - x[12]/60  - x[11];
   S_H = -phi*0.5 + 0.5*sqrt(phi*phi+4*K_w);
   pH = -log10(S_H);
   p_ch4 = x[15]*R*T/16;
   p_co2 = x[16]*R*T/44;
   p_gas = p_ch4 + p_co2 + p_h2o ;
   q_gas = k_p * ( p_gas - p_atm) *p_gas/p_atm;

   // Define inhibtion functions

   inhibition[0] = x[3] / (x[3] + K_I_IN);
   inhibition[1] = pow(10,(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)) / (pow(S_H,(3/(pK_u_ac - pK_l_ac))) + pow(10,(-(3/(pK_u_ac - pK_l_ac))*(pK_l_ac+pK_u_ac)/2)));
   inhibition[2] = K_I_nh3 / (K_I_nh3 + x[14]);

   // Define rate equations

   rate[0] = k_ch * x[5];
   rate[1] = k_pr * x[6];
   rate[2] = k_li * x[7];
   rate[3] = k_m_ac * x[0] / (K_ac + x[0])*x[9] * inhibition[0] * inhibition[1] * inhibition[2];
   rate[4] = k_dec * x[8];
   rate[5] = k_dec * x[9];
   rate[6] = k_AB_ac*(x[12] * (K_a_ac + S_H) - K_a_ac * x[0]);
   rate[7] = k_AB_co2*(x[13] * (K_a_co2 + S_H) - K_a_co2 * x[2]);
   rate[8] = k_AB_IN*(x[14] * (K_a_IN + S_H) - K_a_IN * x[3]);
   rate[9] = k_La*(x[1] -16*(K_H_ch4*p_ch4));
   rate[10] = k_La*(S_co2 - 44*(K_H_co2*p_co2));

   // Define process equations

   process[0] =  0.6555 * rate[0] +  0.9947 * rate[1] +  1.7651 * rate[2] -26.5447 * rate[3];
   process[1] =  0.081837 * rate[0] +  0.069636 * rate[1] +  0.19133 * rate[2] +  6.7367 * rate[3] - rate[9];
   process[2] =  0.2245 * rate[0] +  0.10291 * rate[1] -0.64716 * rate[2] +  18.4808 * rate[3] - rate[10];
   process[3] =  -0.016932 * rate[0] +  0.17456 * rate[1] -0.024406 * rate[2] -0.15056 * rate[3];
   process[4] =  -0.057375 * rate[0] -0.47666 * rate[1] -0.44695 * rate[2] +  0.4778 * rate[3];
   process[5] =  - rate[0] +  0.18 * rate[4] +  0.18 * rate[5];
   process[6] =  - rate[1] +  0.77 * rate[4] +  0.77 * rate[5];
   process[7] =  - rate[2] +  0.05 * rate[4] +  0.05 * rate[5];
   process[8] =  0.11246 * rate[0] +  0.13486 * rate[1] +  0.1621 * rate[2] - rate[4];
   process[9] = rate[3] - rate[5];
   process[10] = 0;
   process[11] = 0;
   process[12] =  - rate[6];
   process[13] =  - rate[7];
   process[14] =  - rate[8];
   process[15] =  (V_liq/V_gas) * rate[9];
   process[16] =  (V_liq/V_gas) * rate[10];

   // Define input equations

   n = ssGetNumInputPorts(S); 

   for (i = 0; i < 17; i++) {
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
      in_out[12] = in_out[12];
      in_out[13] = in_out[13] ;
      in_out[14] = in_out[14] ;
      in_out[15] = in_out[15] - x[15] * q_gas / V_gas;
      in_out[16] = in_out[16] - x[16] * q_gas / V_gas;
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
   dx[15] =  - x[15] * q_gas / V_gas + process[15];
   dx[16] =  - x[16] * q_gas / V_gas + process[16];
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