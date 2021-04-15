# Simplified ADM1
Matlab and Simulink implementations of different simplifications of a mass-based Anaerobic Digestion Model No. 1 (ADM1). 

## Model Description
Various model structures have been developed for practical application on industrial anaerobic digestion plants. Individual model variations greatly differ in their number of implemented process phases, characteristic components and required parameters. Complex models like the mass-based ADM1, ADM1-R1 or ADM1-R2 depict specific degradation pathways and intermediates during acido- and acetogenesis in detail, whereas simplified models (such as the ADM1-R4) combine nutrient degradation and biogas formation based on single first-order sum reactions.

![Simplified ADM1](https://github.com/soerenweinrich/ADM1/blob/main/Documents/Weinrich_2021_Simplified_ADM1.png)

In regard to available measurements on anaerobic digestion plants simplified model structures show clear advantages for practical application, due to the small number of model parameters required and suitable system characteristics. Thus, the ADM1-R4 and ADM1-R3 can be applied as robust estimators to predict gas production rates for plant design, process monitoring and control during full-scale plant operation. Complex model variants, such as the ADM1, ADM1-R1 and ADM1-R2 enable a precise description of characteristic intermediates and allow for a detailed state analysis based on microbial growth conditions (including relevant inhibitors).

## Literature
Detailed description of development, application and evaluation of individual model structures is provided in the following publications: 
Weinrich, S.; Nelles, M., (2021): Systematic model reduction of the Anaerobic Digestion Model No. 1 (ADM1) – Model development and stoichiometric analysis. Bioresource Technology. In press. https://doi.org/10.1016/j.biortech.2021.125124.
Weinrich, S.; Mauky, E.; Schmidt, T.; Krebs, C.; Liebetrau, J.; Nelles, M., (2021): Systematic model reduction of the Anaerobic Digestion Model No. 1 (ADM1) – Laboratory experiments and model application. Bioresource Technology. In press. https://doi.org/10.1016/j.biortech.2021.125104.
