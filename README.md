# Code Description

The code files serve for the simulation and real data application for HTL-GMM algorithm, implemented by R package `glmnet`.

In the folder `simulation`, we show the simulation codes under the settings with $p_{\mathrm{W}} = 150$ or $p_{\mathrm{W}} = 1500$, paired with logistic or linear regressions, respectively. 

[simulate-logistic-pW=150.R](simulation/simulate-logistic-pW=150.R) 

 [simulate-logistic-pW=1500.R](simulation/simulate-logistic-pW=1500.R)  

[simulate-linear-pW=150.R](simulation/simulate-linear-pW=150.R)  

[simulate-linear-pW=1500.R](simulation/simulate-linear-pW=1500.R)  

In the other folder `real_data_application`, we show three parts including 

1. The definition of prevalent, incident cases: 
   [incident_case_definition_cancer.R](real_data_application/incident_case_definition_cancer.R)  
   [incident_case_definition.R](real_data_application/incident_case_definition.R) 

2. The extraction of risk factors: 

   [Postmenopausal_Breast_Cancer_risk_factor_analysis.R](real_data_application/Postmenopausal_Breast_Cancer_risk_factor_analysis.R)  
   [Colorectal_Cancer_risk_factor_analysis.R](real_data_application/Colorectal_Cancer_risk_factor_analysis.R)  
   [CVD_risk_factor_analysis.R](real_data_application/CVD_risk_factor_analysis.R) 
   [Stroke_risk_factor_analysis.R](real_data_application/Stroke_risk_factor_analysis.R) 
   [Asthma_risk_factor_analysis.R](real_data_application/Asthma_risk_factor_analysis.R)  

3. The prediction analysis jointly using proteomics data and risk factors: 
    [UKB_protein_data_analysis.R](real_data_application/UKB_protein_data_analysis.R) 

The prediction performance can demonstrate the superiority of HTL-GMM algorithm under the heterogeneous transfer learning setting. 