# Wavelet-Nonparametric-Quantile-Causality
To run the Quantile Causality Test, ensure that all necessary R packages (quantreg, waveslim, tidyverse, readxl, QCSIS, and wranglR) are installed and loaded. 
The script begins by loading the dataset (DATA.xlsx), where the dependent (DEP) and independent (IND) variables are defined. 
Using Wavelet Multi-Resolution Analysis (MRA), the data is decomposed into different frequency components, capturing short-, medium-, and long-term effects. 
The quantile causality test is then applied using the lrq.causality.test function, which estimates causality at different quantiles (0.05 to 0.95) and across time scales, allowing for the detection of nonlinear, asymmetric, and scale-dependent causal relationships. 
Results are extracted and visualized using level plots, where color gradients represent the strength and statistical significance of causality across different quantiles and periods
