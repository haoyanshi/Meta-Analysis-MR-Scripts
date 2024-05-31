*Input data录入数据
metan OA组样本量 OA组浓度mean OA组浓度sd 对照组样本量 对照组浓度mean 对照组浓度sd, label(namevar=作者, yearvar=年份) random cohen
*Adjust the parameters according to the results of heterogeneity analysis. If the heterogeneity is small, run the following line of code.
metan OA组样本量 OA组浓度mean OA组浓度sd 对照组样本量 对照组浓度mean 对照组浓度sd, label(namevar=作者, yearvar=年份) fixed cohen
*Draw a funnel chart
metafunnel _ES _seES
*Perform a bias test and further check whether the funnel plot is statistically symmetrical.
metabias6 _ES _seES
*Sensitivity analysis was performed, and the choice of parameter random or fixed remained consistent with the main analysis of the meta-analysis.
metainf _ES _seES, random id(作者)//It should be noted that there are often incompatibility issues when performing sensitivity analysis with STATA 17.0 version. You can switch to version 15.0 for this analysis.