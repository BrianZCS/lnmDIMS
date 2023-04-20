# lnmDIMS

To compare potential longitudinal sampling strategies, it is valuable and cost-efficient to simulate hypothetical experimental data and imagine their subsequent analysis results. Packages supporting such power analysis have been published, drawing from various models of longitudinal microbiome data. However, though it would be useful to calibrate these simulators based on available pilot data, there are no existing packages that support this. In addition, specifying model parameters can be difficult, especially in the absence of prior knowledge. However, most of the simulators do not contain functions that allow users to generate data with estimated model parameters obtained from data fitting. To improve the features, we propose a package that applies a logistic multinomial model to simulate longitudinal microbiome count data. Based on the model, we extended our simulations to two practical experiment designs in microbiome studies, one investigating the effect of certain perturbations on the microbiome, and the other modeling the transitions of dominant microbiome clusters. Our simulator allows users to generate data from scratch with default parameters or to estimate reasonable defaults using pilot data. To help with experimental design, we create functions for power analysis and provide visualizations to help determine the optimal experimental conditions.

# Installation
```
#install.packages("devtools")

devtools::install_github("BrianZCS/lnmDIMS")
```

# Website

Welcome to visit our [lnmDIMS website](https://brianzcs.github.io/lnmDIMS/).

# License

Code open source for anyone to use as it's under [MIT License](https://opensource.org/license/mit/).