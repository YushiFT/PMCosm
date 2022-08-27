PMCosm: Probabilistic model for microbial-cosmos
=======

This `PMCosm` package implements statistical inferences for microbial ecology analysis. This package provides methods for 

* categorizing the admixed structure of microbial communities,
* detecting influential microbes in both engineered and natural environments,
* evaluating sub-communities' preferences of life strategies,
* identifying causal relationship between microbes' metabolism functions and fluctuations of surrounding environments,
* predicting the community composition and structure for a given local environment. 

Installation and documentation
------------------------------

To install, open R and type:

```R
install.packages("devtools")
devtools::install_github("YushiFT/PMCosm")
```

Package overview
----------------

### Functions 
* `calc_mle_trio`: Calculate MLE estimates for universal parameter trio for inferences about the admixed community structure. 
* `classify_taxa`: Traditional relative-abundance-based methods to classify microbial taxa into rare or abundant biospheres.
* `is_dispersion`: Hypothesis test for the existence of dispersion.
* `is_overdispersion`: Hypothesis test for the existence of over-dispersion.
* `is_zero_infla`: Check the existence of inflated zeros. 
* `pca_simple`: A simplified principle component analysis for microbial data.
* `plot_trio`: Graphic display of trio estimates to observe the Two-Wing pattern.
* `refine_boundary`: Define a model-driven testable boundary between dispersal vangguards and dispersal laggards.

### Quick start guide

#### Deriving MLE estimates for the admixed structure of environmental microbial communities

```R
library(PMCosm)
# import sample data from microbial communities in hangzhou bay
data(hzmicrobe)
param_trio_bay <- calc_mle_trio(mic_bay, n_sample=10, replicates=3)
param_trio_era <- calc_mle_trio(mic_era, n_sample=12, replicates=3)
```
#### Visualize the Two-Wing admixed structure

```R
library(ggplot2)
library(latex2exp)
# for bay
plot_trio(param_trio_bay, point_size=0.8, a=0.6)
# zoom in for the majority 
plot_trio(param_trio_bay, point_size=0.8, a=0.6, zoom_in=TRUE)
# for era 
plot_trio(param_trio_era, point_size=0.8, a=0.6)
```






