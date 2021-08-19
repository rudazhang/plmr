# Probabilistic Learning on Manifold

R pacakge `plmr` ("plumber") implements methods for Probabilistic Learning on Manifolds in R.

## Installation

Install the developing verion via `devtools`.

``` R
if (!("devtools" %in% installed.packages()[,"Package"])) {
    install.packages(devtools)
}
devtools::install_github("rudazhang/plmr")
```

## Manifest

Source code for package exported objects, `/R`:

- `MParzen.R`, manifold Parzen window [@Vincent2002] density estimation and sampling.
- `SCMS.R`, subspace-constrained mean shift [@Ozertem2011] and variants for ridge estimation.
- `NormalBundleBootstrap.R`, normal-bundle bootstrap (NBB)
  for inference in normal spaces of the density ridge,
  and data augmentation to reduce overfitting.
- `DiffusionSampling.R`, diffusion sampling on density ridge.
- `DiffusionPropagate.R`, heat diffusion on a manifold point cloud for probability propagation.
- `ParameterizedManifolds.R`, some parameterized manifolds and sampling.

Scripts for experiments and figures, `/inst/script`:

- `nbb/`, scripts for NBB.
    - `1-fig-scms-circle.R`, script for SCMS with circle data;
    - `1-fig-scms-parabola.R`, script for SCMS with parabola data;
    - `1-fig-exp-circle-ridge.R`, script for circle;
    - `1-exp-wheel.R`, script for wheel functional data,
    includes a rough implementation of *NBB with smooth frame*;
    - `1-fig-exp-wheel.R`, plot wheel results;
- `misc/`, miscellaneous scripts:
    - `exp-SCMS.R`, ridge estimation by SCMS.
    - `exp-diffusion-propagate.R`, heat diffusion for probability propagation.
    - `exp-earthquake.R`, earthquake data.
    - `exp-misc.R`, miscellaneous toy examples.
    - `exp-wasserstein.R`, optimal transport and Wasserstein distances.

Documents, `/vignettes`:

- `earthquake.md`, earthquake data description.

## References

- Data-driven probability concentration and sampling on manifold.
C Soize, R Ghanem - Journal of Computational Physics, 2016.
https://doi.org/10.1016/j.jcp.2016.05.044
- Normal-bundle Bootstrap.
R Zhang, R Ghanem - arXiv, 2020.
https://arxiv.org/abs/2007.13869

BibTeX citation:
``` bibtex
@Article{ZhangRD2021nbb,
  author        = {Zhang, Ruda and Ghanem, Roger},
  title         = {Normal-Bundle Bootstrap},
  journal       = {SIAM Journal on Mathematics of Data Science},
  year          = {2021},
  volume        = {3},
  number        = {2},
  pages         = {573--592},
  doi           = {10.1137/20m1356002},
}
```
