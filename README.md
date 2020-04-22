# Sampling on Manifold

R pacakge `plmr` ("plumber"), Probabilistic Learning on Manifolds in R.

## Manifest

Source code for package exported objects, `/R`:

- `MParzen.R`, manifold Parzen window [@Vincent2002] density estimation and sampling.
- `SCMS.R`, subspace-constrained mean shift [@Ozertem2011] and variants for ridge estimation.
- `NBB.R`, normal-bundle bootstrap for inference in normal spaces of the density ridge,
  and data augmentation to reduce overfitting.
- `DiffusionSampling.R`, diffusion sampling on density ridge.
- `DiffusionPropagate.R`, heat diffusion on a manifold point cloud for probability propagation.
- `ParameterizedManifolds.R`, some parameterized manifolds and sampling.
- (dismissed) `DMapProjection.R`, Gaussian KDE projected on truncated diffusion maps basis.

Scripts for experiments and figures, `/script`:

- `exp-ridge.R`, ridge estimation by SCMS.
- `exp-diffusion-propagate.R`, heat diffusion for probability propagation.
- `exp-earthquake.R`, earthquake data.
- `exp-misc.R`, miscellaneous toy examples.
- `transport.R`, optimal transport and Wasserstein distances.
- (dismissed) `exp-dmap-project.R`, example script.

Documents, `/vignettes`:

- `earthquake.md`, earthquake data description.
