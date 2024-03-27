# seqimpute 2.0.0

## Breaking changes

* `seqimpute()` now returns an object of class `seqimp`. In particular, the 
`include` and `mice.return` arguments are no longer relevant and have been
removed.

* The `OD` argument has been renamed to `data`. The argument `OD` itself is
deprecated.

* The `CO` argument has been renamed to `covariates`. The argument `CO` itself
is deprecated.

* The `COt` argument has been renamed to `time. covariates`. The argument 
`COt` itself is deprecated.

* The `mi` argument has been renamed to `m`. The argument `mi` itself is
deprecated.

* The dataset provided in the package used to be divided into three parts: 
the trajectories (`OD`), the covariates (`CO`),
and the time-varying covariates (`COt`). They now appear as a 
single dataset, called `gameadd`.

* `seqimpute()` no longer implements linear and ordinal regressions. 

* The default argument of `m` has been set to 5.


## New features

* `seqimpute()` implements the MICT-timing imputation algorithm. The argument
`timing` indicates whether to use this algorithm or the MICT algorithm, and 
`frame.radius` specifies the radius of the time frame.

* The user can now pass a dataset to the `seqimpute()` function and specify 
which columns correspond to the trajectories with the `var` argument, to the 
covariates with the `covariates` argument, and the time-varying covariates with 
the `time.covariates` argument.

* A vignette has been added.

* New `seqmissfplot()` plot function, which plots the most frequent patterns of 
missing data.

* New `seqmissIplot()` plot function, which plots all patterns of missing data.

* New `seqmissimplic()` function for identifying and visualizing the states 
that best characterize sequences with missing data.

* New `fromseqimp()` function, which converts a `seqimp` object into a 
specified format.

* New `addcluster()` function,  which adds a clustering result to a `seqimp` 
object

* New `seqaddNA()` function to simulate missing data.

* New `seqcomplete()` function, which extracts all trajectories without missing 
data.

* New `seqwithmiss()` function,  which extracts all the trajectories with at 
least one missing value.

* `seqimpute()` now returns an object of class `seqimp`. A print, summary, and 
plot functions have been added for this object type.

* `seqTrans()` and `seqQuickLook()` now accept objects of class `stslist` 
built with the `TraMineR` package.

* A `...` argument has been added to `seqimpute()` to pass named arguments to 
the imputation functions.



## Bug fixes

* Fixes issues in `seqimpute()` related to the multinomial model when there is 
only one state in the dependent variable. 

* Fixes issues in `seqimpute()` related to random forest when one state does not appear in the 
dependent variable.

* Fixes bug in `seqimpute()` when a single covariate is specified.

* Fixes bug in `seqimpute()` related to long internal gaps.

* Fixes bug in `seqQuickLook()` induced by datasets without missing data.