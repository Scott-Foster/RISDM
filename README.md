
## RISDM is an R package that provides functionality for estimating, diagnosing, predicting and interpreting Integrated Species Distribution Models. These models take data from disparate sources and utilises their best attributes whilst minimising their worst. See Foster et al. (2024) for more details about the package.

## Summary

At its very base, the package implements the models described in
(Fletcher Jr. et al. 2019, Miller et al. 2019, Isaac et al. 2020). These
models exploit multiple data types into a single species distribution
model. Such models coherently capture uncertainty throughout the entire
estimation and prediction process, unlike most approaches that consist
of multiple analysis stages. The structure of code for the core function
in `RISDM`, namely isdm(), evolved from that specified in Dambly et al.
(2019).

Worked examples are presented in the package’s vignette, as well as in
(Foster et al. 2024).

### Installation

The `RISDM` package can be installed using `devtools` R package.

``` r
install.packages('devtools')
library( devtools)
devtools::install_github( repo="Scott-Foster/RISDM", build_vignettes=FALSE)
```

`RISDM` bases its inference on the `INLA` package. As such, an
installation of `INLA` is required. This can be performed using the
following code.

``` r
library( devtools)
install.packages("INLA",repos=c(getOption("repos"),  
                INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
```

## Funding

This work was part of The National Vertebrate Pests and Weeds
Distribution project, which was partially funded by the Australian
Government Department of Agriculture, Fisheries and Forestry’s
Established Pest Animals and Weeds Management Pipeline Program and
Supporting Communities Manage Pests and Weeds Program.

## Code of Conduct

Please note that the ppmData project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Contributing Software

Fork the `RISDM` repository

Clone your fork using the command:

`git clone https://github.com/<username>/RISDM.git`

Contribute to your forked repository.

Create a pull request.

If your code passes the necessary checks and is documented, your changes
and/or additions will be merged in the main `RISDM` repository.

## Reporting Bugs

If you notice an issue with this repository, please report it using
[Github Issues](https://github.com/Scott-Foster/RISDM/issues). When
reporting an implementation bug, include a small example that helps to
reproduce the error (a non-working minimal example). The issue will be
addressed as quickly as possible.

## Seeking Support

If you have questions or need additional support, please open a [Github
Issues](https://github.com/Scott-Foster/RISDM/issues) or send a direct
email to <scott.foster@data61.csiro.au>.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-dam19" class="csl-entry">

Dambly, L., O’Hara, B. and Golding, N. 2019.
[<span class="nocase">oharar/IM_warbler: Integrated analysis of black-
throated blue warbler data from PA,
USA</span>](https://doi.org/10.5281/zenodo.3363936).

</div>

<div id="ref-fle19" class="csl-entry">

Fletcher Jr., R. J., Hefley, T. J., Robertson, E. P., Zuckerberg, B.,
McCleery, R. A. and Dorazio, R. M. 2019. [A practical guide for
combining data to model species
distributions](https://doi.org/10.1002/ecy.2710). - Ecology 100: e02710.

</div>

<div id="ref-fos24" class="csl-entry">

Foster, S. D., Peel, D., Hosack, G. R., Hoskins, A., Mitchell, D. J.,
Proft, K., Yang, W.-H., Uribe-Rivera, D. E. and Froese, J. G. 2024.
[‘RISDM‘: Species distribution modelling from multiple data sources in
r](https://doi.org/10.1111/ecog.06964). - Ecography n/a: e06964.

</div>

<div id="ref-isa20" class="csl-entry">

Isaac, N. J. B., Jarzyna, M. A., Keil, P., Dambly, L. I., Boersch-Supan,
P. H., Browning, E., Freeman, S. N., Golding, N., Guillera-Arroita, G.,
Henrys, P. A., Jarvis, S., Lahoz-Monfort, J., Pagel, J., Pescott, O. L.,
Schmucki, R., Simmonds, E. G. and O’Hara, R. B. 2020. [Data integration
for large-scale models of species
distributions](https://doi.org/10.1016/j.tree.2019.08.006). - Trends in
Ecology & Evolution 35: 56–67.

</div>

<div id="ref-mil19" class="csl-entry">

Miller, D. A. W., Pacifici, K., Sanderlin, J. S. and Reich, B. J. 2019.
[The recent past and promising future for data integration methods to
estimate species’
distributions](https://doi.org/10.1111/2041-210X.13110). - Methods in
Ecology and Evolution 10: 22–37.

</div>

</div>
