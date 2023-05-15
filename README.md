
## RISDM is an R package that provides functionality for estimating, diagnosing, predicting and interpreting Integrated Species Distribution Models. These models take data from disparate sources and utilises their best attributes whilst minimising their worst.

## Summary

At its very base, the package implements the models described in
(Fletcher Jr. et al. 2019, Miller et al. 2019, Isaac et al. 2020). These
models exploit multiple data types into a single species distribution
model. Such models coherently capture uncertainty throughout the entire
estimation and prediction process, unlike most approaches that consist
of multiple analysis stages.

Worked examples are presented in the package’s vignette, as well as in
(Foster et al. 2023).

### Installation

The `RISDM` package can be installed using `devtools` R package.

``` r
install.packages('devtools')
devtools::install_github( repo="Scott-Foster/RISDM", build_vignettes=FALSE)
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

<div id="ref-fle19" class="csl-entry">

Fletcher Jr., R. J. et al. 2019. [A practical guide for combining data
to model species distributions](https://doi.org/10.1002/ecy.2710). -
Ecology 100: e02710.

</div>

<div id="ref-fos23" class="csl-entry">

Foster, S. D. et al. 2023. RISDM: Species distribution modelling from
multiple data sources in r. - Manuscript in preparation in press.

</div>

<div id="ref-isa20" class="csl-entry">

Isaac, N. J. B. et al. 2020. [Data integration for large-scale models of
species distributions](https://doi.org/10.1016/j.tree.2019.08.006). -
Trends in Ecology & Evolution 35: 56–67.

</div>

<div id="ref-mil19" class="csl-entry">

Miller, D. A. W. et al. 2019. [The recent past and promising future for
data integration methods to estimate species’
distributions](https://doi.org/10.1111/2041-210X.13110). - Methods in
Ecology and Evolution 10: 22–37.

</div>

</div>
