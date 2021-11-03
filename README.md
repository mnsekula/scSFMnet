# scSFMnet

**scSFMnet** is an R package that interfaces with [`rstan`](https://mc-stan.org/users/interfaces/rstan). The `rstan` package should be installed and working properly before installing **scSFMnet**. See [Rstan Getting Started](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for more details on the `rstan` installation process.


## Installation
For most users, this package can be installed from GitHub with:

```{r, warning=FALSE}
# Install from GitHub using devtools
require(devtools)
devtools::install_github("mnsekula/scSFMnet")
```

Note: There may be some `rstan` compiler warnings related to certain dependent packages (e.g., RcppEigen, StanHeaders, etc.). According to the [Brief Guide to Stan’s Warnings](https://mc-stan.org/misc/warnings.html), these warnings can safely be ignored as long as the **scSFMnet** package is installed successfully. The code below should output `TRUE` for a successful installation:

```{r, warning=FALSE}
# Check if scSFMnet was successfully installed
is.element("scSFMnet", installed.packages()[,1])
```


## Getting Started
The package vignette demonstrates how to use the **scSFMnet** package to perform a differential network analysis. This vignette can be viewed online [here](http://htmlpreview.github.io/?https://github.com/mnsekula/scSFMnet/blob/master/scSFMnet.html).

