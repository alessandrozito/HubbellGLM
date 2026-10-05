# Cross-platform distribution

Collaborators on **any operating system can install with one line and no compiler**.

The repository contains clean source code. Platform-specific binaries are produced automatically by hosted build services, so users receive the correct build for their operating system, CPU architecture, and R version.

## How binaries are produced

R packages with compiled code are built for a specific combination of **operating system × CPU architecture × R major version**. This package therefore follows a source-first distribution model:

- `src/*.cpp` contains the compiled source code for the package.
- Hosted build services compile the package separately for Windows, macOS, and Linux.
- Local source installation compiles the package on the target machine through `R CMD INSTALL`.
- Generated build artifacts such as `*.o`, `*.so`, and `*.dll` are excluded from the repository and from R package builds through `.gitignore` and `.Rbuildignore`.

This keeps the repository portable and ensures that each user receives binaries built for their own platform.

## Primary distribution channel. r-universe

r-universe is the distribution channel for the package.

It builds Windows, macOS, and Linux binaries automatically from the GitHub source repository and hosts them in a package repository. Users install the package without setting up a compiler:

```r
install.packages("HubbellGLM",
                 repos = "https://alessandrozito.r-universe.dev")
```

The package is published through a public GitHub repository named `alessandrozito.r-universe.dev`, following the pattern `<your-github-username>.r-universe.dev`.

That repository contains a `packages.json` file at its root. The file `r-universe-packages.json` in this folder provides the package definition and is used as the template for `packages.json`.

Once the r-universe repository is in place, r-universe detects the package source, builds the package, publishes binaries, and creates a pkgdown documentation site.

The r-universe dashboard is available at:

```text
https://r-universe.dev
```

## Developer installation from GitHub source

Developers can install directly from GitHub when they want the latest source version.

This installation path compiles the package locally, so it uses the developer’s local build toolchain:

- macOS: Xcode Command Line Tools
- Windows: Rtools
- Linux: build-essential and R development headers

```r
# install.packages("remotes")
remotes::install_github("alessandrozito/HubbellGLM")
```

This route is mainly used for development and testing before binaries are available through r-universe.

## CRAN distribution

CRAN is the standard public distribution channel once the package is ready for submission.

After CRAN acceptance, users will be able to install the package with the default R command:

```r
install.packages("HubbellGLM")
```

The package is being prepared to pass `R CMD check` cleanly before submission. The remaining work is:

- replace the use of unexported internals in `R/fit.R`, specifically `stats:::C_Cdqrls` and `stats:::qr.lm`
- add a `tests/` directory with package tests

The GitHub Actions workflow at `.github/workflows/R-CMD-check.yaml` provides the cross-platform build check. It verifies that the source package builds successfully on Windows, macOS, and Linux before distribution through r-universe and eventual submission to CRAN.
