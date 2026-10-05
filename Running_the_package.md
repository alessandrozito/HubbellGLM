# Installing & running HubbellGLM (macOS / Windows / Linux)

This fork ships **clean source with no committed compiled objects** — the
package is compiled from `src/*.cpp` on your machine at install time, so it
builds natively on any OS with a C/C++ toolchain. No platform-specific binaries
to go stale, and no `--preclean` workaround needed.

> **Verified vs. reasoned:** the macOS path was reproduced on this machine
> (Apple Silicon arm64, R 4.5.1 and R 4.6.0) — clean install, correct binary,
> correct results. The Windows and Linux paths follow the same source build but
> were **not tested here**.

## 1. Prerequisites

A C/C++ toolchain (one-time), then the R build/install packages.

| OS | Toolchain |
|----|-----------|
| macOS | `xcode-select --install` |
| Windows | Install **Rtools** matching your R version: <https://cran.r-project.org/bin/windows/Rtools/> |
| Linux (Debian/Ubuntu) | `sudo apt-get install r-base-dev build-essential` (Fedora: `sudo dnf install R-core-devel`) |

Then, in R (a fresh R library may have none of these):

```r
install.packages(c("Rcpp", "RcppArmadillo", "remotes"))
```

## 2. Install (identical on all platforms)

```r
remotes::install_github("DNAdiversity/HubbellGLM")
```

or, from a local clone of this fork:

```bash
R CMD INSTALL /path/to/HubbellGLM-fork
```

(No `--preclean` needed — that flag was only required for copies that still
ship committed Linux objects, such as the upstream `alessandrozito/HubbellGLM`.)

## 3. Verify

Checks 2 and 3 are pure R, identical on every OS, with the same expected output.
Only Check 1 (the binary's format) differs per OS.

**1. The compiled library was built for your platform.** A successful
`library()` in Check 2 already implies this, but to confirm explicitly:

```r
library(HubbellGLM)
cat(getLoadedDLLs()[["HubbellGLM"]][["path"]], "\n")
```

then inspect that file:

| OS | command | expected |
|----|---------|----------|
| macOS | `file <path>` | `Mach-O ... arm64` (or `x86_64` on Intel) |
| Linux | `file <path>` | `ELF ... x86-64` |
| Windows | it is a `.dll`; `file <path>` in Git Bash | `PE32+ ... x86-64` — successful `library()` is sufficient proof |

**2. It loads and fits a model:**

```r
data(BarroColorado)
fit <- HubbellGLM(cbind(n, y) ~ EnvHet, family = hubbell(sigma = 0),
                  data = BarroColorado[1:40, ])
cat("converged =", fit$converged, " intercept =", round(coef(fit)[1], 2), "\n")
# expected: converged = TRUE  intercept = 3.55
```

**3. It reproduces the documented result** (correctness, not just "it ran").
The README states `estimate_sigma(...)` on `BarroColorado[1:40, ]` returns
`-0.7852143`:

```r
s <- estimate_sigma("cbind(n, y) ~ EnvHet + Habitat + Stream",
                    data = BarroColorado[1:40, ], verbose = FALSE)
cat("estimate_sigma =", round(s, 4), "\n")
# expected: -0.7852   (README: -0.7852143; last-digit drift is optimizer tolerance)
```

## Included test data

No external files needed — the package ships its data:

- **`BarroColorado`** — 50-row example dataset, `data(BarroColorado)` (used above).
- **`data_regGMTP.tsv`** — the full 2,415-row GMTP dataset behind the paper, at
  `system.file("data", "data_regGMTP.tsv", package = "HubbellGLM")`; read with
  `read.delim()`. Heavier; only needed to reproduce the paper's models.

(There is no installed test suite: `tutorials/test_package.R` is excluded from
the build, and no vignette is installed.)
