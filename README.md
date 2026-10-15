<p align="center">
  <img src="darkin-class-logo.png" alt="darkin-class logo" width="320">
</p>

# darkin-class

**Dark energy–dark matter interactions in CLASS.**

Authors: Elsa M. Teixeira and William Giarè

darkin-class extends the [CLASS Einstein–Boltzmann solver](https://github.com/lesgourg/class_public) to study interacting dark-sector models at the background and linear perturbation levels. It provides scalar-field couplings to fluid dark matter alongside the familiar CLASS workflow for cosmological observables.

[Getting started](#getting-started) · [Coupling guide](doc/input/COUPLINGS.md) · [Examples](#examples-and-documentation) · [Publications](#publications-and-citation) · [Support](#support)

## Features

- **Several interaction sectors:** conformal, disformal, mixed conformal–disformal, entropy, and momentum couplings. The conformal–disformal sector uses a general coupling function `Q` and its perturbation; entropy and momentum interactions have separate source terms.
- **Scalar-field potentials:** exponential and double-exponential forms, with named input parameters or ordered parameter lists. Potential functions can be extended in the source code.
- **Cosmological predictions:** background evolution, CMB angular power spectra, and linear matter power spectra through CLASS.
- **Command-line and Python interfaces:** run `.ini` files or use the `classy` wrapper in scripts and notebooks.
- **Example configurations:** coupling-specific Cobaya YAML files and notebooks to explore the implementation.

The scalar-field-coupled dark matter component is `qcdm`, configured with `omega_qcdm` or `Omega_qcdm`. It is distinct from CLASS's built-in `idm` species. See the [coupling guide](doc/input/COUPLINGS.md) for conventions, initial conditions, and coupling combinations. That guide is currently marked as awaiting review; check the implementation when using it for research.

## Getting started

### Build the command-line solver

You need Git, GNU Make, and C/C++ compilers. Compiler names and flags can be adjusted in the [Makefile](Makefile).

```sh
git clone https://github.com/elsateixeira/darkin-class.git
cd darkin-class
make -j class
mkdir -p output
./class default.ini
```

The supplied `default.ini` requests CMB and matter power spectra. Output files are written under `output/` unless you set a different `root` in your input file. This first run checks the standard CLASS workflow; configure a coupling separately to explore interacting models.

For compiler troubleshooting, including macOS compiler setup, consult the [CLASS installation guide](https://github.com/lesgourg/class_public/wiki/Installation). Use this repository's Makefile as the reference for its current build flags.

### Set up an interacting model

Start with the [coupling guide](doc/input/COUPLINGS.md), which describes:

- the required scalar-field and coupled-dark-matter densities;
- the coupling selector and named parameters for each sector;
- potential parameters and scalar-field initial conditions;
- parameter-list ordering and combined couplings.

Create a separate `.ini` file for your model, then run:

```sh
./class my_model.ini
```

Keep [explanatory.ini](explanatory.ini) as a reference for general CLASS inputs. A precision file can be passed as a second argument, for example `./class my_model.ini cl_permille.pre`.

### Use Python

Install the wrapper in a dedicated environment so it is easy to identify which CLASS implementation you are using:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install setuptools wheel numpy cython
make -j all PYTHON=python
```

The wrapper is imported as `from classy import Class`. Start with the [warmup notebook](notebooks/warmup.ipynb) for the Python interface and the [coupling example](notebooks/Coupling_example.ipynb) for dark-sector inputs. Running notebooks requires Jupyter in the same environment.

## Examples and documentation

| Resource | Use it for |
| --- | --- |
| [Coupling guide](doc/input/COUPLINGS.md) | Interaction parameters, conventions, and input templates; review status noted above |
| [Coupling notebook](notebooks/Coupling_example.ipynb) | Exploring dark-sector inputs in Python |
| [Cobaya configurations](configs/couplings/) | Coupling-specific inference setups; adapt paths, priors, and likelihood requirements to your analysis |
| [Python notebooks](notebooks/) | Examples of backgrounds, spectra, and other CLASS calculations |
| [explanatory.ini](explanatory.ini) | Detailed general CLASS input reference |
| [default.ini](default.ini) | A shorter starting point for command-line runs |
| [Source code](source/) | Background equations, perturbation evolution, and observables |

The inherited plotting utility supports CLASS output files. Run `python CPU.py -h` for usage; a MATLAB counterpart is available in [plot_CLASS_output.m](plot_CLASS_output.m).

## Publications and citation

This repository is associated with work on interacting dark-sector models, including:

- [Dissecting kinetically coupled quintessence: phenomenology and observational tests](https://arxiv.org/abs/2207.13682).
- [Scalar field dark matter and dark energy: a hybrid model for the dark sector](https://arxiv.org/abs/2211.13653).
- [Forecasts on interacting dark energy with standard sirens](https://arxiv.org/abs/2309.06544).
- [Alleviating cosmological tensions with a hybrid dark sector](https://arxiv.org/abs/2412.14139).

When using darkin-class in a publication, cite the work relevant to your model and the underlying CLASS solver, including [CLASS II: Approximation schemes](https://arxiv.org/abs/1104.2933). Record the commit and configuration used for your analysis: the current development version and example files should not be assumed to reproduce every published result without the corresponding historical setup.

## Credits

darkin-class builds on CLASS, developed by Julien Lesgourgues, Thomas Tram, Nils Schöneberg, and many contributors. See the [CLASS project](https://github.com/lesgourg/class_public) for the upstream code, contributors, and documentation.

## Acknowledgements

We thank Vitor Fonseca and Jurgen Mifsud for their help with darkin-class.

## Support

For questions, bug reports, or suggestions specific to darkin-class, [open an issue in this repository](https://github.com/elsateixeira/darkin-class/issues). Include your commit, platform, input file, and full error message so the problem can be reproduced.
