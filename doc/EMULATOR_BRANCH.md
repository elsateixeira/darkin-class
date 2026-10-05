# IDE emulator preparation branch

This branch builds on the small bug fixes in master. It additionally changes:

- Explicit exponential-potential V0 shooting: positive log variables, bracketing,
  damped joint density shooting and continuation recovery. Other target combinations
  and legacy scf_parameters tuning retain the original route.
- Momentum perturbations: reject sampled background trajectories with a nonfinite
  or nonpositive Euler denominator before emitting invalid spectra.
- Accurate lensing: batch angular quadrature with lensing_mu_chunk_size (default
  256; zero retains full tables), and integrate the correction to the unlensed
  correlation to reduce cancellation at high multipoles.

These changes are isolated here for further scientific review. Selected-point
regressions pass; the full proposed emulator parameter domains remain unvalidated.
Some boundary points still fail shooting. A failed numerical search is not a
physical exclusion, and finite outputs alone do not establish stability or accuracy.

## Build and check

```sh
make -j4 class libclass.a
python test/test_ide_numerics.py --executable ./class
cd python
python setup.py build_ext --inplace --force
cd ..
OMP_NUM_THREADS=1 python test/test_runtime_regressions.py
```

Use the Python environment intended for generation (NumPy, Cython and build
dependencies required). Confirm classy.__file__ points to this build.
The IDE regression covers seed independence, density closure, continuation,
lensing batch agreement and the known pole diagnostic. Runtime regressions cover
Python data discovery, CLI agreement and a missing HyRec path.

Use the separately supplied Ivan dataset YAMLs and their numerical profiles.
Do not apply the older handover patch on top of this branch: its code is already
included here, together with review corrections and the newer master cleanup fix.
The branch does not launch training, set the scientific prior, or contain a
nonlinear IDE implementation. Reproduce the checks on the generation host before
bulk production, and resolve failures across the intended training domain.
