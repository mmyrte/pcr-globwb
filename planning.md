# Refactoring & Modernization Plan for pcr-globwb

## Primary Objective

Modernize the `pcr-globwb` hydrological model (specifically targeting a Switzerland-wide domain of roughly 360x230 km, kilometre resolution) into a standard, testable, and highly performant Python package. The ultimate goal is to replace the legacy `pcraster` dependency with a modern, well-vectorized array processing library for CPU usage (such as NumPy/Xarray) that supports compressed I/O formats (NetCDF/Zarr).

## Specific Goals

1. **Codebase Triage & Cleanup:** Identify the primary execution entrypoints and eliminate dead code, stale files, and redundant configuration files (retaining one well-annotated example).
2. **Project Modernization:** Restructure the repository into a standard Python package using `uv` for environment management, `pyproject.toml` for dependencies, and `ruff`/`basedpyright` for linting and type checking.
3. **Robust Test Infrastructure:** Establish an automated testing suite (golden master tests) using a representative slice of the real-world NetCDF input data (`PCRInput_CH`). Floating-point tolerance is acceptable as long as it does not add up to significant differences due to error propagation.
4. **Computational Refactoring:** Decouple the core logic from `pcraster` and migrate to a modern, vectorized framework.

## Working Style

Ask the user for clarification if specifications are ambiguous.
Do not commit changes by yourself; give the user the option to intervene.
Intermediary steps and progress should be written/edited in this file.

---

## Progress Tracker & Execution Plan

### Phase 1: Discovery & Triage

- [x] **T1.1:** Trace the codebase to locate the primary execution entrypoint(s). `model/deterministic_runner.py` identified.
- [x] **T1.2:** Audit the repository for redundant configuration files, scripts, and legacy artifacts.
- [x] **T1.3:** Compile a list of _distinct usages_ and operations of `pcraster` required by the codebase.

**Collected `pcraster` usages (from T1.3):**

```text
pcr.abs
pcr.accufractionstate
pcr.accuthresholdstate
pcr.accutraveltimeflux
pcr.accutraveltimestate
pcr.acos
pcr.aguila
pcr.areaarea
pcr.areaaverage
pcr.areamajority
pcr.areamaximum
pcr.areaminimum
pcr.areaorder
pcr.areatotal
pcr.asin
pcr.atan
pcr.boolean
pcr.catchmenttotal
pcr.cellvalue
pcr.clone
pcr.clump
pcr.cos
pcr.cover
pcr.defined
pcr.downstream
pcr.downstreamdist
pcr.exp
pcr.framework
pcr.ifthen
pcr.ifthenelse
pcr.kinematic
pcr.ldd
pcr.lddcreate
pcr.ldddist
pcr.lddmask
pcr.lddrepair
pcr.ln
pcr.log10
pcr.mapmaximum
pcr.mapminimum
pcr.maptotal
pcr.max
pcr.min
pcr.nominal
pcr.normal
pcr.numpy2pcr
pcr.ordinal
pcr.pcr2numpy
pcr.pcrand
pcr.pcrnot
pcr.pcror
pcr.readmap
pcr.report
pcr.rounddown
pcr.roundoff
pcr.roundup
pcr.scalar
pcr.Scalar
pcr.setclone
pcr.setglobaloption
pcr.sin
pcr.spatial
pcr.sqrt
pcr.subcatchment
pcr.tan
pcr.uniqueid
pcr.upstream
pcr.windowaverage
pcr.windowmajority
pcr.xcoordinate
pcr.ycoordinate
```

### Phase 2: Packaging & Modernization

- [ ] **T2.1:** Use `uv` for local virtual environment management and create a `pyproject.toml`.
- [ ] **T2.2:** Reorganize the code into a standard Python layout (e.g., `src/pcr_globwb/`).
- [ ] **T2.3:** Introduce and enforce basic linting, formatting, and type-checking using `ruff` and `basedpyright`.

### Phase 3: Test Data Preparation & Golden Testing

- [ ] **T3.1:** Inspect `test_data_raw` to identify necessary NetCDF files for a minimum viable run.
- [ ] **T3.2:** Slice out spatial domain using `ncks` (46.3-46.5 N, 7.7-7.9 E).
- [ ] **T3.3:** Verify the test domain touches all model mechanics (glaciers, etc.).
- [ ] **T3.4:** Run `pcraster`-based model on dataset to create Golden Master; write integration tests.

### Phase 4: Abstraction & Refactoring (The Heavy Lifting)

- [ ] **T4.1:** Abstract spatial/I/O operations.
- [ ] **T4.2:** Implement modern backend (NumPy/Xarray).
- [ ] **T4.3:** Run Golden Tests against the new backend.

### Phase 5: I/O and Performance Optimization

- [ ] **T5.1:** Strip out old `pcraster` backend.
- [ ] **T5.2:** Optimize new backend for vectorization and compressed I/O.
