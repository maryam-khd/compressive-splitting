# Compressive Splitting in Brittle Solids

Finite-element implementation accompanying:

**M. Khodadad, F. Barthelat, J. D. Clayton, G. Gazonas, and K. Dayal,  
“Compressive Splitting in Brittle Solids: The Inverse of Wrinkling in Sheets,”  
Physical Review Letters (2026).**

DOI: [10.1103/j8x1-hy4t](https://doi.org/10.1103/j8x1-hy4t)

This repository is intentionally small. It extracts the core pre-fracture
finite-element calculation used to evaluate the **boundary-induced tension
coefficient**

```math
k = \frac{\sigma_{xx,\max}}{\sigma_c}
```

## What the code computes

For a rectangular brittle-solid specimen compressed between rigid, fully
clamped loading faces, the solver computes the plane-strain linear-elastic
stress field and returns `k`.

In the model,

```math
\sigma_{xx,\max}
=
k(L/B,\nu)\,\sigma_c-\sigma_{\mathrm{conf}}
```

so the predicted onset of tensile splitting is

```math
\sigma_{c,\mathrm{crit}}
=
\frac{\sigma_t+\sigma_{\mathrm{conf}}}{k(L/B,\nu)}
```

Here:

- `L/B` is the specimen aspect ratio;
- `nu` is Poisson's ratio;
- `sigma_t` is tensile strength;
- `sigma_conf` is lateral confinement;
- `sigma_c` is the positive magnitude of axial compression.

Because the governing problem is linear elastic, `k` depends on `L/B` and
`nu`, not on the absolute magnitude of `E` or the imposed strain. `E` is kept
as an input so dimensional stress fields remain available and the code mirrors
the finite-element formulation used in the paper.

## Model assumptions

The implementation preserves the assumptions used for the finite-element
calculation in the paper:

- homogeneous, isotropic, linear elasticity;
- small strain;
- 2D **plane strain**;
- compression along the `y` direction;
- top and bottom loading faces fully clamped laterally (`u_x = 0`);
- symmetric prescribed axial displacement on the top and bottom faces.

This code evaluates the **pre-fracture elastic stress state**. It does not
simulate crack propagation or fracture evolution.

The displacement formulation is for `-1 < nu < 0.5`. The exactly
incompressible case `nu = 0.5` requires the mixed displacement-pressure
formulation used separately in the research code.

## Installation

A working FEniCSx / DOLFINx environment is required. The most reliable route is
a conda environment containing DOLFINx, PETSc, mpi4py, NumPy, and UFL.

After activating that environment:

```bash
python -m pip install -e .
```

## Minimal use

```python
from compressive_splitting import CompressiveSplittingModel

model = CompressiveSplittingModel(
    E=1.0,
    nu=0.30,
    aspect_ratio=2.0,
)

result = model.solve()

print(result.k)
```

For the confinement-dependent strength relation:

```python
sigma_c_crit = model.critical_compressive_strength(
    tensile_strength=10.0,
    confinement=5.0,
)
```

For the confinement that exactly suppresses the tensile hotspot at a specified
axial stress:

```python
sigma_conf_crit = model.critical_confinement(
    axial_compressive_stress=100.0
)
```

All stress quantities must use the same units.

A complete runnable example is in `examples/basic_usage.py`.

## Citation

If you use this code or the model in your work, please cite:

> M. Khodadad, F. Barthelat, J. D. Clayton, G. Gazonas, and K. Dayal,  
> “Compressive Splitting in Brittle Solids: The Inverse of Wrinkling in Sheets,”  
> *Physical Review Letters* (2026).  
> https://doi.org/10.1103/j8x1-hy4t

BibTeX:

```bibtex
@article{khodadad2026compressive,
  title   = {Compressive Splitting in Brittle Solids: The Inverse of Wrinkling in Sheets},
  author  = {Khodadad, Maryam and Barthelat, Francois and Clayton, John D. and Gazonas, George and Dayal, Kaushik},
  journal = {Physical Review Letters},
  year    = {2026},
  doi     = {10.1103/j8x1-hy4t}
}
```

GitHub's **Cite this repository** button is also populated by `CITATION.cff`.

## Scope

The goal of this repository is reproducibility of the core coefficient `k`.