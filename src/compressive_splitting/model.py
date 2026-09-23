"""Finite-element evaluation of the boundary-induced tension coefficient k.

This module implements the pre-fracture linear-elastic boundary-value problem
used in:

    Khodadad et al., "Compressive Splitting in Brittle Solids:
    The Inverse of Wrinkling in Sheets", Physical Review Letters (2026).

The implementation intentionally keeps the model assumptions of the paper:
2D plane strain, homogeneous isotropic linear elasticity, small strain,
rectangular geometry, and fully clamped top/bottom loading boundaries.
"""

from dataclasses import dataclass

import dolfinx
import dolfinx.fem.petsc
import numpy as np
import ufl
from dolfinx.mesh import DiagonalType
from mpi4py import MPI


@dataclass(frozen=True)
class SplittingResult:
    """Scalar outputs from one finite-element solve."""

    k: float
    sigma_xx_max: float
    sigma_c: float
    aspect_ratio: float
    poisson_ratio: float
    youngs_modulus: float


class CompressiveSplittingModel:
    """Compute the boundary-induced transverse tension coefficient ``k``.

    Parameters
    ----------
    E
        Young's modulus. Any consistent stress unit may be used.
        In the linear-elastic problem, ``k`` is dimensionless and independent
        of the absolute value of ``E``.
    nu
        Poisson's ratio. The displacement formulation implemented here assumes
        ``-1 < nu < 0.5``. The exactly incompressible case ``nu = 0.5``
        requires the mixed displacement-pressure formulation used separately
        in the research code.
    aspect_ratio
        Specimen aspect ratio ``L/B``, where ``L`` is the loading-direction
        length and ``B`` is the transverse width.
    applied_strain
        Uniform nominal axial compressive strain used to load the model.
        Negative values denote compression. Because the problem is linear,
        ``k`` does not depend on its magnitude.
    nx
        Number of mesh divisions across the specimen width ``B``. The number
        along ``L`` is chosen as ``max(int(nx * L/B), 8)``.

    Notes
    -----
    Boundary conditions match the paper/code:
      * top and bottom faces: u_x = 0;
      * top and bottom move symmetrically in y;
      * side faces are traction free.
    """

    def __init__(
        self,
        E: float,
        nu: float,
        aspect_ratio: float,
        applied_strain: float = -0.01,
        nx: int = 80,
        comm=MPI.COMM_WORLD,
    ):
        if E <= 0:
            raise ValueError("E must be positive.")
        if not (-1.0 < nu < 0.5):
            raise ValueError(
                "This displacement formulation requires -1 < nu < 0.5. "
                "The exactly incompressible nu=0.5 case requires a mixed formulation."
            )
        if aspect_ratio <= 0:
            raise ValueError("aspect_ratio = L/B must be positive.")
        if applied_strain >= 0:
            raise ValueError("applied_strain must be negative for compression.")
        if nx < 4:
            raise ValueError("nx must be at least 4.")

        self.E = float(E)
        self.nu = float(nu)
        self.aspect_ratio = float(aspect_ratio)
        self.applied_strain = float(applied_strain)
        self.nx = int(nx)
        self.comm = comm

        self._result = None
        self._domain = None
        self._uh = None
        self._sigma_xx = None

    def solve(self) -> SplittingResult:
        """Solve the plane-strain elasticity problem and return ``k``."""

        # Geometry convention from the paper: B = width, L = loading length.
        B = 1.0
        L = self.aspect_ratio
        ny = max(int(self.nx * L), 8)

        mu = self.E / (2.0 * (1.0 + self.nu))
        lmbda = (
            self.E
            * self.nu
            / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        )

        sigma_c = abs(self.E * self.applied_strain)
        applied_displacement = L * self.applied_strain

        domain = dolfinx.mesh.create_rectangle(
            self.comm,
            [np.array([-B / 2, -L / 2]), np.array([B / 2, L / 2])],
            [self.nx, ny],
            cell_type=dolfinx.mesh.CellType.triangle,
            diagonal=DiagonalType.crossed,
        )

        V = dolfinx.fem.functionspace(
            domain, ("Lagrange", 1, (domain.geometry.dim,))
        )
        uh = dolfinx.fem.Function(V, name="Displacement")
        u_trial = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)

        def top_bottom_boundary(x):
            return np.isclose(np.abs(x[1]), L / 2)

        facets = dolfinx.mesh.locate_entities_boundary(
            domain, domain.topology.dim - 1, top_bottom_boundary
        )
        boundary_dofs = dolfinx.fem.locate_dofs_topological(
            V, domain.topology.dim - 1, facets
        )

        class BoundaryDisplacement:
            def __call__(self, x):
                values = np.zeros(
                    (domain.geometry.dim, x.shape[1]),
                    dtype=dolfinx.default_scalar_type,
                )
                # Fully clamp lateral Poisson expansion on the loading faces.
                values[0] = 0.0
                # Symmetric axial compression.
                values[1] = (
                    np.sign(x[1]) * applied_displacement / 2.0
                )
                return values

        u_D = dolfinx.fem.Function(V)
        u_D.interpolate(BoundaryDisplacement())
        bc = dolfinx.fem.dirichletbc(u_D, boundary_dofs)

        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        def sigma(u):
            # 3D Lamé constants used in the 2D plane-strain formulation.
            return (
                lmbda * ufl.nabla_div(u) * ufl.Identity(len(u))
                + 2.0 * mu * epsilon(u)
            )

        a = ufl.inner(sigma(u_trial), epsilon(v)) * ufl.dx
        zero_body_force = dolfinx.fem.Constant(
            domain,
            np.array((0.0, 0.0), dtype=dolfinx.default_scalar_type),
        )
        L_form = ufl.inner(zero_body_force, v) * ufl.dx

        problem = dolfinx.fem.petsc.LinearProblem(
            a,
            L_form,
            u=uh,
            bcs=[bc],
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        )
        problem.solve()

        V_scalar = dolfinx.fem.functionspace(domain, ("Lagrange", 1))
        sigma_xx = dolfinx.fem.Function(V_scalar, name="sigma_xx")
        sigma_expr = dolfinx.fem.Expression(
            sigma(uh)[0, 0], V_scalar.element.interpolation_points()
        )
        sigma_xx.interpolate(sigma_expr)

        local_max = float(np.max(sigma_xx.x.array))
        sigma_xx_max = self.comm.allreduce(local_max, op=MPI.MAX)

        k = sigma_xx_max / sigma_c
        if k < 0.0:
            k = 0.0

        result = SplittingResult(
            k=float(k),
            sigma_xx_max=float(sigma_xx_max),
            sigma_c=float(sigma_c),
            aspect_ratio=self.aspect_ratio,
            poisson_ratio=self.nu,
            youngs_modulus=self.E,
        )

        self._result = result
        self._domain = domain
        self._uh = uh
        self._sigma_xx = sigma_xx
        return result

    @property
    def result(self) -> SplittingResult:
        """Return the latest solve result."""
        if self._result is None:
            raise RuntimeError("Call solve() first.")
        return self._result

    def critical_compressive_strength(
        self, tensile_strength: float, confinement: float = 0.0
    ) -> float:
        r"""Return the predicted splitting strength.

        Uses Eq. (4) of the paper:

            sigma_c,crit = (sigma_t + sigma_conf) / k

        ``tensile_strength`` and ``confinement`` must use the same stress unit.
        """
        if self._result is None:
            self.solve()
        if tensile_strength <= 0:
            raise ValueError("tensile_strength must be positive.")
        if confinement < 0:
            raise ValueError("confinement must be non-negative.")
        if self._result.k <= 0:
            raise ZeroDivisionError(
                "k is zero; this model predicts no tensile splitting threshold."
            )
        return (float(tensile_strength) + float(confinement)) / self._result.k

    def critical_confinement(self, axial_compressive_stress: float) -> float:
        r"""Return confinement that suppresses transverse tension.

        From superposition:

            sigma_conf / sigma_c = k

        so ``sigma_conf,crit = k * sigma_c``.
        """
        if self._result is None:
            self.solve()
        if axial_compressive_stress < 0:
            raise ValueError("Pass axial_compressive_stress as a positive magnitude.")
        return self._result.k * float(axial_compressive_stress)
