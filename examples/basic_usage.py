from compressive_splitting import CompressiveSplittingModel

# Example: L/B = 2 and nu = 0.30, matching the representative geometry
# used for the stress-field discussion in the paper.
model = CompressiveSplittingModel(
    E=1.0,
    nu=0.30,
    aspect_ratio=2.0,
    applied_strain=-0.01,
    nx=80,
)

result = model.solve()

if model.comm.rank == 0:
    print(f"k = {result.k:.6f} ({100*result.k:.3f} %)")

    # Optional strength prediction from Eq. (4).
    sigma_t = 10.0      # e.g. MPa
    sigma_conf = 0.0    # same units as sigma_t
    sigma_c_crit = model.critical_compressive_strength(
        tensile_strength=sigma_t,
        confinement=sigma_conf,
    )
    print(f"Predicted splitting strength = {sigma_c_crit:.3f}")

    # Optional confinement required to remove the tensile hotspot
    # at a specified axial compressive stress magnitude.
    sigma_c = 100.0
    sigma_conf_crit = model.critical_confinement(sigma_c)
    print(f"Critical confinement at sigma_c={sigma_c:g}: {sigma_conf_crit:.3f}")
