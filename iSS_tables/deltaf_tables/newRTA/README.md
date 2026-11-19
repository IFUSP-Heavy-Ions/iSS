# Coefficients for the RTA delta f with an energy-dependent relaxation time

Used by `bulk_deltaf_kind = 31`: relaxation time approximation with
tau_R ~ E^gamma, shear and bulk, no conserved charges. gamma is the parameter
`deltaf_newRTA_gamma`; the tables are given at gamma = 0, 0.5, 1 and 1.5 and
interpolated linearly in between.

| folder | file | content |
|---|---|---|
| `shear/` | `beta_pi_gam-EQ-<gamma>-weights-quant.dat` | T [GeV], beta_pi/T^4 [(GeV fm)^-3] |
| `bulk/` | `beta_PI-gam=<gamma>-weights-quant.dat` | T [GeV], beta_Pi/T^4 [(GeV fm)^-3] |
| `bulk/` | `varphi-gam=<gamma>-weights-quant.dat` | T [GeV], varphi |
| `bulk/` | `cs2-cfv-weights-quant.dat` | T [GeV], 1/3 - cs^2 |
| `thermoIntegrals/` | `J_gamma_{BE,FD}.dat`, `J_gammaplus2_{BE,FD}.dat` | m/T, then J at gamma = 0, 0.5, 1, 1.5 |
| `thermoIntegrals/` | `J_2_{BE,FD}.dat` | m/T, J_2 |

The temperature grid is 0.01 to 0.1975 GeV in steps of 0.0025 GeV (76 rows; the
code holds the end values outside it). The m/T grid is 0.5 to 15 in steps of
0.05. `bulk/varphi-gam=0,0-densities-quant.dat` is not read by the code.

## Origin

These are the files of `newRTA.tar.gz` as fetched by
`iSS_tables/deltaf_tables/download_newRTA.sh` on the branch `newRTA` of
`github.com/chunshen1987/iSS` (commit `3af73d1`), downloaded 2026-11-19, md5 of
the tarball `2da50a1134053727f082334131767eb6`. They are kept in the repository
because the download location is not versioned. The thermodynamic integrals can
be regenerated with `utilities/computeThermodynamicJ.py`.
