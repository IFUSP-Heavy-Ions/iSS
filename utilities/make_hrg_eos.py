#!/usr/bin/env python3
"""Hadron-resonance-gas equation of state at mu = 0 from an iSS particle table,
in the layout of iSS_tables/EOS_tables/HRGEOS_PST-*.dat (iSS uses it to turn
the freeze-out energy density back into a temperature, regulateEOS = 1).

    make_hrg_eos.py <pdg-*.dat> <HRGEOS_PST-*.dat> [--Tmax 0.19]

Every species of the table and, for baryons, its antiparticle (iSS adds
those), as ideal gases with their pole masses and spin degeneracies, with
quantum statistics (Bessel series, 20 terms); the photon is left out.  Rows:
e = 0.001 .. 2 GeV/fm^3 in steps of 0.001, columns e, P [GeV/fm^3],
s [1/fm^3], T [GeV], cs^2 = dP/de; above the energy density of T = Tmax the
temperature is held at Tmax and s = (e + P)/T, as in the tables shipped
with iSS.  For pdg-urqmd_v3.3+.dat it reproduces HRGEOS_PST-urqmd_v3.3+.dat
to 3e-5 in P, s and T.
"""
import argparse

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import k0, k1

HBARC = 0.19733  # GeV fm
NTERMS = 20


def read_species(table):
    """(mass, degeneracy, statistics) arrays; statistics -1 Bose, +1 Fermi."""
    m, g, stat = [], [], []
    for line in open(table):
        c = line.split()
        if len(c) < 12 or c[1].lstrip("-").isdigit():
            continue  # decay lines
        pdg, mass, deg, baryon = int(c[0]), float(c[2]), int(c[4]), int(c[5])
        if pdg == 22:
            continue
        copies = 2 if baryon != 0 else 1
        for _ in range(copies):
            m.append(mass)
            g.append(deg)
            stat.append(1 if baryon != 0 else -1)
    return np.array(m), np.array(g, float), np.array(stat, float)


def thermo(T, m, g, stat, chunk=200):
    """P and e [GeV/fm^3] at the temperatures T [GeV] (array)."""
    T = np.atleast_1d(np.asarray(T, float))
    P = np.zeros(T.shape)
    e = np.zeros(T.shape)
    for c in range(0, T.size, chunk):
        t = T[c:c + chunk, None]
        for n in range(1, NTERMS + 1):
            x = n*m[None, :]/t
            sign = (-stat)**(n + 1)
            pref = sign*g/(2*np.pi**2)*m**2/n**2
            ok = x < 300
            xs = np.minimum(x, 300)
            K1 = np.where(ok, k1(xs), 0.0)
            K2 = np.where(ok, k0(xs) + 2*K1/xs, 0.0)
            P[c:c + chunk] += t[:, 0]**2*np.sum(pref*K2, axis=1)
            e[c:c + chunk] += t[:, 0]**2*np.sum(pref*(3*K2 + x*K1), axis=1)
    conv = 1.0/HBARC**3
    return P*conv, e*conv


def temperature(e_target, m, g, stat, T_lo, T_hi):
    """T with e(T) = e_target: a cubic spline of T(e) on a fine T grid,
    then Newton steps on e(T) itself (e increases with T)."""
    T_grid = np.linspace(T_lo, T_hi, 4001)
    e_grid = thermo(T_grid, m, g, stat)[1]
    spline = CubicSpline(e_grid, T_grid)
    T = spline(e_target)
    dTde = spline.derivative()
    for _ in range(3):
        T = T - (thermo(T, m, g, stat)[1] - e_target)*dTde(e_target)
    return T


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--Tmax", type=float, default=0.19)
    a = ap.parse_args()
    m, g, stat = read_species(a.table)
    e_max = thermo(a.Tmax, m, g, stat)[1][0]
    e = np.arange(1, 2001)*1e-3
    T = np.full(e.shape, a.Tmax)
    below = e < e_max
    T[below] = temperature(e[below], m, g, stat, 0.005, a.Tmax)
    P = thermo(T, m, g, stat)[0]
    dT = 1e-5
    P1, e1 = thermo(T - dT, m, g, stat)
    P2, e2 = thermo(T + dT, m, g, stat)
    with open(a.out, "w") as f:
        f.write("# e [GeV/fm^3]  P [GeV/fm^3]  s [1/fm^3]  T [GeV]  cs^2\n")
        for i in range(e.size):
            # s = (e + P)/T with the row's e, also where T is held at Tmax
            f.write("%.6e  %.6e  %.6e  %.6e  %.6e\n"
                    % (e[i], P[i], (e[i] + P[i])/T[i], T[i],
                       (P2[i] - P1[i])/(e2[i] - e1[i])))


if __name__ == "__main__":
    main()
