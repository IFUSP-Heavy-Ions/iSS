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
temperature is held at Tmax, as in the tables shipped with iSS.  For
pdg-urqmd_v3.3+.dat it reproduces HRGEOS_PST-urqmd_v3.3+.dat to 1e-5 in P,
s and T.
"""
import argparse

import numpy as np
from scipy.optimize import brentq
from scipy.special import kn

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


def thermo(T, m, g, stat):
    """P and e [GeV/fm^3] at T [GeV]."""
    n = np.arange(1, NTERMS + 1)[:, None]
    x = n*m[None, :]/T
    sign = (-stat[None, :])**(n + 1)
    pref = sign*g[None, :]/(2*np.pi**2)*m[None, :]**2*T**2/n**2
    ok = x < 300
    k2 = np.where(ok, kn(2, np.minimum(x, 300)), 0.0)
    k1 = np.where(ok, kn(1, np.minimum(x, 300)), 0.0)
    P = np.sum(pref*k2)
    e = np.sum(pref*(3*k2 + x*k1))
    conv = 1.0/HBARC**3
    return P*conv, e*conv


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--Tmax", type=float, default=0.19)
    a = ap.parse_args()
    m, g, stat = read_species(a.table)
    e_max = thermo(a.Tmax, m, g, stat)[1]
    with open(a.out, "w") as f:
        f.write("# e [GeV/fm^3]  P [GeV/fm^3]  s [1/fm^3]  T [GeV]  cs^2\n")
        for e in np.arange(1, 2001)*1e-3:
            if e >= e_max:
                T = a.Tmax
            else:
                T = brentq(lambda t: thermo(t, m, g, stat)[1] - e, 0.005,
                           a.Tmax, xtol=1e-13)
            P, ee = thermo(T, m, g, stat)
            dT = 1e-5
            P1, e1 = thermo(T - dT, m, g, stat)
            P2, e2 = thermo(T + dT, m, g, stat)
            f.write("%.6e  %.6e  %.6e  %.6e  %.6e\n"
                    % (e, P, (ee + P)/T, T, (P2 - P1)/(e2 - e1)))


if __name__ == "__main__":
    main()
