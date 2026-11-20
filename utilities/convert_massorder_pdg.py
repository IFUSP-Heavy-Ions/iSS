#!/usr/bin/env python3
"""Convert a mass-ordered hadron list with decays (the decays_*_massorder.dat
files of github.com/jordissm/PDG2021Plus, also used by Thermal-FIST and CCAKE)
into iSS's particle table and chosen-particle list.

    convert_massorder_pdg.py <decays_*_massorder.dat> <name>

writes pdg-<name>.dat and chosen_particles_<name>.dat in the working
directory.  The input has 12-column species lines (id, name, mass, width,
spin degeneracy, B, S, C, b, isospin I, charge, number of channels), each
followed by 8-column decay lines (id, number of daughters, branching ratio,
five daughter ids), particles and antiparticles alike.  iSS's table is the
same except that

* antibaryons are left out (iSS adds them, with conjugated decays), and
* the isospin column is the degeneracy 2I + 1.

A species is stable for iSS when its first channel has one daughter (itself).
Every daughter of a baryon must appear earlier in the file (iSS builds the
antibaryon channels from the species read so far); the script checks it.
The chosen-particle list (the species iSS samples) is every hadron of the
table and every antibaryon; the photon is not sampled.
"""
import sys


def read_massorder(name):
    rows = [l.split() for l in open(name)
            if l.strip() and not l.lstrip().startswith("#")]
    species = []
    i = 0
    while i < len(rows):
        r = rows[i]
        if len(r) != 12:
            raise ValueError("expected a 12-column species line, got %r" % r)
        n = int(r[11])
        channels = rows[i + 1:i + 1 + n]
        for c in channels:
            if len(c) != 8 or int(c[0]) != int(r[0]):
                raise ValueError("bad decay line for %s: %r" % (r[0], c))
        species.append((r, channels))
        i += 1 + n
    return species


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    source, name = sys.argv[1:3]
    species = read_massorder(source)
    out, chosen, seen = [], [], set()
    for r, channels in species:
        pdg, baryon = int(r[0]), int(r[5])
        if baryon < 0:
            continue  # iSS adds antibaryons itself
        for c in channels if baryon > 0 else []:
            if int(c[1]) == 1:
                continue
            for d in c[3:8]:
                d = int(d)
                if d and d not in seen and -d not in seen:
                    raise ValueError("daughter %d of %d is not listed before it"
                                     % (d, pdg))
        seen.add(pdg)
        isospin_deg = int(round(2*float(r[9]) + 1))
        out.append("%10d  %-22s %10.5f %12.5e %3d %3d %3d %3d %3d %3d %3d %3d"
                   % (pdg, r[1], float(r[2]), float(r[3]), int(r[4]), baryon,
                      int(r[6]), int(r[7]), int(r[8]), isospin_deg, int(r[10]),
                      len(channels)))
        for c in channels:
            out.append("%10d %3d %10.6f %10d %10d %10d %10d %10d"
                       % (pdg, int(c[1]), float(c[2]), *[int(x) for x in c[3:8]]))
        if pdg != 22:
            chosen.append(pdg)
            if baryon > 0:
                chosen.append(-pdg)
    with open("pdg-%s.dat" % name, "w") as f:
        f.write("\n".join(out) + "\n")
    with open("chosen_particles_%s.dat" % name, "w") as f:
        f.write("\n".join(str(p) for p in chosen) + "\n")
    print("pdg-%s.dat: %d species (antibaryons added by iSS); "
          "chosen_particles_%s.dat: %d" % (name, len(seen), name, len(chosen)))


if __name__ == "__main__":
    main()
