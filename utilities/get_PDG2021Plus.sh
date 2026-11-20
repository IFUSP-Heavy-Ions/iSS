#!/usr/bin/env bash
# Build iSS's tables for the PDG2021+ hadron list (particle_list = 21):
#
#   <tables>/pdg-PDG2021Plus.dat               species and decays
#   <tables>/chosen_particles_PDG2021Plus.dat  species sampled
#   <tables>/EOS_tables/HRGEOS_PST-PDG2021Plus.dat  its hadron gas
#   <tables>/PDG2021Plus_SOURCE.txt            where they came from
#
# usage: utilities/get_PDG2021Plus.sh [tables directory, default iSS_tables]
#
# The list is fetched from github.com/jordissm/PDG2021Plus at a pinned
# commit and checked against its sha256.  It is distributed under GPL-3 and
# is not committed here; cite J. Salinas San Martín et al., arXiv:2309.01737,
# when you use it.
set -euo pipefail

COMMIT=002232d34bd3e15df037a2c63d6a6d8d051f5220
FILE=hadron-lists/PDG2021Plus/mass_order/decays_PDG2021Plus-1.0_massorder.dat
SHA256=e815bed12c53aec460d72702f31dc1672a4c6913a08e1fd0d28b6d3c720939dd

here=$(cd "$(dirname "$0")" && pwd)
tables=$(cd "${1:-$here/../iSS_tables}" && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

url=https://raw.githubusercontent.com/jordissm/PDG2021Plus/$COMMIT/$FILE
curl -fsSL --retry 3 "$url" -o "$tmp/massorder.dat"
echo "$SHA256  $tmp/massorder.dat" | sha256sum -c --quiet -

(cd "$tmp" && python3 "$here/convert_massorder_pdg.py" massorder.dat \
     PDG2021Plus)
cp "$tmp/pdg-PDG2021Plus.dat" "$tmp/chosen_particles_PDG2021Plus.dat" \
   "$tables/"
python3 "$here/make_hrg_eos.py" "$tables/pdg-PDG2021Plus.dat" \
    "$tables/EOS_tables/HRGEOS_PST-PDG2021Plus.dat"
cat > "$tables/PDG2021Plus_SOURCE.txt" <<EOF
PDG2021+ hadron list, converted for iSS by utilities/get_PDG2021Plus.sh
source:  https://github.com/jordissm/PDG2021Plus/blob/$COMMIT/$FILE
sha256:  $SHA256
licence: GPL-3.0 (https://github.com/jordissm/PDG2021Plus/blob/$COMMIT/LICENSE)
cite:    J. Salinas San Martín et al., arXiv:2309.01737
EOF
echo "PDG2021Plus tables written to $tables"
