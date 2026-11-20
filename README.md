iSS
===

Monte Carlo sampler for particle distribution from Cooper-Frye freeze-out procedure

Particle lists
--------------

`particle_list` (default 0) selects the hadron list iSS samples, decays and
takes its hadron gas from (`regulateEOS`): 0 is the afterburner's
(`urqmd_v3.3+`, `SMASH`, `s95p-v1`); 21 is PDG2021+ (J. Salinas San Martín et
al., arXiv:2309.01737). Its tables are built by `utilities/get_PDG2021Plus.sh`,
which fetches the list at a pinned commit (GPL-3, so not committed here) and
runs `utilities/convert_massorder_pdg.py` and `utilities/make_hrg_eos.py`;
the same two scripts turn any list in the mass-ordered Thermal-FIST layout
into iSS tables.

Resonance decays (`perform_decays = 1`)
---------------------------------------

A decay channel is drawn with the table's branching ratios, normalised over
the channels the mother can reach. A channel that is closed at the daughters'
pole masses is reached by giving its broad daughters (width above 0.1 MeV)
masses from their Breit-Wigner distributions, each no lower than the mass at
which all of its own channels stay open, so that the daughter keeps the
table's branching ratios; a channel that cannot be reached that way is
skipped and its share goes to the others. Two- and three-body channels use
the original routines and random-number sequence, four- and five-body ones
are uniform in phase space (Raubold-Lynch). A negative number of daughters
(the UrQMD-derived tables mark closed channels with it) is read as its
absolute value. Before this, a mother that drew a four-body or a
negative-count channel vanished without daughters (about 2% of the final
pions with `urqmd_v3.3+`), and a list whose branching ratios do not add up to
one could crash.
