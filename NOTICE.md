# Third-party software

FlexAppeal is MIT licensed: see [`LICENSE`](LICENSE).

It uses, and in two cases redistributes, software under its own licences.
None of them restricts the MIT terms, but their notices are retained here as
those licences require. Kept separate from `LICENSE` so GitHub's licence
detector still reads the repository as MIT.

THIRD-PARTY SOFTWARE

FlexAppeal uses, and in two cases redistributes, software under its own
licences. None of them restricts the MIT terms above, but their notices are
retained here as those licences require.

REDISTRIBUTED IN THIS REPOSITORY
--------------------------------

  Plotly.js 2.35.2                                              MIT
  flexappeal/static/vendor/plotly.min.js
  Copyright (c) 2012-2024 Plotly, Inc.
  https://github.com/plotly/plotly.js

  Mol* 5.11.0                                                   MIT
  flexappeal/static/vendor/molstar.js
  flexappeal/static/vendor/molstar.css
  Copyright (c) 2018-2026 Mol* contributors
  https://github.com/molstar/molstar

REQUIRED AT RUN TIME, NOT REDISTRIBUTED
---------------------------------------

  OpenMM                          MIT and LGPL-3.0-or-later (see its own COPYING)
  PDBFixer                        MIT
  openmmforcefields               MIT
  OpenFF Toolkit                  MIT
  AmberTools                      GPL-3.0 / LGPL-3.0 (component dependent)
  MDTraj                          LGPL-2.1-or-later
  gemmi                           MPL-2.0
  Flask, NumPy, SciPy             BSD-3-Clause

FlexAppeal imports MDTraj and gemmi, and calls the rest either through
generated scripts that run on the user's own machine or not at all. LGPL and
MPL both permit use from differently licensed code where the library is used
unmodified through its published interface, which is the case here: no
third-party source is patched or statically combined.

AmberTools is installed into a generated bundle's own environment only when
AM1-BCC charges are requested, and is invoked as a separate program rather than
linked. It is never distributed as part of FlexAppeal.

FORCE FIELDS AND DATA
---------------------

Force-field parameters (AMBER, CHARMM36, OpenFF, GAFF), chemical definitions
from the RCSB Chemical Component Dictionary, and structures fetched from the
RCSB PDB, OPM or the AlphaFold Database carry their own terms and citation
expectations. Cite the force field and the structures you use.
