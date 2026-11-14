.. _tidal_circ_ftid:

***************
tidal_circ_ftid
***************

.. tags:: binary, binary-evolution, dual-star-evolution, eccentric-orbit, tidal-circularization, low-mass, main-sequence

This test checks that tidal circularisation by star 2 is scaled by
``Ftid_2``.

This test case has 1 part.

* Part 1 (``inlist_project``) evolves a 1.0 Msun and a 0.8 Msun
  zero-age main sequence star in a 1.5 day orbit with an initial
  eccentricity of 0.2. Both stars use ``circ_type = 'Hut_conv'``.
  ``Ftid_1 = 0`` switches off star 1's contribution, so the
  eccentricity can only decrease through star 2's term, scaled by
  ``Ftid_2 = 1``. The run stops at an age of 5e7 years.

How this test passes
====================

``extras_binary_after_evolve`` compares the final eccentricity with
``initial_eccentricity`` and prints
``tidal_circ_ftid: star 2 tides reduced the eccentricity`` if it has
decreased. If star 2's term used ``Ftid_1`` instead, ``edot_tidal``
would be exactly zero and the test would fail. The final model is
``final.mod`` and the restart uses the antepenultimate photo (``auto``).
