.. _restart_tout_wind:

*****************
restart_tout_wind
*****************

.. tags:: binary, mass-transfer, roche-lobe-overflow, point-mass, stellar-winds, hot-wind, tidally-enhanced-wind, restart, intermediate-mass, main-sequence

This test checks that a restart from a photo gives the same result as
an uninterrupted run when the orbit-phase integrals are in use.

This test case has 1 part.

* Part 1 (``inlist_project``) is :ref:`star_plus_point_mass` (a 5 Msun
  donor and a 10 Msun point mass in a 1 day orbit) with a Dutch wind
  and tidally enhanced wind (``do_enhance_wind_1 = .true.``, Tout &
  Eggleton 1988). The enhanced wind is integrated over the orbit using
  ``theta_co`` and ``time_co``. These arrays are not stored in photos,
  so ``binarydata_init`` must rebuild them on restart. The run stops
  when the central hydrogen mass fraction of the donor drops below 1e-2.

How this test passes
====================

The test expects ``termination code: xa_central_lower_limit`` and
checks ``final.mod``. The restart from photo ``x100`` must reproduce the
same ``final.mod`` checksum.
