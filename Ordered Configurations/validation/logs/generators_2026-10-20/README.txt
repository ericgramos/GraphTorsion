Checks of the prime-power and subgraph tools, October 20, 2026, after the code was reorganised around the shared
headers certificates/swiat_cells.h, elim_modq.h and elim_z.h.

Regression (outputs of the reorganised programs compared byte for byte with those of the first versions):
  gens.cmp     certq generator mode on the H-first edge order: identical generator file
  oqA.cmp      orderq image test for the triangle with a leaf: identical output
  solA.cmp     exactsolve: identical cycle y
  yred.cmp     reducecyc: identical coordinates
  c16.*        certq certificate mode for order 16 (same pairing 68275 as before) and its verification
  ytori.*      orderq with the tori of H as relations: y vanishes (order 1, no free coordinates)
Verification of the certificates with validation/verify_cert.py: verify_*.log (order 4, 8, 16, the triangle with a leaf,
and T(1,2,3) as a check that ordinary prime certificates still pass).
tamper_subgraph5.log   the triangle certificate with SUBGRAPH 6 changed to 5 is rejected at check (c').

Second round (October 20, evening), after the spanq error was found and fixed.
check_generator_tools_full_run.log   validation/check_generator_tools.py, all steps: all ten certificates verify; the image
                                     test, the torus span and the two-tori certificate are reproduced; spanq agrees with
                                     orderq on the 14 tori and on 10, 300 and 3000 random tori; the exact check passes.  Two
                                     checks failed because they compared exactsolve's answer byte for byte with the saved one,
                                     and exactsolve now chooses its dropped cells mod a large prime, which gives a different but
                                     equally valid cycle (39,986 instead of 40,839 cells).  The step was changed to check
                                     validity instead (build the certificate and verify it).
check_generator_tools_exact_rerun.log the corrected step: 5 of 5 checks pass (the new y gives a verified certificate).
regression_vs_reference.log          the engine (after comments were added) against oconf26: ALL OK.
