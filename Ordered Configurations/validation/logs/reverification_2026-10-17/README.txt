Re-verification of the main computations, October 16-17, 2026.

Everything here was run from a fresh build of the delivered package (tarball unpacked into a new folder,
then `make`). The summary table is data/accounting/reverification.csv (also a sheet of
data/accounting/accounting.xlsx).

major_cases.jsonl            run_verify.py and run_verify2.py: exact, field and local computations for K4
                             (n = 2..6), the triangle graphs, the doubled K4, the tripled and doubled C4,
                             T(2,2,3), and the two survey cases g4_78 and g4_80. One entry (T(1,3,3) over
                             F_2) failed because another job was using the memory at the same time; it was
                             rerun successfully (see the next file).
certificates_fields_rankchecks_part1.log
                             T(2,2,3) over F_1000003, F_2, F_3; doubled K4 over F_2; tripled C4 over F_3;
                             verification of the T(1,3,3) certificate; the two new certificates (doubled K4,
                             tripled C4) and their verification; T(1,3,3) over F_2; the first 16 rank checks
                             of campaign runs (nonplanar_campaign.py --rational).
rankchecks_part2_and_exact.log
                             two more rank checks, one out-of-memory, and the exact integral homology of
                             Conf_5 of the doubled K4 and of the tripled C4.
regression_vs_v26.log        engine/oconf against history/engine_versions/oconf26.cpp: ALL OK.
analysis_tools.log           smoke tests of the analysis tools on T(2,2,2).
enum_check.py                independent networkx enumeration of the multigraphs of genus 2, 3, 4.
