# Makefile for graph-braid-torsion.
#
# Nothing here is compiled: the package is pure Python.  The Makefile only collects the
# commands you are most likely to want, so that you do not have to remember them.
# Run them from the top folder of the repository, e.g.  `make check`.
#
#   make help        list the targets below
#   make examples    run the four tutorial scripts in examples/ (a few seconds)
#   make check       the fast validation suite (about 5 minutes; see TUTORIAL.md, Lesson 12)
#   make check-full  everything in `check`, plus the full two-route recomputation of the
#                    headline groups and the first million-cell regression case (20+ minutes)
#   make tables      print the summary tables from the stored data (instant)
#   make inventory   rebuild data/inventory.csv from the stored data (seconds; needs networkx)
#   make hunt        resume the large search torsion_hunt.py (first-phase families; hours;
#                    needs numpy and nauty; see TUTORIAL.md, Lesson 11)
#   make logs        run `check` and save every output in logs/ (the logs shipped with the
#                    package were made this way)
#   make clean       delete Python caches and the template's scratch output
#
# Use a different Python with  `make check PYTHON=python3.11`.

PYTHON ?= python3
EXP = experiments

.PHONY: help examples check check-full tables inventory hunt logs clean

help:
	@sed -n '1,20p' Makefile | sed 's/^# \{0,1\}//'

examples:
	$(PYTHON) examples/compute.py bipartite:4,5 3
	$(PYTHON) examples/look_at_the_complex.py --max-cells 3
	$(PYTHON) examples/explicit_torsion_class.py
	$(PYTHON) examples/sweep_template.py --out examples/output/make_examples.json
	rm -rf examples/output

check:
	$(PYTHON) $(EXP)/model_sanity.py
	$(PYTHON) $(EXP)/validate.py
	$(PYTHON) $(EXP)/wilson_block.py
	$(PYTHON) $(EXP)/wilson_cokernel.py
	$(PYTHON) $(EXP)/crosscheck_mod_p.py
	$(PYTHON) $(EXP)/verify_headlines.py --skip-large

check-full: check
	$(PYTHON) $(EXP)/verify_headlines.py
	$(PYTHON) $(EXP)/validate_slow.py 0

tables:
	$(PYTHON) $(EXP)/make_tables.py

inventory:
	$(PYTHON) $(EXP)/inventory.py

hunt:
	$(PYTHON) -u $(EXP)/torsion_hunt.py

logs:
	mkdir -p logs
	$(PYTHON) $(EXP)/model_sanity.py       > logs/model_sanity.txt 2>&1
	$(PYTHON) $(EXP)/validate.py           > logs/validate.txt 2>&1
	$(PYTHON) $(EXP)/wilson_block.py       > logs/wilson_block.txt 2>&1
	$(PYTHON) $(EXP)/wilson_cokernel.py    > logs/wilson_cokernel.txt 2>&1
	$(PYTHON) $(EXP)/crosscheck_mod_p.py   > logs/crosscheck_mod_p.txt 2>&1
	$(PYTHON) $(EXP)/verify_headlines.py   > logs/verify_headlines.txt 2>&1
	$(PYTHON) $(EXP)/inventory.py          > logs/inventory.txt 2>&1
	$(PYTHON) $(EXP)/prime_census.py --skip-compute > logs/prime_census.txt 2>&1
	@echo "logs written to logs/"

clean:
	find . -name __pycache__ -prune -exec rm -rf {} +
	rm -rf examples/output
