# Makefile for the gconf-torsion package.
#
#   make            build the engine, the certificate finders, the exact tools and the analysis tools
#   make engine     build only engine/oconf (this is all that most users need)
#   make check      run the fast validation scripts (a few minutes in total)
#   make clean      remove the compiled binaries
#
# Requirements: a C++17 compiler with OpenMP (g++ 9 or later works; on macOS use g++ from Homebrew or
# clang with libomp), Python 3.8 or later, and the Python package python-flint (pip install python-flint).
# The flag -march=native tunes the code for the machine it is compiled on; remove it if the binaries
# have to run on a different machine.

CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -fopenmp -std=c++17

ANALYSIS = analysis/prod2 analysis/g2tors analysis/g2split analysis/phi analysis/decomp \
           analysis/cutimg3 analysis/vexp
CERTS    = certificates/cert certificates/certp certificates/certq certificates/orderq certificates/exactsolve certificates/reducecyc \
           certificates/reducemod certificates/spanq

all: engine $(CERTS) $(ANALYSIS)

engine: engine/oconf

engine/oconf: engine/oconf.cpp engine/gr42.h
	$(CXX) $(CXXFLAGS) -o $@ engine/oconf.cpp

certificates/%: certificates/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# the prime-power and exact tools share three headers
certificates/certq certificates/orderq certificates/reducemod: certificates/swiat_cells.h certificates/elim_modq.h
certificates/spanq: certificates/elim_modq.h
certificates/exactsolve certificates/reducecyc: certificates/swiat_cells.h certificates/elim_z.h certificates/elim_modq.h

analysis/%: analysis/%.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

check: all
	python3 validation/check_explosion.py
	python3 validation/elimination_vs_flint.py
	python3 validation/abrams_vs_swiatkowski.py
	python3 validation/dense_phase_suite.py
	python3 validation/local_smith_suite.py --quick
	python3 validation/check_quotients.py
	python3 validation/verify_cert.py certificates/cert_T123_n6.txt.gz
	python3 validation/verify_cert.py certificates/cert_K4dbl5_n5_order16.txt.gz
	python3 validation/verify_cert.py certificates/cert_triangle_leaf_into_K4dbl5_order16.txt.gz

clean:
	rm -f engine/oconf $(CERTS) $(ANALYSIS)

.PHONY: all engine check clean
