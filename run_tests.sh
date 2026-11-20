#!/bin/bash
# Rerun every simulation test in order, then redraw the figures.
# Run from the repository root. Each step overwrites its JSON file in results/.
set -e
cd "$(dirname "$0")"
python arm_margin_test.py --tag pilot          # Test 1 pilot (100 runs, gain 60)
python run_all.py                               # Test 1 Parts A and B (pass worker count as argument, default 4)
python run_diagnostic.py                        # Test 1 Part C
python test2_resolution.py                      # Test 2
python test3_realistic.py                       # Test 3
python test3b_followup.py                       # Test 3b
python test4_rule_check.py                      # Test 4
python test5_confirm.py                         # Test 5 (+ traces for Fig. 2)
python test6_settling.py                        # Test 6
python test7_two_link.py                        # Test 7 (+ traces for Fig. 5)
python test8_diagnose.py                        # Test 8
python test9_fix.py                             # Test 9
python make_figures.py                          # Figures 1-5
echo "done: results/ and figures/ updated"
