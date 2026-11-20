# Simulation test 1: observer-aware margin

Run 5 November 2026. Criteria are fixed in arm_margin_test.py. Single-link arm, simulation only, not a proof.

## Part A: fixed test (1,000 runs, observer gain 60, ideal sensor, guaranteed margin)

| Case | Runs breaking velocity limit | 95% CI (%) | Median RMS, t>=2 s (rad) | 95th pct RMS | Peak torque, median (N m) | Median energy |
|---|---|---|---|---|---|---|
| no margin | 0/1000 | 0.00-0.38 | 0.0013 | 0.0019 | 8.00 | 6.74 |
| constant margin | 0/1000 | 0.00-0.38 | 0.1528 | 0.2646 | 8.00 | 5.52 |
| time-varying margin | 0/1000 | 0.00-0.38 | 0.0013 | 0.0094 | 8.00 | 6.57 |

Margin: m0 = 0.554, final = 0.095 rad/s, settling 0.021 s. Guaranteed disturbance bound D = 28.8; largest seen = 8.7.

Verdict: C1 (motive) FAILED. C2 (safety) passed, but trivially. C3 passed on median RMS, although the 95th percentile rose from 0.0019 to 0.0094 rad. C4 passed.

## Part B: exploratory, added after the pilot (300 runs, clearance 1%, guaranteed margin)

| Gain | Encoder | Case | Violating runs | Estimate outside limits | Worst excess (rad/s) | Median RMS | Median saturated time (s of 6) |
|---|---|---|---|---|---|---|---|
| 10 | ideal | no margin | 6 | 0 | 0.047 | 0.0017 | 0.05 |
| 10 | 12-bit | no margin | 6 | 0 | 0.049 | 0.0023 | 0.05 |
| 20 | ideal | no margin | 0 | 0 | 0.000 | 0.0010 | 0.04 |
| 20 | 12-bit | no margin | 0 | 0 | 0.000 | 0.0009 | 0.06 |
| 40 | ideal | no margin | 0 | 0 | 0.000 | 0.0012 | 0.04 |
| 40 | ideal | constant margin | 0 | 0 | 0.000 | 0.1252 | 0.01 |
| 40 | ideal | time-varying margin | 0 | 0 | 0.000 | 0.0013 | 0.03 |
| 40 | 12-bit | no margin | 212 | 211 | 0.847 | 0.1348 | 4.93 |
| 40 | 12-bit | constant margin | 280 | 300 | 0.830 | 0.4313 | 5.26 |
| 40 | 12-bit | time-varying margin | 265 | 299 | 0.721 | 0.2647 | 5.07 |
| 60 | ideal | no margin | 0 | 0 | 0.000 | 0.0013 | 0.05 |
| 60 | ideal | constant margin | 0 | 0 | 0.000 | 0.1215 | 0.02 |
| 60 | ideal | time-varying margin | 0 | 0 | 0.000 | 0.0013 | 0.04 |
| 60 | 12-bit | no margin | 291 | 292 | 1.023 | 0.5196 | 5.81 |
| 60 | 12-bit | constant margin | 207 | 300 | 1.118 | 1.3017 | 5.88 |
| 60 | 12-bit | time-varying margin | 293 | 297 | 1.108 | 0.5908 | 5.84 |

Guaranteed margin did not fit at gains 10 and 20.

## Part C: exploratory diagnostic (300 runs, clearance 1%, margin uses D = 8, set from the pilot maximum of about 7.4 and not guaranteed; the 1,000-run maximum was later found to be 8.7)

| Gain | Encoder | Case | Violating runs | Estimate outside limits | Worst excess | Margin final value | Median RMS | Saturated time | First violation (s) |
|---|---|---|---|---|---|---|---|---|---|
| 10 | ideal | constant margin | 0 | 0 | 0.000 | 0.532 | 0.1571 | 0.01 | - |
| 10 | ideal | no margin | 6 | 0 | 0.047 | 0.532 | 0.0017 | 0.05 | 0.028-0.050 |
| 10 | ideal | time-varying margin | 0 | 0 | 0.000 | 0.532 | 0.1038 | 0.01 | - |
| 10 | 12-bit | constant margin | 0 | 0 | 0.000 | 0.532 | 0.1547 | 0.01 | - |
| 10 | 12-bit | no margin | 6 | 0 | 0.049 | 0.532 | 0.0023 | 0.05 | 0.028-0.050 |
| 10 | 12-bit | time-varying margin | 0 | 0 | 0.000 | 0.532 | 0.1007 | 0.01 | - |
| 15 | ideal | constant margin | 0 | 0 | 0.000 | 0.229 | 0.1226 | 0.01 | - |
| 15 | ideal | no margin | 0 | 0 | 0.000 | 0.229 | 0.0010 | 0.04 | - |
| 15 | ideal | time-varying margin | 0 | 0 | 0.000 | 0.229 | 0.0015 | 0.01 | - |
| 15 | 12-bit | constant margin | 0 | 0 | 0.000 | 0.229 | 0.1185 | 0.01 | - |
| 15 | 12-bit | no margin | 0 | 0 | 0.000 | 0.229 | 0.0012 | 0.04 | - |
| 15 | 12-bit | time-varying margin | 0 | 0 | 0.000 | 0.229 | 0.0019 | 0.01 | - |
| 20 | ideal | constant margin | 0 | 0 | 0.000 | 0.127 | 0.1199 | 0.01 | - |
| 20 | ideal | no margin | 0 | 0 | 0.000 | 0.127 | 0.0010 | 0.04 | - |
| 20 | ideal | time-varying margin | 0 | 0 | 0.000 | 0.127 | 0.0010 | 0.01 | - |
| 20 | 12-bit | constant margin | 0 | 0 | 0.000 | 0.127 | 0.1140 | 0.03 | - |
| 20 | 12-bit | no margin | 0 | 0 | 0.000 | 0.127 | 0.0009 | 0.06 | - |
| 20 | 12-bit | time-varying margin | 0 | 0 | 0.000 | 0.127 | 0.0009 | 0.04 | - |
| 40 | 16-bit | no margin | 0 | 0 | 0.000 | 0.031 | 0.0010 | 0.05 | - |
| 40 | 12-bit | no margin | 212 | 211 | 0.847 | 0.031 | 0.1348 | 4.93 | 0.231-5.494 |
| 60 | 16-bit | no margin | 0 | 0 | 0.000 | 0.014 | 0.0070 | 3.67 | - |
| 60 | 12-bit | no margin | 291 | 292 | 1.023 | 0.014 | 0.5196 | 5.81 | 0.213-3.842 |

## Limits of this test

- One plant, one controller family and one scenario family, all chosen in advance.
- The margin comes from a sampled worst case, not a proved bound.
- Starts close to a limit with large velocity uncertainty were not tested.
- The encoder breakdown may depend on this controller cancelling the observer injection term; filtering was not tested.

# Test 2: resolution test

Rules are fixed in test2_resolution.py. Raw values are in test2.json. 12-bit encoder for Part F (200 runs, one seed, 4 s); Part L used 60 runs per setting.

| Gain | Design | Violating runs | Estimate outside limits | Saturated share | Median RMS |
|---|---|---|---|---|---|
| 20 | base | 0 | 0 | 1.5% | 0.0009 |
| 20 | no cancel | 2 | 1 | 1.7% | 0.0029 |
| 20 | low-pass | 1 | 1 | 1.1% | 0.0009 |
| 20 | dead zone | 1 | 0 | 1.1% | 0.0011 |
| 20 | linear | 29 | 0 | 1.3% | 0.0106 |
| 40 | base | 123 | 122 | 82.2% | 0.0840 |
| 40 | no cancel | 2 | 2 | 1.7% | 0.0019 |
| 40 | low-pass | 0 | 1 | 2.1% | 0.0015 |
| 40 | dead zone | 0 | 0 | 4.7% | 0.0011 |
| 40 | linear | 4 | 0 | 1.2% | 0.0043 |
| 60 | base | 187 | 187 | 96.3% | 0.4779 |
| 60 | no cancel | 3 | 20 | 1.7% | 0.0026 |
| 60 | low-pass | 23 | 93 | 11.2% | 0.0034 |
| 60 | dead zone | 0 | 0 | 5.7% | 0.0012 |
| 60 | linear | 0 | 0 | 1.2% | 0.0024 |

Ceilings (lowest gain with median saturation above 10%): fractional 20, 28, 40, 57 for 10, 12, 14, 16 bits, fitted slope -0.252 (predicted -0.25). Linear 57 and 113 for 10 and 12 bits, none in range above that; two points, so the -0.50 prediction is not properly tested.

Artefact rule: not triggered, narrowly. No cancel missed at gain 60 (1.5% violating) and marginally on RMS at gain 40; low-pass missed at gain 60. The dead zone, not covered by the rule, removed all violations at gains 40 and 60.

Limits: continuous-time style update at 5 kHz, no sensor noise, one seed, coarse gain grid, one controller design.

# Test 3: realistic digital implementation

Rules are fixed in test3_realistic.py. Raw values are in test3.json. 1 kHz control, one-sample delay, 12-bit encoder, 0.5 count noise (S2), seeds 7, 8, 9, 100 runs each.

| Rule | Outcome | Detail |
|---|---|---|
| R1 dead zone robust at gains 40 and 57 | FAILED | Violating runs: one count 14.7% and 48.3%; two counts 15.7% and 41.0% |
| R2 fixed-time beats linear on test seeds | PASSED, with a caveat | One-count dead zone at gain 20 vs linear at gain 57: velocity-estimate RMS 0.0262 vs 0.0484, tracking RMS 0.00099 vs 0.00247. But 4/200 vs 1/200 violating runs, worst excess 0.017 vs 0.042 |
| R3 small motions | FAILED narrowly | Passed at gains 28 and 40; at gain 57 tracking ratio was 1.68 (slow) and 1.57 (hold). Peak-to-peak hold error passed at all three (at most 3.86 counts). Gain 20 was not tested here |
| R4 gain law slope | PASSED | -0.251 from ceilings 14, 20, 28, 40 for 10, 12, 14, 16 bits |

Bias and drift (S3): one-count dead zone at gain 28 gave 18/300 violating runs and tracking RMS 0.0019; bias raises tracking error directly in every design.

## Follow-up (test3b.json), added after seeing Test 3; exploratory

Torque dependence, prediction written before the run: lam_c = sqrt(b0 u_max / q^a2) = 14.3, 20.2, 28.6 for 4, 8, 16 N m. Measured ceilings 14, 20, 28 on a coarse grid. With no fitted constant the same formula gives 14.3, 20.2, 28.6, 40.4 for 10 to 16 bits (measured 14, 20, 28, 40); the c = 1 reading was noticed after Test 3. b0 and a2 were not varied, the 10% threshold is a choice, and Test 2 ceilings without noise or delay were about 1.4x higher.

Clearance 5% vs 1%, violating runs under S2:

| Design and gain | 1% clearance | 5% clearance | Tracking RMS (5%) | Velocity-estimate RMS (5%) | Torque variation (5%) |
|---|---|---|---|---|---|
| Low-pass, 20 | 2.7% | 0.0% | 0.0011 | 0.0331 | 1183 |
| Dead zone 1, 20 | 1.3% | 0.3% | 0.0010 | 0.0257 | 1408 |
| Dead zone 1, 40 | 14.7% | 0.3% | 0.0014 | 0.0372 | 1867 |
| Dead zone 1, 57 | 48.3% | 1.0% | 0.0019 | 0.0643 | 2116 |
| Low-pass, 40 | 55.7% | 0.3% | 0.0012 | 0.0637 | 3810 |
| Linear, 57 | 0.3% | 0.3% | 0.0024 | 0.0477 | 3729 |

These 5% rows pool all seeds and were not selected on separate data. No-cancel design under S2: 1.0% to 2.7% violating at gains 14 to 28, but worst excess above 0.5 rad/s and tracking RMS about three times higher.

Limits: one plant, one controller family, simulated sensors, no real hardware, no proof.

# Test 4: predictions of the derived rule

Predictions and the pass rule are fixed in test4_rule_check.py; raw values are in test4.json. All nine ceilings fell within 15% of prediction. Seven crossings fell within the 8% grid step ending at the prediction; sigma = 1 and 2 counts gave 7% and 14% shortfalls. The grid was centred on each prediction, so agreement is within one step, not exact. The derivation is in Gain_rule_derivation.docx.

# Test 5: confirmation on unseen seeds

Settings were frozen in test5_confirm.py before the run. Seeds 20 to 24, 500 runs per row, 5% clearance, 12-bit encoder. Raw values are in test5.json, and time responses in traces.json. Low-pass at gain 20: 0/500 violations under S2 and S3, tracking RMS 0.00103 (S2) and 0.00198 (S3). Dead zone at gain 20: 0/500 under both. Linear at gain 57: 0/500 under both, tracking RMS 0.00252. Base at gain 40: 49/500 (S2) and 55/500 (S3), 90% saturated.

# Test 6: observer settling time

Undisturbed observer error, 16 directions per size, sizes 0.1 to 1,000,000. Longest settling time approached 13.3 at gain 1 and 0.665 at gain 20. Numerical only.

# Test 7: two-link arm

Rules fixed in test7_two_link.py. W1 passed: ceilings within one 8% step of 23.7 and 28.7. W2 FAILED: proposed design 18/500 (S2) and 20/500 (S3) velocity violations against a 1% target, with 11 angle violations; cause not diagnosed; not retuned. W3: no linear gain met the selection rule, a deviation from plan; the fewest-violation candidate gave 70/500 and higher tracking error. W4: full cancellation at twice the rule failed in 500/500 runs.

# Test 8: diagnosis of two-link violations (test8.json)

All 11 angle violations had a joint 2 reference past, or within 0.02 rad of, a limit: a fault in the Test 7 scenarios. All 18 velocity violations were in start-up (0.15 to 0.44 s), with no saturation at that instant and no pulse nearby. Estimate error at those instants was 0.11 to 0.17 rad/s against a clearance of 0.12 to 0.13 rad/s.

# Test 9: fixes on new seeds 40 to 44 (test9.json)

Rules fixed in test9_fix.py. References kept 0.1 rad inside limits. 5% clearance: 22/500 (S2), 17/500 (S3) velocity violations, 0 angle. Constant formula clearance: 0 violations, joint 2 tracking 0.00748 against 0.0064. Start-up clearance: 0/500 under S2 and S3, tracking 0.00531 and 0.00643. Clearance size and timing use values measured on seeds 30 to 34, so they are not guaranteed bounds.