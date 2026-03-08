# Systematic Audit: pArrayBasicOperations

Target Folder: `D:\OneDrive\Documents\GitHub\phasorArray_Toolbox\Fonctions\pArrayBasicOperations`
Generated on: 2026-03-08

## Results Summary

| Category | Count | Action Recommendation |
|----------|-------|-----------------------|
| 🟢 Essential | 42 | Keep |
| 🟠 Utility / Plot | 12 | Keep |
| 🔴 Redundant | 15 | DELETE |
| ⚪ Placeholder | 10 | DELETE (Folder +shift2pi) |

---

## Root Folder (pArrayBasicOperations)

| File | Size | LOC | Status | Dependencies | Recommendation |
|:-----|:-----|:----|:-------|:-------------|:---------------|
| [`analyzeMonotonicRegions.m`](../Fonctions/pArrayBasicOperations/analyzeMonotonicRegions.m) | 0B | 0 | 🔴 Redundant | findMonotonicRegionAndFirstRevolution.m | **DELETE** |
| [`array2BToepliz.m`](../Fonctions/pArrayBasicOperations/array2BToepliz.m) | 3.3KB | 101 | 🟢 Essential | PhasorArray.m | Keep |
| [`Array2TBHankel.m`](../Fonctions/pArrayBasicOperations/Array2TBHankel.m) | 3.6KB | 116 | 🟢 Essential | PhasorArray.m | Keep |
| [`array2TBlocks.m`](../Fonctions/pArrayBasicOperations/array2TBlocks.m) | 3.5KB | 117 | 🟢 Essential | PhasorArray.m, LyapHarmonic.m, PR_v1.m, PR_v2.m, PR_v3.m... | Keep |
| [`array2TBlocks2.m`](../Fonctions/pArrayBasicOperations/array2TBlocks2.m) | 2.1KB | 75 | 🔴 Redundant | None | **DELETE** |
| [`array2TFTB.m`](../Fonctions/pArrayBasicOperations/array2TFTB.m) | 1.9KB | 60 | 🟢 Essential | PhasorArrayTimes2.m | Keep |
| [`better_lyap.m`](../Fonctions/pArrayBasicOperations/better_lyap.m) | 0B | 0 | 🔴 Redundant | None | **DELETE** |
| [`BT2array.m`](../Fonctions/pArrayBasicOperations/BT2array.m) | 1.9KB | 60 | 🟢 Essential | None | Keep |
| [`BToeplitz.m`](../Fonctions/pArrayBasicOperations/BToeplitz.m) | 2.3KB | 57 | 🟢 Essential | None | Keep |
| [`BT_2_TB.m`](../Fonctions/pArrayBasicOperations/BT_2_TB.m) | 992B | 28 | 🟢 Essential | shuffle_matrix.m | Keep |
| [`compare_shift2pi_versions.m`](../Fonctions/pArrayBasicOperations/compare_shift2pi_versions.m) | 0B | 0 | 🔴 Redundant | None | **DELETE** |
| [`Copy_of_PhasorArrayAdd.m`](../Fonctions/pArrayBasicOperations/Copy_of_PhasorArrayAdd.m) | 2.8KB | 88 | 🔴 Redundant | None | **DELETE** |
| [`Copy_of_PhasorUnif.m`](../Fonctions/pArrayBasicOperations/Copy_of_PhasorUnif.m) | 1.9KB | 53 | 🔴 Redundant | None | **DELETE** |
| [`dephase.m`](../Fonctions/pArrayBasicOperations/dephase.m) | 560B | 20 | 🟢 Essential | GettingStarted.m, BasicToolbox.m, PhasorArray.m | Keep |
| [`detLeibnizHmc.m`](../Fonctions/pArrayBasicOperations/detLeibnizHmc.m) | 4.6KB | 149 | 🟢 Essential | PhasorArrayTest.m, test_PhasorArray_advanced.m, test_PhasorArray_basic.m | Keep |
| [`F_tb_2_PhasorArray.m`](../Fonctions/pArrayBasicOperations/F_tb_2_PhasorArray.m) | 1.8KB | 41 | 🟢 Essential | PhasorArray.m, PhasorArrayTimes2.m, RicHarmonicKlein.m | Keep |
| [`hmq_sim.m`](../Fonctions/pArrayBasicOperations/hmq_sim.m) | 12.8KB | 385 | 🟢 Essential | Exemple_Toolbox_LMI.m, PhasorArray.m, FloquetDec.m, TransitionMatrixOverT.m | Keep |
| [`isphasor.m`](../Fonctions/pArrayBasicOperations/isphasor.m) | 878B | 27 | 🟢 Essential | N_tb.m, PhasorArrayPad.m, PhasorArrayPad2.m | Keep |
| [`isrealp.m`](../Fonctions/pArrayBasicOperations/isrealp.m) | 1.8KB | 50 | 🟢 Essential | PhasorArray2time.m | Keep |
| [`isToeplitz.m`](../Fonctions/pArrayBasicOperations/isToeplitz.m) | 1.2KB | 34 | 🟠 Utility | None | Keep |
| [`LyapHarmonic.m`](../Fonctions/pArrayBasicOperations/LyapHarmonic.m) | 1.8KB | 62 | 🟢 Essential | None | Keep |
| [`N_bt.m`](../Fonctions/pArrayBasicOperations/N_bt.m) | 388B | 22 | 🟢 Essential | GettingStarted.m, BasicToolbox.m, array2BToepliz.m, spN_bt.m | Keep |
| [`N_tb.m`](../Fonctions/pArrayBasicOperations/N_tb.m) | 986B | 46 | 🟢 Essential | GettingStarted.m, BasicToolbox.m, ECC_ex.m, Exemple_Toolbox_LMI.m, SPMSM_template.m... | Keep |
| [`pagekron.m`](../Fonctions/pArrayBasicOperations/pagekron.m) | 676B | 27 | 🟢 Essential | None | Keep |
| [`Phasor2CosSin.m`](../Fonctions/pArrayBasicOperations/Phasor2CosSin.m) | 2.8KB | 82 | 🟢 Essential | None | Keep |
| [`PhasorArray2time.m`](../Fonctions/pArrayBasicOperations/PhasorArray2time.m) | 26.3KB | 748 | 🟢 Essential | expm.m, PhasorArray.m, hmq_sim.m, PhasorDet.m, PhasorInv.m... | Keep |
| [`PhasorArrayAdd.m`](../Fonctions/pArrayBasicOperations/PhasorArrayAdd.m) | 2.1KB | 53 | 🟢 Essential | PhasorArray.m, Copy_of_PhasorArrayAdd.m, PhasorArrayBilAdd.m | Keep |
| [`PhasorArrayBilAdd.m`](../Fonctions/pArrayBasicOperations/PhasorArrayBilAdd.m) | 1.1KB | 54 | 🔴 Redundant | None | **DELETE** |
| [`PhasorArrayKron.m`](../Fonctions/pArrayBasicOperations/PhasorArrayKron.m) | 388B | 32 | 🟢 Essential | PhasorArray.m | Keep |
| [`PhasorArrayOplus.m`](../Fonctions/pArrayBasicOperations/PhasorArrayOplus.m) | 1.9KB | 84 | 🟢 Essential | PhasorArray.m | Keep |
| [`PhasorArrayPad.m`](../Fonctions/pArrayBasicOperations/PhasorArrayPad.m) | 1.5KB | 60 | 🔴 Redundant | PhasorArray.m, array2TFTB.m, PhasorArrayPad2.m | **DELETE** |
| [`PhasorArrayPad2.m`](../Fonctions/pArrayBasicOperations/PhasorArrayPad2.m) | 1.5KB | 47 | 🔴 Redundant | None | **DELETE** |
| [`PhasorArrayTimes.m`](../Fonctions/pArrayBasicOperations/PhasorArrayTimes.m) | 2.6KB | 110 | 🟢 Essential | PhasorArray.m, PhasorArrayKron.m, PhasorArrayTimes2.m, PhasorInv.m | Keep |
| [`PhasorArrayTimes2.m`](../Fonctions/pArrayBasicOperations/PhasorArrayTimes2.m) | 4.4KB | 168 | 🟢 Essential | PhasorArrayTimes.m | Keep |
| [`PhasorDet.m`](../Fonctions/pArrayBasicOperations/PhasorDet.m) | 5.4KB | 165 | 🟢 Essential | PhasorArray.m, detLeibnizHmc.m | Keep |
| [`PhasorInv.m`](../Fonctions/pArrayBasicOperations/PhasorInv.m) | 8.2KB | 222 | 🟢 Essential | PhasorArray.m | Keep |
| [`phasorPad.m`](../Fonctions/pArrayBasicOperations/phasorPad.m) | 6.6KB | 212 | 🟢 Essential | quick_test.m, run_tests.m, run_validation_mcp.m, PhasorArray.m, array2TBlocks2.m... | Keep |
| [`PhasorPow.m`](../Fonctions/pArrayBasicOperations/PhasorPow.m) | 2.4KB | 83 | 🟢 Essential | PhasorArray.m | Keep |
| [`PhasorUnif.m`](../Fonctions/pArrayBasicOperations/PhasorUnif.m) | 2.1KB | 50 | 🟢 Essential | PhasorArray.m, Copy_of_PhasorUnif.m, PhasorArrayAdd.m | Keep |
| [`PhasorVec2plot.m`](../Fonctions/pArrayBasicOperations/PhasorVec2plot.m) | 217B | 7 | 🔴 Redundant | None | **DELETE** |
| [`PosPart2PhasorArray.m`](../Fonctions/pArrayBasicOperations/PosPart2PhasorArray.m) | 1.3KB | 54 | 🟢 Essential | GettingStarted.m, BasicToolbox.m, PhasorArray.m | Keep |
| [`pr.m`](../Fonctions/pArrayBasicOperations/pr.m) | 1.0KB | 38 | 🟢 Essential | BasicToolbox.m, Exemple_Toolbox_LMI.m, PR_v1.m, PR_v2.m | Keep |
| [`PR_In.m`](../Fonctions/pArrayBasicOperations/PR_In.m) | 1.1KB | 40 | 🟢 Essential | LyapHarmonic.m, SylvHarmonic.m | Keep |
| [`PR_v1.m`](../Fonctions/pArrayBasicOperations/PR_v1.m) | 920B | 33 | 🔴 Redundant | None | **DELETE** |
| [`PR_v2.m`](../Fonctions/pArrayBasicOperations/PR_v2.m) | 903B | 33 | 🔴 Redundant | None | **DELETE** |
| [`PR_v3.m`](../Fonctions/pArrayBasicOperations/PR_v3.m) | 951B | 33 | 🔴 Redundant | None | **DELETE** |
| [`pvalue.m`](../Fonctions/pArrayBasicOperations/pvalue.m) | 430B | 15 | 🟢 Essential | PhasorArray.m | Keep |
| [`rand_phasor.m`](../Fonctions/pArrayBasicOperations/rand_phasor.m) | 5.4KB | 129 | 🟢 Essential | GettingStarted.m, BasicToolbox.m, Exemple_Toolbox_LMI.m, PhasorArray.m | Keep |
| [`ReduceArray.m`](../Fonctions/pArrayBasicOperations/ReduceArray.m) | 3.4KB | 91 | 🟢 Essential | expm.m, PhasorArray.m, PhasorArrayTimes.m, PhasorArrayTimes2.m, PhasorDet.m... | Keep |
| [`RicHarmonicKlein.m`](../Fonctions/pArrayBasicOperations/RicHarmonicKlein.m) | 2.8KB | 104 | 🟢 Essential | ECC_ex.m, truncationAdvisor.m | Keep |
| [`ScalarPhasorArray.m`](../Fonctions/pArrayBasicOperations/ScalarPhasorArray.m) | 2.3KB | 59 | 🟢 Essential | PhasorArrayTest.m, GettingStarted.m, BasicToolbox.m, SPMSM_template.m, PhasorArray.m... | Keep |
| [`shift2pi.m`](../Fonctions/pArrayBasicOperations/shift2pi.m) | 21.5KB | 543 | 🟢 Essential | angularsft.m, angularsft_v2.m, dSpaceDataExplorer.m, shift2pi_old.m | Keep |
| [`shift2pi_modular.m`](../Fonctions/pArrayBasicOperations/shift2pi_modular.m) | 0B | 0 | 🔴 Redundant | None | **DELETE** |
| [`shift2pi_old.m`](../Fonctions/pArrayBasicOperations/shift2pi_old.m) | 21.4KB | 542 | 🔴 Redundant | None | **DELETE** |
| [`shuffle_matrix.m`](../Fonctions/pArrayBasicOperations/shuffle_matrix.m) | 1.5KB | 32 | 🟢 Essential | BT_2_TB.m, TB_2_BT.m | Keep |
| [`spArray2TBHankel.m`](../Fonctions/pArrayBasicOperations/spArray2TBHankel.m) | 2.7KB | 102 | 🟢 Essential | PhasorArray.m, Array2TBHankel.m | Keep |
| [`sparray2TBlocks.m`](../Fonctions/pArrayBasicOperations/sparray2TBlocks.m) | 2.2KB | 80 | 🟢 Essential | PhasorArray.m, PhasorArrayTimes2.m, pr.m, PR_In.m, spPR_In.m | Keep |
| [`spN_bt.m`](../Fonctions/pArrayBasicOperations/spN_bt.m) | 412B | 22 | 🟢 Essential | GettingStarted.m, BasicToolbox.m | Keep |
| [`spN_tb.m`](../Fonctions/pArrayBasicOperations/spN_tb.m) | 493B | 22 | 🟢 Essential | GettingStarted.m, BasicToolbox.m | Keep |
| [`spPR_In.m`](../Fonctions/pArrayBasicOperations/spPR_In.m) | 1.2KB | 43 | 🟢 Essential | None | Keep |
| [`stemPhasor.m`](../Fonctions/pArrayBasicOperations/stemPhasor.m) | 9.3KB | 248 | 🟠 Utility | PhasorArray.m, PhasorArray2time.m | Keep |
| [`SylvHarmonic.m`](../Fonctions/pArrayBasicOperations/SylvHarmonic.m) | 2.4KB | 74 | 🟢 Essential | GettingStarted.m, BasicToolbox.m, PhasorArray.m, testCOne.m | Keep |
| [`S_bt.m`](../Fonctions/pArrayBasicOperations/S_bt.m) | 265B | 12 | 🟢 Essential | None | Keep |
| [`S_tb.m`](../Fonctions/pArrayBasicOperations/S_tb.m) | 240B | 8 | 🟢 Essential | PhasorArrayTest.m, test_PhasorArray_advanced.m | Keep |
| [`TB2array.m`](../Fonctions/pArrayBasicOperations/TB2array.m) | 2.9KB | 79 | 🟢 Essential | GettingStarted.m, BasicToolbox.m, PhasorArray.m, rand_phasor.m | Keep |
| [`TB_2_BT.m`](../Fonctions/pArrayBasicOperations/TB_2_BT.m) | 1.8KB | 37 | 🟢 Essential | shuffle_matrix.m | Keep |
| [`TFTB_2_array.m`](../Fonctions/pArrayBasicOperations/TFTB_2_array.m) | 2.9KB | 79 | 🟢 Essential | PhasorArray.m, testCOne.m | Keep |
| [`TimeArray2Phasors.m`](../Fonctions/pArrayBasicOperations/TimeArray2Phasors.m) | 4.3KB | 108 | 🟢 Essential | expm.m, PhasorArray.m, PhasorArray2time.m, PhasorDet.m, PhasorInv.m... | Keep |
| [`TimeArray2plot.m`](../Fonctions/pArrayBasicOperations/TimeArray2plot.m) | 2.7KB | 75 | 🟠 Utility | None | Keep |

---

## Folder: +shift2pi

| File | Status | Dependencies | Recommendation |
|:-----|:-------|:-------------|:---------------|
| [`clampToValidRange.m`](../Fonctions/pArrayBasicOperations/+shift2pi/clampToValidRange.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |
| [`computeInterpolationCoefficients.m`](../Fonctions/pArrayBasicOperations/+shift2pi/computeInterpolationCoefficients.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |
| [`computeShiftIndexMapping.m`](../Fonctions/pArrayBasicOperations/+shift2pi/computeShiftIndexMapping.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |
| [`findBracketingInterval.m`](../Fonctions/pArrayBasicOperations/+shift2pi/findBracketingInterval.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |
| [`findMonotonicRegionStart.m`](../Fonctions/pArrayBasicOperations/+shift2pi/findMonotonicRegionStart.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |
| [`generateShiftedOutputs.m`](../Fonctions/pArrayBasicOperations/+shift2pi/generateShiftedOutputs.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |
| [`main.m`](../Fonctions/pArrayBasicOperations/+shift2pi/main.m) | ⚪ Placeholder | PhasorArray.m, dSpaceDataExplorer.m... | **DELETE** |
| [`preprocessInputs.m`](../Fonctions/pArrayBasicOperations/+shift2pi/preprocessInputs.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |
| [`README.md`](../Fonctions/pArrayBasicOperations/+shift2pi/README.md) | ⚪ Placeholder | — | **DELETE** |
| [`validateInputDimensions.m`](../Fonctions/pArrayBasicOperations/+shift2pi/validateInputDimensions.m) | ⚪ Placeholder | shift2pi.m, shift2pi_old.m | **DELETE** |

