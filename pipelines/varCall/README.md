This folder stores major program code used for WEScall varCall.

`${WK_DIR}/varCall/out/aux/union` stores the variants discovery results accumulated across all samples.

## Notes on `vt genotyping`:

Genotype likelihoods (`tmp_pls`) are calculated at method `vt/joint_genotyping_record.cpp:void JointGenotypingRecord::add_allele(..)`.

Afterwards, in `void JointGenotypingRecord::flush_sample`, `tmp_pls` is normalized and capped at `255`, the maximal value of phred score, stored as `p_pls`. `imax` is the genotype (coded as 0,1, or 2) with largest likelihood.

Genotype probabilities are calculated at method `vt/joint_genotyping_record.cpp:bcf1_t* JointGenotypingRecord::flush_variant` as `gp`.

$hello$ 

Potential bug: `vt/joint_genotyping_record.cpp:427`:

`p_ads[i] = ((tmp_ads[i] > 255) ? 255 : (uint8_t)tmp_ads[i]);`


