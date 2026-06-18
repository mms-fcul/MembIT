# Future performance optimization notes

This document records the performance state after adding and validating XTC trajectory support in MembIT. It is intended as a handoff note for the next developer who works on MembIT performance.

## Current status

MembIT now supports a production-style reduced-XTC workflow:

```text
treated centered trajectory
  -> reduced Protein + phosphate-marker XTC
  -> matching reduced GRO
  -> matching MembIT index
  -> validation against trusted reduced-PDB output
  -> profile timing
```

The exploratory direct full-system XTC route was not kept as the recommended workflow. The recommended path is to prepare an analysis-specific trajectory and matching index.

## Most recent profile summary

Profile timing was run for reduced PDB and reduced XTC inputs using:

```bash
FRAME_COUNTS="5 100 500" PROFILE_FORMATS="pdb xtc" ./15_profile_reduced_xtc.sh
```

The key results were:

| Frames | PDB total wall clock | XTC total wall clock | Interpretation |
| ---: | ---: | ---: | --- |
| 5 | 0.653 s | 2.356 s | XTC startup dominates small runs |
| 100 | 12.236 s | 13.183 s | PDB and XTC are near parity |
| 500 | 110.908 s | 64.954 s | XTC is substantially faster |

For 500 frames, the main XTC profile sections were:

| Section | Seconds | Percent |
| --- | ---: | ---: |
| total_wall_clock | 64.954 | 100.000 |
| frame_analysis_total | 46.656 | 71.830 |
| thickness_total | 46.655 | 71.830 |
| thickness_calc_half_z | 29.695 | 45.720 |
| trajectory_reader_frame | 18.273 | 28.130 |
| mda_populate_membit_atoms | 14.644 | 22.550 |
| thickness_top_calculation | 14.595 | 22.470 |
| trajectory_index_match_checks | 1.781 | 2.740 |

## Main conclusion

The XTC reader optimization worked. For production-scale reduced trajectories, XTC is now faster than PDB. At 500 frames, PDB spent 63.118 s in `trajectory_reader_frame`, while XTC spent 18.273 s.

The performance bottleneck has shifted away from file-format parsing and toward the actual thickness/deformation calculation.

## Highest-priority future target

The next optimization effort should focus on:

```text
thickness_calc_half_z
```

This section took 29.695 s for the 500-frame reduced-XTC run, or about 45.7% of the total runtime.

This is the most important target because improving it would benefit both PDB and XTC workflows.

## Secondary target

The next XTC-specific target is:

```text
mda_populate_membit_atoms
```

This took 14.644 s for 500 XTC frames, or about 22.5% of total runtime.

Possible directions:

- reduce Python-level per-atom update overhead;
- batch coordinate updates where possible;
- cache atom update records more aggressively;
- avoid repeated branching between Protein, CoI, and membrane atoms.

## Validation overhead target

The section:

```text
trajectory_index_match_checks
```

took 1.781 s for 500 frames. This is not the main bottleneck, but it is likely safe to reduce after validation.

Possible direction:

- run index consistency checks once at startup or only under a debug flag;
- keep full checks under `--diagnose-index`;
- avoid repeated per-frame validation in production mode.

## Important asymmetry to investigate

The profile showed a large top/bottom calculation asymmetry:

| Section | 500-frame XTC seconds |
| --- | ---: |
| thickness_top_calculation | 14.595 |
| thickness_bottom_calculation | 1.211 |

The same asymmetry appears in the PDB run, so it is not caused by XTC or MDAnalysis. It is probably in the analysis algorithm or in the geometry/group assignment of this system.

Future profiling should count, per frame and leaflet:

- number of CoI atoms assigned to each leaflet;
- number of membrane marker atoms considered per leaflet;
- number of radial windows with occurrences;
- number of atom-pair or distance checks;
- number of atoms contributing to top and bottom calculations.

This should reveal whether the top/bottom asymmetry is biological/geometric or algorithmic.

## Suggested next profiling patch

Add an optional deeper profiler, for example:

```bash
--profile-thickness-detail
```

or extend `--profile-timing` to include lower-level counters inside `calcHalfMembraneZ`.

Useful counters:

```text
n_frames
n_coi_top
n_coi_bottom
n_membrane_top
n_membrane_bottom
n_distance_checks_top
n_distance_checks_bottom
n_nonempty_windows_top
n_nonempty_windows_bottom
time_half_z_top
time_half_z_bottom
time_window_assignment_top
time_window_assignment_bottom
```

## Suggested optimization ideas for `thickness_calc_half_z`

Potential directions to evaluate:

1. Precompute static group mappings when possible.
2. Replace repeated nested loops with vectorized distance calculations.
3. Use spatial cutoffs or neighbor-search structures for radial shell selection.
4. Avoid recalculating values that are constant over frames.
5. Separate algorithmic cost from system-geometry cost using detailed counters.
6. Check whether the top/bottom asymmetry reflects unbalanced CoI assignment or inefficient loop structure.

## Suggested acceptance criteria for future performance work

Any future optimization should preserve numerical output.

Minimum validation:

```bash
FRAME_COUNTS="5 25 50 100 250 500" FORCE_REPROCESS=1 ./13_prepare_reduced_xtc_benchmark_inputs.sh
FRAME_COUNTS="5 25 50 100 250 500" FORCE_REPROCESS=1 ./14_validate_reduced_xtc_against_pdb_reference.sh
FRAME_COUNTS="5 100 500" PROFILE_FORMATS="pdb xtc" FORCE_REPROCESS=1 ./15_profile_reduced_xtc.sh
```

Recommended checks:

- all PDB-vs-XTC comparisons pass;
- `*_thicknessTop.xvg` and `*_thicknessBottom.xvg` remain identical or within an explicitly justified tolerance;
- no change in documented sign convention;
- total runtime and section-level runtime are recorded before and after the change.

## Current recommended benchmark scripts

The current local benchmark scripts are:

```text
13_prepare_reduced_xtc_benchmark_inputs.sh
14_validate_reduced_xtc_against_pdb_reference.sh
15_profile_reduced_xtc.sh
```

These scripts live outside the MembIT source repository in the local analysis script folder. They are project-specific but serve as the current validation and profiling protocol for this system.
