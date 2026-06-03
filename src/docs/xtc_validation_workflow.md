# MembIT XTC validation workflow

This note documents the validation convention used while adding direct XTC
support to MembIT.

## Atom-number convention

MembIT has two trajectory readers after this change:

1. **PDB reader**
   - Keeps the historical behavior.
   - Index numbers refer to atom serials in the PDB trajectory being read.

2. **MDAnalysis reader for XTC/TRR/DCD/NC**
   - Requires `-s/--structure`.
   - Index numbers refer to 1-based atom positions in the structure/topology
     supplied with `-s`, matching GROMACS `.ndx` numbering.

Do not use a reduced-PDB MembIT index with a full-system XTC unless the atom
numbers have been remapped to the full system.

## Recommended testing order

1. Re-run the existing 5-frame PDB benchmark and compare against the known-good
   XVG outputs. This checks that PDB behavior was not changed.
2. Create a reduced 5-frame XTC/GRO using the same Protein+Phos atom set as the
   PDB benchmark. This isolates reader differences from atom-number remapping.
3. Create a full-system 5-frame XTC plus a full-system MembIT index. This tests
   the production path intended to avoid huge PDB trajectories.
4. Benchmark PDB, reduced-XTC, and full-XTC runs with `tools/benchmark_membit.py`.

## Performance metrics to collect

For each path, record:

- wall time;
- user CPU time;
- system CPU time;
- peak resident memory where available;
- number of frames;
- trajectory size on disk;
- whether the run used PDB, reduced XTC, or full-system XTC.

The goal is not only to know whether MDAnalysis is faster or slower per frame,
but whether the complete workflow avoids the expensive PDB conversion and 8+ GB
intermediate file.
