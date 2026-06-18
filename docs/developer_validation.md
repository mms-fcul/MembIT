# Developer validation

This document describes validation checks useful when modifying MembIT.

## PDB vs XTC equivalence

For reader changes, compare a trusted PDB workflow to an XTC workflow generated
from the same processed atom set.

Recommended checks:

1. Use the same frames.
2. Use the same structure/index atom numbering.
3. Run the same MembIT command.
4. Compare `*_thicknessTop.xvg` and `*_thicknessBottom.xvg`.
5. Require exact or near-exact equality depending on expected floating-point
   behavior.

## Profile timing

Use:

```bash
--profile-timing
```

to separate:

- trajectory reader time;
- MDAnalysis universe initialization;
- atom mapping/population;
- frame analysis time;
- thickness/deformation calculation time.

## Useful milestones

Before committing reader or diagnostic changes:

```bash
python -m py_compile membit.py
python membit.py --help
python membit.py ... --diagnose-index
```

Then run at least one short PDB/XTC validation set.
