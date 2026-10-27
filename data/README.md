# EOS data and distribution warning

Six `*_standardized.csv` files were supplied for the user's thesis. The source,
derivation, data rights, and physical validity of each original EOS table have
**not been independently established**. Before publishing this repository,
confirm the original sources, required attributions and redistribution terms.
The project's MIT software license **does not grant permission for these data**.
Remove the `data/*_standardized.csv` files from a public repository if any
licenses or provenance are uncertain; end users can supply their own tables.

Required columns: `rho_g_cm-3`, `epsilon_MeV_fm-3`, `P_MeV_fm-3`.
The `rho_g_cm-3` column equals total epsilon/c^2, not independently tabulated
baryonic density. Pressure and energy density are mapped into km^-2 with
G=c=1. We use monotonic log-log PCHIP (NOT the supplied eighth-order fits).
Some supplied lower-density crust tables have small pressure reversals. The
loader warns and filters non-invertible points; the metadata records counts.
APR contains some superluminal tabulated local finite-difference slopes.
All physics conclusions therefore require an independent EOS provenance audit.

The stellar surface is approximated by the *lowest tabulated pressure*.
This is not mathematically identical to p(R)=0; investigate its effect on radii.
