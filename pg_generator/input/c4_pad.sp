*** C4 Pad Subcircuit ***
.SUBCKT C4_PAD n_prt
*L1 n_vdd 1 1e-11
*C1 1 0 1e-12
R1 n_prt n_vdd 0.05
*C2 2 0 1e-12
*L2 n_vdd 2 1e-11
*C3 3 0 1e-12
*R2 2 3 0.01
.ends
