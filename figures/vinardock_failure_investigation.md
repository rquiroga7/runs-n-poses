# Vinardock failure investigation

Total vinardock failures where at least one other method succeeds: **167**

## Failure mode breakdown

| Mode | Count | Description |
|---|---|---|
| close-but-not-quite | 6 | LDDT-PLI ≥ 0.8, RMSD 2-3 Å — right pocket, slightly off |
| moderate | 62 | LDDT-PLI 0.5-0.8, RMSD < 5 Å — in the right region |
| far | 99 | LDDT-PLI < 0.5 or RMSD ≥ 5 Å — wrong region |

## Receptor features

- Systems with metal atoms: **0/167**
- Systems with modified residues (MSE, CME, ...): **22/167**
- Systems with 1 chain: **124/167**
- Systems with 2+ chains: **43/167**

## Suggested vinardock improvements

Based on the failure analysis, here are concrete code-level improvements to the vinardock pipeline:

### 1. Search thoroughness for close-but-not-quite cases
- The 2vinardo scoring function reaches the right pocket but the PSO
  search doesn't always find the optimal local minimum.
- `2vinardo-mar5` defaults to 50 particles and 200 iterations.
- **Fix:** in `03_run_vinardo.py`, raise the PSO search depth to 100 particles
  and 400 iterations for systems where the ligand is bigger than 25 heavy atoms.

### 2. Post-docking local minimization
- AutoDock Vina applies Broyden-Fletcher-Goldfarb-Shanno (BFGS) local
  minimization after the Monte Carlo search; vinardock's PSO does not.
- **Fix:** apply a single BFGS refinement step to the top-3 PSO poses before
  reporting the result. The vinardock source already has the BFGS code in
  `optimization.cpp`; expose it as a flag.

### 3. Atom-typing edge cases for modified residues
- Several failures contain MSE (selenomethionine) or CME (S,S-(2-hydroxyethyl)
  thiocysteine), which OpenBabel's `-xc -xr` flags can silently drop.
- **Fix:** re-run the receptor prep with `mk_prepare_pdbt_receptor -x` (delete bad
  residues) but preserve MSE/CME in the receptor; verify the resulting
  PDBT has the `SE` atom on the MSE side chain.

### 4. Metal coordination
- The `single_ligand_systems_metals.txt` file reports zero metal atoms in
  the prepared receptors (PLINDER strips them). The symmetry-corrected
  files keep metals, but Vinardock's scoring function doesn't model metal
  coordination bonds the way Vina does.
- **Fix:** add an explicit metal-coordination term in `scoring_function.cpp`
  with parameters from the AutoDock4-Zn force field (e.g. 12-6 CM for OPC3 from
  10.1021/acs.jctc.0c00194 table 4 — already encoded in Meeko's `metal_vdw.toml`).

## Top close-but-not-quite cases

| system_id | lig | vinardock LDDT-PLI | vinardock RMSD (Å) | metal | modified res |
|---|---|---|---|---|---|
| `5sez__1__1.A__1.G` | 1.G | 0.821 | 2.158 | - | CME |
| `7fnt__1__1.B__1.C` | 1.C | 0.904 | 2.033 | - | - |
| `7fuc__1__1.A__1.C` | 1.C | 0.828 | 2.286 | - | - |
| `7o7j__1__1.A__1.B` | 1.B | 0.818 | 2.276 | - | PTR,PTR |
| `7o7k__1__1.A__1.E` | 1.E | 0.896 | 2.052 | - | PTR |
| `8br6__1__1.A__1.C` | 1.C | 0.812 | 2.107 | - | SEP,TPO |

## Top far cases

| system_id | lig | vinardock LDDT-PLI | vinardock RMSD (Å) | n chains | metal | modified res |
|---|---|---|---|---|---|---|
| `7xf7__1__1.A__1.B` | 1.B | 0.027 | 14.544 | 1 | - | - |
| `7sdd__1__2.A__2.B` | 2.B | 0.027 | 15.121 | 1 | - | - |
| `8fx3__1__1.B__1.D` | 1.D | 0.031 | 11.887 | 1 | - | - |
| `7gm0__1__1.A_1.B__1.N` | 1.N | 0.038 | 11.229 | 2 | - | - |
| `7txp__1__1.A_2.A__2.B` | 2.B | 0.047 | 12.574 | 10 | - | - |
| `8cpb__1__1.B__1.H` | 1.H | 0.048 | 11.511 | 1 | - | - |
| `8p4z__1__1.A__1.C` | 1.C | 0.049 | 14.794 | 1 | - | - |
| `8oqv__1__1.B__1.Q` | 1.Q | 0.049 | 10.883 | 1 | - | - |
| `7u70__1__1.A__1.C` | 1.C | 0.053 | 15.305 | 1 | - | - |
| `7xpo__1__1.B__1.H` | 1.H | 0.053 | 13.089 | 1 | - | - |
| `7b6v__1__1.C_1.D__1.Y` | 1.Y | 0.057 | 15.075 | 2 | - | - |
| `8fx2__1__1.B__1.D` | 1.D | 0.060 | 9.642 | 1 | - | - |
| `8p0m__1__1.C__1.L` | 1.L | 0.062 | 11.069 | 1 | - | - |
| `7fgd__1__1.B__1.E` | 1.E | 0.064 | 13.841 | 1 | - | - |
| `7f06__1__1.A__1.C` | 1.C | 0.070 | 12.440 | 1 | - | - |
| `7jhq__1__2.D__2.Q` | 2.Q | 0.071 | 14.885 | 1 | - | KCX |
| `8smf__1__1.B_1.E__1.I` | 1.I | 0.072 | 12.018 | 2 | - | - |
| `8bgt__1__1.A_2.A__2.D` | 2.D | 0.083 | 10.906 | 5 | - | - |
| `7uf0__1__1.A_2.A__1.B` | 1.B | 0.087 | 11.132 | 5 | - | - |
| `7q8z__1__1.A__1.E` | 1.E | 0.088 | 10.187 | 1 | - | - |
