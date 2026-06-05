# darkin-class Coupling Guide (AI-generated, needs to be checked and modified)

This guide explains how to run the scalar field + interacting-CDM (`qcdm`) couplings in `darkin-class` in a user-friendly way.

## 1) Coupling types

For one coupling, the old shortcut still works:

```ini
scf_coupling_type = none | conformal | disformal | mixed | entropy | momentum
```

Current behavior:

- `conformal`, `disformal`, `mixed` use the `Q`-type coupling sector.
- `entropy` activates entropy-source terms.
- `momentum` activates Type-3 momentum coupling terms.
- `mixed` means conformal + disformal, as before.
- To combine the `Q` sector with entropy and/or momentum, use the explicit flags below.
- Vanilla CLASS `idm` is left for its built-in interactions (`idm_b`, `idm_g`, `idm_dr`).
- Scalar-field couplings use the separate `qcdm` component, so use `omega_qcdm` or `Omega_qcdm` for the dark matter that couples to the scalar field.

Multiple coupling sectors:

```ini
scf_allow_multiple_couplings = yes
scf_use_conformal = yes
scf_use_disformal = no
scf_use_entropy = yes
scf_use_momentum = no
```

Notes:

- `scf_use_*` flags can be used without `scf_coupling_type`.
- If both are present, `scf_coupling_type` first sets a default selection and explicit `scf_use_* = yes/no` entries override or add to it.
- `scf_allow_multiple_couplings = yes` is required when more than one sector is active, counting conformal/disformal together as the single `Q` sector.
- The combined equations add the active source terms in perturbations. For momentum + `Q`, the code uses the momentum-modified KG denominator with the replacement `dV -> dV-Q`.

## 2) Minimal requirements

A coupling is only active if both sectors are present:

- non-zero scalar-field-coupled CDM density (`omega_qcdm` / `Omega_qcdm`)
- non-zero scalar field density (`Omega_scf` / `Omega0_scf`)

Also set:

- `scf_potential = exp` or `double_exp`
- either `scf_coupling_type = ...` or explicit `scf_use_*` flags

Important runtime behavior:

- If `scf_coupling_type != none` and `Omega_cdm = 0` in synchronous gauge, the code sets a tiny floor for `Omega_cdm` automatically to avoid gauge singularities.

## 3) Two input styles

You can configure SCF parameters in either style.

### A) Explicit parameters (`scf_*`)

You provide named entries like `scf_lambda`, `scf_V0`, `scf_beta`, etc.

### B) Parameter list (`scf_parameters`)

You provide one ordered list. This is compact but easier to mis-order.

If explicit potential/coupling parameters are present, they take precedence over `scf_parameters`.

## 4) Core parameters by coupling

Always needed for potential:

- `exp`: `scf_lambda`, `scf_V0`
- `double_exp`: `scf_lambda`, `scf_V0`, `scf_lambda_2`, `scf_V0_2`

Coupling-specific parameters:

- `conformal`: typically `scf_C0`, `scf_beta`
- `disformal`: typically `scf_D0`, optionally `scf_alpha`
- `mixed`: combine conformal + disformal parameters
- `entropy`: `scf_g0`, `scf_h0`, `scf_As`, `scf_ns`, `scf_kp`, `scf_kc`, `scf_pc`
- `momentum`: `scf_gamma0` (or `scf_log10minus_gamma0`, but not both). For the quadratic model used in the momentum-transfer papers with positive parameter `beta`, use `scf_gamma0 = -beta`.

Momentum safety condition:

- avoid `scf_gamma0 = 0.5` (`1 - 2*scf_gamma0` singular denominator)
- more generally (for custom non-quadratic `gamma_scf(Z)`), avoid backgrounds where `1 - d²gamma_scf/dZ² = 0`

## 5) `scf_parameters` ordering

For `scf_potential = exp`:

```text
[lambda, V0, C0, beta, alpha, D0, (optional gamma for momentum), (optional entropy block...)]
```

For `scf_potential = double_exp`:

```text
[lambda, V0, lambda_2, V0_2, C0, beta, alpha, D0, (optional gamma for momentum), (optional entropy block...)]
```

Optional entropy block order:

```text
[g0, h0, As, ns, kp, kc, pc]
```

If `attractor_ic_scf = no` and you use `scf_parameters`, the code reads the last two entries of the list as:

```text
[... , scf_phi_ini, scf_phi_prime_ini]
```

So keep that in mind when building long lists.

## 6) Initial conditions for the scalar field

You can use:

- `attractor_ic_scf = yes` (automatic attractor ICs)
- `attractor_ic_scf = no` with either:
- explicit `scf_phi_ini`, `scf_phi_prime_ini`, or
- `scf_parameters` where the final two entries are interpreted as ICs

## 7) Minimal templates (INI-style)

### Conformal

```ini
Omega_scf = -1
omega_qcdm = 0.12
omega_cdm = 1e-10

scf_potential = exp
scf_coupling_type = conformal
scf_lambda = 0.1
scf_V0 = 1.0
scf_C0 = 1.0
scf_beta = 0.05
attractor_ic_scf = no
scf_phi_ini = 1.0
scf_phi_prime_ini = -1e-8
```

### Disformal

```ini
Omega_scf = -1
omega_qcdm = 0.12
omega_cdm = 1e-10

scf_potential = exp
scf_coupling_type = disformal
scf_lambda = 0.1
scf_V0 = 1.0
scf_D0 = 1.0
scf_alpha = 0.0
```

### Mixed

```ini
scf_coupling_type = mixed
scf_C0 = 1.0
scf_beta = 0.05
scf_D0 = 1.0
scf_alpha = 0.0
```

### Entropy

```ini
scf_coupling_type = entropy
scf_g0 = 0.0
scf_h0 = 0.0
scf_As = 1.0
scf_ns = 0.0
scf_kp = 0.05
scf_kc = 1.0
scf_pc = 2.0
```

Implementation note:

- entropy helper functions are user-editable in `source/background.c`: `g_scf`, `dg_scf`, `ddg_scf`, `h_scf`, `dh_scf`, `ddh_scf` (default `dh_scf=ddh_scf=0`)
- perturbations with `scf_coupling_type = entropy` are currently forced to the RK evolver internally for stability (NDF15 can hit linearisation failures in this branch)

### Momentum

```ini
scf_coupling_type = momentum
scf_gamma0 = -0.1
```

Implementation note:

- the code is written for a generic `gamma_scf(Z)` through `gamma_scf`, `dgamma_scf`, `ddgamma_scf` in `source/background.c`
- default model is quadratic: `gamma_scf(Z)=scf_gamma0*Z^2`
- with `Z=-phi'/a`, the positive paper convention `beta>0` corresponds to `gamma_scf(Z)=-beta*Z^2`, hence `scf_gamma0=-beta`
- background outputs now include cached momentum quantities: `scf_mom`, `gamma_scf`, `dgamma_scf`, `ddgamma_scf` (when momentum coupling is active)

Or:

```ini
scf_coupling_type = momentum
scf_log10minus_gamma0 = -1.0
```

`scf_log10minus_gamma0 = x` sets `scf_gamma0 = -10^x`, so it is the safest way to scan the positive paper parameter `beta=10^x`. `scf_gamma0` and `scf_log10minus_gamma0` are mutually exclusive.

### Combined Q + entropy

```ini
scf_allow_multiple_couplings = yes
scf_use_conformal = yes
scf_use_entropy = yes

scf_C0 = 1.0
scf_beta = 0.05
scf_g0 = 0.0
scf_h0 = 0.0
scf_As = 1.0
scf_ns = 0.0
scf_kp = 0.05
scf_kc = 1.0
scf_pc = 2.0
```

### Combined Q + momentum

```ini
scf_allow_multiple_couplings = yes
scf_use_conformal = yes
scf_use_momentum = yes

scf_C0 = 1.0
scf_beta = 0.05
scf_gamma0 = -0.1
```

## 8) Cobaya YAML usage

If you run through Cobaya, put the same keys under:

```yaml
theory:
  classy:
    extra_args:
      scf_coupling_type: "entropy"
      ...
```

Reference files in this repo:

- `conformal.yaml`
- `disformal.yaml`
- `mixed.yaml`
- `entropy_algebraic.yaml`
- `entropy_derivative.yaml`
- `momentum.yaml`

## 9) Practical troubleshooting

- If coupling appears inactive, first verify `Omega_scf != 0` and `omega_qcdm != 0`.
- Use `omega_idm` only when you intentionally want the built-in CLASS IDM species; it is not the scalar-field-coupled component.
- If momentum runs fail, check `scf_gamma0 != 0.5`.
- If you select `scf_coupling_type = entropy`, RK is used internally even if the global evolver is set to NDF15.
- If using `scf_parameters`, check ordering first; most setup issues come from list order.
- If using non-attractor ICs with `scf_parameters`, ensure the final two entries are intended to be `phi_ini` and `phi_prime_ini`.
