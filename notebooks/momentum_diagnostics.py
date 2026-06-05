#!/usr/bin/env python3
"""
Momentum-coupling diagnostics for darkin-class vs class_IDE.

What this checks:
1) Whether your Python `classy` import is the build you expect.
2) Coupled/uncoupled P(k) ratio in each code for matched inputs.
3) Background coupling-strength diagnostics:
   R(z) = Z * gamma_Z / (3 rho_DM), gamma/rho_DM.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--classide-bin",
        default="/private/tmp/class_IDE_local/class",
        help="Path to runnable class_IDE binary",
    )
    parser.add_argument(
        "--darkin-bin",
        default="class",
        help="Path to darkin-class binary",
    )
    parser.add_argument(
        "--out-dir",
        default="/private/tmp/momentum_diagnostics",
        help="Output directory for temporary runs",
    )
    parser.add_argument(
        "--mode",
        choices=["classide_like", "paper2603_like"],
        default="classide_like",
        help="Benchmark parameter preset",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=-0.8,
        help="Momentum coupling value: scf_veta in class_IDE, scf_gamma0 in darkin",
    )
    parser.add_argument(
        "--scan",
        action="store_true",
        help="Scan a preset range of momentum couplings in darkin-class only",
    )
    return parser.parse_args()


def header_cols(path: Path) -> List[str]:
    line = None
    for raw in path.read_text().splitlines()[:60]:
        if raw.startswith("#") and "1:z" in raw:
            line = raw[1:]
            break
    if line is None:
        raise RuntimeError(f"Cannot parse header columns in {path}")
    marks = list(re.finditer(r"\s(\d+):", line))
    cols: List[str] = []
    for i, mark in enumerate(marks):
        start = mark.end()
        end = marks[i + 1].start() if i + 1 < len(marks) else len(line)
        cols.append(line[start:end].strip())
    return cols


def write_ini(path: Path, entries: Dict[str, object]) -> None:
    path.write_text("\n".join(f"{k} = {v}" for k, v in entries.items()) + "\n")


def newest_match(out_dir: Path, prefix: str, suffix: str) -> Path:
    cands = sorted(out_dir.glob(f"{prefix}_*_{suffix}.dat")) + sorted(
        out_dir.glob(f"{prefix}_{suffix}.dat")
    )
    if not cands:
        raise FileNotFoundError(f"No {suffix}.dat for prefix {prefix} in {out_dir}")
    return cands[-1]


def run_class(exe: Path, out_dir: Path, tag: str, params: Dict[str, object]) -> None:
    ini = out_dir / f"{tag}.ini"
    p = dict(params)
    p["root"] = str(out_dir / f"{tag}_")
    write_ini(ini, p)
    proc = subprocess.run([str(exe), str(ini)], capture_output=True, text=True)
    (out_dir / f"{tag}.log").write_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    if proc.returncode != 0:
        tail = "\n".join(((proc.stdout or "") + "\n" + (proc.stderr or "")).splitlines()[-40:])
        raise RuntimeError(f"Run failed for {tag}:\n{tail}")


def load_pk(out_dir: Path, prefix: str) -> Tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(newest_match(out_dir, prefix, "pk"))
    return arr[:, 0], arr[:, 1]


@dataclass
class DiagPoint:
    z_target: float
    z_actual: float
    z_var: float
    gamma_over_rho_dm: float
    r_transfer: float


def background_diag_ide(
    bg: np.ndarray, cols: List[str], beta: float, z_points: List[float]
) -> List[DiagPoint]:
    c = {name: i for i, name in enumerate(cols)}
    z = bg[:, c["z"]]
    a = 1.0 / (1.0 + z)
    rho_cdm = bg[:, c["(.)rho_cdm"]]
    phi_prime = bg[:, c["phi'_scf"]]
    z_var = -phi_prime / a
    gamma = beta * z_var * z_var
    gamma_z = 2.0 * beta * z_var
    r_transfer = z_var * gamma_z / (3.0 * rho_cdm)
    out: List[DiagPoint] = []
    for zt in z_points:
        i = int(np.argmin(np.abs(z - zt)))
        out.append(
            DiagPoint(
                z_target=zt,
                z_actual=float(z[i]),
                z_var=float(z_var[i]),
                gamma_over_rho_dm=float(gamma[i] / rho_cdm[i]),
                r_transfer=float(r_transfer[i]),
            )
        )
    return out


def background_diag_darkin(bg: np.ndarray, cols: List[str], z_points: List[float]) -> List[DiagPoint]:
    c = {name: i for i, name in enumerate(cols)}
    z = bg[:, c["z"]]
    rho_idm = bg[:, c["(.)rho_idm"]]
    z_var = bg[:, c["scf_mom"]]
    gamma = bg[:, c["gamma_scf"]]
    gamma_z = bg[:, c["dgamma_scf"]]
    r_transfer = z_var * gamma_z / (3.0 * rho_idm)
    out: List[DiagPoint] = []
    for zt in z_points:
        i = int(np.argmin(np.abs(z - zt)))
        out.append(
            DiagPoint(
                z_target=zt,
                z_actual=float(z[i]),
                z_var=float(z_var[i]),
                gamma_over_rho_dm=float(gamma[i] / rho_idm[i]),
                r_transfer=float(r_transfer[i]),
            )
        )
    return out


def maybe_print_classy_probe() -> None:
    try:
        import classy
        from classy import Class  # noqa: F401
    except Exception as err:
        print(f"[classy probe] import failed: {err}")
        return

    print(f"[classy probe] imported from: {classy.__file__}")
    test_params = {
        "omega_b": 0.022,
        "omega_cdm": 0.12,
        "h": 0.67,
        "A_s": 2.1e-9,
        "n_s": 0.965,
        "tau_reio": 0.054,
        "output": "mPk",
        "P_k_max_h/Mpc": 1.0,
        "z_pk": 0.0,
        "Omega_scf": -1,
        "Omega_Lambda": 0.0,
        "Omega_fld": 0.0,
        "scf_potential": "exp",
        "scf_coupling_type": "momentum",
        "scf_lambda": 0.5,
        "scf_V0": 1.0,
        "attractor_ic_scf": "no",
        "scf_phi_ini": 1e-4,
        "scf_phi_prime_ini": 0.0,
        "scf_shooting_target": "scf_V0",
    }
    for key, value in [("scf_gamma0", -0.8), ("scf_gamma", -0.8), ("scf_veta", -0.8)]:
        from classy import Class

        c = Class()
        p = dict(test_params)
        p[key] = value
        try:
            c.set(p)
            c.compute()
            print(f"[classy probe] accepts `{key}` (sigma8={c.sigma8():.6f})")
        except Exception as err:
            print(f"[classy probe] rejects `{key}` ({type(err).__name__})")
        finally:
            try:
                c.struct_cleanup()
                c.empty()
            except Exception:
                pass


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    beta = float(args.beta)

    class_ide = Path(args.classide_bin).resolve()
    darkin = Path(args.darkin_bin).resolve()
    if not class_ide.exists():
        raise FileNotFoundError(f"class_IDE binary not found: {class_ide}")
    if not darkin.exists():
        raise FileNotFoundError(f"darkin binary not found: {darkin}")

    maybe_print_classy_probe()

    common = {
        "h": 0.673,
        "A_s": 2.215e-9,
        "n_s": 0.9619,
        "tau_reio": 0.05564022,
        "N_ncdm": 1,
        "N_ur": 2.0328,
        "m_ncdm": 0.06,
        "Omega_Lambda": 0.0,
        "Omega_fld": 0.0,
        "Omega_scf": -1,
        "attractor_ic_scf": "no",
        "output": "mPk",
        "P_k_max_h/Mpc": 10.0,
        "z_pk": 0,
        "gauge": "synchronous",
        "format": "class",
        "write background": "yes",
    }

    if args.mode == "classide_like":
        scf_lambda = 1.2247
        ide_scf_parameters = "1.2247, 0.0, 0.0, 0.0, 1e-4, 0.0"
    else:
        scf_lambda = 0.5
        ide_scf_parameters = "0.5, 0.0, 0.0, 0.0, 1e-4, 0.0"

    ide_base = dict(common)
    ide_base.update(
        {
            "omega_b": 0.022032,
            "omega_cdm": 0.12038,
            "scf_parameters": ide_scf_parameters,
            "scf_tuning_index": 2,
            "a_ini_over_a_today_default": 1e-9,
        }
    )

    dk_base = dict(common)
    dk_base.update(
        {
            "omega_b": 0.022032,
            "omega_cdm": 1e-10,
            "omega_idm": 0.12038,
            "scf_potential": "exp",
            "scf_coupling_type": "momentum",
            "scf_lambda": scf_lambda,
            "scf_V0": 1.0,
            "scf_phi_ini": 1e-4,
            "scf_phi_prime_ini": 0.0,
            "scf_shooting_target": "scf_V0",
            "a_ini_over_a_today_default": 8e-10,
        }
    )

    if args.scan:
        gammas = [-0.01, -0.05, -0.1, -0.3, -0.8, -1.0, -3.0, -10.0, -100.0]
        run_class(darkin, out_dir, "scan_ref", {**dk_base, "scf_gamma0": 0.0})
        k_ref, p_ref = load_pk(out_dir, "scan_ref")
        print("\nDarkin momentum scan:")
        print("gamma0       min P/P0   max P/P0   P/P0(k~1e-2)   P/P0(k~1)   P/P0(k~5)")
        for gamma0 in gammas:
            tag = f"scan_{str(gamma0).replace('-', 'm').replace('.', 'p')}"
            try:
                run_class(darkin, out_dir, tag, {**dk_base, "scf_gamma0": gamma0})
                k_g, p_g = load_pk(out_dir, tag)
                r = np.interp(k_ref, k_g, p_g) / p_ref
                vals = []
                for k_test in [1e-2, 1.0, 5.0]:
                    i = int(np.argmin(np.abs(k_ref - k_test)))
                    vals.append(r[i])
                print(
                    f"{gamma0:9.3g}  {np.nanmin(r):9.6f}  {np.nanmax(r):9.6f}"
                    f"  {vals[0]:13.6f}  {vals[1]:10.6f}  {vals[2]:10.6f}"
                )
            except Exception as err:
                print(f"{gamma0:9.3g}  failed: {err}")
        print(f"\nScan outputs written to: {out_dir}")
        return

    run_class(class_ide, out_dir, "ide_unc", {**ide_base, "scf_veta": 0.0})
    run_class(class_ide, out_dir, "ide_cpl", {**ide_base, "scf_veta": beta})
    run_class(darkin, out_dir, "dk_unc", {**dk_base, "scf_gamma0": 0.0})
    run_class(darkin, out_dir, "dk_cpl", {**dk_base, "scf_gamma0": beta})

    k_iu, p_iu = load_pk(out_dir, "ide_unc")
    k_ic, p_ic = load_pk(out_dir, "ide_cpl")
    k_du, p_du = load_pk(out_dir, "dk_unc")
    k_dc, p_dc = load_pk(out_dir, "dk_cpl")

    r_ide = np.interp(k_iu, k_ic, p_ic) / p_iu
    r_dk = np.interp(k_du, k_dc, p_dc) / p_du
    r_ide_on_dk = np.interp(k_du, k_iu, r_ide)

    print("\nP(k) coupled/uncoupled:")
    print(f"class_IDE min/max = {r_ide.min():.6f} / {r_ide.max():.6f}")
    print(f"darkin    min/max = {r_dk.min():.6f} / {r_dk.max():.6f}")
    print(
        "darkin / class_IDE (ratio of ratios) min/max = "
        f"{(r_dk / r_ide_on_dk).min():.6f} / {(r_dk / r_ide_on_dk).max():.6f}"
    )

    for k_test in [1e-3, 1e-2, 1e-1, 1.0]:
        i = int(np.argmin(np.abs(k_du - k_test)))
        print(
            f"k={k_du[i]:.4e}: "
            f"IDE={np.interp(k_du[i], k_iu, r_ide):.6f}, "
            f"darkin={r_dk[i]:.6f}"
        )

    bg_ide_path = newest_match(out_dir, "ide_cpl", "background")
    bg_dk_path = newest_match(out_dir, "dk_cpl", "background")
    bg_ide = np.loadtxt(bg_ide_path)
    bg_dk = np.loadtxt(bg_dk_path)
    ide_cols = header_cols(bg_ide_path)
    dk_cols = header_cols(bg_dk_path)

    z_points = [1000.0, 100.0, 10.0, 1.0, 0.0]
    ide_diag = background_diag_ide(bg_ide, cols=ide_cols, beta=beta, z_points=z_points)
    dk_diag = background_diag_darkin(bg_dk, cols=dk_cols, z_points=z_points)

    print("\nBackground coupling diagnostics:")
    for i, zt in enumerate(z_points):
        a = ide_diag[i]
        b = dk_diag[i]
        print(
            f"z~{zt:.0f} IDE: Z={a.z_var:.3e}, gamma/rho={a.gamma_over_rho_dm:.3e}, R={a.r_transfer:.3e}"
        )
        print(
            f"      DK : Z={b.z_var:.3e}, gamma/rho={b.gamma_over_rho_dm:.3e}, R={b.r_transfer:.3e}"
        )

    print(f"\nLogs and outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
