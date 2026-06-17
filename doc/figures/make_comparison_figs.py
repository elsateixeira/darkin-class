#!/usr/bin/env python3
"""Regenerate the implementation-note comparison figures from CLASS background output.

Reads the *_background.dat files produced by
    ./class fig3_2505_matter.ini
    ./class notes_fig1_model1.ini
    ./class notes_fig1_model2.ini
and writes
    matter_coupling_fig3_comparison.{pdf,svg}
    radiation_coupling_fig1_comparison.{pdf,svg}

Self-contained: no project dependencies beyond numpy/matplotlib.
"""
import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "output")


def load(tag):
    fn = os.path.join(OUT, f"{tag}_background.dat")
    hdr = None
    with open(fn) as f:
        for line in f:
            if line.startswith("#") and re.search(r"\d+:z", line):
                hdr = line
    cols = re.findall(r"\d+:([^\s].*?)(?=\s{2,}\d+:|\s*$)", hdr.strip().lstrip("#"))
    cols = [c.strip() for c in cols]
    data = np.loadtxt(fn)
    idx = {c: i for i, c in enumerate(cols)}
    return data, idx


def col(data, idx, name):
    return data[:, idx[name]]


def matter_figure():
    d, i = load("fig3_2505_matter")
    z = col(d, i, "z")
    N = np.log(1.0 / (1.0 + z))
    rho_tot = col(d, i, "(.)rho_tot")
    rho_de = col(d, i, "(.)rho_DE_eff")
    w_de = col(d, i, "w_DE_eff")
    Om_r = col(d, i, "Omega_r(z)")
    A_m = np.sqrt(col(d, i, "C_scf"))
    rho_m_bar = col(d, i, "(.)rho_qcdm") / A_m
    Om_m_bar = rho_m_bar / rho_tot
    Om_de = rho_de / rho_tot
    # z=0 reference for f_DE
    r0 = np.argmin(np.abs(z))
    f_de = rho_de / rho_de[r0]
    order = np.argsort(N)

    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.3))
    fig.suptitle("Matter coupling: reproduced with Fig. 3 axes from arXiv:2505.10410v2")

    a = ax[0, 0]
    a.plot(N[order], Om_r[order], "C3", label=r"$\Omega_r$")
    a.plot(N[order], Om_m_bar[order], "C0", label=r"$\bar\Omega_m$")
    a.plot(N[order], Om_de[order], "C2", label=r"$\Omega_{\rm DE}$")
    a.set_xlabel(r"$N=\log a$"); a.set_ylabel(r"$\Omega_i$")
    a.set_xlim(-21, 0.2); a.set_title("(a) Density fractions, CLASS"); a.legend(loc="center left")

    zsel = z <= 4.0
    b = ax[0, 1]
    b.plot(z[zsel], f_de[zsel], "C2")
    b.axhline(0, color="k", lw=0.8); b.set_xlabel("z"); b.set_ylabel(r"$f_{\rm DE}$")
    b.set_xlim(0, 4); b.set_title("(b) Recent DE density, CLASS")

    c = ax[1, 0]
    c.plot(N[order], w_de[order], "C2")
    c.axhline(-1, color="k", ls="--", lw=0.8); c.set_xlabel(r"$N=\log a$"); c.set_ylabel(r"$w_{\rm DE}$")
    c.set_xlim(-21, 0.2); c.set_title("(c) DE equation of state, CLASS")

    dd = ax[1, 1]
    dd.plot(z[zsel], w_de[zsel], "C2")
    dd.axhline(-1, color="k", ls="--", lw=0.8); dd.set_xlabel("z"); dd.set_ylabel(r"$w_{\rm DE}$")
    dd.set_xlim(0, 4); dd.set_title("(d) Recent DE equation of state, CLASS")

    for row in ax:
        for a_ in row:
            a_.grid(alpha=0.3)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(HERE, f"matter_coupling_fig3_comparison.{ext}"))
    plt.close(fig)


def radiation_panel(axcol, tag, title):
    d, i = load(tag)
    z = col(d, i, "z")
    N = np.log(1.0 / (1.0 + z))
    rho_tot = col(d, i, "(.)rho_tot")
    A_r = col(d, i, "A_r_scf")
    rho_dr_bar = col(d, i, "(.)rho_idr") / A_r
    A_m = np.sqrt(col(d, i, "C_scf"))
    rho_m_bar = col(d, i, "(.)rho_qcdm") / A_m
    rho_phi = col(d, i, "(.)rho_scf")
    rho_de = col(d, i, "(.)rho_DE_eff")
    w_de = col(d, i, "w_DE_eff")
    order = np.argsort(N)
    # CLASS density units are such that rho_tot = H^2; normalise by H0^2 = rho_tot(z=0).
    H02 = rho_tot[np.argmin(np.abs(z))]

    a = axcol[0]
    a.plot(N[order], (rho_dr_bar / rho_tot)[order], "C3", label=r"$\bar\Omega_{\rm DR}$")
    a.plot(N[order], (rho_m_bar / rho_tot)[order], "C0", label=r"$\bar\Omega_m$")
    a.plot(N[order], (rho_de / rho_tot)[order], "C2", label=r"$\Omega_{\rm DE}$")
    a.set_xlabel(r"$N=\log a$"); a.set_ylabel(r"$\Omega_i$"); a.set_xlim(-30, 0.5)
    a.set_title(title); a.legend(loc="center left")

    zz = z > 0
    lz = np.log10(z[zz])
    b = axcol[1]
    b.plot(lz, np.log10(rho_dr_bar[zz] / H02), "C3", label=r"$\bar\rho_{\rm DR}$")
    b.plot(lz, np.log10(rho_m_bar[zz] / H02), "C0", label=r"$\bar\rho_m$")
    b.plot(lz, np.log10(rho_phi[zz] / H02), "C2", label=r"$\rho_\phi$")
    b.set_xlabel(r"$\log_{10}z$"); b.set_ylabel(r"$\log_{10}(\rho/H_0^2)$")
    b.set_xlim(0, 7); b.legend(loc="upper right")

    zsel = z <= 4.0
    c = axcol[2]
    c.plot(z[zsel], w_de[zsel], "C2")
    c.axhline(-1, color="k", ls="--", lw=0.8)
    c.set_xlabel("z"); c.set_ylabel(r"$w_{\rm DE}$"); c.set_xlim(0, 4)
    for a_ in axcol:
        a_.grid(alpha=0.3)


def radiation_figure():
    fig, ax = plt.subplots(3, 2, figsize=(11.5, 12.5))
    fig.suptitle("Radiation coupling: reproduced with Fig. 1 axes from NotesRadiationCoupling")
    radiation_panel(ax[:, 0], "notes_fig1_model1", "Model 1: CLASS, Notes Fig. 1 axes")
    radiation_panel(ax[:, 1], "notes_fig1_model2", "Model 2: CLASS, Notes Fig. 1 axes")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(HERE, f"radiation_coupling_fig1_comparison.{ext}"))
    plt.close(fig)


if __name__ == "__main__":
    matter_figure()
    radiation_figure()
    print("wrote matter_coupling_fig3_comparison and radiation_coupling_fig1_comparison")
