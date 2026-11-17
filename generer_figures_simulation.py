#!/usr/bin/env python3
"""Laboratoire 4 (GEL-7072) : simulations LTspice et figures du rapport.

Toutes les courbes de simulation proviennent des fichiers .raw produits par
LTspice a partir des netlists du dossier ltspice/ (modeles TI de l'INA823 et
du TL082). Ce script ne calcule aucune reponse en frequence lui-meme : il
lance LTspice, lit les resultats et les met en forme.

Usage :
    python generer_figures_simulation.py           # lance LTspice puis trace
    python generer_figures_simulation.py --no-run  # trace a partir des .raw existants

Le chemin de LTspice peut etre impose avec la variable d'environnement LTSPICE.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ICI = Path(__file__).resolve().parent
DOSSIER_SPICE = ICI / "ltspice"
DOSSIER_FIG = ICI / "figures"
NETLISTS = ["ina823_offset_dc", "ina823_bode", "filtres_bode",
            "bioampli_ac", "drl_ac", "bioampli_tran"]

# Valeurs nominales (doivent correspondre a ltspice/bioampli.sub)
VDD, VREF = 3.3, 1.65
G1_THEO = 1 + 100e3 / (2 * 5.6e3)          # 9,93 V/V
G2_THEO = 1 + 100e3 / 1e3                  # 101 V/V

# Palette categorielle validee (slots 1 a 3) et encres de texte
BLEU, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GRIS = "#8a8984"
ENCRE, ENCRE2 = "#0b0b0b", "#52514e"

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.labelsize": 10.5,
    "axes.titlesize": 10.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 8.5, "lines.linewidth": 1.5, "axes.edgecolor": ENCRE2,
    "axes.labelcolor": ENCRE, "xtick.color": ENCRE2, "ytick.color": ENCRE2,
    "grid.color": "#d9d8d4", "grid.linewidth": 0.6, "grid.linestyle": "-",
    "legend.frameon": True, "legend.framealpha": 0.92, "legend.edgecolor": "#d9d8d4",
    "savefig.dpi": 300, "savefig.bbox": "tight",
})


# ---------------------------------------------------------------------------
# ECG synthetique pour la simulation temporelle (fichier PWL de LTspice)
# ---------------------------------------------------------------------------
def ecg_synthetique(t, fc=68.0):
    """ECG de derivation I simplifie (en volts), onde R de 1,2 mV."""
    periode = 60.0 / fc
    v = np.zeros_like(t)
    for t0 in np.arange(0.2, t[-1] + periode, periode):
        v += 0.15 * np.exp(-((t - (t0 - 0.16)) / 0.025) ** 2)   # P
        v -= 0.15 * np.exp(-((t - (t0 - 0.04)) / 0.010) ** 2)   # Q
        v += 1.20 * np.exp(-((t - t0) / 0.015) ** 2)            # R
        v -= 0.30 * np.exp(-((t - (t0 + 0.035)) / 0.012) ** 2)  # S
        v += 0.35 * np.exp(-((t - (t0 + 0.22)) / 0.050) ** 2)   # T
    return v * 1e-3


def ecrire_pwl():
    t = np.arange(0.0, 4.0 + 1e-9, 0.5e-3)
    v = ecg_synthetique(t)
    with open(DOSSIER_SPICE / "ecg_pwl.txt", "w") as f:
        for ti, vi in zip(t, v):
            f.write(f"{ti:.4f} {vi:.6e}\n")


# ---------------------------------------------------------------------------
# Execution de LTspice
# ---------------------------------------------------------------------------
def trouver_ltspice():
    if os.environ.get("LTSPICE"):
        return os.environ["LTSPICE"]
    candidats = [
        "/Applications/LTspice.app/Contents/MacOS/LTspice",
        r"C:\Program Files\ADI\LTspice\LTspice.exe",
        os.path.expandvars(r"%LOCALAPPDATA%\Programs\ADI\LTspice\LTspice.exe"),
        r"C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe",
    ]
    for c in candidats:
        if Path(c).exists():
            return c
    return shutil.which("LTspice") or shutil.which("ltspice")


def lancer_ltspice(exe):
    for nom in NETLISTS:
        cir = DOSSIER_SPICE / f"{nom}.cir"
        raw = DOSSIER_SPICE / f"{nom}.raw"
        if raw.exists():
            raw.unlink()
        print(f"LTspice : {cir.name} ...", flush=True)
        subprocess.run([exe, "-b", str(cir)], cwd=DOSSIER_SPICE, check=False, timeout=900)
        if not raw.exists():
            log = DOSSIER_SPICE / f"{nom}.log"
            texte = log.read_text(errors="replace") if log.exists() else "(pas de .log)"
            sys.exit(f"Echec de la simulation {cir.name} :\n{texte}")


# ---------------------------------------------------------------------------
# Lecture des fichiers .raw de LTspice (binaire ou ASCII, UTF-16 ou UTF-8)
# ---------------------------------------------------------------------------
def lire_raw(chemin):
    donnees = Path(chemin).read_bytes()
    codec = "utf-16-le" if donnees[1:2] == b"\x00" else "latin-1"
    marqueurs = {m: m.encode(codec) for m in ("Binary:\n", "Values:\n")}
    pos, mode = -1, None
    for m, octets in marqueurs.items():
        pos = donnees.find(octets)
        if pos >= 0:
            mode = m
            break
    if pos < 0:
        raise ValueError(f"{chemin} : en-tete non reconnu")
    entete = donnees[:pos].decode(codec)
    debut = pos + len(marqueurs[mode])

    lignes = entete.splitlines()
    info, noms = {}, []
    for i, ligne in enumerate(lignes):
        if ligne.startswith("Variables:"):
            for v in lignes[i + 1:]:
                champs = v.split()
                if len(champs) >= 2 and champs[0].isdigit():
                    noms.append(champs[1])
            break
        if ":" in ligne:
            cle, val = ligne.split(":", 1)
            info[cle.strip().lower()] = val.strip()
    n, m = int(info["no. variables"]), int(info["no. points"])
    drapeaux = info.get("flags", "").lower()
    complexe = "complex" in drapeaux

    if mode == "Binary:\n":
        reste = len(donnees) - debut
        if complexe:
            tab = np.frombuffer(donnees, dtype="<c16", count=n * m, offset=debut).reshape(m, n)
            cols = [tab[:, k] for k in range(n)]
        elif "double" in drapeaux or reste == n * m * 8:
            tab = np.frombuffer(donnees, dtype="<f8", count=n * m, offset=debut).reshape(m, n)
            cols = [tab[:, k] for k in range(n)]
        else:
            dt = np.dtype([("x", "<f8")] + [(f"v{k}", "<f4") for k in range(1, n)])
            tab = np.frombuffer(donnees, dtype=dt, count=m, offset=debut)
            cols = [tab["x"]] + [tab[f"v{k}"].astype(float) for k in range(1, n)]
    else:
        jetons = donnees[debut:].decode(codec).split()
        valeurs = []
        k = 0
        for _ in range(m):
            k += 1  # indice du point
            for _ in range(n):
                j = jetons[k]
                k += 1
                if complexe:
                    re, im = j.split(",")
                    valeurs.append(complex(float(re), float(im)))
                else:
                    valeurs.append(float(j))
        tab = np.array(valeurs).reshape(m, n)
        cols = [tab[:, k] for k in range(n)]

    res = {}
    for nom, col in zip(noms, cols):
        res[nom.lower()] = col
    x = np.real(cols[0])  # temps, frequence ou valeur balayee
    if "transient" in info.get("plotname", "").lower():
        x = np.abs(x)  # LTspice marque certains points de temps par un signe negatif
    res["_x"] = x
    return res


def v(raw, noeud):
    return raw[f"v({noeud.lower()})"]


# ---------------------------------------------------------------------------
# Outils d'analyse
# ---------------------------------------------------------------------------
def db(h):
    return 20 * np.log10(np.maximum(np.abs(h), 1e-30))


def valeur_a(f, y, f0):
    """Interpolation en echelle logarithmique de frequence."""
    return np.interp(np.log10(f0), np.log10(f), y)


def _croise(f, y, k, seuil):
    lf = np.log10(f[k:k + 2])
    return 10 ** np.interp(seuil, sorted(y[k:k + 2]), lf if y[k] < y[k + 1] else lf[::-1])


def coupures(f, mag_db):
    """Gain max et frequences a -3 dB (basse et haute) d'une courbe en dB."""
    i = int(np.argmax(mag_db))
    gmax = mag_db[i]
    seuil = gmax - 3.0
    fl = fh = None
    bas = np.where(mag_db[:i] < seuil)[0]
    if bas.size:
        fl = _croise(f, mag_db, bas[-1], seuil)
    haut = np.where(mag_db[i:] < seuil)[0]
    if haut.size:
        fh = _croise(f, mag_db, i + haut[0] - 1, seuil)
    return gmax, fl, fh


def fmt_hz(x):
    if x is None:
        return "n/a"
    if x >= 1e3:
        return f"{x / 1e3:.3g} kHz".replace(".", ",")
    return f"{x:.3g} Hz".replace(".", ",")


def virgule(x, nd=2):
    return f"{x:.{nd}f}".replace(".", ",")


def axes_freq(ax):
    ax.set_xscale("log")
    ax.grid(True, which="major")
    ax.grid(True, which="minor", alpha=0.35)


def marquer(ax, x, y, texte, dx=1.6, dy=0.0, ha="left"):
    ax.plot([x], [y], "o", ms=5, mfc="white", mec=ENCRE, mew=1.1, zorder=6)
    ax.annotate(texte, xy=(x, y), xytext=(x * dx, y + dy), color=ENCRE, fontsize=8.5,
                ha=ha, va="center",
                arrowprops=dict(arrowstyle="-", color=ENCRE2, lw=0.7))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_offset(res):
    raw = lire_raw(DOSSIER_SPICE / "ina823_offset_dc.raw")
    voff = raw["_x"]
    vout = np.real(v(raw, "out"))
    gain = np.gradient(vout, voff)
    g0 = float(np.interp(0.0, voff, gain))
    ok = np.abs(gain - g0) <= 0.05 * g0
    i0 = int(np.argmin(np.abs(voff)))
    g = i0
    while g > 0 and ok[g - 1]:
        g -= 1
    d = i0
    while d < len(ok) - 1 and ok[d + 1]:
        d += 1
    vmin, vmax = voff[g], voff[d]
    res["offset"] = {"gain_0": g0, "voff_min_mV": vmin * 1e3, "voff_max_mV": vmax * 1e3,
                     "vout_min": float(vout[g]), "vout_max": float(vout[d])}

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 4.6), sharex=True,
                                 gridspec_kw={"height_ratios": [1.4, 1]})
    a1.axvspan(vmin * 1e3, vmax * 1e3, color=BLEU, alpha=0.08, lw=0)
    a1.plot(voff * 1e3, VREF + G1_THEO * voff, color=GRIS, ls="--", lw=1.2,
            label=f"Idéal : $V_{{ref}} + G_1 V_{{off}}$ ($G_1 = {virgule(G1_THEO)}$)")
    a1.plot(voff * 1e3, vout, color=BLEU, label="INA823 (modèle TI), sortie")
    for yv in (0, VDD):
        a1.axhline(yv, color=ENCRE2, lw=0.8, ls=":")
    a1.set_ylim(-0.6, 3.9)
    a1.set_ylabel("$V_{OUT}$ (V)")
    a1.legend(loc="upper left")
    a1.grid(True)
    a1.set_title("Sortie de l'INA823 en fonction de l'offset continu entre les électrodes")
    a2.axvspan(vmin * 1e3, vmax * 1e3, color=BLEU, alpha=0.08, lw=0)
    a2.plot(voff * 1e3, gain, color=BLEU)
    a2.set_ylim(-0.5, 11.5)
    a2.set_ylabel("Gain petit signal (V/V)")
    a2.set_xlabel("Offset différentiel continu $V_{off}$ (mV)")
    a2.grid(True)
    a2.text((vmin + vmax) / 2 * 1e3, 4.2,
            f"Zone linéaire (gain à ±5 %)\n{virgule(vmin * 1e3, 0).replace('-', '−')} à +{virgule(vmax * 1e3, 0)} mV",
            ha="center", va="center", color=ENCRE, fontsize=9)
    fig.tight_layout()
    fig.savefig(DOSSIER_FIG / "fig_6_1_offset_ina823.png")
    plt.close(fig)


def fig_bode_ina(res):
    raw = lire_raw(DOSSIER_SPICE / "ina823_bode.raw")
    f = raw["_x"]
    ad = v(raw, "out_d") / (v(raw, "pos_d") - v(raw, "neg_d"))
    acm = v(raw, "out_c") / ((v(raw, "pos_c") + v(raw, "neg_c")) / 2)
    ad_db, acm_db = db(ad), db(acm)
    trmc = ad_db - acm_db
    gmax, _, fh = coupures(f, ad_db)
    r = {"Ad_dB_1Hz": valeur_a(f, ad_db, 1.0), "Ad_dB_60Hz": valeur_a(f, ad_db, 60.0),
         "fH_INA_Hz": fh, "Acm_dB_60Hz": valeur_a(f, acm_db, 60.0),
         "TRMC_dB_1Hz": valeur_a(f, trmc, 1.0), "TRMC_dB_60Hz": valeur_a(f, trmc, 60.0),
         "TRMC_dB_1kHz": valeur_a(f, trmc, 1e3)}
    res["ina823"] = r

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 5.0), sharex=True)
    a1.plot(f, ad_db, color=BLEU, label=r"Mode différentiel $V(OUT)/(V(POS)-V(NEG))$")
    a1.plot(f, acm_db, color=ORANGE, label=r"Mode commun $V(OUT)/\frac{V(POS)+V(NEG)}{2}$")
    a1.axvline(60, color=ENCRE2, lw=0.8, ls=":")
    marquer(a1, 60, r["Ad_dB_60Hz"], f"{virgule(r['Ad_dB_60Hz'])} dB", dx=2.2, dy=-14)
    marquer(a1, 60, r["Acm_dB_60Hz"], f"{virgule(r['Acm_dB_60Hz'], 1).replace('-', '−')} dB", dx=2.2, dy=-14)
    a1.set_ylabel("Gain (dB)")
    a1.set_ylim(-130, 40)
    axes_freq(a1)
    a1.legend(loc="center left")
    a1.set_title("INA823, $G_1 = 9{,}93$ V/V, alimentation 3,3 V (modèle TI, LTspice)")
    a2.plot(f, trmc, color=AQUA, label=r"TRMC $= 20\log_{10}|A_d/A_{cm}|$")
    a2.axvline(60, color=ENCRE2, lw=0.8, ls=":")
    marquer(a2, 60, r["TRMC_dB_60Hz"], f"{virgule(r['TRMC_dB_60Hz'], 1)} dB à 60 Hz", dx=2.2, dy=-12)
    a2.set_ylabel("TRMC (dB)")
    a2.set_xlabel("Fréquence (Hz)")
    axes_freq(a2)
    a2.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(DOSSIER_FIG / "fig_6_2a_bode_ina823.png")
    plt.close(fig)


def fig_bode_filtres(res):
    raw = lire_raw(DOSSIER_SPICE / "filtres_bode.raw")
    f = raw["_x"]

    def rapport(sortie_p, sortie_n, entree_p, entree_n):
        return (v(raw, sortie_p) - v(raw, sortie_n)) / (v(raw, entree_p) - v(raw, entree_n))

    rfi_a = db(rapport("p1a", "n1a", "posa", "nega"))
    rfi_b = db(rapport("p1b", "n1b", "posb", "negb"))
    rfi_c = db(rapport("p1c", "n1c", "posc", "negc"))
    hp = db(v(raw, "hout") / v(raw, "hin"))
    lp5 = db(v(raw, "lout5") / v(raw, "lin5"))
    lp3 = db(v(raw, "lout3") / v(raw, "lin3"))
    rc = db(v(raw, "oout") / v(raw, "oin"))

    r = {}
    r["RFI_fc_Hz"] = coupures(f, rfi_a)[2]
    r["RFI_elec_fc_Hz"] = coupures(f, rfi_b)[2]
    r["RFI_prec_elec_fc_Hz"] = coupures(f, rfi_c)[2]
    r["RFI_elec_dB_150Hz"] = valeur_a(f, rfi_b, 150.0)
    r["RFI_prec_elec_dB_150Hz"] = valeur_a(f, rfi_c, 150.0)
    r["HP_fL_Hz"] = coupures(f, hp)[1]
    g5, _, fh5 = coupures(f, lp5)
    r["LP5_gain_dB"], r["LP5_fH_Hz"] = g5, fh5
    r["LP5_plancher_dB"] = valeur_a(f, lp5, 1e5)
    r["LP3_gain_dB_10Hz"] = valeur_a(f, lp3, 10.0)
    r["RC_fc_Hz"] = coupures(f, rc)[2]
    res["filtres"] = r

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 6.0), sharex=True)
    a1.plot(f, rfi_a, color=BLEU, label=f"RFI retenu, source idéale ($f_c$ = {fmt_hz(r['RFI_fc_Hz'])})")
    a1.plot(f, rfi_b, color=ORANGE, label=f"RFI retenu, à travers les électrodes ($f_c$ = {fmt_hz(r['RFI_elec_fc_Hz'])})")
    a1.plot(f, rfi_c, color=GRIS, ls="--",
            label=f"RFI précédent (1 nF / 10 nF), électrodes ($f_c$ = {fmt_hz(r['RFI_prec_elec_fc_Hz'])})")
    a1.axhline(-3, color=ENCRE2, lw=0.8, ls=":")
    a1.axvspan(0.5, 150, color=AQUA, alpha=0.07, lw=0)
    a1.text(9, -9, "bande ECG", color=ENCRE2, fontsize=8.5, ha="center")
    a1.set_ylim(-30, 3)
    a1.set_ylabel("Gain (dB)")
    a1.set_title("(a) Filtre RFI différentiel")
    axes_freq(a1)
    a1.legend(loc="lower left")

    a2.plot(f, hp, color=BLEU, label=f"Passe-haut $R_H C_H$ ($f_L$ = {fmt_hz(r['HP_fL_Hz'])})")
    a2.plot(f, lp5 - db(G2_THEO), color=ORANGE,
            label=f"Passe-bas et gain, TL082 ±5 V, normalisé ($f_H$ = {fmt_hz(fh5)})")
    a2.plot(f, rc, color=AQUA, label=f"RC de sortie 1 kΩ / 100 nF ($f_c$ = {fmt_hz(r['RC_fc_Hz'])})")
    a2.plot(f, lp3 - db(G2_THEO), color=GRIS, ls="--", label="Passe-bas et gain, TL082 sur 3,3 V, normalisé")
    a2.axhline(-3, color=ENCRE2, lw=0.8, ls=":")
    a2.set_ylim(-84, 6)
    a2.set_ylabel("Gain normalisé (dB)")
    a2.set_xlabel("Fréquence (Hz)")
    a2.set_title("(b) Passe-haut, passe-bas et gain, RC de sortie")
    axes_freq(a2)
    a2.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(DOSSIER_FIG / "fig_6_2b_bode_filtres.png")
    plt.close(fig)


def fig_bode_global(res):
    raw = lire_raw(DOSSIER_SPICE / "bioampli_ac.raw")
    f = raw["_x"]
    ha = v(raw, "outa") / (v(raw, "inpa") - v(raw, "inna"))
    hpeau = v(raw, "outa") / v(raw, "vecg")
    hb = v(raw, "outb") / (v(raw, "inpb") - v(raw, "innb"))
    ha_db, hpeau_db, hb_db = db(ha), db(hpeau), db(hb)
    gmax, fl, fh = coupures(f, ha_db)
    gp, flp, fhp = coupures(f, hpeau_db)
    phase = np.degrees(np.angle(ha))
    r = {"gain_max_dB": gmax, "gain_max_VV": 10 ** (gmax / 20), "fL_Hz": fl, "fH_Hz": fh,
         "peau_gain_max_dB": gp, "peau_fL_Hz": flp, "peau_fH_Hz": fhp,
         "phase_fL": valeur_a(f, phase, fl), "phase_fH": valeur_a(f, phase, fh),
         "phase_20Hz": valeur_a(f, phase, 20.0),
         "B33_gain_dB_10Hz": valeur_a(f, hb_db, 10.0),
         "gain_dB_1kHz": valeur_a(f, ha_db, 1e3), "gain_dB_10kHz": valeur_a(f, ha_db, 1e4)}
    res["global"] = r

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 5.4), sharex=True,
                                 gridspec_kw={"height_ratios": [1.6, 1]})
    a1.plot(f, ha_db, color=BLEU, label=r"TL082 en ±5 V : $V(OUT)/(V(POS)-V(NEG))$")
    a1.plot(f, hpeau_db, color=ORANGE, ls="--", label=r"TL082 en ±5 V, vu de la peau : $V(OUT)/V_{ECG}$")
    a1.plot(f, hb_db, color=GRIS, label="Tout sur 3,3 V (TL082 hors plage)")
    a1.axhline(gmax - 3, color=ENCRE2, lw=0.8, ls=":")
    marquer(a1, fl, gmax - 3, f"$f_L$ = {fmt_hz(fl)}", dx=0.12, dy=-9, ha="right")
    marquer(a1, fh, gmax - 3, f"$f_H$ = {fmt_hz(fh)}", dx=4.0, dy=-9)
    a1.text(3.0, gmax + 3.5, f"Gain max {virgule(gmax)} dB ({10 ** (gmax / 20):.0f} V/V)",
            color=ENCRE, fontsize=8.5, ha="center")
    a1.set_ylim(-40, 72)
    a1.set_ylabel("Gain (dB)")
    a1.set_title("Bioamplificateur complet (modèles TI, LTspice)")
    axes_freq(a1)
    a1.legend(loc="lower center")
    a2.plot(f, phase, color=BLEU)
    a2.axvline(fl, color=ENCRE2, lw=0.8, ls=":")
    a2.axvline(fh, color=ENCRE2, lw=0.8, ls=":")
    a2.set_ylim(-200, 200)
    a2.set_yticks([-180, -90, 0, 90, 180])
    a2.set_ylabel("Phase (°)")
    a2.set_xlabel("Fréquence (Hz)")
    axes_freq(a2)
    fig.tight_layout()
    fig.savefig(DOSSIER_FIG / "fig_6_2c_bode_global.png")
    plt.close(fig)


def marge_drl(f, z_drl, ro=100e3, cbody=200e-12, ccm=100e-12, rrfi=1e3):
    """Gain de boucle du DRL extrait de la simulation LTspice.

    z_drl est l'impedance du corps (V_corps / i_d) avec le DRL actif. Avec l'AOP
    du DRL inactif (sortie a la masse AC), le corps ne voit que des elements
    passifs (memes valeurs que PATIENT et BIOAMPLI) : Z_p. On a alors
    Z_drl = Z_p / (1 + L), d'ou L = Z_p / Z_drl - 1.
    """
    s = 2j * np.pi * f

    def ze(rd, cd=47e-9, rs=200.0):
        return rs + rd / (1 + s * rd * cd)

    y = (1 / (ro + ze(51e3)) + s * cbody
         + 1 / (ze(51e3) + rrfi + 1 / (s * ccm)) + 1 / (ze(40.8e3) + rrfi + 1 / (s * ccm)))
    boucle = (1 / y) / z_drl - 1
    k = int(np.where(np.abs(boucle) < 1)[0][0])
    return {"L_dB_60Hz": valeur_a(f, db(boucle), 60.0), "f_unite_Hz": float(f[k]),
            "marge_phase_deg": float(180 + np.degrees(np.angle(boucle[k]))),
            "S_max": float(np.max(np.abs(1 / (1 + boucle))))}


def fig_drl(res):
    raw = lire_raw(DOSSIER_SPICE / "drl_ac.raw")
    f = raw["_x"]
    # V(VID) = 0,2 V correspond a i_d = 0,2 uA ; la reponse AC est deja a l'echelle
    vcm_a = np.abs(v(raw, "bodya")) * 1e3
    vcm_b = np.abs(v(raw, "bodyb")) * 1e3
    out_a = np.abs(v(raw, "outa")) * 1e3
    out_b = np.abs(v(raw, "outb")) * 1e3
    out_c = np.abs(v(raw, "outc")) * 1e3
    r = {"Vcm_DRL_mV_60": valeur_a(f, vcm_a, 60.0), "Vcm_sans_mV_60": valeur_a(f, vcm_b, 60.0),
         "Vout_DRL_mV_60": valeur_a(f, out_a, 60.0), "Vout_sans_mV_60": valeur_a(f, out_b, 60.0),
         "Vout_DRL_RFIprec_mV_60": valeur_a(f, out_c, 60.0),
         "pic_Vcm_DRL_mV": float(vcm_a.max()), "f_pic_Vcm_DRL_Hz": float(f[np.argmax(vcm_a)])}
    r["attenuation_Vcm_dB_60"] = 20 * np.log10(r["Vcm_sans_mV_60"] / r["Vcm_DRL_mV_60"])
    r["attenuation_Vout_dB_60"] = 20 * np.log10(r["Vout_sans_mV_60"] / r["Vout_DRL_mV_60"])
    r.update(marge_drl(f, v(raw, "bodya") / 0.2e-6))
    res["drl"] = r

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 5.4), sharex=True)
    a1.plot(f, vcm_b, color=ORANGE, label="Sans DRL (jambe droite reliée à REF)")
    a1.plot(f, vcm_a, color=BLEU, label="Avec DRL")
    a1.axvline(60, color=ENCRE2, lw=0.8, ls=":")
    marquer(a1, 60, r["Vcm_sans_mV_60"], f"{virgule(r['Vcm_sans_mV_60'])} mV", dx=2.2)
    marquer(a1, 60, r["Vcm_DRL_mV_60"], f"{virgule(r['Vcm_DRL_mV_60'])} mV", dx=2.2)
    a1.set_yscale("log")
    a1.set_ylabel("$|V_{cm}|$ du corps (mV)")
    a1.set_title("Courant de déplacement $i_d$ = 0,2 µA, électrodes déséquilibrées de 20 %")
    axes_freq(a1)
    a1.legend(loc="lower left")
    a2.plot(f, out_b, color=ORANGE, label="Sans DRL")
    a2.plot(f, out_a, color=BLEU, label="Avec DRL")
    a2.plot(f, out_c, color=GRIS, ls="--", label="Avec DRL, filtre RFI précédent (1 nF / 10 nF)")
    a2.axvline(60, color=ENCRE2, lw=0.8, ls=":")
    marquer(a2, 60, r["Vout_sans_mV_60"], f"{virgule(r['Vout_sans_mV_60'])} mV", dx=2.2)
    marquer(a2, 60, r["Vout_DRL_mV_60"], f"{virgule(r['Vout_DRL_mV_60'])} mV", dx=2.2)
    a2.set_yscale("log")
    a2.set_ylabel("$|V_{OUT}|$ parasite (mV)")
    a2.set_xlabel("Fréquence (Hz)")
    axes_freq(a2)
    a2.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(DOSSIER_FIG / "fig_6_2d_drl.png")
    plt.close(fig)


def fig_temporel(res):
    raw = lire_raw(DOSSIER_SPICE / "bioampli_tran.raw")
    t = raw["_x"]
    m = (t >= 1.0) & (t <= 4.0)
    tt = t[m] - 1.0
    vecg = np.real(v(raw, "vecg"))[m] * 1e3
    ina = np.real(v(raw, "inaa"))[m]
    outa = np.real(v(raw, "outa"))[m]
    outb = np.real(v(raw, "outb"))[m]
    body = np.real(v(raw, "bodya"))[m]
    r = {"INA_moy_V": float(ina.mean()), "INA_min_V": float(ina.min()), "INA_max_V": float(ina.max()),
         "OUTA_min_V": float(outa.min()), "OUTA_max_V": float(outa.max()), "OUTA_moy_V": float(outa.mean()),
         "OUTB_min_V": float(outb.min()), "OUTB_max_V": float(outb.max()),
         "BODYA_moy_V": float(body.mean())}
    res["temporel"] = r

    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(7.0, 6.4), sharex=True,
                                      gridspec_kw={"height_ratios": [1, 1, 1.4]})
    a1.plot(tt, vecg, color=ENCRE, lw=1.0)
    a1.set_ylabel("$V_{ECG}$ (mV)")
    a1.set_title("Entrée différentielle à la peau : ECG + offset de 50 mV + respiration")
    a1.grid(True)
    a2.plot(tt, ina, color=BLEU, lw=1.0)
    a2.set_ylim(ina.mean() - 0.04, ina.mean() + 0.04)
    a2.text(0.01, 0.12, f"Moyenne {virgule(ina.mean())} V = 1,65 V + 9,93 × 50 mV",
            transform=a2.transAxes, color=ENCRE, fontsize=8.5)
    a2.set_ylabel("$V_{INA}$ (V)")
    a2.set_title("Sortie de l'INA823 ($G_1$ = 9,93) : offset amplifié, sans saturation")
    a2.grid(True)
    a3.plot(tt, outb, color=GRIS, lw=1.2, label="Tout sur 3,3 V : TL082 bloqué")
    a3.plot(tt, outa, color=BLEU, lw=1.2, label="TL082 en ±5 V (montage proposé)")
    for yv in (0, VDD):
        a3.axhline(yv, color=ORANGE, lw=0.9, ls=":")
    a3.axhline(VREF, color=ENCRE2, lw=0.8, ls="--")
    a3.text(2.98, VDD - 0.05, "Pleine échelle ADC 3,3 V", color=ENCRE2, fontsize=8, ha="right", va="top")
    a3.set_ylim(-0.2, 3.6)
    a3.set_ylabel("$V_{OUT}$ (V)")
    a3.set_xlabel("Temps (s), après 1 s de simulation")
    a3.set_title("Sortie vers l'ADC du BioBLU")
    a3.grid(True)
    a3.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(DOSSIER_FIG / "fig_6_2e_temporel.png")
    plt.close(fig)


def fig_schema():
    fig, ax = plt.subplots(figsize=(9.0, 4.0))
    ax.axis("off")

    def bloc(x, y, w, h, titre, sous_titre, fond, bord):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.03",
                                    fc=fond, ec=bord, lw=1.3))
        ax.text(x + w / 2, y + h * 0.74, titre, ha="center", va="center", weight="bold",
                fontsize=8.4, color=ENCRE)
        ax.text(x + w / 2, y + h * 0.30, sous_titre, ha="center", va="center", fontsize=7.2,
                color=ENCRE2)

    def fleche(x1, y1, x2, y2, style="-|>", ls="-"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle=style, color=ENCRE2, lw=1.2, ls=ls))

    bw, bh, y = 0.15, 0.34, 0.52
    xs = [0.03, 0.225, 0.42, 0.615, 0.81]
    fond_bio, bord_bio = "#e8f1fb", BLEU
    fond_op, bord_op = "#fdeee7", ORANGE
    bloc(xs[0], y, bw, bh, "1. Filtre RFI", "1 kΩ, 100 pF / 1 nF\n$f_{c,diff}$ ≈ 76 kHz", fond_bio, bord_bio)
    bloc(xs[1], y, bw, bh, "2. INA823", "$R_G$ = 2 × 5,6 kΩ\n$G_1$ = 9,93 V/V", fond_bio, bord_bio)
    bloc(xs[2], y, bw, bh, "3. Passe-haut", "330 kΩ, 1 µF\n$f_L$ = 0,48 Hz", fond_bio, bord_bio)
    bloc(xs[3], y, bw, bh, "4. Passe-bas\net gain", "TL082, $G_2$ = 101\n$f_H$ = 159 Hz", fond_op, bord_op)
    bloc(xs[4], y, bw, bh, "5. Sortie ADC", "1 kΩ, 100 nF\ndiodes vers 0 / 3,3 V", fond_bio, bord_bio)
    bloc(0.20, 0.06, 0.42, 0.27, "6. Driven-Right-Leg (TL082)",
         "point milieu de $R_G$, inverseur $A$ = 100 (pôle 48 Hz)\n$R_o$ = 100 kΩ vers l'électrode de jambe droite",
         fond_op, bord_op)
    yc = y + bh / 2
    fleche(-0.03, yc, xs[0], yc)
    ax.text(-0.03, yc + 0.04, "POS\nNEG", fontsize=7.5, color=ENCRE, ha="left", va="bottom")
    for a, b in zip(xs[:-1], xs[1:]):
        fleche(a + bw, yc, b, yc)
    fleche(xs[4] + bw, yc, 1.02, yc)
    ax.text(1.0, yc + 0.04, "ADC", fontsize=7.5, color=ENCRE, ha="right", va="bottom")
    fleche(xs[1] + bw / 2, y, xs[1] + bw / 2, 0.33, ls="--")
    fleche(0.20, 0.195, -0.03, 0.195, ls="--")
    ax.text(0.085, 0.215, "jambe droite", fontsize=7.5, color=ENCRE2, ha="center")
    ax.text(0.83, 0.20, "Bleu : BioBLU, AVCC = 3,3 V, REF = 1,65 V\nOrange : TL082 en ±5 V (labo)",
            fontsize=7.6, color=ENCRE, ha="center", va="center")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(0, 0.95)
    fig.savefig(DOSSIER_FIG / "fig_6_schema_bioampli.png")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--no-run", action="store_true", help="ne pas relancer LTspice")
    args = p.parse_args()
    DOSSIER_FIG.mkdir(exist_ok=True)
    ecrire_pwl()
    if not args.no_run:
        exe = trouver_ltspice()
        if not exe:
            sys.exit("LTspice introuvable : installez-le ou definissez LTSPICE=/chemin/vers/LTspice")
        lancer_ltspice(exe)
    res = {}
    fig_schema()
    fig_offset(res)
    fig_bode_ina(res)
    fig_bode_filtres(res)
    fig_bode_global(res)
    fig_drl(res)
    fig_temporel(res)
    with open(DOSSIER_SPICE / "resultats_simulation.json", "w") as f:
        json.dump(res, f, indent=2, ensure_ascii=False, default=float)
    print(json.dumps(res, indent=2, ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()
