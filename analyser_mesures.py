#!/usr/bin/env python3
"""Laboratoire 4 (GEL-7072) : traitement des mesures de la question 7.

Lit les captures binaires de l'oscilloscope (Data_et_notes/SDS0000x.bin),
extrait le gain et la phase de chaque balayage par demodulation synchrone,
calcule les frequences de coupure et trace les figures du rapport en
superposant les simulations LTspice de la question 6.

Usage :
    python analyser_mesures.py

Sorties : figures/fig_7_*.png et Data_et_notes/resultats_mesures.json
Dependances : numpy, matplotlib (et generer_figures_simulation.py pour la
lecture des .raw de LTspice, qui doivent exister dans ltspice/).
"""
import json
import struct
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import generer_figures_simulation as sim   # palette, style et lecture des .raw

ICI = Path(__file__).resolve().parent
DONNEES = ICI / "Data_et_notes"
FIG = ICI / "figures"
BLEU, ORANGE, AQUA, GRIS = sim.BLEU, sim.ORANGE, sim.AQUA, sim.GRIS
ENCRE, ENCRE2 = sim.ENCRE, sim.ENCRE2

# ---------------------------------------------------------------------------
# Lecture des fichiers binaires Siglent (8 bits, en-tete de 0x800 octets)
# ---------------------------------------------------------------------------
CODES_PAR_DIV = 25.0
DIV_HORIZ = 14


def lire_siglent(chemin):
    entete = open(chemin, "rb").read(0x800)
    actives = struct.unpack("<4i", entete[0:16])
    vdiv = [struct.unpack("<d", entete[0x10 + 16 * i:0x18 + 16 * i])[0] for i in range(4)]
    decal = [struct.unpack("<d", entete[0x50 + 16 * i:0x58 + 16 * i])[0] for i in range(4)]
    sonde = [struct.unpack("<d", entete[0x11c + 8 * i:0x124 + 8 * i])[0] for i in range(4)]
    tdiv = struct.unpack("<d", entete[0xd4:0xdc])[0]
    retard = struct.unpack("<d", entete[0xe4:0xec])[0]
    n = struct.unpack("<i", entete[0xf4:0xf8])[0]
    fe = struct.unpack("<d", entete[0xf8:0x100])[0]
    brut = np.memmap(chemin, dtype=np.uint8, mode="r", offset=0x800)
    voies, k = {}, 0
    for i in range(4):
        if actives[i]:
            code = brut[k * n:(k + 1) * n]
            volts = ((code.astype(np.float32) - 128.0) / CODES_PAR_DIV * vdiv[i] - decal[i]) * sonde[i]
            sature = (code <= 2) | (code >= 253)
            voies[i + 1] = (volts, sature)
            k += 1
    t = np.arange(n) / fe - tdiv * DIV_HORIZ / 2 - retard
    return dict(t=t, voies=voies, fe=fe, vdiv=[v * s for v, s in zip(vdiv, sonde)], sonde=sonde)


def decimer(d, q):
    """Moyenne par blocs de q echantillons (filtre anti-repliement simple)."""
    t, fe = d["t"], d["fe"]
    x, sat = {}, {}
    for c, (v, s) in d["voies"].items():
        v = np.asarray(v, dtype=np.float64)
        if q > 1:
            m = len(v) // q
            v = v[:m * q].reshape(m, q).mean(1)
            s = s[:m * q].reshape(m, q).any(1)
        x[c], sat[c] = v, np.asarray(s)
    if q > 1:
        m = len(t) // q
        t, fe = t[:m * q].reshape(m, q).mean(1), fe / q
    return t, x, sat, fe


# ---------------------------------------------------------------------------
# Demodulation synchrone le long d'un balayage logarithmique
# ---------------------------------------------------------------------------
def demoduler(t, x, sat, ref, f1, f2, T, t0, fcs, ncyc, nmin, dmax, fe, segs=(0,), recherche=0.03):
    """Pour chaque frequence fc : fenetre centree sur l'instant ou le generateur
    passe par fc (loi f = f1 (f2/f1)^((t-t0)/T)), correction locale de la
    frequence sur la voie ref, puis ajustement par moindres carres de
    a cos(phi) + b sin(phi) + c + d t sur chaque voie. Le modele etant evalue
    aux instants d'echantillonnage, la methode reste valable au-dela de la
    frequence de Nyquist (repliement)."""
    s = np.log(f2 / f1) / T
    res = []
    for k in segs:
        for fc in fcs:
            tc = t0 + k * T + np.log(fc / f1) / s
            Tw = max(ncyc / fc, nmin / fe)
            if Tw * s > dmax:
                Tw = max(dmax / s, min(Tw, 1.2 / fc))
            a, b = tc - Tw / 2, tc + Tw / 2
            if a < max(t[0], t0 + k * T) or b > min(t[-1], t0 + (k + 1) * T):
                continue
            i0, i1 = np.searchsorted(t, [a, b])
            N = i1 - i0
            if N < 20:
                continue
            tt = t[i0:i1] - tc
            phi0 = 2 * np.pi * (fc * tt + 0.5 * s * fc * tt ** 2)
            df = 0.0
            if fc * Tw >= 10:
                xr = x[ref][i0:i1] - x[ref][i0:i1].mean()
                y = xr * np.exp(-1j * phi0) * np.hanning(N)
                nfft = 1 << int(np.ceil(np.log2(N * 16)))
                Y = np.fft.fft(y, nfft)
                fr = np.fft.fftfreq(nfft, 1 / fe)
                m = np.abs(fr) <= recherche * fc
                df = fr[m][np.argmax(np.abs(Y[m]))]
            phi = phi0 + 2 * np.pi * df * tt
            M = np.vstack([np.cos(phi), np.sin(phi), np.ones_like(tt), tt / Tw]).T
            r = dict(f=fc + df, k=k, n=N)
            for c, v in x.items():
                coef, *_ = np.linalg.lstsq(M, v[i0:i1], rcond=None)
                r[c] = coef[0] - 1j * coef[1]
                r[f"bruit{c}"] = float(np.std(v[i0:i1] - M @ coef))
                r[f"sat{c}"] = float(sat[c][i0:i1].mean())
            res.append(r)
    return res


def trouver_t0(t, x, sat, ref, f1, f2, T, t0_approx, bande, fe, ncyc, nmin, dmax):
    """Debut du balayage : maximise l'amplitude demodulee de la voie ref."""
    fcs = np.geomspace(bande[0], bande[1], 25)
    meilleur = (None, t0_approx)
    for largeur, npts in ((0.03 * T, 61), (0.002 * T, 41), (0.0002 * T, 41)):
        centre = meilleur[1]
        for t0 in np.linspace(centre - largeur, centre + largeur, npts):
            r = demoduler(t, {ref: x[ref]}, {ref: sat[ref]}, ref, f1, f2, T, t0, fcs,
                          ncyc, nmin, dmax, fe, recherche=0)
            score = sum(abs(e[ref]) for e in r)
            if meilleur[0] is None or score > meilleur[0]:
                meilleur = (score, t0)
    return meilleur[1]


# ---------------------------------------------------------------------------
# Description des mesures (notes de laboratoire du 2026-11-17)
# ---------------------------------------------------------------------------
# entree/sortie : voie de l'oscilloscope ; ref : voie qui sert de reference de phase
MESURES = {
    "rfi_diff": dict(fichier="SDS00002.bin", entree=2, sortie=1, ref=1, f1=100, f2=1e6, T=14,
                     t0=-6.79, q=1, ncyc=6, nmin=2000, dmax=0.05, segs=(0,), bande=(1e3, 1e5),
                     bruit_max=1.0, fmin=800),
    "rfi_mc": dict(fichier="SDS00003.bin", entree=2, sortie=1, ref=2, f1=1e4, f2=5e6, T=14,
                   t0=-6.64, q=1, ncyc=6, nmin=2000, dmax=0.05, segs=(0,), bande=(1.2e4, 2e5),
                   bruit_max=5.0, fmin=1.2e4),
    "passe_haut": dict(fichier="SDS00004.bin", entree=2, sortie=1, ref=2, f1=0.1, f2=10, T=200,
                       t0=-139.1, q=25, ncyc=2, nmin=50, dmax=0.3, segs=(0, 1), bande=(1, 9),
                       bruit_max=1.0, fmin=0.1),
    "gain_pb": dict(fichier="SDS00005.bin", entree=2, sortie=1, ref=1, f1=10, f2=5e4, T=20,
                    t0=-13.4, q=1, ncyc=6, nmin=2000, dmax=0.05, segs=(0,), bande=(30, 1e3),
                    bruit_max=1.0, fmin=30),
    "rc_sortie": dict(fichier="SDS00006.bin", entree=2, sortie=1, ref=2, f1=100, f2=5e4, T=2,
                      t0=-1.03, q=4, ncyc=6, nmin=500, dmax=0.05, segs=(-1, 0), bande=(120, 1.2e3),
                      bruit_max=0.6, fmin=100, fmax=2.5e3),
    "chaine_bas": dict(fichier="SDS00007.bin", entree=2, sortie=1, ref=1, f1=0.1, f2=20, T=200,
                       t0=-139.7, q=25, ncyc=2, nmin=50, dmax=0.3, segs=(0, 1), bande=(1, 15),
                       bruit_max=3.0, fmin=0.1),
    "chaine_haut": dict(fichier="SDS00008.bin", entree=2, sortie=1, ref=1, f1=10, f2=2e3, T=20,
                        t0=-14.0, q=10, ncyc=6, nmin=500, dmax=0.05, segs=(0,), bande=(15, 500),
                        bruit_max=1.0, fmin=14),
}


def mesurer(nom, p):
    d = lire_siglent(DONNEES / p["fichier"])
    t, x, sat, fe = decimer(d, p["q"])
    t0 = trouver_t0(t, x, sat, p["ref"], p["f1"], p["f2"], p["T"], p["t0"], p["bande"], fe,
                    max(p["ncyc"], 8), p["nmin"], p["dmax"])
    fcs = np.geomspace(p["f1"], p["f2"], int(40 * np.log10(p["f2"] / p["f1"])) + 1)
    r = demoduler(t, x, sat, p["ref"], p["f1"], p["f2"], p["T"], t0, fcs, p["ncyc"], p["nmin"],
                  p["dmax"], fe, segs=p["segs"])
    e, s = p["entree"], p["sortie"]
    f = np.array([q["f"] for q in r])
    ai = np.array([q[e] for q in r])
    ao = np.array([q[s] for q in r])
    bruit_in = np.array([q[f"bruit{e}"] for q in r]) / np.abs(ai)
    satur = np.array([max(q[f"sat{e}"], q[f"sat{s}"]) for q in r])
    ok = (bruit_in < p["bruit_max"]) & (satur < 1e-3) & (f >= p["fmin"]) & (f <= p.get("fmax", np.inf))
    ordre = np.argsort(f)
    return dict(nom=nom, t0=t0, f=f[ordre], H=(ao / ai)[ordre], vin=2 * np.abs(ai)[ordre],
                vout=2 * np.abs(ao)[ordre], ok=ok[ordre], seg=np.array([q["k"] for q in r])[ordre],
                ai=ai[ordre], ao=ao[ordre])


# ---------------------------------------------------------------------------
# Outils : lissage, coupures, ajustements
# ---------------------------------------------------------------------------
def lisser(y, n=5):
    if len(y) < n:
        return y
    pad = n // 2
    yp = np.pad(y, pad, mode="edge")
    return np.array([np.median(yp[i:i + n]) for i in range(len(y))])


def croisement(f, g, seuil, sens, depart=None):
    """Premiere frequence ou g passe sous le seuil (interpolation en log f),
    en partant de l'indice depart (par defaut le maximum de g).
    sens = +1 : vers les hautes frequences ; -1 : vers les basses."""
    i_max = int(np.argmax(g)) if depart is None else depart
    idx = range(i_max, len(f) - 1) if sens > 0 else range(i_max, 0, -1)
    for i in idx:
        j = i + sens
        if (g[i] - seuil) * (g[j] - seuil) <= 0 and g[i] != g[j]:
            u = (seuil - g[i]) / (g[j] - g[i])
            return float(np.exp(np.log(f[i]) + u * (np.log(f[j]) - np.log(f[i]))))
    return None


def ajuster_complexe(f, H, modele, grilles, gain_libre=True):
    """Ajuste un modele complexe (module et phase) : minimise |ln(H / (g m))|^2.
    grilles : liste de tableaux de parametres ; renvoie (g, parametres)."""
    meilleur = None
    for params in np.array(np.meshgrid(*grilles, indexing="ij")).reshape(len(grilles), -1).T:
        m = modele(f, *params)
        z = np.log(H / m)
        lg = np.mean(z.real) if gain_libre else 0.0
        err = np.sum((z.real - lg) ** 2 + np.angle(H / m) ** 2)
        if meilleur is None or err < meilleur[0]:
            meilleur = (err, np.exp(lg), params)
    return meilleur[1], meilleur[2]


def ajuster_pole(f, H, type_, g0=None):
    """Ajuste un filtre du 1er ordre (passe-bas ou passe-haut) : renvoie (g0, fc)."""
    meilleur = None
    for fc in np.geomspace(f.min() / 20, f.max() * 20, 4000):
        x = f / fc
        mod = 1 / np.sqrt(1 + x ** 2) if type_ == "pb" else x / np.sqrt(1 + x ** 2)
        g = g0 if g0 is not None else np.exp(np.mean(np.log(np.abs(H) / mod)))
        err = np.sum((np.log(np.abs(H)) - np.log(g * mod)) ** 2)
        if meilleur is None or err < meilleur[0]:
            meilleur = (err, g, fc)
    return meilleur[1], meilleur[2]


def ajuster_gain_pb(f, H):
    """G(f) = 1 + K / (1 + j f / fH) (non-inverseur avec Rf || Cf)."""
    meilleur = None
    for fh in np.geomspace(50, 500, 600):
        for K in np.linspace(30, 110, 401):
            mod = np.abs(1 + K / (1 + 1j * f / fh))
            err = np.sum((np.log(np.abs(H)) - np.log(mod)) ** 2)
            if meilleur is None or err < meilleur[0]:
                meilleur = (err, K, fh)
    return meilleur[1], meilleur[2]


def db(h):
    return 20 * np.log10(np.abs(h))


def deg(h):
    return np.degrees(np.angle(h))


# ---------------------------------------------------------------------------
# Simulations LTspice (question 6) pour comparaison
# ---------------------------------------------------------------------------
def simulations():
    rs = {n: sim.lire_raw(sim.DOSSIER_SPICE / f"{n}.raw") for n in ("ina823_bode", "filtres_bode", "bioampli_ac")}
    v = sim.v
    out = {}
    r = rs["ina823_bode"]
    out["f_ina"] = np.real(r["frequency"])
    out["ina"] = v(r, "out_d") / (v(r, "pos_d") - v(r, "neg_d"))
    r = rs["filtres_bode"]
    out["f_filt"] = np.real(r["frequency"])
    out["rfi"] = (v(r, "p1a") - v(r, "n1a")) / (v(r, "posa") - v(r, "nega"))
    out["hp"] = v(r, "hout") / v(r, "hin")
    out["pb"] = v(r, "lout5") / v(r, "lin5")
    out["rc"] = v(r, "oout") / v(r, "oin")
    r = rs["bioampli_ac"]
    out["f_ch"] = np.real(r["frequency"])
    out["chaine"] = v(r, "outa") / (v(r, "inpa") - v(r, "inna"))
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def axes_bode(ax_m, ax_p, titre):
    ax_m.set_title(titre, loc="left", color=ENCRE)
    for ax in (ax_m, ax_p):
        ax.set_xscale("log")
        ax.grid(True, which="major")
        ax.grid(True, which="minor", alpha=0.35)
    ax_m.set_ylabel("Gain (dB)")
    ax_p.set_ylabel("Phase (°)")
    ax_p.set_xlabel("Fréquence (Hz)")


def points(ax, f, y, ok, label="Mesure"):
    ax.plot(f[ok], y[ok], "o", ms=3.2, color=BLEU, mec="white", mew=0.4, label=label, zorder=3)


def marquer(ax, x, y, texte, dx=1.15, dy=0.0, ha="left"):
    ax.plot([x], [y], "o", ms=6, mfc="white", mec=ENCRE, mew=1.0, zorder=4)
    ax.annotate(texte, (x, y), xytext=(x * dx, y + dy), fontsize=8.5, color=ENCRE, ha=ha, va="center")


def fmt(x):
    if x is None:
        return "n/a"
    if x >= 1e6:
        return f"{x / 1e6:.3g} MHz".replace(".", ",")
    return sim.fmt_hz(x)


def phase_deroulee(h):
    return np.degrees(np.unwrap(np.angle(h)))


def main():
    S = simulations()
    M = {nom: mesurer(nom, p) for nom, p in MESURES.items()}
    R = {}

    # ---------------- RFI differentiel + INA823 (SDS00002)
    m = M["rfi_diff"]
    ok = m["ok"]
    bande = ok & (m["f"] > 1e3) & (m["f"] < 3e4)
    g1 = float(np.median(np.abs(m["H"][bande])))
    gl = lisser(np.abs(m["H"][ok]))
    fc_diff = croisement(m["f"][ok], gl, g1 / np.sqrt(2), +1)
    h_sim_tot = S["rfi"] * S["ina"]
    _, _, fc_sim_tot = sim.coupures(S["f_ina"], db(h_sim_tot))
    _, _, fc_ina = sim.coupures(S["f_ina"], db(S["ina"]))
    # pole equivalent du RFI une fois la reponse de l'INA823 (simulee) retiree
    hina = np.interp(m["f"][ok], S["f_ina"], np.abs(S["ina"]))
    sel = m["f"][ok] > 5e4
    _, fc_rfi_eq = ajuster_pole(m["f"][ok][sel], (np.abs(m["H"][ok]) / hina * sim.G1_THEO)[sel] / g1, "pb", g0=1.0)
    def gain_a(mm, f0):
        sel = mm["ok"]
        return float(np.interp(np.log(f0), np.log(mm["f"][sel]), np.abs(mm["H"][sel])))
    R["rfi_diff_points"] = {f"{int(f0)}": gain_a(m, f0) for f0 in (7.5e4, 1e5, 2e5, 4e5)}
    R["rfi_diff"] = dict(G1=g1, G1_dB=20 * np.log10(g1), fc_mesure=fc_diff, fc_sim=fc_sim_tot,
                         fc_ina_sim=fc_ina, fc_rfi_equivalent=fc_rfi_eq,
                         CD_equivalent_pF=1e12 / (2 * np.pi * 2e3 * fc_rfi_eq))

    # ---------------- RFI mode commun (SDS00003)
    m2 = M["rfi_mc"]
    ok2 = m2["ok"]
    ref = ok2 & (m2["f"] < 1e5)
    g_mc = float(np.median(np.abs(m2["H"][ref])))
    hn = np.abs(m2["H"]) / g_mc
    fc_mc = croisement(m2["f"][ok2], lisser(hn[ok2]), 1 / np.sqrt(2), +1)
    _, fc_mc_fit = ajuster_pole(m2["f"][ok2], hn[ok2], "pb", g0=1.0)
    R["rfi_mc"] = dict(gain_bande=g_mc, fc_mesure=fc_mc, fc_ajuste=fc_mc_fit,
                       C_effective_pF=1e12 / (2 * np.pi * 1e3 * (fc_mc or fc_mc_fit)),
                       fc_theorique=1 / (2 * np.pi * 1e3 * 100e-12))

    # ---------------- Passe-haut (SDS00004)
    m3 = M["passe_haut"]
    ok3 = m3["ok"]
    g_hp = float(np.max(lisser(np.abs(m3["H"][ok3]))))
    fl = croisement(m3["f"][ok3], lisser(np.abs(m3["H"][ok3])), g_hp / np.sqrt(2), -1)
    _, (fl_fit,) = ajuster_complexe(m3["f"][ok3], m3["H"][ok3],
                                    lambda f, fc: 1j * f / fc / (1 + 1j * f / fc), [np.geomspace(0.2, 1, 801)])
    _, fl_sim, _ = sim.coupures(S["f_filt"], db(S["hp"]))
    # couplage AC de l'oscilloscope : l'amplitude de l'entree (generateur, 1 Vpp constant) chute sous 1 Hz
    a0_ac, (f_ac,) = ajuster_complexe(m3["f"][ok3], m3["vin"][ok3].astype(complex),
                                      lambda f, fc: np.abs(1j * f / fc / (1 + 1j * f / fc)) + 0j,
                                      [np.geomspace(0.2, 3, 801)])
    R["passe_haut"] = dict(gain_max=g_hp, fL_mesure=fl, fL_ajuste=fl_fit, fL_sim=fl_sim,
                           phase_a_fL=float(np.interp(np.log(fl), np.log(m3["f"][ok3]), deg(m3["H"][ok3]))),
                           couplage_AC_fc=f_ac, couplage_AC_vin_pp=a0_ac)

    # ---------------- Gain et passe-bas (SDS00005)
    m4 = M["gain_pb"]
    ok4 = m4["ok"]
    g2max = float(np.max(lisser(np.abs(m4["H"][ok4]))))
    fh = croisement(m4["f"][ok4], lisser(np.abs(m4["H"][ok4])), g2max / np.sqrt(2), +1)
    sel4 = ok4 & (m4["f"] < 2e3)
    _, (K, fh_fit) = ajuster_complexe(m4["f"][sel4], m4["H"][sel4], lambda f, K, fh: 1 + K / (1 + 1j * f / fh),
                                      [np.linspace(40, 100, 241), np.geomspace(100, 300, 241)], gain_libre=False)
    _, _, fh_sim = sim.coupures(S["f_filt"], db(S["pb"]))
    R["gain_pb"] = dict(G2_max_mesure=g2max, G2_dB=20 * np.log10(g2max), fH_mesure=fh,
                        K_ajuste=K, G2_dc_ajuste=1 + K, fH_ajuste=fh_fit, fH_sim=fh_sim,
                        R1_effective_ohm=100e3 / K, pot_estime_ohm=100e3 / K - 1e3,
                        gain_haute_freq=float(np.abs(m4["H"][ok4][-5:]).mean()))

    # ---------------- RC de sortie (SDS00006)
    m5 = M["rc_sortie"]
    ok5 = m5["ok"]
    g_rc = float(np.max(lisser(np.abs(m5["H"][ok5]))))
    fc_rc = croisement(m5["f"][ok5], lisser(np.abs(m5["H"][ok5])), g_rc / np.sqrt(2), +1)
    g_rc_fit, (fc_rc_fit,) = ajuster_complexe(m5["f"][ok5], m5["H"][ok5], lambda f, fc: 1 / (1 + 1j * f / fc),
                                              [np.geomspace(1e3, 3e3, 801)])
    _, _, fc_rc_sim = sim.coupures(S["f_filt"], db(S["rc"]))
    R["rc_sortie"] = dict(gain_max=g_rc, fc_mesure=fc_rc, fc_ajuste=fc_rc_fit, gain_ajuste=g_rc_fit,
                          fc_sim=fc_rc_sim, C_effective_nF=1e9 / (2 * np.pi * 1e3 * fc_rc_fit),
                          f_max_exploitable=float(m5["f"][ok5].max()))

    # ---------------- Chaine complete (SDS00007 + SDS00008 + SDS00001)
    # L'entree de 4 mVpp ne couvre que 5 codes de l'oscilloscope (20 mV/div) : son amplitude
    # mesuree fenetre par fenetre est bruitee. Le generateur fournissant une amplitude constante,
    # on prend son amplitude mediane (balayage 10 Hz - 2 kHz), corrigee sous 10 Hz par le
    # couplage AC de l'oscilloscope mesure avec le passe-haut (meme sonde 10x). La phase reste
    # celle du rapport sortie / entree (le couplage AC agit de la meme facon sur les deux voies).
    mb, mh = M["chaine_bas"], M["chaine_haut"]
    vin0 = float(np.median(mh["vin"][mh["ok"]]))
    def h_ac(f):
        return np.abs(1j * f / f_ac / (1 + 1j * f / f_ac))
    for mm in (mb, mh):
        mm["Hr"] = mm["H"].copy()
        mm["H"] = mm["vout"] / (vin0 * h_ac(mm["f"])) * np.exp(1j * np.angle(mm["Hr"]))
    vin_bas_mesure = mb["vin"] / h_ac(mb["f"])
    fb = np.concatenate([mb["f"][mb["ok"] & (mb["f"] < 14)], mh["f"][mh["ok"]]])
    hb = np.concatenate([mb["H"][mb["ok"] & (mb["f"] < 14)], mh["H"][mh["ok"]]])
    o = np.argsort(fb)
    fb, hb = fb[o], hb[o]
    gl = lisser(np.abs(hb), 5)
    zone = (fb > 1) & (fb < 100)
    i_max = int(np.flatnonzero(zone)[np.argmax(gl[zone])])
    gmax = float(gl[i_max])
    fL_ch = croisement(fb, gl, gmax / np.sqrt(2), -1, depart=i_max)
    fH_ch = croisement(fb, gl, gmax / np.sqrt(2), +1, depart=i_max)
    selb = fb < 3
    _, fL_fit = ajuster_pole(fb[selb], np.abs(hb[selb]), "ph", g0=gmax)
    gs, fls, fhs = sim.coupures(S["f_ch"], db(S["chaine"]))
    # point fixe a 20 Hz (SDS00001)
    d1 = lire_siglent(DONNEES / "SDS00001.bin")
    t1, x1, _, _ = decimer(d1, 250)
    meilleur = None
    for f0 in np.linspace(19.5, 20.5, 1001):
        Mx = np.vstack([np.cos(2 * np.pi * f0 * t1), np.sin(2 * np.pi * f0 * t1), np.ones_like(t1)]).T
        c, *_ = np.linalg.lstsq(Mx, x1[2], rcond=None)
        if meilleur is None or np.hypot(c[0], c[1]) > meilleur[0]:
            meilleur = (np.hypot(c[0], c[1]), f0)
    f0 = meilleur[1]
    Mx = np.vstack([np.cos(2 * np.pi * f0 * t1), np.sin(2 * np.pi * f0 * t1), np.ones_like(t1)]).T
    a = {}
    for c in (1, 2):
        co, *_ = np.linalg.lstsq(Mx, x1[c], rcond=None)
        a[c] = co[0] - 1j * co[1]
    Mh = np.vstack([np.cos(4 * np.pi * f0 * t1), np.sin(4 * np.pi * f0 * t1)]).T
    ch2, *_ = np.linalg.lstsq(Mh, x1[1] - x1[1].mean(), rcond=None)
    h20 = a[1] / a[2]
    R["chaine"] = dict(vin0_pp=vin0, vin_bas_corrige_median=float(np.median(vin_bas_mesure[mb["f"] > 0.3])),
                       G_max=gmax, G_max_dB=20 * np.log10(gmax), f_G_max=float(fb[i_max]),
                       fL_mesure=fL_ch, fL_ajuste=fL_fit,
                       fH_mesure=fH_ch, largeur_bande=(fH_ch - fL_ch) if fL_ch and fH_ch else None,
                       G_sim_dB=gs, fL_sim=fls, fH_sim=fhs,
                       phase_bande_deg=float(np.median(deg(mh["H"][mh["ok"] & (mh["f"] < 30)]))),
                       point_20Hz=dict(f=f0, vin_pp=2 * abs(a[2]), vout_pp=2 * abs(a[1]), G=abs(h20),
                                       G_dB=20 * np.log10(abs(h20)), phase=float(np.degrees(np.angle(h20))),
                                       h2_pourcent=100 * np.hypot(*ch2) / abs(a[1]),
                                       vout_min=float(x1[1].min()), vout_max=float(x1[1].max())))

    # ---------------- ECG (SDS00009)
    d9 = lire_siglent(DONNEES / "SDS00009.bin")
    t9, x9, _, fe9 = decimer(d9, 1000)
    e = x9[1] - np.median(x9[1])
    seuil = 0.5 * e.min()
    pics = []
    for i in range(1, len(e) - 1):
        if e[i] < seuil and e[i] == e[max(0, i - 500):i + 500].min():
            if not pics or t9[i] - t9[pics[-1]] > 0.3:
                pics.append(i)
    rr = np.diff(t9[pics])
    X = np.fft.rfft(e * np.hanning(len(e)))
    fr = np.fft.rfftfreq(len(e), 1 / fe9)
    Mx = np.vstack([np.cos(2 * np.pi * 60 * t9), np.sin(2 * np.pi * 60 * t9), np.ones_like(t9)]).T
    c60, *_ = np.linalg.lstsq(Mx, e, rcond=None)
    # bruit rapide : ecart a une moyenne glissante de 10 ms, hors complexes QRS (+-60 ms)
    w = int(0.010 * fe9)
    lisse = np.convolve(e, np.ones(w) / w, mode="same")
    hors_qrs = np.ones(len(e), bool)
    for i in pics:
        hors_qrs[max(0, i - int(0.06 * fe9)):i + int(0.06 * fe9)] = False
    hors_qrs[:w] = hors_qrs[-w:] = False
    bruit = float(np.std((e - lisse)[hors_qrs]))
    R["ecg"] = dict(n_battements=len(pics), rr_s=rr.tolist(), fc_bpm=float(60 / rr.mean()) if len(rr) else None,
                    amplitude_R_V=float(-np.mean(e[pics])), crete_a_crete_V=float(e.max() - e.min()),
                    v60_pp_V=float(2 * np.hypot(c60[0], c60[1])), bruit_rapide_rms_V=bruit,
                    gain_suppose=gmax, amplitude_R_entree_mV=float(-np.mean(e[pics]) / gmax * 1e3),
                    bruit_entree_uV_rms=bruit / gmax * 1e6, v60_entree_uVpp=float(2 * np.hypot(c60[0], c60[1]) / gmax * 1e6))

    # ======================= Figures =======================
    # Figure 7.1a : filtre RFI
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 5.6), sharex="col",
                            gridspec_kw=dict(height_ratios=[1.35, 1], hspace=0.08, wspace=0.28))
    axes_bode(axs[0, 0], axs[1, 0], "(a) Mode différentiel : filtre RFI + INA823")
    axs[0, 0].plot(S["f_ina"], db(h_sim_tot), color=ORANGE, label="Simulation (conception)")
    axs[0, 0].plot(S["f_ina"], db(S["ina"]), color=GRIS, ls="--", lw=1.2, label="Simulation, INA823 seul")
    points(axs[0, 0], m["f"], db(m["H"]), ok)
    axs[0, 0].set_xlim(500, 1.2e6)
    axs[0, 0].set_ylim(4, 24)
    axs[1, 0].plot(S["f_ina"], phase_deroulee(h_sim_tot), color=ORANGE)
    axs[1, 0].plot(S["f_ina"], phase_deroulee(S["ina"]), color=GRIS, ls="--", lw=1.2)
    points(axs[1, 0], m["f"], deg(m["H"]), ok)
    axs[1, 0].set_ylim(-230, 20)
    marquer(axs[0, 0], fc_diff, 20 * np.log10(g1 / np.sqrt(2)), f"−3 dB : {fmt(fc_diff)}", dx=0.6, dy=-1.2, ha="right")
    axs[0, 0].legend(loc="lower left")

    axes_bode(axs[0, 1], axs[1, 1], "(b) Mode commun : sortie du filtre (+IN de l'INA823)")
    ff = np.geomspace(1e4, 6e6, 400)
    for C, sty, lab in ((100e-12, "-", "Théorie, C$_C$ = 100 pF"), (115e-12, ":", "Théorie, 100 pF + sonde 15 pF")):
        h = 1 / (1 + 2j * np.pi * ff * 1e3 * C)
        axs[0, 1].plot(ff, db(h), color=ORANGE, ls=sty, label=lab)
        axs[1, 1].plot(ff, deg(h), color=ORANGE, ls=sty)
    points(axs[0, 1], m2["f"], db(m2["H"] / g_mc), ok2, label="Mesure (normalisée)")
    points(axs[1, 1], m2["f"], deg(m2["H"]), ok2)
    axs[0, 1].set_xlim(1e4, 6e6)
    axs[0, 1].set_ylim(-12, 2)
    axs[1, 1].set_ylim(-80, 5)
    if fc_mc:
        marquer(axs[0, 1], fc_mc, -3, f"−3 dB : {fmt(fc_mc)}", dx=0.6, dy=-1.6, ha="right")
    axs[0, 1].legend(loc="lower left")
    fig.savefig(FIG / "fig_7_1a_bode_rfi.png")
    plt.close(fig)

    # Figure 7.1b : passe-haut, gain et passe-bas, RC de sortie
    fig, axs = plt.subplots(2, 3, figsize=(13.5, 5.6), sharex="col",
                            gridspec_kw=dict(height_ratios=[1.35, 1], hspace=0.08, wspace=0.3))
    axes_bode(axs[0, 0], axs[1, 0], "(a) Passe-haut")
    axs[0, 0].plot(S["f_filt"], db(S["hp"]), color=ORANGE, label="Simulation")
    for seg, mk in ((0, "o"), (1, "s")):
        sel = ok3 & (m3["seg"] == seg)
        axs[0, 0].plot(m3["f"][sel], db(m3["H"][sel]), mk, ms=3.2, color=BLEU, mec="white", mew=0.4,
                       label="Mesure, balayage 1" if seg == 0 else "Mesure, balayage 2")
        axs[1, 0].plot(m3["f"][sel], deg(m3["H"][sel]), mk, ms=3.2, color=BLEU, mec="white", mew=0.4)
    axs[1, 0].plot(S["f_filt"], deg(S["hp"]), color=ORANGE)
    axs[0, 0].set_xlim(0.08, 12)
    axs[0, 0].set_ylim(-16, 2)
    axs[1, 0].set_ylim(-5, 95)
    marquer(axs[0, 0], fl, 20 * np.log10(g_hp / np.sqrt(2)), f"−3 dB : {fmt(fl)}", dx=1.3)
    axs[0, 0].legend(loc="lower right")

    axes_bode(axs[0, 1], axs[1, 1], "(b) Gain et passe-bas (U2A)")
    axs[0, 1].plot(S["f_filt"], db(S["pb"]), color=ORANGE, label="Simulation (G$_2$ = 101)")
    hm = 1 + K / (1 + 1j * S["f_filt"] / fh_fit)
    axs[0, 1].plot(S["f_filt"], db(hm), color=GRIS, ls="--", lw=1.2,
                   label="Modèle ajusté")
    axs[1, 1].plot([], [])
    points(axs[0, 1], m4["f"], db(m4["H"]), ok4)
    axs[1, 1].plot(S["f_filt"], deg(S["pb"]), color=ORANGE)
    axs[1, 1].plot(S["f_filt"], deg(hm), color=GRIS, ls="--", lw=1.2)
    points(axs[1, 1], m4["f"], deg(m4["H"]), ok4)
    axs[0, 1].set_xlim(20, 6e4)
    axs[0, 1].set_ylim(-2, 44)
    axs[1, 1].set_ylim(-95, 5)
    marquer(axs[0, 1], fh, 20 * np.log10(g2max / np.sqrt(2)), f"−3 dB : {fmt(fh)}", dx=1.35)
    axs[0, 1].legend(loc="lower left")

    axes_bode(axs[0, 2], axs[1, 2], "(c) RC de sortie")
    axs[0, 2].plot(S["f_filt"], db(S["rc"]), color=ORANGE, label="Simulation")
    hr = 1 / (1 + 1j * S["f_filt"] / fc_rc_fit)
    axs[0, 2].plot(S["f_filt"], db(hr), color=GRIS, ls="--", lw=1.2, label=f"1$^{{er}}$ ordre ajusté ({fmt(fc_rc_fit)})")
    points(axs[0, 2], m5["f"], db(m5["H"]), ok5)
    axs[1, 2].plot(S["f_filt"], deg(S["rc"]), color=ORANGE)
    axs[1, 2].plot(S["f_filt"], deg(hr), color=GRIS, ls="--", lw=1.2)
    points(axs[1, 2], m5["f"], deg(m5["H"]), ok5)
    marquer(axs[0, 2], fc_rc, 20 * np.log10(g_rc / np.sqrt(2)), f"−3 dB : {fmt(fc_rc)}", dx=1.3)
    axs[0, 2].set_xlim(80, 2e4)
    axs[0, 2].set_ylim(-22, 2)
    axs[1, 2].set_ylim(-95, 5)
    axs[0, 2].legend(loc="lower left")
    fig.savefig(FIG / "fig_7_1b_bode_filtres.png")
    plt.close(fig)

    # Figure 7.1c : chaine complete
    fig, (am, ap) = plt.subplots(2, 1, figsize=(9.5, 5.8), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1.35, 1], hspace=0.08))
    axes_bode(am, ap, "")
    am.plot(S["f_ch"], db(S["chaine"]), color=ORANGE, label="Simulation (pot. à 0 Ω)")
    am.plot(S["f_ch"], db(S["chaine"]) - (gs - 20 * np.log10(gmax)), color=GRIS, ls="--", lw=1.2,
            label="Simulation ramenée au gain mesuré")
    for mm, mk, lab in ((mb, "o", "Mesure, 0,1 à 20 Hz"), (mh, "s", "Mesure, 10 Hz à 2 kHz")):
        sel = mm["ok"]
        am.plot(mm["f"][sel], db(mm["H"][sel]), mk, ms=3.2, color=BLEU, mec="white", mew=0.4, label=lab)
        ap.plot(mm["f"][sel], deg(-mm["H"][sel]), mk, ms=3.2, color=BLEU, mec="white", mew=0.4)
    am.plot([f0], [20 * np.log10(abs(h20))], "D", ms=6, color=AQUA, mec="white", mew=0.5,
            label="Point à 20 Hz, autre réglage du pot.")
    ap.plot(S["f_ch"], phase_deroulee(S["chaine"]), color=ORANGE)
    am.set_xlim(0.08, 3e3)
    am.set_ylim(28, 64)
    ap.set_ylim(-110, 100)
    marquer(am, fL_ch, 20 * np.log10(gmax / np.sqrt(2)), f"f$_L$ = {fmt(fL_ch)}", dx=0.75, ha="right")
    marquer(am, fH_ch, 20 * np.log10(gmax / np.sqrt(2)), f"f$_H$ = {fmt(fH_ch)}", dx=1.35)
    am.legend(loc="lower center", fontsize=8)
    ap.set_ylabel("Phase (°)\n(+180° retirés)")
    fig.savefig(FIG / "fig_7_1c_bode_chaine.png")
    plt.close(fig)

    # Figure 7.9 : ECG
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw=dict(width_ratios=[2.2, 1], wspace=0.25))
    a1.plot(t9 - t9[0], e, color=BLEU, lw=0.7)
    a1.plot(t9[pics] - t9[0], e[pics], "v", color=ENCRE, ms=5)
    a1.set_xlabel("Temps (s)")
    a1.set_ylabel("Sortie (V, couplage AC)")
    a1.set_title(f"(a) ECG à la sortie, DRL actif ({R['ecg']['fc_bpm']:.0f} battements/min)", loc="left")
    a1.grid(True)
    a1.set_xlim(0, t9[-1] - t9[0])
    spec = np.convolve(np.abs(X), np.ones(3) / 3, mode="same")
    a2.plot(fr, 20 * np.log10(spec / spec.max()), color=BLEU, lw=0.8)
    for fh_ in (60, 120, 180):
        a2.axvline(fh_, color=GRIS, ls=":", lw=1)
        a2.text(fh_ + 3, -3, f"{fh_} Hz", fontsize=7.5, color=ENCRE2, va="top")
    a2.set_xlim(0, 250)
    a2.set_ylim(-60, 2)
    a2.set_xlabel("Fréquence (Hz)")
    a2.set_ylabel("Amplitude (dB, relative)")
    a2.set_title("(b) Spectre", loc="left")
    a2.grid(True, which="major")
    fig.savefig(FIG / "fig_7_9_ecg.png")
    plt.close(fig)

    # ======================= Resultats =======================
    def nettoyer(o):
        if isinstance(o, dict):
            return {k: nettoyer(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [nettoyer(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        return o
    R["debut_balayages"] = {nom: float(m_["t0"]) for nom, m_ in M.items()}
    with open(DONNEES / "resultats_mesures.json", "w") as fh_:
        json.dump(nettoyer(R), fh_, indent=2, ensure_ascii=False)
    print(json.dumps(nettoyer(R), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
