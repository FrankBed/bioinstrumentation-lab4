#!/usr/bin/env python3
"""Laboratoire 4 (GEL-7072) : schema electrique du montage a realiser au laboratoire.

Les valeurs sont celles de ltspice/bioampli.sub (conception simulee a la question 6).
Sortie : figures/schema_montage.png et figures/schema_montage.pdf
Dependance : schemdraw (pip install schemdraw)
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import schemdraw
import schemdraw.elements as elm

ICI = Path(__file__).resolve().parent
SORTIE = ICI / "figures" / "schema_montage"

BLEU, ORANGE, GRIS = "#2a78d6", "#eb6834", "#52514e"
PIN = dict(fontsize=7.5, color=GRIS)


def tag(d, point, texte, sens="right", couleur="black"):
    """Etiquette de reseau (meme nom = meme noeud)."""
    sens = {"up": "right", "down": "right"}.get(sens, sens)  # etiquettes toujours horizontales
    t = getattr(elm.Tag(), sens)().at(point).label(texte, fontsize=8.5).color(couleur)
    d.add(t)
    return t


def num(d, point, texte, dx=0.0, dy=0.25):
    """Numero de broche."""
    d.add(elm.Label().at((point[0] + dx, point[1] + dy)).label(texte, **PIN))


def pt(a):
    return (float(a[0]), float(a[1]))


with schemdraw.Drawing(show=False) as d:
    d.config(unit=2.5, fontsize=9, lw=1.1, font="serif")

    yP, yN = 9.0, 3.0

    # ------------------------------------------------------------------ Electrodes et filtre RFI
    for y, nom in ((yP, "POS\n(bras gauche)"), (yN, "NEG\n(bras droit)")):
        d.add(elm.Dot(open=True, radius=0.18).at((0, y)).label(nom, loc="left", fontsize=8.5))
    d.add(elm.Resistor().at((0.2, yP)).to((3.2, yP)).label("1 kΩ", loc="top"))
    d.add(elm.Resistor().at((0.2, yN)).to((3.2, yN)).label("1 kΩ", loc="bottom"))
    d.add(elm.Line().at((3.2, yP)).to((7.6, yP)))
    d.add(elm.Line().at((3.2, yN)).to((7.6, yN)))
    for y in (yP, yN):
        d.add(elm.Dot().at((4.3, y)))
        d.add(elm.Capacitor().at((4.3, y)).down().length(1.8).label("100 pF", loc="bottom", ofst=0.15))
        d.add(elm.Ground())
    d.add(elm.Dot().at((6.0, yP)))
    d.add(elm.Dot().at((6.0, yN)))
    d.add(elm.Capacitor().endpoints((6.0, yP), (6.0, yN)).label("1 nF", loc="bottom", ofst=0.15))
    d.add(elm.Label().at((3.6, yP + 1.7)).label("Filtre RFI", fontsize=8.5, color=BLEU))

    # ------------------------------------------------------------------ INA823 (U1)
    ina = elm.Ic(size=(6.6, 5.0), pins=[
        elm.IcPin(name="+IN", pin="3", side="L", slot="3/4", anchorname="inp"),
        elm.IcPin(name="−IN", pin="2", side="L", slot="2/4", anchorname="inn"),
        elm.IcPin(name="+VS", pin="7", side="T", slot="1/1", anchorname="vp"),
        elm.IcPin(name="OUT", pin="6", side="R", slot="3/4", anchorname="out"),
        elm.IcPin(name="REF", pin="5", side="R", slot="2/4", anchorname="ref"),
        elm.IcPin(name="−VS", pin="4", side="B", slot="2/5", anchorname="vm"),
        elm.IcPin(name="RG", pin="1", side="B", slot="3/5", anchorname="rg1"),
        elm.IcPin(name="RG", pin="8", side="B", slot="4/5", anchorname="rg8"),
    ], edgepadding=1.3, leadlen=0.6).right().label("U1\nINA823", fontsize=10)
    d.add(ina.anchor("inp").at((7.6, yP)))
    d.add(elm.Line().at((7.6, yN)).to(pt(ina.inn)))
    d.add(elm.Line().at(pt(ina.vp)).up().length(0.5))
    tag(d, d.here, "AVCC 3,3 V", "up")
    d.add(elm.Line().at(pt(ina.ref)).right().length(0.6))
    tag(d, d.here, "REF", "right")
    d.add(elm.Line().at(pt(ina.vm)).down().length(0.4))
    d.add(elm.Ground())

    # RG = 2 x 5,6 k, point milieu vers le DRL
    rg1, rg8 = pt(ina.rg1), pt(ina.rg8)
    d.add(elm.Resistor().at(rg1).down().length(2.6).label("5,6 kΩ", loc="top", fontsize=8))
    b1 = d.here
    d.add(elm.Resistor().at(rg8).down().length(2.6).label("5,6 kΩ", loc="bottom", fontsize=8))
    b8 = d.here
    d.add(elm.Line().at(b1).to(b8))
    mid = ((b1[0] + b8[0]) / 2, b1[1])
    d.add(elm.Dot().at(mid))
    d.add(elm.Label().at((mid[0] - 1.3, mid[1] - 0.45)).label("point milieu de RG", fontsize=7.5, color=GRIS))

    # ------------------------------------------------------------------ Passe-haut
    y_out = pt(ina.out)[1]
    d.add(elm.Capacitor().at(pt(ina.out)).right().length(2.6).label("1 µF", loc="top"))
    hp = d.here
    d.add(elm.Dot().at(hp))
    d.add(elm.Resistor().at(hp).down().length(2.6).label("330 kΩ", loc="bottom"))
    tag(d, d.here, "REF", "down")
    d.add(elm.Label().at((hp[0] - 1.3, y_out + 1.5)).label("Passe-haut", fontsize=8.5, color=BLEU))

    # ------------------------------------------------------------------ Gain et passe-bas (U2A)
    xa = hp[0] + 3.4
    d.add(elm.Line().at(hp).to((xa, y_out)))
    u2a = elm.Opamp(leads=True).right().anchor("in2").at((xa, y_out)).label("U2A\nTL082", loc="center", ofst=(-0.35, 0), fontsize=8)
    d.add(u2a)
    num(d, pt(u2a.in2), "3", dx=0.45, dy=-0.55)
    num(d, pt(u2a.in1), "2", dx=0.45, dy=0.22)
    num(d, pt(u2a.out), "1", dx=-0.35, dy=0.22)
    d.add(elm.Line().at(pt(u2a.vd)).up().length(0.45))
    tag(d, d.here, "+5 V", "right")
    num(d, pt(u2a.vd), "8", dx=-0.3, dy=0.15)
    d.add(elm.Line().at(pt(u2a.vs)).down().length(0.45))
    tag(d, d.here, "−5 V", "right")
    num(d, pt(u2a.vs), "4", dx=-0.3, dy=-0.45)
    nfb = pt(u2a.in1)
    yfb = nfb[1] + 2.6
    d.add(elm.Line().at(nfb).to((nfb[0], yfb + 1.3)))
    d.add(elm.Dot().at((nfb[0], yfb)))
    d.add(elm.Resistor().at((nfb[0], yfb)).to((pt(u2a.out)[0] + 0.6, yfb)).label("100 kΩ", loc="bottom"))
    d.add(elm.Capacitor().at((nfb[0], yfb + 1.3)).to((pt(u2a.out)[0] + 0.6, yfb + 1.3)).label("10 nF", loc="top"))
    d.add(elm.Line().at((pt(u2a.out)[0] + 0.6, yfb + 1.3)).to((pt(u2a.out)[0] + 0.6, pt(u2a.out)[1])))
    d.add(elm.Dot().at((pt(u2a.out)[0] + 0.6, yfb)))
    d.add(elm.Line().at(pt(u2a.out)).right().length(0.6))
    nout = d.here
    d.add(elm.Dot().at(nout))
    # R1 + potentiometre vers REF
    d.add(elm.Resistor().at((nfb[0], yfb)).left().length(2.6).label("1 kΩ", loc="bottom"))
    d.add(elm.Potentiometer().left().length(2.8).label("pot. 1 kΩ (rhéostat)", loc="bottom", fontsize=8))
    tag(d, d.here, "REF", "left")
    d.add(elm.Label().at((nfb[0] - 2.8, yfb + 1.1)).label("Gain et passe-bas\nG2 = 51 à 101, fH = 159 Hz", fontsize=8.5, color=ORANGE))

    # ------------------------------------------------------------------ Sortie vers l'ADC
    d.add(elm.Resistor().at(nout).right().length(2.6).label("1 kΩ", loc="top"))
    adc = d.here
    d.add(elm.Dot().at(adc))
    d.add(elm.Schottky().at(adc).up().length(2.0).label("BAT54", loc="bottom", fontsize=7.5))
    tag(d, d.here, "AVCC", "up")
    d.add(elm.Capacitor().at(adc).down().length(2.0).label("100 nF", loc="bottom"))
    d.add(elm.Ground())
    d.add(elm.Line().at(adc).right().length(1.6))
    nd = d.here
    d.add(elm.Dot().at(nd))
    d.add(elm.Schottky().at((nd[0], nd[1] - 2.0)).up().length(2.0).label("BAT54", loc="bottom", fontsize=7.5))
    d.add(elm.Ground().at((nd[0], nd[1] - 2.0)))
    d.add(elm.Line().at(nd).right().length(1.0))
    tag(d, d.here, "ADC0", "right", BLEU)
    d.add(elm.Label().at((adc[0] + 0.6, adc[1] + 3.6)).label("Sortie ADC", fontsize=8.5, color=BLEU))

    # ------------------------------------------------------------------ Driven-Right-Leg (U2B)
    d.add(elm.Line().at(mid).down().length(2.6))
    d.add(elm.Resistor().right().length(3.0).label("10 kΩ", loc="top"))
    din = d.here
    d.add(elm.Dot().at(din))
    u2b = elm.Opamp(leads=True).right().anchor("in1").at((din[0] + 0.0, din[1])).label("U2B\nTL082", loc="center", ofst=(-0.35, 0), fontsize=8)
    d.add(u2b)
    num(d, pt(u2b.in1), "6", dx=0.45, dy=0.22)
    num(d, pt(u2b.in2), "5", dx=0.45, dy=-0.55)
    num(d, pt(u2b.out), "7", dx=-0.35, dy=0.22)
    d.add(elm.Line().at(pt(u2b.in2)).left().length(0.6))
    tag(d, d.here, "REF", "left")
    yfb2 = din[1] + 2.0
    d.add(elm.Line().at(din).to((din[0], yfb2 + 1.3)))
    d.add(elm.Dot().at((din[0], yfb2)))
    xo2 = pt(u2b.out)[0] + 0.6
    d.add(elm.Resistor().at((din[0], yfb2)).to((xo2, yfb2)).label("1 MΩ", loc="bottom"))
    d.add(elm.Capacitor().at((din[0], yfb2 + 1.3)).to((xo2, yfb2 + 1.3)).label("3,3 nF", loc="top"))
    d.add(elm.Line().at((xo2, yfb2 + 1.3)).to((xo2, pt(u2b.out)[1])))
    d.add(elm.Dot().at((xo2, yfb2)))
    d.add(elm.Line().at(pt(u2b.out)).to((xo2, pt(u2b.out)[1])))
    d.add(elm.Dot().at((xo2, pt(u2b.out)[1])))
    yrl = pt(u2b.out)[1] - 3.2
    d.add(elm.Line().at((xo2, pt(u2b.out)[1])).to((xo2, yrl)))
    d.add(elm.Resistor().at((xo2, yrl)).left().length(3.0).label("100 kΩ\n(protection)", loc="bottom", fontsize=8))
    d.add(elm.Line().to((0.2, yrl)))
    d.add(elm.Dot(open=True, radius=0.18).at((0, yrl)).label("RL\n(jambe droite)", loc="left", fontsize=8.5))
    d.add(elm.Label().at((xo2 + 0.5, pt(u2b.out)[1] - 0.9)).label("Driven-Right-Leg\nA = 100, pôle 48 Hz", fontsize=8.5, color=ORANGE, halign="left"))
    d.add(elm.Label().at((5.2, yrl - 1.5)).label(
        "Sans DRL (question 7.7) : débrancher la résistance de 100 kΩ et relier l'électrode RL à REF.",
        fontsize=8, color=GRIS))

    # ------------------------------------------------------------------ Alimentations
    xa0, ya0 = 0.0, 16.0
    bio = elm.Ic(size=(3.0, 4), pins=[
        elm.IcPin(name="AVCC", side="R", slot="4/4", anchorname="avcc"),
        elm.IcPin(name="REF", side="R", slot="3/4", anchorname="ref"),
        elm.IcPin(name="ADC0", side="R", slot="2/4", anchorname="adc"),
        elm.IcPin(name="AGND", side="R", slot="1/4", anchorname="agnd"),
    ], leadlen=0.5).right().label("BioBLU", fontsize=9, loc="top").color(BLEU)
    d.add(bio.anchor("avcc").at((xa0 + 2.6, ya0 + 1.2)))
    d.add(elm.Line().at(pt(bio.avcc)).right().length(1.2))
    n_av = d.here
    tag(d, d.here, "AVCC 3,3 V", "right")
    d.add(elm.Line().at(pt(bio.ref)).right().length(3.2))
    tag(d, d.here, "REF 1,65 V", "right")
    d.add(elm.Line().at(pt(bio.agnd)).right().length(0.6))
    d.add(elm.Ground())
    d.add(elm.Line().at(pt(bio.adc)).right().length(0.6))
    tag(d, d.here, "ADC0", "right", BLEU)

    lab = elm.Ic(size=(3.0, 3), pins=[
        elm.IcPin(name="+5 V", side="R", slot="3/3", anchorname="p5"),
        elm.IcPin(name="−5 V", side="R", slot="2/3", anchorname="m5"),
        elm.IcPin(name="COM", side="R", slot="1/3", anchorname="com"),
    ], leadlen=0.5).right().label("Alimentation labo", fontsize=9, loc="top").color(ORANGE)
    d.add(lab.anchor("p5").at((xa0 + 13.0, ya0 + 1.2)))
    d.add(elm.Line().at(pt(lab.p5)).right().length(3.2))
    tag(d, d.here, "+5 V", "right")
    d.add(elm.Line().at(pt(lab.com)).right().length(0.6))
    d.add(elm.Ground())
    d.add(elm.Line().at(pt(lab.m5)).right().length(3.2))
    tag(d, d.here, "−5 V", "right")

    # Notes de decouplage et REF
    d.add(elm.Label().at((xa0 + 23.0, ya0 + 1.0)).label(
        "• 100 kΩ entre REF et AGND (stabilise REF, Lab 2)\n"
        "• Découplage 100 nF + 10 µF près de chaque broche :\n"
        "   U1 : broche 7 (AVCC) vers AGND\n"
        "   U2 : broche 8 (+5 V) et broche 4 (−5 V) vers COM\n"
        "• COM de l'alimentation du labo relié à AGND du BioBLU\n"
        "• Les noeuds de même nom sont reliés\n"
        "• U2A et U2B : même boîtier TL082 (U2, DIP-8)\n"
        "• Brochage INA823 = boîtier SOIC-8 : vérifier l'adaptateur du labo",
        fontsize=8, halign="left"))

    d.save(str(SORTIE) + ".png", dpi=220)
    d.save(str(SORTIE) + ".pdf")
print("ok")
