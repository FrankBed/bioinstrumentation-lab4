========================================================================
GEL-7072 / GEL-4072 : Bio-instrumentation et microsystèmes biomédicaux
Laboratoire 4 : Bioélectronique – filtrage et interfaces
Fichiers de simulation LTspice pour l'Étape 6
========================================================================

Ce répertoire contient les fichiers de simulation SPICE / LTspice conçus pour
valider l'Étape 6 du Laboratoire 4 :

1. `ina823_bode.cir` :
   Simulation AC selon l'« approche systémique » demandée à la question 6.2 :
   - Mode différentiel : Tracez l'expression :
     V(OUT_D)/(V(POS_D) - V(NEG_D))
     Gain mesuré : 20 dB (G = 10 V/V).
   - Mode commun : Tracez l'expression :
     V(OUT_CM)/((V(POS_CM) + V(NEG_CM))/2)
     Gain mesuré : -95 dB à DC, -85 dB à 60 Hz.
   - TRMC (CMRR en dB) : Tracez l'expression :
     20*log10(abs((V(OUT_D)/(V(POS_D)-V(NEG_D))) / (V(OUT_CM)/((V(POS_CM)+V(NEG_CM))/2))))
     Valeur mesurée : 115 dB à DC, > 104 dB à 60 Hz.

2. `bioampli_complet.cir` :
   Simulation complète de toute la chaîne analogique :
   - Filtre RFI en entrée (R = 1k, C_diff = 10nF, C_cm = 1nF, fc = 7.6 kHz)
   - Préamplificateur INA823 (G1 = 10 V/V, RG = 11.111 kOhm, REF = 1.65 V)
   - Filtre Passe-Haut du 1er ordre (CH = 1 uF, RH = 330 kOhm, fL = 0.48 Hz)
   - Filtre Passe-Bas & Étage de gain TL082 (R1 = 1k, Rf = 100k, Cf = 10nF, G2 = 101, fH = 159 Hz)
   - Circuit Driven-Right-Leg (DRL) avec AOP TL082 (Ra = 47k, Rf = 470k, Ro = 100k)
   - Alimentation unipolaire 3.3 V et masse virtuelle 1.65 V

Pour ouvrir dans LTspice :
- Lancez LTspice.
- Allez dans File -> Open, sélectionnez les fichiers `.cir`.
- Cliquez sur "Run" (icône du petit bonhomme qui court).
- Pour tracer un signal ou un ratio, faites un clic droit sur le graphique et choisissez "Add Traces", puis entrez les formules ci-dessus.
========================================================================

