========================================================================
GEL-7072 : Bio-instrumentation et microsystemes biomedicaux
Laboratoire 4 : Bioelectronique - filtrage et interfaces
Simulations LTspice (questions 6.1 et 6.2)
========================================================================

Modeles (dossier modeles/, telecharges de ti.com, non modifies)
  INA823.LIB : INA823 PSpice Model rev. B (SBOMBU8B)
  TL082.301  : TL082 PSpice Model (SLOJ070), concu pour +/-15 V

Fichiers
  bioampli.sub          Sous-circuits partages : PATIENT (corps, courant de
                        deplacement, 3 electrodes) et BIOAMPLI (RFI, INA823,
                        passe-haut, passe-bas et gain TL082, sortie ADC, DRL)
  ina823_offset_dc.cir  Balayage DC de l'offset entre electrodes (plage lineaire)
  ina823_bode.cir       Bode de l'INA823, modes differentiel et commun, TRMC
  filtres_bode.cir      Bode de chaque filtre pris isolement
  bioampli_ac.cir       Bode de la chaine complete (TL082 en +/-5 V et sur 3,3 V)
  drl_ac.cir            Effet du DRL sur le 60 Hz (avec / sans DRL)
  bioampli_tran.cir     Simulation temporelle avec un ECG synthetique
  ecg_pwl.txt           ECG synthetique (genere par generer_figures_simulation.py)

Lancer les simulations
  macOS : LTspice -> File -> Open, choisir le type "Netlists (*.cir)",
          ouvrir chaque .cir puis Run. (Sous macOS, LTspice tourne dans Wine :
          le dossier personnel est le lecteur Y:.) Le .raw est ecrit a cote
          du .cir.
  Windows : python generer_figures_simulation.py lance LTspice en mode batch.

Tracer les figures
  python generer_figures_simulation.py --no-run
  (lit les .raw, ecrit figures/*.png et ltspice/resultats_simulation.json)

Expressions pour tracer dans LTspice (approche systemique, Add Trace)
  ina823_bode   Ad   : V(OUT_D)/(V(POS_D)-V(NEG_D))
                Acm  : V(OUT_C)/((V(POS_C)+V(NEG_C))/2)
                TRMC : 20*log10(abs((V(OUT_D)/(V(POS_D)-V(NEG_D)))/(V(OUT_C)/((V(POS_C)+V(NEG_C))/2))))
  bioampli_ac   V(OUTA)/(V(INPA)-V(INNA))
  filtres_bode  (V(P1A)-V(N1A))/(V(POSA)-V(NEGA)), V(HOUT)/V(HIN), V(LOUT5)/V(LIN5), V(OOUT)/V(OIN)
========================================================================
