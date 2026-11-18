# Laboratoire 4 : Bioélectronique – filtrage et interfaces

**Cours :** GEL-7072 : Bio-instrumentation et microsystèmes biomédicaux  
**Auteurs :** François Bédard & Émile Tremblay-Antoine  
**Programme :** Génie Physique (GPH)  
**Professeur :** Benoit Gosselin | **Auxiliaire :** Karim Bouzid  
**Département de génie électrique et génie informatique – Université Laval**

---

## Contenu du dépôt

Ce dépôt contient le code source en LaTeX du rapport du Laboratoire 4, les scripts de simulation numérique et les schémas/modèles de simulation LTspice :

* **`main.tex`** : Code source LaTeX du rapport (sections 5, 6 et 7 ; les questions GEL/GIF sont ignorées).
* **`main.pdf`** / **`TREMBLAY-ANTOINE_BEDARD_LAB4.pdf`** : Rendu PDF prêt pour révision et conforme aux exigences de soumission.
* **`generer_figures_simulation.py`** : Lit les résultats LTspice (`ltspice/*.raw`) et trace les figures du rapport. Avec `--no-run`, il n'appelle pas LTspice (sous macOS, lancer les `.cir` dans LTspice d'abord).
* **`analyser_mesures.py`** : Traite les captures de l'oscilloscope (`Data_et_notes/SDS0000x.bin`) : gain et phase de chaque balayage par démodulation synchrone, fréquences de coupure, figures `fig_7_*` et `Data_et_notes/resultats_mesures.json`. Nécessite les `.raw` de LTspice pour superposer les simulations.
* **`generer_schema_montage.py`** : Schéma de câblage du montage (`figures/schema_montage.pdf`).
* **`Data_et_notes/`** : Captures de l'oscilloscope du 2026-11-17 et notes de laboratoire (`2026-11-17.md`).
* **`figures/`** : Figures du rapport :
  * `fig_6_schema_bioampli.png` : Schéma fonctionnel de la chaîne.
  * `fig_6_1_offset_ina823.png` : Plage d'offset continu admissible par l'INA823 (balayage DC).
  * `fig_6_2a_bode_ina823.png` : Bode systémique de l'INA823 (modes différentiel et commun, TRMC).
  * `fig_6_2b_bode_filtres.png` : Bode de chaque filtre.
  * `fig_6_2c_bode_global.png` : Bode de la chaîne complète.
  * `fig_6_2d_drl.png` : Effet du Driven-Right-Leg sur le 60 Hz.
  * `fig_6_2e_temporel.png` : Simulation temporelle avec un ECG synthétique.
  * `fig_7_1a_bode_rfi.png`, `fig_7_1b_bode_filtres.png`, `fig_7_1c_bode_chaine.png` : Bode mesurés (question 7.1), simulations superposées.
  * `fig_7_9_ecg.png` : ECG mesuré, DRL actif (question 7.9).
  * `schema_montage.png` / `.pdf` : Schéma de câblage pour le laboratoire.
* **`ltspice/`** : Netlists LTspice et modèles TI (INA823, TL082) ; voir `ltspice/README_ltspice.txt`.

---

## Comment modifier et compiler le rapport

### Option A : Sur Overleaf (Recommandé pour travailler à deux)
1. Téléversez le fichier `main.tex`, le dossier `figures/` et `logo_ul.png` dans un projet Overleaf.
2. Définissez le compilateur sur **pdfLaTeX**.
3. Recompilez en direct et éditez à plusieurs simultanément.

### Option B : En local avec Tectonic ou LaTeX
Si vous avez Tectonic sur votre machine, la compilation se fait en une ligne :
```bash
tectonic main.tex
```
Ou avec un compilateur LaTeX standard :
```bash
pdflatex main.tex
```

---

## Synchronisation avec GitHub

Pour envoyer les modifications vers votre dépôt GitHub partagé :
```bash
git add .
git commit -m "Mise à jour du rapport"
git push origin main
```

