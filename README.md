# Laboratoire 4 : Bioélectronique – filtrage et interfaces

**Cours :** GEL-7072 : Bio-instrumentation et microsystèmes biomédicaux  
**Auteurs :** François Bédard & Émile Tremblay-Antoine  
**Programme :** Génie Physique (GPH)  
**Professeur :** Benoit Gosselin | **Auxiliaire :** Karim Bouzid  
**Département de génie électrique et génie informatique – Université Laval**

---

## Contenu du dépôt

Ce dépôt contient le code source en LaTeX du rapport du Laboratoire 4, les scripts de simulation numérique et les schémas/modèles de simulation LTspice :

* **`main.tex`** : Code source LaTeX du rapport (Étapes 5 et 6 complétées selon le gabarit du Laboratoire 3).
* **`main.pdf`** / **`TREMBLAY-ANTOINE_BEDARD_LAB4.pdf`** : Rendu PDF prêt pour révision et conforme aux exigences de soumission.
* **`generer_figures_simulation.py`** : Script Python autonome générant l'ensemble des figures, calculs analytiques et diagrammes de Bode.
* **`figures/`** : Répertoire contenant toutes les figures vectorielles/haute résolution insérées dans le rapport :
  * `fig_6_schema_bioampli.png` : Schéma fonctionnel complet de la chaîne de bioamplification.
  * `fig_6_1_bode_ina823.png` : Diagramme de Bode systémique de l'INA823 (modes différentiel et commun, TRMC).
  * `fig_6_2_bode_filtres.png` : Diagrammes de Bode individuels des étages de filtrage.
  * `fig_6_3_bode_global.png` : Réponse en fréquence globale (magnitude et phase).
  * `fig_6_4_drl_attenuation.png` : Simulation de l'atténuation du mode commun 60 Hz par le DRL.
  * `fig_6_5_simulation_temporelle.png` : Analyse temporelle transitoire (Transient) sur tracé ECG réaliste.
* **`ltspice/`** : Modèles de simulation SPICE conformes aux exigences du protocole :
  * `ina823_bode.cir` : Analyse AC systémique de l'INA823.
  * `bioampli_complet.cir` : Netlist SPICE de la chaîne complète (AC + Transitoire).
  * `README_ltspice.txt` : Guide d'exécution dans LTspice.

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

