import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

# Set plotting style matching the Lab 3 LaTeX document
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'lines.linewidth': 1.6,
    'grid.alpha': 0.5,
    'grid.linestyle': '--'
})

os.makedirs('figures', exist_ok=True)

# -----------------------------------------------------------------------------
# 1. PARAMETERS & TRANSFER FUNCTIONS
# -----------------------------------------------------------------------------
# RFI Filter
R_rfi = 1e3        # 1 kOhm
C_diff = 10e-9     # 10 nF
C_cm = 1e-9        # 1 nF
tau_rfi_diff = 2 * R_rfi * (C_diff + C_cm / 2.0)
f_c_rfi_diff = 1.0 / (2.0 * np.pi * tau_rfi_diff)  # ~7.58 kHz
f_c_rfi_cm = 1.0 / (2.0 * np.pi * R_rfi * C_cm)     # ~159 kHz

# INA823 Pre-amplifier
# G1 = 1 + 100k / RG = 10 -> RG = 100k / 9 = 11.111 kOhm
G1 = 10.0
G1_dB = 20 * np.log10(G1)
f_c_ina = 380e3 # Bandwidth at G=10 is 380 kHz for INA823
cmrr_dc_dB = 115.0 # Typical CMRR at DC for G=10
cmrr_60hz_dB = 105.0

# High-Pass Filter (1st order)
# f_L = 0.482 Hz approx 0.5 Hz
C_H = 1.0e-6       # 1 uF
R_H = 330e3        # 330 kOhm
f_L = 1.0 / (2.0 * np.pi * R_H * C_H) # 0.482 Hz
w_L = 2.0 * np.pi * f_L

# Low-Pass Filter & Gain Stage (TL082)
# G2 = 1 + Rf / R1 = 1 + 100k / 1k = 101 V/V ~ 40.09 dB
# Cf = 10 nF in parallel with Rf -> f_H = 1 / (2*pi*Rf*Cf) = 159.15 Hz ~ 150 Hz
R1_gain = 1e3      # 1 kOhm
Rf_gain = 100e3    # 100 kOhm
Cf_gain = 10e-9    # 10 nF
G2 = 1.0 + Rf_gain / R1_gain # 101 V/V
G2_dB = 20 * np.log10(G2)
f_H = 1.0 / (2.0 * np.pi * Rf_gain * Cf_gain) # 159.15 Hz
w_H = 2.0 * np.pi * f_H

# Total Bioamplifier Gain & Bandwidth
G_tot = G1 * G2
G_tot_dB = 20 * np.log10(G_tot) # ~60.09 dB
BW = f_H - f_L

print(f"Bioamplifier Design Summary:")
print(f"RFI diff cutoff: {f_c_rfi_diff:.2f} Hz, cm cutoff: {f_c_rfi_cm:.2f} Hz")
print(f"INA823 Gain G1: {G1:.2f} ({G1_dB:.2f} dB), BW: {f_c_ina/1e3:.1f} kHz")
print(f"High-pass cutoff f_L: {f_L:.3f} Hz")
print(f"Low-pass cutoff f_H: {f_H:.2f} Hz")
print(f"Gain second stage G2: {G2:.2f} ({G2_dB:.2f} dB)")
print(f"Total Gain G_tot: {G_tot:.2f} ({G_tot_dB:.2f} dB)")
print(f"Bandwidth: {BW:.2f} Hz")

# Frequency vector for Bode plots
freqs = np.logspace(-2, 6, 2000) # 0.01 Hz to 1 MHz
w = 2.0 * np.pi * freqs
s = 1j * w

# -----------------------------------------------------------------------------
# FIGURE 6.1: DIAGRAMME DE BODE INA823 (APPROCHE SYSTÉMIQUE)
# -----------------------------------------------------------------------------
# Differential transfer function: Ad(s) = G1 / (1 + s / (2*pi*f_c_ina))
Ad = G1 / (1.0 + s / (2.0 * np.pi * f_c_ina))

# Common Mode Rejection Ratio CMRR(s) modeling:
# CMRR is 115 dB at DC, single pole roll-off around 1 kHz
w_cmrr_pole = 2.0 * np.pi * 1000.0
cmrr_linear_dc = 10**(cmrr_dc_dB / 20.0)
cmrr_s = cmrr_linear_dc / np.sqrt(1.0 + (w / w_cmrr_pole)**2)
cmrr_s_dB = 20.0 * np.log10(cmrr_s)

# Common Mode Gain: Acm(s) = Ad(s) / cmrr(s)
Acm = np.abs(Ad) / cmrr_s
Acm_dB = 20.0 * np.log10(Acm)
Ad_dB = 20.0 * np.log10(np.abs(Ad))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)

ax1.semilogx(freqs, Ad_dB, 'b-', label=r'Mode différentiel $A_d = \frac{V(OUT)}{V(POS)-V(NEG)}$ ($G_1 = 10$)')
ax1.semilogx(freqs, Acm_dB, 'r--', label=r'Mode commun $A_{cm} = \frac{V(OUT)}{(V(POS)+V(NEG))/2}$')
ax1.axvline(60, color='gray', linestyle=':', alpha=0.8, label=r'Fréquence secteur (60 Hz)')
ax1.scatter([60], [20*np.log10(np.interp(60, freqs, Acm))], color='red', zorder=5)
ax1.annotate(f'Acm(60 Hz) = {20*np.log10(np.interp(60, freqs, Acm)):.1f} dB', 
             xy=(60, 20*np.log10(np.interp(60, freqs, Acm))), xytext=(120, -75),
             arrowprops=dict(arrowstyle='->', color='red'))

ax1.set_ylabel('Gain (dB)')
ax1.set_ylim(-110, 30)
ax1.grid(True, which='both')
ax1.legend(loc='upper right')
ax1.set_title("Diagramme de Bode de l'étage préamplificateur INA823 (approche systémique)")

ax2.semilogx(freqs, cmrr_s_dB, 'g-', label=r'TRMC (CMRR) $= 20\,\log_{10}(|A_d / A_{cm}|)$')
ax2.axvline(60, color='gray', linestyle=':', alpha=0.8)
cmrr_at_60 = np.interp(60, freqs, cmrr_s_dB)
ax2.scatter([60], [cmrr_at_60], color='green', zorder=5)
ax2.annotate(f'TRMC(60 Hz) = {cmrr_at_60:.1f} dB\n(Fiche technique : $\\geq 104$ dB)', 
             xy=(60, cmrr_at_60), xytext=(150, 95),
             arrowprops=dict(arrowstyle='->', color='green'))

ax2.set_xlabel('Fréquence (Hz)')
ax2.set_ylabel('TRMC (dB)')
ax2.set_ylim(40, 125)
ax2.grid(True, which='both')
ax2.legend(loc='lower left')

plt.tight_layout()
fig.savefig('figures/fig_6_1_bode_ina823.png', dpi=300)
fig.savefig('fig_6_1_bode_ina823.png', dpi=300)
plt.close(fig)
print("Figure 6.1 saved.")

# -----------------------------------------------------------------------------
# FIGURE 6.2: DIAGRAMME DE BODE DES ÉTAGES DE FILTRAGE
# -----------------------------------------------------------------------------
# High-pass filter
H_HP = (s / w_L) / (1.0 + s / w_L)
H_HP_dB = 20.0 * np.log10(np.abs(H_HP))

# Low-pass filter (active TL082 non-inverting with Cf)
# H_LP_gain(s) = (1 + Rf/R1) * (1 + s*tau_z) / (1 + s*tau_p)
tau_p = Rf_gain * Cf_gain
tau_z = (R1_gain * Rf_gain * Cf_gain) / (R1_gain + Rf_gain)
# To provide true low-pass attenuation beyond zero, cascaded with op-amp GBW (TL082 GBW ~ 3 MHz)
w_opamp = 2.0 * np.pi * 3e6
H_opamp = 1.0 / (1.0 + (G2 * s) / w_opamp)
H_LP_gain = ((1.0 + s * tau_z) / (1.0 + s * tau_p)) * G2 * H_opamp
H_LP_dB = 20.0 * np.log10(np.abs(H_LP_gain))

# Normalized Low-pass filter alone (0 dB in passband)
H_LP_norm_dB = H_LP_dB - G2_dB

# RFI differential filter
H_RFI = 1.0 / (1.0 + s * tau_rfi_diff)
H_RFI_dB = 20.0 * np.log10(np.abs(H_RFI))

fig, ax = plt.subplots(figsize=(7.2, 4.2))
ax.semilogx(freqs, H_HP_dB, 'b-', label=f'Passe-haut 1er ordre ($R_H=330\\,\\mathrm{{k\\Omega}}, C_H=1\\,\\mu\\mathrm{{F}}$, $f_L={f_L:.2f}\\,\\mathrm{{Hz}}$)')
ax.semilogx(freqs, H_LP_norm_dB, 'r-', label=f'Passe-bas 1er ordre ($R_f=100\\,\\mathrm{{k\\Omega}}, C_f=10\\,\\mathrm{{nF}}$, $f_H={f_H:.1f}\\,\\mathrm{{Hz}}$)')
ax.semilogx(freqs, H_RFI_dB, 'm:', label=f'Filtre RFI diff. ($R=1\\,\\mathrm{{k\\Omega}}, C=10\\,\\mathrm{{nF}}$, $f_c={f_c_rfi_diff/1e3:.2f}\\,\\mathrm{{kHz}}$)')

# Cutoff markers
ax.axhline(-3, color='k', linestyle='--', alpha=0.6, label='Atténuation de -3 dB')
ax.scatter([f_L, f_H], [-3, -3], color='black', zorder=5)
ax.annotate(f'$f_L = {f_L:.2f}\\,\\mathrm{{Hz}}$\n(-3 dB)', xy=(f_L, -3), xytext=(0.04, -18),
            arrowprops=dict(arrowstyle='->', color='blue'))
ax.annotate(f'$f_H = {f_H:.1f}\\,\\mathrm{{Hz}}$\n(-3 dB)', xy=(f_H, -3), xytext=(400, -18),
            arrowprops=dict(arrowstyle='->', color='red'))

# Slopes annotations
ax.text(0.02, -25, '+20 dB/décade', color='blue', fontsize=9, rotation=28)
ax.text(1200, -25, '-20 dB/décade', color='red', fontsize=9, rotation=-28)

ax.set_xlim(0.01, 100e3)
ax.set_ylim(-40, 5)
ax.set_xlabel('Fréquence (Hz)')
ax.set_ylabel('Gain normalisé (dB)')
ax.set_title('Diagrammes de Bode des étages de filtrage individuels')
ax.grid(True, which='both')
ax.legend(loc='lower center', framealpha=0.9)

plt.tight_layout()
fig.savefig('figures/fig_6_2_bode_filtres.png', dpi=300)
fig.savefig('fig_6_2_bode_filtres.png', dpi=300)
plt.close(fig)
print("Figure 6.2 saved.")

# -----------------------------------------------------------------------------
# FIGURE 6.3: DIAGRAMME DE BODE GLOBAL DU BIOAMPLIFICATEUR
# -----------------------------------------------------------------------------
# Total Bioamplifier Transfer Function
H_tot = H_RFI * Ad * H_HP * H_LP_gain
H_tot_mag = np.abs(H_tot)
H_tot_dB = 20.0 * np.log10(H_tot_mag)
H_tot_phase = np.angle(H_tot, deg=True)

# Find exact -3dB points numerically from total response
max_gain_dB = np.max(H_tot_dB)
threshold_3dB = max_gain_dB - 3.0
passband_indices = np.where(H_tot_dB >= threshold_3dB)[0]
f_low_sim = freqs[passband_indices[0]]
f_high_sim = freqs[passband_indices[-1]]

fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(7.2, 5.2), sharex=True)

# Magnitude Plot
ax_mag.semilogx(freqs, H_tot_dB, 'b-', label='Gain total $A_v(f)$')
ax_mag.axhline(max_gain_dB, color='darkblue', linestyle=':', alpha=0.7, label=f'Gain max = {max_gain_dB:.2f} dB ({G_tot:.0f} V/V)')
ax_mag.axhline(threshold_3dB, color='gray', linestyle='--', alpha=0.7, label=f'Niveau -3 dB = {threshold_3dB:.2f} dB')

ax_mag.scatter([f_low_sim, f_high_sim], [threshold_3dB, threshold_3dB], color='red', zorder=5)
ax_mag.annotate(f'$f_L = {f_low_sim:.2f}\\,\\mathrm{{Hz}}$', xy=(f_low_sim, threshold_3dB), xytext=(0.04, 45),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.2))
ax_mag.annotate(f'$f_H = {f_high_sim:.1f}\\,\\mathrm{{Hz}}$', xy=(f_high_sim, threshold_3dB), xytext=(400, 45),
                arrowprops=dict(arrowstyle='->', color='red', lw=1.2))

# Highlight ECG diagnostic passband
ax_mag.axvspan(f_low_sim, f_high_sim, color='green', alpha=0.1, label=f'Bande passante ECG ({f_low_sim:.2f} - {f_high_sim:.0f} Hz)')

ax_mag.text(0.02, 20, '+20 dB/décade\n(Rejet dérive DC)', color='blue', fontsize=8)
ax_mag.text(2500, 20, '-20 dB/décade\n(Rejet bruit HF/EMG)', color='blue', fontsize=8)

ax_mag.set_ylabel('Gain en tension (dB)')
ax_mag.set_ylim(-20, 70)
ax_mag.grid(True, which='both')
ax_mag.legend(loc='lower left', fontsize=8.5)
ax_mag.set_title("Diagramme de Bode complet du bioamplificateur ECG ($A_{v,\\mathrm{total}}$)")

# Phase Plot
ax_phase.semilogx(freqs, H_tot_phase, 'purple', label='Déphasage (degrés)')
ax_phase.axhline(0, color='gray', linestyle=':')
ax_phase.axvline(f_low_sim, color='red', linestyle=':', alpha=0.5)
ax_phase.axvline(f_high_sim, color='red', linestyle=':', alpha=0.5)

phase_fl = np.interp(f_low_sim, freqs, H_tot_phase)
phase_fh = np.interp(f_high_sim, freqs, H_tot_phase)
ax_phase.scatter([f_low_sim, f_high_sim], [phase_fl, phase_fh], color='red', zorder=5)
ax_phase.annotate(f'{phase_fl:.1f}°', xy=(f_low_sim, phase_fl), xytext=(0.05, 50),
                  arrowprops=dict(arrowstyle='->', color='purple'))
ax_phase.annotate(f'{phase_fh:.1f}°', xy=(f_high_sim, phase_fh), xytext=(350, -50),
                  arrowprops=dict(arrowstyle='->', color='purple'))

ax_phase.set_xlabel('Fréquence (Hz)')
ax_phase.set_ylabel('Phase (degrés)')
ax_phase.set_ylim(-180, 180)
ax_phase.set_yticks([-180, -90, -45, 0, 45, 90, 180])
ax_phase.grid(True, which='both')
ax_phase.legend(loc='lower left', fontsize=8.5)

plt.tight_layout()
fig.savefig('figures/fig_6_3_bode_global.png', dpi=300)
fig.savefig('fig_6_3_bode_global.png', dpi=300)
plt.close(fig)
print(f"Figure 6.3 saved: f_L = {f_low_sim:.3f} Hz, f_H = {f_high_sim:.2f} Hz, Gain = {max_gain_dB:.2f} dB.")

# -----------------------------------------------------------------------------
# FIGURE 6.4: ATTÉNUATION DU MODE COMMUN AVEC ET SANS DRL
# -----------------------------------------------------------------------------
f_drl = np.logspace(0, 4, 500)
A_drl_dc = 20.0
f_c_drl = 338.0 # compensation capacitor 1 nF with 470k
A_drl_f = A_drl_dc / np.sqrt(1.0 + (f_drl / f_c_drl)**2)

# Vcm without DRL
vcm_raw = 10e-3 # 10 mV at 60 Hz
vcm_without = vcm_raw * np.ones_like(f_drl)
vcm_with = vcm_raw / (1.0 + A_drl_f)

attenuation_dB = 20 * np.log10(vcm_without / vcm_with)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 4.8), sharex=True)

ax1.semilogx(f_drl, vcm_without * 1e3, 'r--', label='Mode commun corporel sans DRL ($V_{cm}$ brut)')
ax1.semilogx(f_drl, vcm_with * 1e3, 'g-', label='Mode commun corporel avec DRL')
ax1.axvline(60, color='gray', linestyle=':')
vcm_60_no = np.interp(60, f_drl, vcm_without * 1e3)
vcm_60_with = np.interp(60, f_drl, vcm_with * 1e3)
ax1.scatter([60, 60], [vcm_60_no, vcm_60_with], color=['red', 'green'], zorder=5)
ax1.annotate(f'{vcm_60_no:.1f} mV', xy=(60, vcm_60_no), xytext=(12, 11.5),
             arrowprops=dict(arrowstyle='->', color='red'))
ax1.annotate(f'{vcm_60_with:.2f} mV (-26,4 dB)', xy=(60, vcm_60_with), xytext=(120, 2),
             arrowprops=dict(arrowstyle='->', color='green'))

ax1.set_ylabel(r'$V_{cm}$ corporel ($\mathrm{mV}$)')
ax1.set_ylim(-0.5, 14)
ax1.grid(True, which='both')
ax1.legend(loc='center right')
ax1.set_title("Effet du circuit Driven-Right-Leg (DRL) sur la réjection du mode commun 60 Hz")

ax2.semilogx(f_drl, attenuation_dB, 'b-', label=r'Facteur d’atténuation du DRL : $20\,\log_{10}(1 + A_{\mathrm{DRL}})$')
ax2.axvline(60, color='gray', linestyle=':')
att_60 = np.interp(60, f_drl, attenuation_dB)
ax2.scatter([60], [att_60], color='blue', zorder=5)
ax2.annotate(f'{att_60:.1f} dB à 60 Hz\n(atténuation d’un facteur {10**(att_60/20):.1f})',
             xy=(60, att_60), xytext=(120, 18),
             arrowprops=dict(arrowstyle='->', color='blue'))

ax2.set_xlabel('Fréquence (Hz)')
ax2.set_ylabel('Atténuation (dB)')
ax2.set_ylim(0, 32)
ax2.grid(True, which='both')
ax2.legend(loc='lower left')

plt.tight_layout()
fig.savefig('figures/fig_6_4_drl_attenuation.png', dpi=300)
fig.savefig('fig_6_4_drl_attenuation.png', dpi=300)
plt.close(fig)
print("Figure 6.4 saved.")

# -----------------------------------------------------------------------------
# FIGURE 6.5: SIMULATION TEMPORELLE (TRANSIENT ANALYSIS) DU SIGNAL ECG
# -----------------------------------------------------------------------------
# Generate synthetic realistic ECG signal over 5 seconds, plotted after settling (t = 1.0 to 4.0 s)
np.random.seed(42)
fs = 2000.0 # 2 kHz sampling rate
t_full = np.arange(0, 5.0, 1.0/fs)

def ecg_wave(t_vec, hr=68):
    period = 60.0 / hr
    ecg = np.zeros_like(t_vec)
    for t0 in np.arange(0.2, t_vec[-1] + period, period):
        # P-wave: Gaussian peak around t0 - 0.16
        ecg += 0.15 * np.exp(-((t_vec - (t0 - 0.16)) / 0.025)**2)
        # Q-wave: negative dip around t0 - 0.04
        ecg -= 0.15 * np.exp(-((t_vec - (t0 - 0.04)) / 0.01)**2)
        # R-wave: sharp positive peak at t0
        ecg += 1.2 * np.exp(-((t_vec - t0) / 0.015)**2)
        # S-wave: negative dip around t0 + 0.035
        ecg -= 0.3 * np.exp(-((t_vec - (t0 + 0.035)) / 0.012)**2)
        # T-wave: wider positive peak around t0 + 0.22
        ecg += 0.35 * np.exp(-((t_vec - (t0 + 0.22)) / 0.05)**2)
    return ecg

ecg_clean_full = ecg_wave(t_full, hr=68) # ~1.2 mV R-wave peak

# Perturbations:
# 1) DC electrode half-cell offset: 120 mV
dc_offset = 120.0 # mV
# 2) Respiration baseline drift (0.25 Hz): 0.5 mV
drift_full = 0.5 * np.sin(2 * np.pi * 0.25 * t_full)
# 3) Residual 60 Hz after DRL & TRMC (common mode 10 mV -> DRL 0.48 mV -> TRMC 105 dB gives < 0.01 mV diff, plus 0.04 mV impedance mismatch)
noise_60hz_full = 0.05 * np.sin(2 * np.pi * 60.0 * t_full)
# 4) Wideband thermal noise
noise_hf_full = 0.015 * np.random.randn(len(t_full))

# Raw input bio-signal (mV)
ecg_raw_full = ecg_clean_full + dc_offset + drift_full + noise_60hz_full + noise_hf_full

# Step 1: INA823 output (Gain 10, offset amplifié, centré sur Vref = 1.65 V)
v_ina_full = 1.65 + G1 * (ecg_raw_full / 1000.0)

# Step 2: High-pass filtering (1st order Butterworth, f_L = 0.48 Hz)
sos_hp = signal.butter(1, f_L, 'highpass', fs=fs, output='sos')
ecg_hp_full = signal.sosfilt(sos_hp, ecg_raw_full - dc_offset)

# Step 3: Low-pass filtering (1st order Butterworth, f_H = 159 Hz) & Gain G2 = 101
sos_lp = signal.butter(1, f_H, 'lowpass', fs=fs, output='sos')
ecg_lp_full = signal.sosfilt(sos_lp, ecg_hp_full)

# Final output: amplified by G2 = 101, centered at Vref = 1.65 V
v_final_full = 1.65 + (G_tot * ecg_lp_full / 1000.0)

# Extract settled window for visualization: t from 1.0 to 4.0 s (re-zero time to 0 .. 3 s)
mask = (t_full >= 1.0) & (t_full <= 4.0)
t = t_full[mask] - 1.0
ecg_raw_input = ecg_raw_full[mask]
v_ina_out = v_ina_full[mask]
v_final_out = v_final_full[mask]

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7.2, 5.8), sharex=True)

# Panneau 1: Signal d'entrée brut
ax1.plot(t, ecg_raw_input, 'k-', lw=1.0, label=r'Biosignal d’entrée brut $V_{in}$ (mV)')
ax1.axhline(dc_offset, color='r', linestyle='--', alpha=0.7, label=f'Décalage DC demi-cellule ({dc_offset:.0f} mV)')
ax1.set_ylabel(r'$V_{in}$ ($\mathrm{mV}$)')
ax1.set_ylim(dc_offset - 3, dc_offset + 5)
ax1.set_title('Simulation temporelle (Transient Analysis) du bioamplificateur conçu')
ax1.grid(True)
ax1.legend(loc='upper right', fontsize=8)

# Panneau 2: Sortie préamplificateur INA823
ax2.plot(t, v_ina_out, 'darkblue', lw=1.1, label=r'Sortie INA823 $V_{\mathrm{INA}}$ (V), $G_1 = 10$ ($V_{ref}=1,65\,\mathrm{V}$)')
ax2.axhline(1.65 + G1*(dc_offset/1e3), color='r', linestyle=':', label=f'Offset amplifié ({1.65 + G1*(dc_offset/1e3):.2f} V, sous saturation)')
ax2.set_ylabel(r'$V_{\mathrm{INA}}$ ($\mathrm{V}$)')
ax2.set_ylim(2.7, 3.1)
ax2.grid(True)
ax2.legend(loc='upper right', fontsize=8)

# Panneau 3: Sortie finale filtrée et amplifiée
ax3.plot(t, v_final_out, 'g-', lw=1.3, label=r'Sortie finale $V_{out}$ (V), $A_v = 1010$ ($V_{ref}=1,65\,\mathrm{V}$)')
ax3.axhline(1.65, color='gray', linestyle='--', label=r'Masse virtuelle $V_{ref} = 1,65\,\mathrm{V}$')
ax3.axhline(3.3, color='r', linestyle=':', label=r'Pleine échelle ADC ($V_{DD} = 3,3\,\mathrm{V}$)')
ax3.axhline(0.0, color='r', linestyle=':', label=r'GND ($0\,\mathrm{V}$)')

# Annotate R wave peak
r_idx = np.argmax(v_final_out[int(0.5*fs):int(1.5*fs)]) + int(0.5*fs)
ax3.annotate(f'Onde R amplifiée\n$V_{{peak}} = {v_final_out[r_idx]:.2f}\\,\\mathrm{{V}}$\n($\\Delta V = {v_final_out[r_idx]-1.65:.2f}\\,\\mathrm{{V}}$)',
             xy=(t[r_idx], v_final_out[r_idx]), xytext=(t[r_idx]+0.16, 2.05),
             arrowprops=dict(arrowstyle='->', color='green', lw=1.2))

ax3.set_xlabel('Temps stabilisé (s)')
ax3.set_ylabel(r'$V_{out}$ ($\mathrm{V}$)')
ax3.set_ylim(-0.2, 3.6)
ax3.grid(True)
ax3.legend(loc='lower right', fontsize=8)

plt.tight_layout()
fig.savefig('figures/fig_6_5_simulation_temporelle.png', dpi=300)
fig.savefig('fig_6_5_simulation_temporelle.png', dpi=300)
plt.close(fig)
print("Figure 6.5 saved.")

# -----------------------------------------------------------------------------
# FIGURE 6.0: SCHÉMA BLOCS / ARCHITECTURE DU BIOAMPLIFICATEUR (CORRIGÉ & AÉRÉ)
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.0, 3.8))
ax.axis('off')

from matplotlib.patches import FancyBboxPatch

def draw_block(ax, x, y, w, h, title, subtitle, color='#EBF3FB', edgecolor='#2B5B84'):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.04",
                       fc=color, ec=edgecolor, lw=1.5)
    ax.add_patch(p)
    ax.text(x + w/2.0, y + h*0.73, title, ha='center', va='center', weight='bold', fontsize=8.6, color='#1A364A')
    ax.text(x + w/2.0, y + h*0.27, subtitle, ha='center', va='center', fontsize=7.4, color='#333333')

def draw_arrow(ax, x1, y1, x2, y2, label=""):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color="#1A364A", lw=1.5))
    if label:
        ax.text((x1+x2)/2.0, (y1+y2)/2.0 + 0.05, label, ha='center', va='bottom', fontsize=8, weight='bold')

# Parameters for clear spacing
bw = 0.16
bh = 0.36
y_top = 0.53

# Blocks coordinates
x1 = 0.06
x2 = 0.30
x3 = 0.54
x4 = 0.78

draw_block(ax, x1, y_top, bw, bh, "1. Filtre RFI", "Passif EMI/RFI\n$f_{c,diff}=7,6\\,\\mathrm{kHz}$", color='#FEF3D6', edgecolor='#B8860B')
draw_block(ax, x2, y_top, bw, bh, "2. Préampli.\nINA823", "$G_1 = 10\\,\\mathrm{V/V}$ (20 dB)\nTRMC $\\geq 104\\,\\mathrm{dB}$", color='#DCEEFB', edgecolor='#1F77B4')
draw_block(ax, x3, y_top, bw, bh, "3. Passe-Haut", "1er ordre $R_H C_H$\n$f_L = 0,48\\,\\mathrm{Hz}$", color='#E2F0D9', edgecolor='#385723')
draw_block(ax, x4, y_top, bw, bh, "4. Passe-Bas\n& Gain TL082", "$G_2 = 101\\,\\mathrm{V/V}$ (40 dB)\n$f_H = 159\\,\\mathrm{Hz}$", color='#FCE4D6', edgecolor='#C65911')

# DRL block at bottom
draw_block(ax, 0.30, 0.06, 0.40, 0.32, "5. Circuit Driven-Right-Leg (DRL)", 
           "Inverseur TL082 ($A_{\\mathrm{DRL}}=20$) + protection $R_o = 100\\,\\mathrm{k\\Omega}$\nRéjection du mode commun 60 Hz (-26,4 dB)", 
           color='#EDEDED', edgecolor='#595959')

# Arrows main chain
y_arrow = y_top + bh/2.0
draw_arrow(ax, -0.04, y_arrow, x1, y_arrow, "ECG")
draw_arrow(ax, x1 + bw, y_arrow, x2, y_arrow)
draw_arrow(ax, x2 + bw, y_arrow, x3, y_arrow)
draw_arrow(ax, x3 + bw, y_arrow, x4, y_arrow)
draw_arrow(ax, x4 + bw, y_arrow, x4 + bw + 0.07, y_arrow)
ax.text(x4 + bw + 0.08, y_arrow, "Vers ADC\n(0 - 3,3 V)", va='center', ha='left', fontsize=8, weight='bold', color='#1A364A')

# DRL connections
ax.annotate("", xy=(x2 + 0.08, 0.38), xytext=(x2 + 0.08, y_top),
            arrowprops=dict(arrowstyle="->", color="#595959", lw=1.2, ls="--"))
ax.text(x2 + 0.10, 0.45, "Mode commun\n$V_{cm}$ ($R_a$)", fontsize=7.5, color='#595959')

ax.annotate("", xy=(-0.04, 0.22), xytext=(0.30, 0.22),
            arrowprops=dict(arrowstyle="->", color="#595959", lw=1.2, ls="--"))
ax.text(0.13, 0.24, "Jambe droite (REF)", fontsize=7.5, color='#595959', ha='center')

# Power supply banner
ax.text(0.5, 0.97, "Alimentation unipolaire : $V_{DD} = 3,3\\,\\mathrm{V}$, $\\mathrm{GND} = 0\\,\\mathrm{V}$, polarisation $V_{ref} = 1,65\\,\\mathrm{V}$ (BioBLU)",
        ha='center', va='top', fontsize=8.5, weight='bold', color='#1A364A',
        bbox=dict(boxstyle='square,pad=0.2', fc='#F2F2F2', ec='none'))

ax.set_xlim(-0.06, 1.06)
ax.set_ylim(0, 1.0)
plt.tight_layout()
fig.savefig('figures/fig_6_schema_bioampli.png', dpi=300)
fig.savefig('fig_6_schema_bioampli.png', dpi=300)
plt.close(fig)
print("Figure 6.0 updated cleanly.")

