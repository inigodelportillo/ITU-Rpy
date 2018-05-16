# -*- coding: utf-8 -*-
"""
This example reproduces the graphs plotted in ITU-R P.676 recommendation.
"""
import itur
import itur.models.itu676 as itu676
import itur.models.itu835 as itu835
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
#                     Comparison of wet and dry atmospheres                   #
###############################################################################

# Define atmospheric parameters
rho_wet = 7.5 * itur.u.g / itur.u.m**3
rho_dry = 0 * itur.u.g / itur.u.m**3
P = 1013.25 * itur.u.hPa
T = 15 * itur.u.deg_C

# Define frequency logspace parameters
N_freq = 1000
fs = np.linspace(0, 1000, N_freq)

# Initialize result vectors
atts_wet = []
atts_dry = []

# Loop over frequencies and compute values
for i, f in enumerate(fs):

    att_wet = itu676.gamma_exact(f, P, rho_wet, T)
    att_dry = itu676.gamma_exact(f, P, rho_dry, T)

    atts_wet.append(att_wet.value)
    atts_dry.append(att_dry.value)
    print('Computing frequency {0:3.0f} GHz  |  {1:3d}/{2}'
          '  |  [Wet] {3:3.2f} dB [Dry] {4:3.2f} dB'
          .format(f, i, N_freq, att_wet.value, att_dry.value))

# Plot the results
plt.figure()
plt.plot(fs, atts_wet, 'b--', label='Wet atmosphere')
plt.plot(fs, atts_dry, 'r', label='Dry atmosphere')
plt.xlabel('Frequency [GHz]')
plt.ylabel('Specific attenuation [dB]')
plt.yscale('log')
plt.xscale('linear')
plt.xlim(0, 1000)
plt.ylim(1e-3, 1e5)
plt.legend()
plt.grid(which='both', linestyle=':', color='gray', linewidth=0.3, alpha=0.5)
plt.grid(which='major', linestyle=':', color='black')
plt.title('FIGURE 1. - Specific attenuation due to atmospheric gases,\n '
          'calculated at 1 GHz intervals, including line centres')
plt.tight_layout()

###############################################################################
#                  Specific attenuation at different altitudes                #
###############################################################################

# Define atmospheric parameters
hs = np.array([0, 5, 10, 15, 20]) * itur.u.km

# Define frequency logspace parameters
N_freq = 2001
fs = np.linspace(50, 70, N_freq)

# Plot the results
plt.figure()

# Loop over frequencies and compute values
for h in hs:
    # Initialize result vectors
    atts = []
    rho = itu835.standard_water_vapour_density(h)
    P = itu835.standard_pressure(h)
    T = itu835.standard_temperature(h)
    for i, f in enumerate(fs):

        att = itu676.gamma_exact(f * itur.u.GHz, P, rho, T)

        atts.append(att.value)
        print('Computing frequency {0:3.0f} GHz  |  {1:3d}/{2}'
              '  |  {3:3.2f} dB ({4:.2f} g/m3, {5:.2f} hPa, {6:.2f} K)'
              .format(f, i, N_freq, att.value, rho.value, P.value, T.value))

    plt.plot(fs, atts, label='Altitude {0} km'.format(h.value))

plt.xlabel('Frequency [GHz]')
plt.ylabel('Specific attenuation [dB]')
plt.yscale('log')
plt.xscale('linear')
plt.xlim(50, 70)
plt.ylim(1e-3, 1e2)
plt.legend()
plt.grid(which='both', linestyle=':', color='gray', linewidth=0.3, alpha=0.5)
plt.grid(which='major', linestyle=':', color='black')
plt.title('FIGURE 2. - Specific attenuation in the range 50-70 GHz at the\n'
          ' altitudes indicated, calculated at intervals of 10 MHz\n'
          ' including line centers (0, 5, 10 15, 20) km')
plt.tight_layout()

###############################################################################
#               Comparison of line-by-line and approximate method             #
###############################################################################
# Define atmospheric parameters
el = 90
rho = 7.5 * itur.u.g / itur.u.m**3
P = 1013.25 * itur.u.hPa
T = 15 * itur.u.deg_C

# Define frequency logspace parameters
N_freq = 350
fs = np.linspace(0, 350, N_freq)

# Initialize result vectors
atts_approx = []
atts_exact = []

# Loop over frequencies and compute values
for i, f in enumerate(fs):

    att_approx = itu676.gaseous_attenuation_slant_path(
        f, el, rho, P, T, mode='approx')

    att_exact = itu676.gaseous_attenuation_slant_path(
        f, el, rho, P, T, mode='exact')

    atts_approx.append(att_approx.value)
    atts_exact.append(att_exact.value)
    print('Computing frequency {0:3.0f} GHz  |  {1:3d}/{2}'
          '  |  [Approx] {3:3.2f} dB [Exact] {4:3.2f} dB'
          .format(f, i, N_freq, att_approx.value, att_exact.value))

# Plot the results
plt.figure()
plt.plot(fs, atts_approx, 'b--', label='Approximate method Annex 2')
plt.plot(fs, atts_exact, 'r', label='Exact line-by-line method')
plt.xlabel('Frequency [GHz]')
plt.ylabel('Attenuation [dB]')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.grid(which='both', linestyle=':', color='gray', linewidth=0.3, alpha=0.5)
plt.grid(which='major', linestyle=':', color='black')
plt.title('Comparison of line-by-line method to approximate method')
plt.tight_layout()
