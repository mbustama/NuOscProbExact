# -*- coding: utf-8 -*-
r"""The scenes of the NuOscProbExact demonstration reel.

Each scene builder returns ``(figure, update, n_frames)``, where ``update``
takes a frame index and redraws.  Every scene uses the same figure size, so
that the rendered clips can be concatenated into one reel without rescaling.

The scenes differ in an instructive way, and the difference is the reason
there are four of them rather than one.

* `scene_cp` and `scene_sterile` recompute a whole 2D map every frame, tens
  of thousands of probabilities, in a *single* broadcast call.  These are
  the cheap ones --- a few milliseconds of physics against tens of
  milliseconds of drawing.
* `scene_earth` cannot do that.  `earth.probabilities_3nu_earth` takes a
  scalar energy, because the chord geometry and the PREM slabbing change
  with it, so a curve of a hundred energies is a hundred calls.  That is
  about 70 ms a frame, still fine for an animation, and it is honest about
  where the library batches and where it does not.
* `scene_slabs` animates the one approximation the library makes --- how
  finely a continuous profile is cut --- and shows it converging.

Physics conventions worth stating, because each is easy to get wrong:

* Antineutrinos need **both** changes, conjugating the vacuum Hamiltonian
  *and* reversing the matter potential.
* Probabilities come back with the initial flavor varying slowest, so with
  ``n`` flavors ``P[n*alpha + beta]`` is P(nu_alpha -> nu_beta).
* The colour scale of an animated map is fixed over the whole sweep.  Taking
  it from the first frame lets later frames clip silently, and clipping
  looks like structure rather than like saturation.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import os
import sys

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

import earth                                                      # noqa: E402
import hamiltonians3nu                                            # noqa: E402
import hamiltonians4nu                                            # noqa: E402
import oscprob3nu                                                 # noqa: E402
import oscprob4nu                                                 # noqa: E402
import slabs                                                      # noqa: E402
from globaldefs import (CONV_KM_TO_INV_EV, D21_NO_BF, D31_NO_BF,  # noqa: E402
                        DCP_NO_BF, EARTH_RADIUS, S12_NO_BF, S13_NO_BF,
                        S23_NO_BF, VCC_EARTH_CRUST, VNC_EARTH_CRUST)

FIGSIZE = (11.0, 4.6)
DPI = 120

# Indices into the flavor-ordered probability tuples.
P3_MUE = 3          # three flavors: P(nu_mu -> nu_e)
P3_MUMU = 4         # three flavors: P(nu_mu -> nu_mu)
P4_MUS = 7          # four flavors:  P(nu_mu -> nu_s)

CAPTION = '#334155'
ACCENT = '#1d4ed8'
MARK = '#dc2626'


def _vacuum_3nu(dcp=None):
    r"""Returns the energy-independent three-flavor vacuum Hamiltonian."""
    return hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        S12_NO_BF, S23_NO_BF, S13_NO_BF,
        DCP_NO_BF if dcp is None else dcp, D21_NO_BF, D31_NO_BF)


def _vacuum_4nu(d41):
    r"""Returns the 3+1 vacuum Hamiltonian for a mass-squared splitting."""
    return hamiltonians4nu.hamiltonian_4nu_vacuum_energy_independent(
        S12_NO_BF, S23_NO_BF, S13_NO_BF,
        np.sqrt(0.10), np.sqrt(0.10), 0.0,
        DCP_NO_BF, D21_NO_BF, D31_NO_BF, d41)


def _frame(fig, note):
    r"""Adds the standing caption and returns the mutable heading."""
    fig.patch.set_facecolor('white')
    fig.text(0.5, 0.028, note, ha='center', fontsize=9.5, color=CAPTION)
    heading = fig.suptitle('', fontsize=13)
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    return heading


def scene_cp(frames=240, grid=200):
    r"""The CP phase sweeping through 2*pi.

    Left, an oscillogram of P_mue in matter, recomputed each frame.  Right,
    the bi-probability ellipse --- which is the locus traced *by* the phase,
    so it is drawn once and a marker runs around it, saying where on the
    ellipse the left panel currently sits.
    """
    energy_ell, baseline_ell = 0.8e9, 1300.0
    phases = np.linspace(0.0, 2.0*np.pi, frames, endpoint=False)
    energies = np.logspace(-1.0, 1.0, grid)*1.e9
    baselines = np.linspace(50.0, 12000.0, grid)*CONV_KM_TO_INV_EV

    def oscillogram(dcp):
        stack = hamiltonians3nu.hamiltonian_3nu_matter(
            _vacuum_3nu(dcp), energies, VCC_EARTH_CRUST)
        return oscprob3nu.probabilities_3nu(
            stack[:, None, :, :], baselines[None, :])[:, :, P3_MUE]

    def ellipse_point(dcp):
        length = baseline_ell*CONV_KM_TO_INV_EV
        h_nu = hamiltonians3nu.hamiltonian_3nu_matter(
            _vacuum_3nu(dcp), energy_ell, VCC_EARTH_CRUST)
        h_bar = hamiltonians3nu.hamiltonian_3nu_matter(
            np.conj(_vacuum_3nu(dcp)), energy_ell, -VCC_EARTH_CRUST)
        return (oscprob3nu.probabilities_3nu(h_nu, length)[P3_MUE],
                oscprob3nu.probabilities_3nu(h_bar, length)[P3_MUE])

    locus = np.array([ellipse_point(d)
                      for d in np.linspace(0.0, 2.0*np.pi, 240)])
    ceiling = max(float(oscillogram(d).max())
                  for d in np.linspace(0.0, 2.0*np.pi, 12, endpoint=False))

    fig, (ax_map, ax_ell) = plt.subplots(
        1, 2, figsize=FIGSIZE, dpi=DPI,
        gridspec_kw={'width_ratios': [1.35, 1.0]})

    image = ax_map.imshow(
        oscillogram(phases[0]), origin='lower', aspect='auto', cmap='viridis',
        vmin=0.0, vmax=ceiling,
        extent=[50.0, 12000.0, -1.0, 1.0])
    ax_map.set_xlabel('Baseline [km]')
    ax_map.set_ylabel('Energy [GeV]')
    ax_map.set_yticks([-1, 0, 1])
    ax_map.set_yticklabels(['0.1', '1', '10'])
    ax_map.set_title(r'$P_{\mu e}$ in matter, %d$\,\times\,$%d'
                     % (grid, grid), fontsize=11)
    fig.colorbar(image, ax=ax_map, pad=0.02).set_label(r'$P_{\mu e}$')

    ax_ell.plot(locus[:, 0], locus[:, 1], color=ACCENT, lw=1.6)
    lo, hi = float(min(locus.min(), 0.0)), float(locus.max())*1.10
    ax_ell.plot([lo, hi], [lo, hi], ls=':', color='#64748b', lw=1.0,
                label='CP symmetric')
    marker, = ax_ell.plot([], [], 'o', ms=9, color=MARK, zorder=5)
    ax_ell.set_xlim(lo, hi)
    ax_ell.set_ylim(lo, hi)
    ax_ell.set_xlabel(r'$P(\nu_\mu \to \nu_e)$')
    ax_ell.set_ylabel(r'$P(\bar{\nu}_\mu \to \bar{\nu}_e)$')
    ax_ell.set_title(r'Bi-probability, %.1f GeV, %d km'
                     % (energy_ell/1.e9, baseline_ell), fontsize=11)
    ax_ell.legend(loc='upper left', fontsize=9, frameon=False)
    ax_ell.set_aspect('equal', adjustable='box')

    heading = _frame(fig, 'One broadcast call a frame  ---  %s probabilities'
                     % format(grid*grid, ','))

    def update(index):
        dcp = phases[index]
        image.set_data(oscillogram(dcp))
        marker.set_data(*[[v] for v in ellipse_point(dcp)])
        heading.set_text(r'The CP phase:  $\delta_{CP} = %.2f\pi$'
                         % (dcp/np.pi))
        return image, marker

    return fig, update, len(phases)


def scene_sterile(frames=180, grid=180):
    r"""A sterile state appearing as its mass-squared splitting sweeps.

    Four flavors, in matter of constant density so that the whole map is
    one call.  The sterile state feels neither potential, so the
    neutral-current term stops cancelling and sits on the sterile entry ---
    which is what places the resonance that moves across the frame.
    """
    splittings = np.logspace(-1.0, 1.0, frames)          # Dm41^2 [eV^2]
    energies = np.logspace(0.0, 2.0, grid)*1.e9
    baselines = np.linspace(100.0, 12742.0, grid)*CONV_KM_TO_INV_EV
    slice_at = grid - 1                                  # the longest chord

    def sterile_map(d41):
        stack = hamiltonians4nu.hamiltonian_4nu_matter(
            _vacuum_4nu(d41), energies, VCC_EARTH_CRUST, VNC_EARTH_CRUST)
        return oscprob4nu.probabilities_4nu(
            stack[:, None, :, :], baselines[None, :])[:, :, P4_MUS]

    ceiling = max(float(sterile_map(d).max())
                  for d in np.logspace(-1.0, 1.0, 8))

    fig, (ax_map, ax_cut) = plt.subplots(
        1, 2, figsize=FIGSIZE, dpi=DPI,
        gridspec_kw={'width_ratios': [1.35, 1.0]})

    image = ax_map.imshow(
        sterile_map(splittings[0]), origin='lower', aspect='auto',
        cmap='magma', vmin=0.0, vmax=ceiling,
        extent=[100.0, 12742.0, 0.0, 2.0])
    ax_map.set_xlabel('Baseline [km]')
    ax_map.set_ylabel('Energy [GeV]')
    ax_map.set_yticks([0, 1, 2])
    ax_map.set_yticklabels(['1', '10', '100'])
    ax_map.set_title(r'$P(\nu_\mu \to \nu_s)$, %d$\,\times\,$%d'
                     % (grid, grid), fontsize=11)
    fig.colorbar(image, ax=ax_map, pad=0.02).set_label(r'$P(\nu_\mu\to\nu_s)$')
    ax_map.axvline(12742.0, color='white', lw=1.0, ls='--', alpha=0.7)

    curve, = ax_cut.plot([], [], color=ACCENT, lw=1.8)
    ax_cut.set_xlim(0.0, 2.0)
    ax_cut.set_ylim(0.0, ceiling*1.05)
    ax_cut.set_xticks([0, 1, 2])
    ax_cut.set_xticklabels(['1', '10', '100'])
    ax_cut.set_xlabel('Energy [GeV]')
    ax_cut.set_ylabel(r'$P(\nu_\mu \to \nu_s)$')
    ax_cut.set_title('Through the diameter, 12742 km', fontsize=11)
    ax_cut.grid(alpha=0.25)

    heading = _frame(
        fig, 'Four flavors, constant density  ---  the sterile state feels '
             'neither potential')
    axis = np.log10(energies/1.e9)

    def update(index):
        d41 = splittings[index]
        grid_now = sterile_map(d41)
        image.set_data(grid_now)
        curve.set_data(axis, grid_now[:, slice_at])
        heading.set_text(r'A sterile state:  $\Delta m^2_{41} = %.2f$ eV$^2$'
                         % d41)
        return image, curve

    return fig, update, len(splittings)


def scene_earth(frames=150, n_energies=100):
    r"""A chord through the Earth, as the zenith angle swings.

    Left, the PREM density in cross-section with the neutrino's path drawn
    on it; right, the survival probability along that path.  Unlike the
    other scenes this one cannot be a single call: the chord, and the slabs
    cut from it, change with the angle, so each energy is its own call.
    """
    cosines = np.linspace(-1.0, -0.05, frames)
    energies = np.logspace(0.0, 2.0, n_energies)*1.e9
    vacuum = _vacuum_3nu()

    fig, (ax_earth, ax_prob) = plt.subplots(
        1, 2, figsize=FIGSIZE, dpi=DPI,
        gridspec_kw={'width_ratios': [1.0, 1.35]})

    # The density disk is static, so it is drawn once: concentric circles
    # shaded by PREM, outermost first so the inner ones sit on top.
    radii = np.linspace(EARTH_RADIUS, 0.0, 240)
    densities = earth.density_prem(radii)
    colours = plt.cm.YlOrBr(0.15 + 0.75*densities/densities.max())
    for radius, colour in zip(radii, colours):
        ax_earth.add_patch(Circle((0.0, 0.0), radius, facecolor=colour,
                                  edgecolor='none', zorder=1))
    ax_earth.add_patch(Circle((0.0, 0.0), EARTH_RADIUS, fill=False,
                              edgecolor='#334155', lw=1.2, zorder=3))
    span = EARTH_RADIUS*1.08
    ax_earth.set_xlim(-span, span)
    ax_earth.set_ylim(-span, span)
    ax_earth.set_aspect('equal')
    ax_earth.axis('off')
    ax_earth.set_title('PREM, and the path through it', fontsize=11)
    ax_earth.plot([0.0], [EARTH_RADIUS], 'v', color='#0f172a', ms=9,
                  zorder=5)
    path, = ax_earth.plot([], [], color=MARK, lw=2.2, zorder=4)

    curve, = ax_prob.plot([], [], color=ACCENT, lw=1.8)
    ax_prob.set_xlim(0.0, 2.0)
    ax_prob.set_ylim(0.0, 1.02)
    ax_prob.set_xticks([0, 1, 2])
    ax_prob.set_xticklabels(['1', '10', '100'])
    ax_prob.set_xlabel('Energy [GeV]')
    ax_prob.set_ylabel(r'$P(\nu_\mu \to \nu_\mu)$')
    ax_prob.set_title('Survival along the chord', fontsize=11)
    ax_prob.grid(alpha=0.25)

    heading = _frame(
        fig, 'The chord is re-cut into PREM slabs at every angle  ---  '
             '%d energies, %d calls' % (n_energies, n_energies))
    axis = np.log10(energies/1.e9)

    def update(index):
        costhz = cosines[index]
        length = earth.distance_traveled_inside_earth(costhz)
        sin_thz = np.sqrt(max(0.0, 1.0 - costhz*costhz))
        entry = (length*sin_thz, EARTH_RADIUS + length*costhz)
        path.set_data([entry[0], 0.0], [entry[1], EARTH_RADIUS])
        curve.set_data(axis, [earth.probabilities_3nu_earth(
            vacuum, float(e), costhz)[P3_MUMU] for e in energies])
        heading.set_text(r'Through the Earth:  $\cos\theta_z = %+.2f$, '
                         r'%d km' % (costhz, round(length)))
        return path, curve

    return fig, update, len(cosines)


def scene_slabs(frames=120, baseline=4000.0, energy=1.0e9):
    r"""A continuous profile cut ever more finely, converging.

    This is the one approximation the library makes.  Within a slab nothing
    is approximated; the only question is how finely a profile that really
    varies is sliced, and that is the caller's to answer.  The reference is
    the same calculation at 600 slabs.
    """
    # Duplicates are kept rather than collapsed: at small slab counts the
    # ramp repeats a value for several frames, which reads as a hold on the
    # part worth looking at, and de-duplicating cut the scene to a third of
    # its length.
    counts = np.round(np.logspace(0.0, np.log10(80.0), frames)).astype(int)
    vacuum = _vacuum_3nu()
    total = baseline*CONV_KM_TO_INV_EV

    def profile(x):
        r"""A smooth, deliberately awkward density profile [g/cm^3].

        The amplitudes sum to more than the offset, so the positivity of
        this is a property of the phase offset and the frequency ratio
        rather than of the coefficients: the two sines do not reach their
        minima together, and the measured minimum is 0.28 g/cm^3.  A
        negative density would quietly reverse the sign of the matter
        potential, which is the antineutrino case and not what this scene
        is about --- so re-check the minimum before changing any of the
        five numbers below.
        """
        return 3.0 + 2.2*np.sin(2.0*np.pi*x/baseline) \
            + 1.1*np.sin(6.0*np.pi*x/baseline + 0.7)

    def probability(n):
        edges = np.linspace(0.0, baseline, n + 1)
        mid = 0.5*(edges[:-1] + edges[1:])
        h = hamiltonians3nu.hamiltonian_3nu_matter(
            vacuum, energy, earth.matter_potential(profile(mid)))
        return slabs.probabilities_3nu_slabs(
            h, np.full(n, total/n))[P3_MUE]

    reference = probability(600)
    values = [probability(int(n)) for n in counts]

    fig, (ax_rho, ax_conv) = plt.subplots(
        1, 2, figsize=FIGSIZE, dpi=DPI,
        gridspec_kw={'width_ratios': [1.15, 1.0]})

    fine = np.linspace(0.0, baseline, 600)
    ax_rho.plot(fine, profile(fine), color='#94a3b8', lw=1.6,
                label='the true profile')
    stair, = ax_rho.step([], [], where='post', color=ACCENT, lw=1.8,
                         label='what is solved')
    ax_rho.set_xlabel('Distance along the trajectory [km]')
    ax_rho.set_ylabel(r'Density [g cm$^{-3}$]')
    ax_rho.set_title('Cutting a profile into slabs', fontsize=11)
    ax_rho.legend(loc='upper right', fontsize=9, frameon=False)
    ax_rho.grid(alpha=0.25)

    ax_conv.axhline(reference, color='#94a3b8', ls='--', lw=1.2,
                    label='600 slabs')
    trace, = ax_conv.plot([], [], color=ACCENT, lw=1.6)
    dot, = ax_conv.plot([], [], 'o', ms=7, color=MARK, zorder=5)
    ax_conv.set_xscale('log')
    ax_conv.set_xlim(1, 90)
    lo, hi = min(min(values), reference), max(max(values), reference)
    pad = 0.12*(hi - lo) if hi > lo else 0.01
    ax_conv.set_ylim(lo - pad, hi + pad)
    ax_conv.set_xlabel('Number of slabs')
    ax_conv.set_ylabel(r'$P_{\mu e}$')
    ax_conv.set_title('Converging on the answer', fontsize=11)
    ax_conv.legend(loc='lower right', fontsize=9, frameon=False)
    ax_conv.grid(alpha=0.25)

    heading = _frame(
        fig, 'Each slab is solved exactly  ---  the only approximation is '
             'how finely the profile is cut')

    def update(index):
        n = int(counts[index])
        edges = np.linspace(0.0, baseline, n + 1)
        mid = 0.5*(edges[:-1] + edges[1:])
        stair.set_data(np.append(edges[:-1], baseline),
                       np.append(profile(mid), profile(mid[-1])))
        trace.set_data(counts[:index + 1], values[:index + 1])
        dot.set_data([n], [values[index]])
        heading.set_text(r'Layered matter:  %d slab%s,  $P_{\mu e} = %.5f$'
                         % (n, '' if n == 1 else 's', values[index]))
        return stair, trace, dot

    return fig, update, len(counts)


SCENES = {
    'cp': scene_cp,
    'sterile': scene_sterile,
    'earth': scene_earth,
    'slabs': scene_slabs,
}
