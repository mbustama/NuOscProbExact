# -*- coding: utf-8 -*-
r"""Renders a short animation of the CP phase sweeping through 2*pi.

The point of the animation is the point of the library: every frame is a
*single* broadcast call.  A 200x200 oscillogram is 40 000 probabilities,
and computing one costs about 3 ms against roughly 25 ms to draw it --- so
the physics is around a tenth of the work and the animation is essentially
free.  Anything that had to loop in Python could not be animated this way.

Two panels, sharing one value of :math:`\delta_{CP}` per frame:

* **left**, an oscillogram of :math:`P_{\mu e}` in matter of constant
  density, over energy and baseline, recomputed for each phase;
* **right**, the bi-probability ellipse at a fixed energy and baseline.
  The ellipse itself is the locus traced *by* :math:`\delta_{CP}`, so it is
  drawn once and a marker runs around it --- which is what ties the two
  panels together: the marker says where on the ellipse the oscillogram
  currently is.

Antineutrinos need both changes, not one: the vacuum Hamiltonian is
conjugated *and* the sign of the matter potential is flipped.  Doing only
the first is the classic way to get a bi-probability plot subtly wrong.

This is a demonstration script, not part of the package.  It is deliberately
outside ``src/``: it imports matplotlib, which the library never does, and
writing a video needs ffmpeg, which is not a dependency of anything here.

Usage
-----

::

    python tools/make_demo_video.py --out demo.mp4
    python tools/make_demo_video.py --out demo.gif --frames 120 --grid 140

Writing ``.gif`` uses Pillow and needs no ffmpeg, at a considerable cost in
file size; ``.mp4`` needs ffmpeg on the path.

.. warning::

   On this machine ``ffmpeg`` is a **snap**, so it has a private ``/tmp``
   and cannot write an ``--out`` under the real one --- it fails with
   ``No such file or directory`` naming a path that plainly exists, which
   is a confusing way to be told about confinement.  Write somewhere under
   ``$HOME``.  ``gh`` in this repository has the same constraint, for the
   same reason.

   Pillow is a Python library and is not confined, so a ``.gif`` output is
   unaffected --- which makes the failure look, misleadingly, like a codec
   problem rather than a filesystem one.

A ten-second 200x200 sweep is about 400 KB as H.264 and takes roughly half
a minute to render.  The same animation as a GIF is twenty times larger, so
for anywhere that accepts video, use the mp4; where a GIF is unavoidable,
generating it through an ffmpeg palette is far smaller than writing one
directly with Pillow::

    ffmpeg -i demo.mp4 -vf "fps=12,scale=860:-1:flags=lanczos,\
palettegen=max_colors=128" palette.png
    ffmpeg -i demo.mp4 -i palette.png \
-lavfi "fps=12,scale=860:-1:flags=lanczos[x];[x][1:v]paletteuse" demo.gif
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import argparse
import os
import sys
import time

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib.animation import FuncAnimation                    # noqa: E402

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

import hamiltonians3nu                                            # noqa: E402
import oscprob3nu                                                 # noqa: E402
from globaldefs import (CONV_KM_TO_INV_EV, D21_NO_BF, D31_NO_BF,  # noqa: E402
                        S12_NO_BF, S13_NO_BF, S23_NO_BF,
                        VCC_EARTH_CRUST)

# P_mue is index 3 of the nine, the initial flavor varying slowest.
P_MUE = 3

ELLIPSE_ENERGY = 0.8e9        # [eV], near the first oscillation maximum
ELLIPSE_BASELINE = 1300.0     # [km], a DUNE-like baseline


def vacuum(dcp):
    r"""Returns the energy-independent vacuum Hamiltonian for a phase."""
    return hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent(
        S12_NO_BF, S23_NO_BF, S13_NO_BF, dcp, D21_NO_BF, D31_NO_BF)


def oscillogram(dcp, energies, baselines):
    r"""Returns P_mue over the energy-baseline grid, in one call."""
    h_stack = hamiltonians3nu.hamiltonian_3nu_matter(
        vacuum(dcp), energies, VCC_EARTH_CRUST)
    return oscprob3nu.probabilities_3nu(
        h_stack[:, None, :, :], baselines[None, :])[:, :, P_MUE]


def bi_probability(dcp):
    r"""Returns (P for neutrinos, P for antineutrinos) at one phase.

    The antineutrino Hamiltonian conjugates the vacuum term *and* reverses
    the potential; either alone is wrong.
    """
    length = ELLIPSE_BASELINE*CONV_KM_TO_INV_EV
    h_nu = hamiltonians3nu.hamiltonian_3nu_matter(
        vacuum(dcp), ELLIPSE_ENERGY, VCC_EARTH_CRUST)
    h_nubar = hamiltonians3nu.hamiltonian_3nu_matter(
        np.conj(vacuum(dcp)), ELLIPSE_ENERGY, -VCC_EARTH_CRUST)
    return (oscprob3nu.probabilities_3nu(h_nu, length)[P_MUE],
            oscprob3nu.probabilities_3nu(h_nubar, length)[P_MUE])


def build(frames, grid):
    r"""Returns the figure and its per-frame update function."""
    phases = np.linspace(0.0, 2.0*np.pi, frames, endpoint=False)
    energies = np.logspace(-1.0, 1.0, grid)*1.e9
    baselines = np.linspace(50.0, 12000.0, grid)*CONV_KM_TO_INV_EV

    # The ellipse is the locus over the whole sweep, so it is computed once
    # and drawn as a fixed curve; only the marker moves.
    ellipse = np.array([bi_probability(d)
                        for d in np.linspace(0.0, 2.0*np.pi, 240)])

    fig, (ax_map, ax_ell) = plt.subplots(
        1, 2, figsize=(11.0, 4.6), dpi=120,
        gridspec_kw={'width_ratios': [1.35, 1.0]})
    fig.patch.set_facecolor('white')

    first = oscillogram(phases[0], energies, baselines)

    # The colour scale is fixed for the whole sweep, so that a change in
    # the picture is a change in the physics rather than in the mapping.
    # It has to be the maximum over the *sweep*, not over the first frame:
    # the phase moves the maximum, so scaling to frame zero would clip
    # later frames silently, and the clipping would look like structure.
    ceiling = max(float(oscillogram(d, energies, baselines).max())
                  for d in np.linspace(0.0, 2.0*np.pi, 12, endpoint=False))

    image = ax_map.imshow(
        first, origin='lower', aspect='auto', cmap='viridis',
        vmin=0.0, vmax=ceiling,
        extent=[baselines[0]/CONV_KM_TO_INV_EV,
                baselines[-1]/CONV_KM_TO_INV_EV,
                np.log10(energies[0]/1.e9), np.log10(energies[-1]/1.e9)])
    ax_map.set_xlabel('Baseline [km]')
    ax_map.set_ylabel('Energy [GeV]')
    ax_map.set_yticks([-1, 0, 1])
    ax_map.set_yticklabels(['0.1', '1', '10'])
    ax_map.set_title(r'$P_{\mu e}$ in matter, %d$\,\times\,$%d grid'
                     % (grid, grid), fontsize=11)
    bar = fig.colorbar(image, ax=ax_map, pad=0.02)
    bar.set_label(r'$P_{\mu e}$')

    ax_ell.plot(ellipse[:, 0], ellipse[:, 1], color='#1d4ed8', lw=1.6)
    lo = float(min(ellipse.min(), 0.0))
    hi = float(ellipse.max())*1.08
    ax_ell.plot([lo, hi], [lo, hi], ls=':', color='#64748b', lw=1.0,
                label='CP symmetric')
    marker, = ax_ell.plot([], [], 'o', ms=9, color='#dc2626', zorder=5)
    ax_ell.set_xlim(lo, hi)
    ax_ell.set_ylim(lo, hi)
    ax_ell.set_xlabel(r'$P(\nu_\mu \to \nu_e)$')
    ax_ell.set_ylabel(r'$P(\bar{\nu}_\mu \to \bar{\nu}_e)$')
    ax_ell.set_title(r'Bi-probability, %.1f GeV, %d km'
                     % (ELLIPSE_ENERGY/1.e9, ELLIPSE_BASELINE), fontsize=11)
    ax_ell.legend(loc='lower right', fontsize=9, frameon=False)
    ax_ell.set_aspect('equal', adjustable='box')

    caption = fig.text(
        0.5, 0.028,
        'NuOscProbExact  ---  every frame is one broadcast call, '
        '%s probabilities' % format(grid*grid, ','),
        ha='center', fontsize=9.5, color='#334155')
    heading = fig.suptitle('', fontsize=13)
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])

    def update(index):
        dcp = phases[index]
        image.set_data(oscillogram(dcp, energies, baselines))
        marker.set_data(*[[value] for value in bi_probability(dcp)])
        heading.set_text(r'$\delta_{CP} = %.2f\pi$' % (dcp/np.pi))
        return image, marker, heading, caption

    return fig, update, len(phases)


def main():
    r"""Parses the arguments and writes the animation."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', default='nuoscprobexact_demo.mp4',
                        help='output file; .mp4 needs ffmpeg, .gif does not')
    parser.add_argument('--frames', type=int, default=240,
                        help='number of frames in the sweep')
    parser.add_argument('--fps', type=int, default=24, help='frames a second')
    parser.add_argument('--grid', type=int, default=200,
                        help='oscillogram is grid x grid')
    args = parser.parse_args()

    started = time.perf_counter()
    fig, update, frames = build(args.frames, args.grid)
    animation = FuncAnimation(fig, update, frames=frames, blit=False)

    writer = 'pillow' if args.out.lower().endswith('.gif') else 'ffmpeg'
    animation.save(args.out, writer=writer, fps=args.fps,
                   savefig_kwargs={'facecolor': 'white'})
    plt.close(fig)

    elapsed = time.perf_counter() - started
    size = os.path.getsize(args.out)/1024.0
    print('wrote %s  --  %d frames, %.1f s of video, %.0f KB, %.1f s to make'
          % (args.out, frames, frames/float(args.fps), size, elapsed))


if __name__ == '__main__':
    main()
