# -*- coding: utf-8 -*-
r"""Renders the NuOscProbExact demonstration reel, or one scene of it.

The scenes live in :mod:`demo_scenes`; this is the command line around
them.  ``--scene all`` renders each in turn and, if ffmpeg is available,
concatenates them into a single reel.

Why an animation is a reasonable thing to build here at all: for the two
map scenes, every frame is a *single* broadcast call, so the physics is
around a tenth of the cost of drawing it.  A library that had to loop in
Python could not be animated this way, which is the property being shown.
The Earth scene is the deliberate exception --- see `demo_scenes`.

Usage
-----

::

    python tools/make_demo_video.py --scene cp      --out ~/demo/cp.mp4
    python tools/make_demo_video.py --scene all     --out ~/demo/reel.mp4
    python tools/make_demo_video.py --scene slabs   --out ~/demo/slabs.gif

Writing ``.gif`` goes through Pillow and needs no ffmpeg, at roughly twenty
times the file size; ``.mp4`` needs ffmpeg on the path.

.. warning::

   On this machine ``ffmpeg`` is a **snap**, so it has a private ``/tmp``
   and cannot write an ``--out`` under the real one.  It fails with
   ``No such file or directory`` naming a path that plainly exists, which
   is a confusing way to be told about confinement.  Write somewhere under
   ``$HOME``.  ``gh`` in this repository has the same constraint.

   Pillow is not confined, so a ``.gif`` from the same run succeeds ---
   which makes the failure look like a codec problem rather than a
   filesystem one.  It is neither.

Where a GIF is unavoidable, going through an ffmpeg palette is far smaller
than writing one directly::

    ffmpeg -i reel.mp4 -vf "fps=12,scale=860:-1:flags=lanczos,\
palettegen=max_colors=128" palette.png
    ffmpeg -i reel.mp4 -i palette.png \
-lavfi "fps=12,scale=860:-1:flags=lanczos[x];[x][1:v]paletteuse" reel.gif
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import argparse
import os
import shutil
import subprocess
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
from matplotlib.animation import FuncAnimation                    # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from demo_scenes import SCENES                                    # noqa: E402

# The order the reel plays in: what it computes, then what it computes it
# on, then what it costs to be honest about a continuous profile.
REEL = ['cp', 'sterile', 'earth', 'slabs']


def render(name, path, fps):
    r"""Renders one scene to `path` and returns how long it took."""
    started = time.perf_counter()
    fig, update, frames = SCENES[name]()
    animation = FuncAnimation(fig, update, frames=frames, blit=False)
    writer = 'pillow' if path.lower().endswith('.gif') else 'ffmpeg'
    animation.save(path, writer=writer, fps=fps,
                   savefig_kwargs={'facecolor': 'white'})
    plt.close(fig)
    elapsed = time.perf_counter() - started
    print('  %-8s %3d frames  %6.1f s of video  %7.0f KB  %5.1f s to render'
          % (name, frames, frames/float(fps),
             os.path.getsize(path)/1024.0, elapsed))
    return elapsed


def concatenate(parts, path):
    r"""Joins rendered clips with ffmpeg's concat demuxer.

    Stream copy rather than re-encode: every scene is rendered at the same
    size and frame rate by construction, so there is nothing to reconcile
    and no reason to lose a generation of quality.

    The listing is written *beside the output*, not to a temporary
    directory.  `tempfile` would put it under ``/tmp``, which the snap
    ffmpeg cannot read --- the same confinement the module docstring warns
    about for ``--out``, and which this function fell foul of first time
    round.  If ffmpeg can write the output it can read a file next to it.
    """
    listing = os.path.join(os.path.dirname(os.path.abspath(path)),
                           '.concat_parts.txt')
    with open(listing, 'w') as handle:
        for part in parts:
            handle.write("file '%s'\n" % os.path.abspath(part))
    subprocess.run(
        ['ffmpeg', '-hide_banner', '-loglevel', 'error', '-y',
         '-f', 'concat', '-safe', '0', '-i', listing, '-c', 'copy', path],
        check=True)
    os.remove(listing)


def main():
    r"""Parses the arguments and writes the animation."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--scene', default='cp',
                        choices=sorted(SCENES) + ['all'],
                        help='which scene to render, or all of them')
    parser.add_argument('--out', default='nuoscprobexact_demo.mp4',
                        help='output file; .mp4 needs ffmpeg, .gif does not')
    parser.add_argument('--fps', type=int, default=24, help='frames a second')
    args = parser.parse_args()

    if os.path.dirname(args.out):
        os.makedirs(os.path.dirname(args.out), exist_ok=True)

    started = time.perf_counter()

    if args.scene != 'all':
        render(args.scene, args.out, args.fps)
    else:
        if args.out.lower().endswith('.gif'):
            parser.error('--scene all writes one reel, which needs ffmpeg; '
                         'use an .mp4 output, or render scenes separately')
        stem, _ = os.path.splitext(args.out)
        parts = []
        for name in REEL:
            part = '%s_%s.mp4' % (stem, name)
            render(name, part, args.fps)
            parts.append(part)
        concatenate(parts, args.out)
        for part in parts:
            os.remove(part)
        print('  joined %d scenes -> %s (%.0f KB)'
              % (len(parts), args.out, os.path.getsize(args.out)/1024.0))

    print('total %.1f s' % (time.perf_counter() - started))


if __name__ == '__main__':
    if shutil.which('ffmpeg') is None:
        print('note: ffmpeg is not on the path; only .gif output will work',
              file=sys.stderr)
    main()
