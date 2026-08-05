# -*- coding: utf-8 -*-
r"""Join and shrink the demonstration clips.

The scenes themselves live in ``notebooks/19_animations.ipynb``, which
draws them as stills and, with ``RENDER = True``, writes each as a GIF
into ``img/``.  This is the step after that: joining those clips into one
reel, and getting a GIF down to a size worth putting in a README.

Splitting it that way is deliberate.  The notebook owns the physics and is
regenerated and executed by CI, so it cannot drift; this owns the
encoding, which is where the useful knowledge is and which has nothing to
do with neutrinos.  An earlier version of this script carried its own copy
of the four scenes, and that copy had already drifted --- it still
explained that an Earth curve of a hundred energies was a hundred calls,
which stopped being true in 1.12.0.

Usage
-----

::

    # 1. Render the scenes: set RENDER = True in notebook 19.
    # 2. Join them into one reel:
    python tools/make_demo_video.py --join img/anim_*.gif --out ~/reel.mp4

    # 3. Or shrink one clip, which is what makes a GIF publishable:
    python tools/make_demo_video.py --shrink ~/reel.mp4 --out ~/reel.gif

Why a palette, and why it matters
---------------------------------

Writing a GIF straight out of matplotlib's Pillow writer gives every
frame its own colour table.  Going through a *shared* palette --- one
pass to compute it, one to apply it --- removes that duplication, and on
its own is worth about 1.6x, measured on a thirty-frame clip::

    ffmpeg -i reel.mp4 -vf "fps=12,scale=860:-1:flags=lanczos,\
palettegen=max_colors=128" palette.png
    ffmpeg -i reel.mp4 -i palette.png \
-lavfi "fps=12,scale=860:-1:flags=lanczos[x];[x][1:v]paletteuse" reel.gif

``--shrink`` runs exactly those two commands.  The palette is only part
of it, and worth saying plainly: the large factors --- the committed
clips are around a twentieth of what Pillow wrote --- come from
combining it with a lower frame rate and a smaller width, not from the
palette alone.  The three knobs are the frame rate, the width and the
number of colours, and they trade size against fidelity in that order
of effect.  A dense, high-frequency image
--- the sterile-state map is the worst of the four --- loses visibly to a
small palette where a smooth one does not, so check the scene you care
about rather than trusting a default.

.. warning::

   On this machine ``ffmpeg`` is a **snap**, so it has a private ``/tmp``
   and cannot read or write under the real one.  It fails with ``No such
   file or directory`` naming a path that plainly exists, which is a
   confusing way to be told about confinement.  Work somewhere under
   ``$HOME``.  ``gh`` in this repository has the same constraint.

   Pillow is not confined, so a GIF written by the notebook succeeds from
   the same directory that fails here --- which makes the failure look
   like a codec problem rather than a filesystem one.  It is neither.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import argparse
import os
import shutil
import subprocess
import sys


def _run(command):
    r"""Runs ffmpeg quietly, and says something useful when it fails."""
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        sys.exit('ffmpeg failed:\n  %s\n%s'
                 % (' '.join(command), result.stderr.strip()[-800:]))


def _width_of(path):
    r"""Returns a clip's pixel width, so that nothing is upsampled.

    Read with Pillow rather than ``ffprobe``: the snap ffmpeg does not
    ship one, and a dependency that is present for the GIFs the notebook
    writes is better than one that is not.
    """
    try:
        from PIL import Image
        with Image.open(path) as image:
            return image.size[0]
    except Exception:                                       # noqa: BLE001
        return 860


def shrink(source, out, fps, width, colors):
    r"""Writes `source` as a GIF through a shared palette.

    Two passes: one to choose a palette from the whole clip, one to apply
    it.  A palette per frame --- which is what writing a GIF directly
    gives --- is the reason those files are so large.
    """
    palette = os.path.join(os.path.dirname(os.path.abspath(out)),
                           '.palette.png')
    chain = 'fps=%d,scale=%d:-1:flags=lanczos' % (fps, width)
    _run(['ffmpeg', '-hide_banner', '-loglevel', 'error', '-y',
          '-i', source, '-vf', '%s,palettegen=max_colors=%d' % (chain, colors),
          palette])
    _run(['ffmpeg', '-hide_banner', '-loglevel', 'error', '-y',
          '-i', source, '-i', palette,
          '-lavfi', '%s[x];[x][1:v]paletteuse' % chain, out])
    os.remove(palette)


def join(parts, out, fps, width):
    r"""Concatenates clips into one, at one frame rate and one width.

    The concat *filter* rather than the demuxer, because the parts may be
    GIFs written by the notebook rather than streams that can be copied.

    Both ``-r`` and the scale are explicit, and neither is optional in
    practice.  A GIF stores its delays in hundredths of a second, so
    ffmpeg reads a 90 ms frame as a stream of about 100 fps and, left to
    itself, writes every frame nine times over --- which turned a
    four-hundred-frame reel into ninety megabytes and climbing before it
    was noticed.  Setting the output rate is what stops that.  Scaling
    defaults to the width of the first clip, so the common case of
    joining clips that already match neither upsamples nor resamples
    them.
    """
    inputs = []
    for part in parts:
        inputs += ['-i', part]
    scale = ''.join(
        '[%d:v]fps=%d,scale=%d:-2:flags=lanczos,setsar=1[v%d];'
        % (index, fps, width, index) for index in range(len(parts)))
    chain = ''.join('[v%d]' % index for index in range(len(parts)))
    _run(['ffmpeg', '-hide_banner', '-loglevel', 'error', '-y'] + inputs
         + ['-filter_complex',
            '%s%sconcat=n=%d:v=1[out]' % (scale, chain, len(parts)),
            '-map', '[out]', '-r', str(fps), '-crf', '20',
            '-pix_fmt', 'yuv420p', out])


def main():
    r"""Parses the arguments and does the one thing asked for."""
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        epilog='The scenes are rendered by notebooks/19_animations.ipynb; '
               'this joins and shrinks what it writes.')
    parser.add_argument('--join', nargs='+', metavar='CLIP',
                        help='clips to concatenate, in order')
    parser.add_argument('--shrink', metavar='CLIP',
                        help='one clip to rewrite as a small GIF')
    parser.add_argument('--out', required=True, help='the file to write')
    parser.add_argument('--fps', type=int, default=12,
                        help='frames a second, when shrinking (default 12)')
    parser.add_argument('--width', type=int, default=0,
                        help='width in pixels; 0 keeps the first clip\'s')
    parser.add_argument('--colors', type=int, default=128,
                        help='palette size, when shrinking (default 128)')
    args = parser.parse_args()

    if bool(args.join) == bool(args.shrink):
        parser.error('give exactly one of --join and --shrink')
    if shutil.which('ffmpeg') is None:
        sys.exit('ffmpeg is not on the path; both of these need it')
    if os.path.dirname(args.out):
        os.makedirs(os.path.dirname(args.out), exist_ok=True)

    width = args.width or _width_of(args.join[0] if args.join
                                    else args.shrink)
    if args.join:
        join(args.join, args.out, args.fps, width)
        print('joined %d clips -> %s (%.1f MB)'
              % (len(args.join), args.out,
                 os.path.getsize(args.out)/1024.0/1024.0))
    else:
        before = os.path.getsize(args.shrink)/1024.0/1024.0
        shrink(args.shrink, args.out, args.fps, width, args.colors)
        after = os.path.getsize(args.out)/1024.0/1024.0
        print('%s -> %s   %.1f MB -> %.1f MB  (%.1fx smaller)'
              % (args.shrink, args.out, before, after,
                 before/after if after else float('inf')))


if __name__ == '__main__':
    main()
