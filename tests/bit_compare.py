"""Compare two captures for exact bit equality.

Run as:  python bit_compare.py <baseline.npz> <candidate.npz>

Exit status 0 only if every value in every array is bitwise identical.
Anything else prints which arrays moved and by how much, in ulps as well
as absolutely, so that a failure is a measurement rather than a verdict.
"""
import sys

import numpy as np


def ulp_distance(a, b):
    r"""Returns the largest gap between two float arrays, in ulps."""
    a = np.asarray(a)
    b = np.asarray(b)
    if np.iscomplexobj(a):
        return max(ulp_distance(a.real, b.real), ulp_distance(a.imag, b.imag))
    diff = np.abs(a - b)
    spacing = np.spacing(np.maximum(np.abs(a), np.abs(b)))
    return float(np.max(np.where(diff == 0.0, 0.0, diff/spacing)))


baseline = np.load(sys.argv[1])
candidate = np.load(sys.argv[2])

missing = sorted(set(baseline.files) - set(candidate.files))
added = sorted(set(candidate.files) - set(baseline.files))
if missing or added:
    print('FAIL: the two captures do not cover the same ground')
    for name in missing[:10]:
        print('  only in baseline :', name)
    for name in added[:10]:
        print('  only in candidate:', name)
    sys.exit(2)

moved = []
values = 0
for name in baseline.files:
    a, b = baseline[name], candidate[name]
    values += a.size
    if a.shape != b.shape:
        moved.append((name, float('inf'), float('inf'), 'shape changed'))
    elif not np.array_equal(a, b):
        moved.append((name, float(np.max(np.abs(a - b))),
                      ulp_distance(a, b), ''))

print('arrays compared : %d' % len(baseline.files))
print('values compared : %d' % values)

if not moved:
    print('VERDICT: bit-identical')
    sys.exit(0)

moved.sort(key=lambda row: -row[1])
print('VERDICT: NOT bit-identical -- %d of %d arrays moved'
      % (len(moved), len(baseline.files)))
print()
print('%-52s %12s %8s' % ('array', 'max |diff|', 'max ulp'))
for name, absolute, ulps, note in moved[:25]:
    print('%-52s %12.3e %8.1f %s' % (name, absolute, ulps, note))
if len(moved) > 25:
    print('... and %d more' % (len(moved) - 25))
sys.exit(1)
