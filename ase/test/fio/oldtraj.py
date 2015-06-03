"""Check that we can read old version 1 Trajectories."""
import functools
import pickle
import sys
from io import BytesIO

from ase import Atoms
from ase.constraints import FixAtoms
from ase.io.pickletrajectory import PickleTrajectory
from ase.test import NotAvailable

if sys.version_info[0] == 3:
    raise NotAvailable


Trajectory = functools.partial(PickleTrajectory, _warn=False)

a = Atoms('FOO')


def v1(a):
    """Create old version-1 trajectory."""
    fd = BytesIO()
    fd.write(b'PickleTrajectory')
    d = {'pbc': a.pbc,
         'numbers': a.numbers,
         'tags': None,
         'masses': None,
         'constraints': a.constraints}
    pickle.dump(d, fd, protocol=-1)
    d = {'positions': a.positions,
         'cell': a.cell,
         'momenta': None}
    pickle.dump(d, fd, protocol=-1)
    return BytesIO(fd.getvalue())


def v2(a):
    """Create new version-2 trajectory."""
    fd = BytesIO()
    t = Trajectory(fd, 'w')
    t.write(a)
    return BytesIO(fd.getvalue())


class MyFixAtoms(FixAtoms):
    pass


a.constraints = FixAtoms(indices=[2])
t1 = v1(a)
t2 = v2(a)

# Read old trajectory:
c1 = Trajectory(t1)[0].constraints
assert c1[0].index[0] == 2

# Read new trajectory:
c2 = Trajectory(t2)[0].constraints
assert c2[0].index[0] == 2

a.constraints = MyFixAtoms(indices=[1])
t3 = v2(a)

# Read new trajectory with missing constraint class.  This should
# ignore the constraint and issue a warning:
del MyFixAtoms, a
import warnings
warnings.filterwarnings('error')
try:
    c3 = Trajectory(t3)[0].constraints
except UserWarning:
    pass
else:
    assert False
