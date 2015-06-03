import sys

from ase.test import NotAvailable, must_raise

if sys.platform in ['win32']:
    raise NotAvailable('Fails on Windows '
                       'https://trac.fysik.dtu.dk/projects/ase/ticket/62')

import os
from ase import Atom, Atoms
from ase.io import Trajectory, read

co = Atoms([Atom('C', (0, 0, 0)),
            Atom('O', (0, 0, 1.2))])
traj = Trajectory('1.traj', 'w', co)
for i in range(5):
    co.positions[:, 2] += 0.1
    traj.write()

traj = Trajectory('1.traj', 'a')
co = read('1.traj')
print(co.positions)
co.positions[:] += 1
traj.write(co)

for a in Trajectory('1.traj'):
    print(1, a.positions[-1, 2])
co.positions[:] += 1
t = Trajectory('1.traj', 'a')
t.write(co)
assert len(t) == 7

co[0].number = 1
with must_raise(ValueError):
    t.write(co)

co[0].number = 6
co.pbc = True
with must_raise(ValueError):
    t.write(co)

co.pbc = False
o = co.pop(1)
with must_raise(ValueError):
    t.write(co)

co.append(o)
t.write(co)

# append to a nonexisting file:
fname = '2.traj'
if os.path.isfile(fname):
    os.remove(fname)
t = Trajectory(fname, 'a', co)
os.remove(fname)

t = Trajectory('empty.traj', 'w')
t.close()
assert os.path.getsize('empty.traj') == 0
