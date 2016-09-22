#  Copyright (C) 2016
#      Max Planck Institute for Polymer Research & JGU Mainz
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""
*********************************************
**H5MD** - IO Object
*********************************************

* `write()`
  write configuration to trajectory H5MD file. By default filename is "out.xyz",
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.xyz"

* `iomode`
  Iomode: 0 serial, 1 Nto1, 2 NtoN

* `data_to_store`
  List containing the particle properties one wants to save: ['all'] saves the
  all particle data struct, i.e. (pid, pos, vel, lambda, force, charge, mass,
  etc.) or ['pid', 'mass', 'charge'] means only the pid,
  mass and charge of each particle need to be saved. Default - ['all']

* `unfolded`
  False if coordinates are folded, True if unfolded. By default - False

* `append`
  True if new trajectory data is appended to existing trajectory file. By default - True

* `length_factor`
  If length dimension in current system is nm, and unit is 0.23 nm, for example, then
  length_factor should be 0.23

* `length_unit`
  It is length unit. Can be LJ, nm or A. By default - LJ

usage:

writing down trajectory

>>> dump_conf_h5md = espressopp.io.H5MD(system, integrator, filename='trajectory.h5', iomode=1)
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_h5md.write()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_h5md = espressopp.io.H5MD(system, integrator, filename='trajectory.h5', iomode=1, data_to_store=['all'])
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_conf_h5md, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both examples will give the same result: 200 configurations in trajectory .h5 file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]`

>>> dump_conf_h5md = espressopp.io.H5MD(system, integrator, filename='trj.h5', iomode=1, data_to_store=['pid', 'mass', 'velocity'], unfolded=False, length_factor=0.34, length_unit='nm', append=True)

will produce trj.h5 with  in nanometers // Federico P. comment: what in nanometers? It's clear coordinate but please don't leave sentence hanging!

.. function:: espressopp.io.H5MD(system, integrator, filename, iomode, data_to_store, unfolded, length_factor, length_unit, append)

        :param system:
        :param integrator:
        :param filename: (default: 'out.h5')
        :param iomode: (default: 1)
        :param data_to_store: (default: ['all'])
        :param unfolded: (default: False)
        :param length_factor: (default: 1.0)
        :param length_unit: (default: 'LJ')
        :param append: (default: True)
        :type system:
        :type integrator:
        :type filename:
        :type iomode:
        :type data_to_store:
        :type unfolded:
        :type length_factor: real
        :type length_unit:
        :type append:

.. function:: espressopp.io.H5MD.write()

        :rtype:
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.ParticleAccess import *
from _espressopp import io_H5MD

class H5MDLocal(ParticleAccessLocal, io_H5MD):

  def __init__(self, system, integrator, filenam e='out.h5', iomode=1, data_to_store=['all'], unfolded=False, length_factor=1.0, length_unit='LJ', append=True):
    cxxinit(self, io_H5MD, system, integrator, filename, iomode, data_to_store, unfolded, length_factor, length_unit, append)

  def write(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.write(self)
      
  def read(self):
    if not pmi._PMIComm or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
      self.cxxclass.read(self)


if pmi.isController :
  class H5MD(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.io.H5MDLocal',
      pmicall = [ 'write', 'read' ],
      pmiproperty = ['filename', 'iomode', 'data_to_store', 'unfolded', 'length_factor', 'length_unit', 'append']
    )
