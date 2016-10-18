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
**H5MDReader** - IO Object
*********************************************

* `write()`
  write configuration to H5MD file. By default filename is "out.h5",
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.h5"

usage:

writing down trajectory

>>> read_conf_h5md = espressopp.io.H5MDReader(filename='trajectory.h5')
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_h5md.write()

reading trajectory

>>> dump_conf_h5md = espressopp.io.H5MDReader(filename='trajectory.h5')
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_conf_h5md, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both examples will give the same result: 200 configurations in trajectory .h5 file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]`

>>> dump_conf_h5md = espressopp.io.H5MDReader(system, integrator, filename='trj.h5', iomode=1, data_to_store=['pid', 'mass', 'velocity'], unfolded=False, length_factor=0.34, length_unit='nm', append=True)

will produce trj.h5 with  in nanometers // Federico P. comment: what in nanometers? It's clear coordinate but please don't leave sentence hanging!

.. function:: espressopp.io.H5MDReader(filename)

        :param filename: (default: 'out.h5')

        :type filename: string

.. function:: espressopp.io.H5MDReader.read()

        :rtype:
"""

import os
from espressopp.esutil import cxxinit
from espressopp import pmi
import espressopp

#from espressopp.ParticleAccess import *
from _espressopp import io_H5MDReader

class H5MDReaderLocal(io_H5MDReader):

  def __init__(self, filename):
    cxxinit(self, io_H5MDReader, filename)

  def read(self):
    if pmi.workerIsActive():
      self.cxxclass.read(self)


if pmi.isController :
  class H5MDReader(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.io.H5MDReaderLocal',
      pmicall = [ 'read' ],
      pmiproperty = ['filename']
    )
