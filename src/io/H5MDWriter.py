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
**H5MDWriter** - IO Object
*********************************************

* `dump()`
  write configuration to H5MD file. By default filename is "out.h5",
  coordinates are folded.
  
* `close()`
  write configuration to H5MD file. By default filename is "out.h5",
  coordinates are folded.

  Properties

* `filename`
  Name of trajectory file. By default trajectory file name is "out.h5"

* `author`
  Name of the user that is running the program.

* `author_email`
  Email of the person that is running the program.

* `creator`
  Program that created the file. Default: EspResso++

* `creator_version`
  Version of the program that created the file.

* `iomode`
  Iomode: 0 serial, 1 Nto1, 2 NtoN. Right now N-to-1 is activated by default.

* `data_to_store`
  List containing the particle properties that one wants to save: ['all'] saves the
  all particle data struct, i.e. (pid, position, velocity, type, force, charge, mass,
  state, drift, lambda, lambdaDeriv) or ['pid', 'mass', 'charge'] means only the pid,
  mass and charge of each particle need to be saved. Default - ['all']
  Example:
      data_to_store = ['all']  --> saves the all particle struct
      data_to_store = ['pid', 'position', 'type'] --> save pid, position and type
      for each particle
  Possible fields per particle to be saved:
  pid, type, position, velocity, force, image, charge, mass, state, drift, lambda, lambdaDeriv.

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

>>> dump_conf_h5md = espressopp.io.H5MDWriter(system, integrator, filename='trajectory.h5', iomode=1)
>>> for i in range (200):
>>>   integrator.run(10)
>>>   dump_conf_h5md.dump()

writing down trajectory using ExtAnalyze extension

>>> dump_conf_h5md = espressopp.io.H5MDWriter(system, integrator, filename='trajectory.h5', iomode=1, data_to_store=['all'])
>>> ext_analyze = espressopp.integrator.ExtAnalyze(dump_conf_h5md, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

Both examples will give the same result: 200 configurations in trajectory .h5 file.

setting up length scale

For example, the Lennard-Jones model for liquid argon with :math:`\sigma=0.34 [nm]`

>>> dump_conf_h5md = espressopp.io.H5MDWriter(system, integrator, filename='trj.h5', iomode=1, data_to_store=['pid', 'mass', 'velocity'], unfolded=False, length_factor=0.34, length_unit='nm', append=True)

will produce trj.h5 with  in nanometers // Federico P. comment: what in nanometers? It's clear coordinate but please don't leave sentence hanging!

.. function:: espressopp.io.H5MDWriter(system, integrator, filename, iomode, data_to_store, unfolded, length_factor, length_unit, append)

        :param system:
        :param integrator:
        :param filename: (default: 'out.h5')
        :param author: (default: 'out.h5')
        :param author_email: (default: 'NULL')
        :param creator: (default: 'EsPResso++')
        :param creator_version: (default: xxxx)

        :param iomode: (default: 1)
        :param data_to_store: (default: ['all'])
        :param unfolded: (default: False)
        :param length_factor: (default: 1.0)
        :param length_unit: (default: 'LJ')
        :param append: (default: True)
        :param writers: (default: 0)
        :type system:
        :type integrator:
        :type filename: string
        :type author: string
        :type author_email: string
        :type creator: string
        :type creator_version: string
        :type iomode: int
        :type data_to_store: list of strings
        :type unfolded: bool
        :type length_factor: real
        :type length_unit: real
        :type append: bool
        :type writers: int

.. function:: espressopp.io.H5MDWriter.write()

        :rtype:
"""

import os
from espressopp.esutil import cxxinit
from espressopp import pmi
import espressopp

from espressopp.ParticleAccess import *
from _espressopp import io_H5MDWriter

class H5MDWriterLocal(ParticleAccessLocal, io_H5MDWriter):

  def __init__(self, system, integrator, filename=os.path.join(os.getcwd(), 'out.h5'), author=os.environ['USER'],
               author_email="username@domain.state", creator = espressopp.VersionLocal().name,
               #creator_version = ".".join([str(espressopp.VersionLocal().major),
               #                           str(espressopp.VersionLocal().minor),
               #                           str(espressopp.VersionLocal().patchlevel)]),
               creator_version = espressopp.VersionLocal().info(),
                iomode=1, data_to_store=['all'], unfolded=False, length_factor=1.0,
               length_unit='LJ', sort_pids=False, append=True, writers=1):
    cxxinit(self, io_H5MDWriter, system, integrator, filename, author,
            author_email, creator, creator_version, iomode, data_to_store,
            unfolded, length_factor, length_unit, sort_pids, append, writers)

  def write(self):
    if pmi.workerIsActive():
      self.cxxclass.write(self)


if pmi.isController :
  class H5MDWriter(ParticleAccess):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls =  'espressopp.io.H5MDWriterLocal',
      pmicall = [ 'dump', 'close' ],
      pmiproperty = ['filename', 'author', 'author_email', 'creator', 'creator_version', 'iomode', 'data_to_store', 'unfolded', 'length_factor', 'length_unit', 'sort_pids', 'append', 'writers']
    )
