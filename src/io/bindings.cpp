/*
  Copyright (C) 2014,2015,2016
      Max Planck Institute for Polymer Research & Johannes Gutenberg-Universität Mainz
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "bindings.hpp"
#include "DumpXYZ.hpp"
#include "DumpGRO.hpp"
#include "DumpGROAdress.hpp"
#ifdef HDF5_LAYER
    #include "H5MDFile.hpp"
#endif
#include "FileBackup.hpp"

namespace espressopp {
  namespace io{
    void registerPython() {
      DumpXYZ::registerPython();
      DumpGRO::registerPython();
      DumpGROAdress::registerPython();
#ifdef HDF5_LAYER    
      H5MDFile::registerPython();
#endif
    }
  }
}
