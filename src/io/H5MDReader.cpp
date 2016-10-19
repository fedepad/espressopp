/*
  Copyright (C) 2016
      JGU Mainz

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


// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "H5MDReader.hpp"  // keep python.hpp on top
#include "storage/Storage.hpp"
//#include "System.hpp"
#include "storage/DomainDecomposition.hpp"
#include "bc/BC.hpp"

#include "iterator/CellListIterator.hpp"

#include "boost/mpi.hpp"
#include "boost/mpi/communicator.hpp"
#include "mpi.h"

#include <fstream>
#include <sstream>

//#ifdef HDF5_LAYER
    //#include "hdf5.h"
    //#include "hdf5_hl.h"
//#endif

//#ifdef H5MD_LAYER
//    #include "ch5md.h"
//#endif

using namespace espressopp;
using namespace std;
using namespace boost::python;


namespace espressopp {
  namespace io {



    // Python wrapping
    void H5MDReader::registerPython() {

      using namespace espressopp::python;

      class_<H5MDReader, boost::noncopyable >
      ("io_H5MDReader", init< std::string>())
//        .add_property("filename", &H5MDReader::getFilename,
//                                  &H5MDReader::setFilename)
//	    .add_property("author", &H5MDReader::getFilename,
//								  &H5MDReader::setFilename)
//	    .add_property("author_email", &H5MDReader::getFilename,
//								  &H5MDReader::setFilename)
//	    .add_property("creator", &H5MDReader::getFilename,
//								  &H5MDReader::setFilename)
//		.add_property("creator_version", &H5MDReader::getFilename,
//								  &H5MDReader::setFilename)
//	    .add_property("iomode", &H5MDReader::getIomode,
//								&H5MDReader::setIomode)
//		.add_property("data_to_store", &H5MDReader::getDataToStore,
//										&H5MDReader::set_init_table)
//
//        .add_property("unfolded", &H5MDReader::getUnfolded,
//                                  &H5MDReader::setUnfolded)
//        .add_property("length_factor", &H5MDReader::getLengthFactor,
//                                       &H5MDReader::setLengthFactor)
//        .add_property("length_unit", &H5MDReader::getLengthUnit,
//                                     &H5MDReader::setLengthUnit)
//        .add_property("append", &H5MDReader::getAppend,
//                                  &H5MDReader::setAppend)
//        .def("write", &H5MDReader::write)
      ;
    } // end registerPython
  } // end namespace io
} // end namespace espressopp

