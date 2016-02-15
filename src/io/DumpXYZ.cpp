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

// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "DumpXYZ.hpp"  // keep python.hpp on top
#include "storage/Storage.hpp"
#include "Particle.hpp"
#include "bc/BC.hpp"
#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"

#include <fstream>

using namespace espressopp;
using namespace espressopp::analysis;
using namespace std;

namespace espressopp {
  namespace io {
      
    void DumpXYZ::dump(){
      shared_ptr<System> system = getSystem();
      ConfigurationsExt conf( system );
      conf.setUnfolded(unfolded);
      conf.gather_modified();
      if( system->comm->rank()==0 ){
        ConfigurationExtPtr conf_real = conf.back();
        
        size_t num_of_particles = conf_real->getSize();

        char *ch_f_name = new char[file_name.length() + 1];
        strcpy(ch_f_name, file_name.c_str());
        ofstream myfile (ch_f_name, ios::out | ios::app);
        if (myfile.is_open()){
          
          myfile << num_of_particles << endl;
          
          Real3D Li = system->bc->getBoxL();
          

          //python::list pids_part = system->storage->getRealParticleIDs();

          // for noncubic simulation boxes
          //TODO: let's decide if you want these hanging 0.0s as in the default repo version or get rid!
          myfile << Li[0] * length_factor << "  0.0  0.0  0.0  "<< 
                  Li[1] * length_factor << "  0.0  0.0  0.0  "<< Li[2] * length_factor;
          
          // most probable fix to the above line...
          //myfile << Li[0] * length_factor << " "<< Li[1] * length_factor << " "<< Li[2] * length_factor;
          
          // additional info to comment line
          myfile << "  currentStep " << integrator->getStep() << "  lengthUnit "<< length_unit << endl;
          
          ConfigurationExtIterator cei = conf_real-> getIterator();
          
          for (size_t i = 0; i < num_of_particles; i++) {
          	  size_t idi = cei.Id();
        	  size_t tipolo = conf_real->getType(idi);// better without querying pid. don't care now
        	  myfile << tipolo << " " << length_factor * cei.Properties() << endl; // strictly XYZ/V: type XYZ VXVYVZ
        	  cei.nextParticle();
          }

          myfile.close();
        }
        else cout << "Unable to open file: "<< file_name <<endl;

        delete [] ch_f_name;
      }
    }
      
    // Python wrapping
    void DumpXYZ::registerPython() {

      using namespace espressopp::python;

      class_<DumpXYZ, bases<ParticleAccess>, boost::noncopyable >
      ("io_DumpXYZ", init< shared_ptr< System >, 
                           shared_ptr< integrator::MDIntegrator >, 
                           std::string, 
                           bool,
                           real,
                           std::string , 
                           bool>())
        .add_property("filename", &DumpXYZ::getFilename, 
                                  &DumpXYZ::setFilename)
        .add_property("unfolded", &DumpXYZ::getUnfolded, 
                                  &DumpXYZ::setUnfolded)
        .add_property("length_factor", &DumpXYZ::getLengthFactor, 
                                       &DumpXYZ::setLengthFactor)
        .add_property("length_unit", &DumpXYZ::getLengthUnit, 
                                     &DumpXYZ::setLengthUnit)
        .add_property("append", &DumpXYZ::getAppend, 
                                  &DumpXYZ::setAppend)
        .def("dump", &DumpXYZ::dump)
      ;
    }
  }
}
