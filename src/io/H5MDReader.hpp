/*
  Copyright (C) 2016
      Max Planck Institute for Polymer Research & JGU Mainz

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

// ESPP_CLASS
#ifndef _IO_H5MDReader_HPP
#define _IO_H5MDReader_HPP

// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "ParticleAccess.hpp"  // keep python.hpp on top
#include "System.hpp"
#include "integrator/MDIntegrator.hpp"
#include "esutil/Error.hpp"
#include "Version.hpp"
#include "template_py_to_cpp_containers.hpp"
#include <string>
#include <set>

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/args.hpp>

#include "ch5md.h"

namespace espressopp {
  namespace io{

    struct info_table_ini{ // what would be the advantage to have it as bools?
        int all;
        int pid;
        int type;
        int mass;
        int charge;
        int lambda;
        int drift;
        int lambdaDeriv;
        int state;
        int position;
        int velocity;
        int force;
        info_table_ini() : all(), pid(), type(), mass(), charge(), lambda(), drift(), lambdaDeriv(), state(), position(), velocity(), force() {}
    };



    class H5MDReader {

    public:

    	H5MDReader(std::string _file_name):
//    	                          ParticleAccess(system),
//    	                          integrator(_integrator),
    	                          file_name( _file_name)
//    							  author(_author),
//    							  author_email(_author_email),
//    							  creator(_creator),
//    							  creator_version(_creator_version),
//    							  iomode(_iomode),
//    	                          data_to_store(_data_to_store),
//    	                          unfolded(_unfolded),
//    	                          length_factor(_length_factor),
//    							  length_unit(_length_unit),
//    	                          append(_append)
    {

    		  	  //setLengthUnit(_length_unit);
    	          //info_table_ini table_store;
    			  //printf("datas before setting: pid %d\n", datas.pid);
 //   	          set_init_table(&datas, data_to_store);
    	          //printf("datas after setting: pid %d\n", datas.pid);
    	          // we can create the file here, but we need to create with MPI in parallel...


//    	          if (iomode == 1 || iomode == 0) {
//
//    	  			if (system->comm->rank() == 0  && !append) {
//    	  				FileBackup backup(file_name);
//    	  			}
//    	          } else if (iomode == 2) {
//
//    	          	if (!append) {
//    	  				int rank = system->comm->rank();
//    	  				size_t filename_length = file_name.length();
//    	  				std::string suffix = file_name.substr(filename_length-3, 3);
//    	  				std::string base_filename = file_name.substr(0,filename_length-3);
//    	  				std::string rankstring = static_cast<std::ostringstream*>( &(std::ostringstream() << rank) )->str();
//    	  				std::string final_name = base_filename + "_" + rankstring + suffix;
//    	  				FileBackup backup(final_name);
//    	          	}
//
//    	          }
//


    		// open the H5MD file
    		hid_t opened_file = h5md_open_file (file_name.c_str());
    		//opened_file.





    }


      ~H5MDReader() {}


      /*
       * essential helper routines
       *
       */

      void printAuthorName();
      void getAuthorEmail();
      void printCreator();
      void printCreatorVersion();
      void printH5MDVersion();
      void printParameters();
      void printNumberParameters();
      void printNumberParticles();
      void printNumberTimesteps();
      void printTimestepsValues();
      void printTimesteps();
      void printPositionsforGivenStep()

      std::string getFilename(){return file_name;}


//      std::string getAuthor(){return author;}
//      std::string getAuthorEmail(){return author_email;}
//      std::string getCreator(){return creator;}
//      void setCreator(std::string v){creator = v;}
//      std::string getCreatorVersion(){return creator;}
//      void setCreatorVersion(std::string v){creator_version = v;}
//
//      inline int getIomode(){return iomode;}
//      inline void setIomode(int v){iomode = v;}
//      // set and get data_to_store
//
//
//      bool getUnfolded(){return unfolded;}
//      void setUnfolded(bool v){unfolded = v;}
//      bool getAppend(){return append;}
//      void setAppend(bool v){append = v;}

	  std::set<std::string> getSetfrompythonlist(boost::python::list data_to_store) {

		  return python_list_to_set<std::string>(data_to_store);

	  };

//      void set_init_table(info_table_ini* table_init, boost::python::list initial_list){
//    	  //printf("Are we called ?\n");
//       std::set<std::string> ciao = getSetfrompythonlist(initial_list);
//
//        if (ciao.find("all") != ciao.end())
//        {
//        	//printf("Should NOT be inside ALLLLLLL!!!!!\n");
//            table_init->all = 1;
//        } else {
//
//            if (ciao.find("pid") != ciao.end()) {table_init->pid = 1;}
//            if (ciao.find("type") != ciao.end()) {table_init->type = 1;}
//            if (ciao.find("mass") != ciao.end()) {table_init->mass = 1;}
//            if (ciao.find("charge") != ciao.end()) {table_init->charge = 1;}
//            if (ciao.find("lambda") != ciao.end()) {table_init->lambda = 1;}
//            if (ciao.find("drift") != ciao.end()) {table_init->drift = 1;}
//            if (ciao.find("lambdaDeriv") != ciao.end()) {table_init->lambdaDeriv = 1;}
//            if (ciao.find("state") != ciao.end()) {table_init->state = 1;}
//            if (ciao.find("position") != ciao.end()) {table_init->position = 1;}
//            if (ciao.find("velocity") != ciao.end()) {table_init->velocity = 1;}
//            if (ciao.find("force") != ciao.end()) {table_init->force = 1;}
//
//        }
//
//        //return std::set<std::string>
//
//      };


//      boost::python::list getDataToStore() {return data_to_store;}


//      std::string getLengthUnit(){return length_unit;}
//      void setLengthUnit(std::string v){
//        esutil::Error err( getSystem()->comm );
//        if( v != "LJ" && v != "nm" && v != "A" ){
//          std::stringstream msg;
//          msg<<"Wrong unit length: "<< v << "  It should be string: LJ, nm or A" <<"\n";
//          err.setException( msg.str() );
//          err.checkException();
//        }
//
//        length_unit = v;
//      }
//      real getLengthFactor(){return length_factor;}
//      void setLengthFactor(real v){length_factor = v;}

      static void registerPython();

    protected:

      //static LOG4ESPP_DECL_LOGGER(logger);

    private:

      // integrator we need to know an integration step
//      shared_ptr<integrator::MDIntegrator> integrator;

      std::string file_name;

//      std::string author;
//      std::string author_email;
//      std::string creator;
//      std::string creator_version;
//
//      int iomode; // 0: serial, 1: N-to-1, 2: N-to-N; real 0 not there now
//      boost::python::list data_to_store;  // python list: can pass either 'all' and all particle data structure info
//      // is stored or ['pid', 'mass', 'position'] and only pid, mass and position of the particles will be stored
//      // Default is: 'all'
//      info_table_ini datas;
//
//      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
//      bool append; //append to existing trajectory file or create a new one
//      real length_factor;  // for example
//      std::string length_unit; // length unit: {could be LJ, nm, A} it is just for user info
    };
  }
}

#endif

