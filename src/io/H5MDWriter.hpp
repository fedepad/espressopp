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
#ifndef _IO_H5MDWRITER_HPP
#define _IO_H5MDWRITER_HPP

#define BOOST_PYTHON_MAX_ARITY 16

// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "ParticleAccess.hpp"  // keep python.hpp on top
#include "System.hpp"
#include "integrator/MDIntegrator.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"
#include "io/FileBackup.hpp"
#include "esutil/Error.hpp"
#include "Version.hpp"
#include "template_py_to_cpp_containers.hpp"
#include "boost/mpi.hpp"
#include "boost/mpi/communicator.hpp"
//#include <string>
//#include <set>

//#include <boost/python/def.hpp>
//#include <boost/python/module.hpp>
//#include <boost/python/args.hpp>
#include "ch5md.h"
#include "esutil/RNG.hpp"


#define ALL              1  //  2^0, bit 0
#define POSITION         2  //  2^1, bit 1
#define VELOCITY         4  //  2^2, bit 2
#define FORCE            8  //  2^3, bit 3
#define PID             16  //  2^4, bit 4
#define TYPE            32  //  2^5, bit 5
#define IMAGE           64  //  2^6, bit 6
#define MASS           128  //  2^7, bit 7
#define CHARGE         256  //  2^8, bit 8
#define LAMBDA         512  //  2^9, bit 9
#define LAMBDADERIV   1024  // 2^10, bit 10
#define DRIFT         2048  // 2^11, bit 11
#define STATE         4096  // 2^12, bit 12


namespace espressopp {
  namespace io{

    class H5MDWriter : public ParticleAccess {

    public:
    	H5MDWriter(const shared_ptr<System>& system,
			const shared_ptr<integrator::MDIntegrator>& _integrator,
			const std::string& _file_name,
			const std::string& _author,
			const std::string& _author_email,
			const std::string& _creator,
			const std::string& _creator_version,
			const int& _iomode,
			const boost::python::list& _data_to_store,
			const bool& _unfolded,
			const real& _length_factor,
			const std::string& _length_unit,
			const bool& _sort_pids,
			const bool& _append,
			const int& _writers):
				  ParticleAccess(system),
				  integrator(_integrator),
				  file_name( _file_name),
				  author(_author),
				  author_email(_author_email),
				  creator(_creator),
				  creator_version(_creator_version),
				  iomode(_iomode),
				  data_to_store(_data_to_store),
				  unfolded(_unfolded),
				  length_factor(_length_factor),
				  length_unit(_length_unit),
				  sort_pids(_sort_pids),
				  append(_append),
				  writers(_writers){

    			  STORE_OPTIONS = 0; // sets to 0 all the bits
    	          set_init_table(STORE_OPTIONS, data_to_store);
    	          // we can create the file here, but we need to create with MPI in parallel...
    	          int mpi_ranks = system->comm->size();
    	          int rank = system->comm->rank();
    	          std::string aggregators_string;
    	          if (writers == 0) {
    	        	aggregators_string = static_cast<std::ostringstream*>( &(std::ostringstream() << mpi_ranks) )->str();
    	          } else {
    	          	aggregators_string = static_cast<std::ostringstream*>( &(std::ostringstream() << writers) )->str();
    	          }
    	          if (iomode == 1 || iomode == 0) {

    	  		if (!append && (rank == 0)) {
    	  			FileBackup backup(file_name);
    	  		}
        	        ilfile = h5md_create_file (file_name.c_str(), author.c_str(), author_email.c_str(), creator.c_str(), creator_version.c_str(), 1, aggregators_string.c_str());
        	        //ilfile = h5md_create_file (file_name.c_str(), author.c_str(), author_email.c_str(), creator.c_str(), creator_version.c_str(), 1);
    	          } else if (iomode == 2) {

    	        	  	size_t filename_length = file_name.length();
						std::string file_extension = file_name.substr(filename_length-3, 3);
						std::string base_filename = file_name.substr(0,filename_length-3);
						std::string rankstring = static_cast<std::ostringstream*>( &(std::ostringstream() << rank) )->str();
						std::string final_name = base_filename + "_" + rankstring + file_extension;

    	          	if (!append) {
    	          		FileBackup backup(final_name);
    	          	}

    	          		ilfile = h5md_create_file (final_name.c_str(), author.c_str(), author_email.c_str(), creator.c_str(), creator_version.c_str(), 0, aggregators_string.c_str());
			//ilfile = h5md_create_file (final_name.c_str(), author.c_str(), author_email.c_str(), creator.c_str(), creator_version.c_str(), 0);

    	          }


    	          int ierr;
    	          int myN = system->storage->getNRealParticles();  // could avoid this call by putting a counter in the Cells loop below...
    	          	//int maxN;   // maximal number of particles one processor has
    	          int totalN; // total number of particles all processors have

    	          	//boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
    	          boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace

    	          	//int* array_nparticles = new int [mpi_ranks];   // to write contiguous in file

    	          	//boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguous writing

    	          	//ilfile = h5md_create_file (file_name.c_str(), author.c_str(), author_email.c_str(), creator.c_str(), creator_version.c_str(), 1);

    	          const char *boundary[] = {"periodic", "periodic", "periodic"};

    	          Real3D L = system->bc->getBoxL();
    	          double box_edges[3];
    	          box_edges[0] = L[0];
    	          box_edges[1] = L[1];
    	          box_edges[2] = L[2];
    	          long long step = integrator->getStep();
    	          real timestep_unit = integrator->getTimeStep();  // can change over time...now fixed...
    	          //real time_ = step * timestep_unit;

    	          real maxcut = system->maxCutoff;

    	          part_group = h5md_create_particles_group(ilfile, "atoms");

    	          double the_skin = system->getSkin();
    	          long the_seed = system->rng->get_seed();

				  hid_t aid2  = H5Screate(H5S_SCALAR);
				  int attr2 = H5Acreate2(ilfile.parameters, "max_cutoff", H5T_NATIVE_DOUBLE, aid2,
											  H5P_DEFAULT, H5P_DEFAULT);
				  int ret = H5Awrite(attr2, H5T_NATIVE_DOUBLE, &maxcut);
				  ret = H5Sclose(aid2);
				  ret = H5Aclose(attr2);

				  aid2  = H5Screate(H5S_SCALAR);
				  attr2 = H5Acreate2(ilfile.parameters, "skin", H5T_NATIVE_DOUBLE, aid2,
									H5P_DEFAULT, H5P_DEFAULT);
				  ret = H5Awrite(attr2, H5T_NATIVE_DOUBLE, &the_skin);
				  ret = H5Sclose(aid2);
				  ret = H5Aclose(attr2);
				  aid2  = H5Screate(H5S_SCALAR);
				  attr2 = H5Acreate2(ilfile.parameters, "rng_seed", H5T_NATIVE_LONG, aid2,
						  H5P_DEFAULT, H5P_DEFAULT);
				  ret = H5Awrite(attr2, H5T_NATIVE_LONG, &the_seed);
				  ret = H5Sclose(aid2);
				  ret = H5Aclose(attr2);
				  aid2  = H5Screate(H5S_SCALAR);
				  attr2 = H5Acreate2(ilfile.parameters, "length_factor", H5T_NATIVE_DOUBLE, aid2,
					H5P_DEFAULT, H5P_DEFAULT);
				  ret = H5Awrite(attr2, H5T_NATIVE_DOUBLE, &length_factor);

				  ret = H5Sclose(aid2);
				  ret = H5Aclose(attr2);

				  h5md_write_string_attribute(ilfile.id, "parameters", "length_unit", length_unit.c_str());
    	          h5md_create_box(&part_group, 3, boundary, false, box_edges, NULL);

    	          int RANK = 2;
    	          int dims[RANK];
    	          dims[0] = totalN;
    	          dims[1] = 3;
    	          int dims_1d[RANK];
    	          dims_1d[0] = totalN;
    	          dims_1d[1] = 1;

    	          if (iomode == 2) {
    	        	dims[0] = myN;
    	        	dims_1d[0] = myN;
    	          }

    	    // mind the types: "errors" happen if one sets H5T_NATIVE_INT as datatype for pids or types and size_t for the array when writing!
    	    part_group.id = h5md_create_time_data(part_group.group, "id", RANK, dims_1d, H5T_NATIVE_INT, NULL);
			if (STORE_OPTIONS & (IMAGE | POSITION | ALL) )  part_group.position = h5md_create_time_data(part_group.group, "position", RANK, dims, H5T_NATIVE_DOUBLE, NULL);
			if (STORE_OPTIONS & (IMAGE | ALL)) 	  			part_group.image = h5md_create_time_data(part_group.group, "image", RANK, dims, H5T_NATIVE_INT, NULL);
			//if (STORE_OPTIONS & (PID | ALL))      	  		part_group.id = h5md_create_time_data(part_group.group, "id", RANK, dims_1d, H5T_NATIVE_INT, NULL);
			if (STORE_OPTIONS & (TYPE | ALL))     	  		part_group.species = h5md_create_time_data(part_group.group, "species", RANK, dims_1d, H5T_NATIVE_INT, NULL);
			if (STORE_OPTIONS & (VELOCITY | ALL)) 	  		part_group.velocity = h5md_create_time_data(part_group.group, "velocity", RANK, dims, H5T_NATIVE_DOUBLE, NULL);
			if (STORE_OPTIONS & (FORCE | ALL))    	  		part_group.force = h5md_create_time_data(part_group.group, "force", RANK, dims, H5T_NATIVE_DOUBLE, NULL);
			if (STORE_OPTIONS & (MASS | ALL))     	  		part_group.mass = h5md_create_time_data(part_group.group, "mass", RANK, dims_1d, H5T_NATIVE_DOUBLE, NULL);
			if (STORE_OPTIONS & (CHARGE | ALL))   	  		part_group.charge = h5md_create_time_data(part_group.group, "charge", RANK, dims_1d, H5T_NATIVE_DOUBLE, NULL);
			if (STORE_OPTIONS & (STATE | ALL))       		part_group.state = h5md_create_time_data(part_group.group, "state", RANK, dims_1d, H5T_NATIVE_INT, NULL);
			if (STORE_OPTIONS & (DRIFT | ALL))       		part_group.drift = h5md_create_time_data(part_group.group, "drift", RANK, dims_1d, H5T_NATIVE_DOUBLE, NULL);
			if (STORE_OPTIONS & (LAMBDA | ALL))      		part_group.lambda = h5md_create_time_data(part_group.group, "lambda", RANK, dims_1d, H5T_NATIVE_DOUBLE, NULL);
			if (STORE_OPTIONS & (LAMBDADERIV | ALL)) 		part_group.lambdaDeriv = h5md_create_time_data(part_group.group, "lambdaDeriv", RANK, dims_1d, H5T_NATIVE_DOUBLE, NULL);


    	}


      ~H5MDWriter() {}

      void perform_action(){
        write();
      }

      /*
       * essential helper routines
       *
       */
      //void open();
      void sort_by_pid();
      void flush_file_stable_storage();
      void close();
      /*
       * write helper routines
       *
       */
      void write_n_to_1();
      void write_n_to_n();
      void write();

      std::string getFilename(){return file_name;}
      void setFilename(std::string v){file_name = v;}


      std::string getAuthor(){return author;}
      void setAuthor(std::string v){author = v;}
      std::string getAuthorEmail(){return author_email;}
      void setAuthorEmail(std::string v){author_email = v;}
      std::string getCreator(){return creator;}
      void setCreator(std::string v){creator = v;}
      std::string getCreatorVersion(){return creator_version;}
      void setCreatorVersion(std::string v){creator_version = v;}

      inline int getIomode(){return iomode;}
      inline void setIomode(int v){iomode = v;}
      // set and get data_to_store


      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}
      bool getSortPids(){return sort_pids;}
      void setSortPids(bool v){sort_pids = v;}
      bool getAppend(){return append;}
      void setAppend(bool v){append = v;}

      int getAggregators(){return writers;}
      void setAggregators(int v){writers = v;}

      std::set<std::string> getSetfrompythonlist(const boost::python::list& data_to_store) {

	  return python_list_to_set<std::string>(data_to_store);

	};

/*
      void set_init_table(uint16_t& options, const boost::python::list& initial_list){

       std::set<std::string> ciao = getSetfrompythonlist(initial_list);

        if (ciao.find("all") != ciao.end())
        {
            options |= ALL;
        } else {

            options |= PID;
            if (ciao.find("type") != ciao.end()) {options |= TYPE;}
            if (ciao.find("mass") != ciao.end()) {options |= MASS;}
            if (ciao.find("charge") != ciao.end()) {options |= CHARGE;}
            if (ciao.find("lambda") != ciao.end()) {options |= LAMBDA;}
            if (ciao.find("drift") != ciao.end()) {options |= DRIFT;}
            if (ciao.find("lambdaDeriv") != ciao.end()) {options |= LAMBDADERIV;}
            if (ciao.find("state") != ciao.end()) {options |= STATE;}
            if ((ciao.find("image") != ciao.end()) || (ciao.find("position") != ciao.end())) {options |= POSITION;} // if image is desired then position must be there according to the H5MD paper standard
            if (ciao.find("velocity") != ciao.end()) {options |= VELOCITY;}
            if (ciao.find("force") != ciao.end()) {options |= FORCE;}
            if (ciao.find("image") != ciao.end()) {options |= IMAGE;}

        }

      };
*/

      void set_init_table(uint16_t& options, const boost::python::list& initial_list){

		std::string current;
		for (int i = 0; i < boost::python::len(initial_list); ++i) {

			current = boost::python::extract<std::string>(initial_list[i]);
			if (current.compare("all") == 0) {
				options |= ALL;
				break;
			}
			// no check on pid, is always saved!

			else if (current.compare("type") == 0) { options |= TYPE;}
			else if (current.compare("mass") == 0) { options |= MASS;}
			else if (current.compare("charge") == 0) { options |= CHARGE;}
			else if (current.compare("lambda") == 0) { options |= LAMBDA;}
			else if (current.compare("drift") == 0) { options |= DRIFT;}
			else if (current.compare("lambdaDeriv") == 0) { options |= LAMBDADERIV;}
			else if (current.compare("state") == 0) { options |= STATE;}
			else if (current.compare("position") == 0) { options |= POSITION;}
			else if (current.compare("velocity") == 0) { options |= VELOCITY;}
			else if (current.compare("force") == 0) { options |= FORCE;}
			else if (current.compare("image") == 0) { options |= IMAGE;}
			options |= PID;
		}
      };

      boost::python::list getDataToStore() {return data_to_store;}


      std::string getLengthUnit(){return length_unit;}
      void setLengthUnit(std::string v){
        esutil::Error err( getSystem()->comm );
        if( v != "LJ" && v != "nm" && v != "A" ){
          std::stringstream msg;
          msg<<"Wrong unit length: "<< v << "  It should be string: LJ, nm or A" <<"\n";
          err.setException( msg.str() );
          err.checkException();
        }

        length_unit = v;
      }
      real getLengthFactor(){return length_factor;}
      void setLengthFactor(real v){length_factor = v;}

      static void registerPython();

    protected:

      //static LOG4ESPP_DECL_LOGGER(logger);

    private:

      shared_ptr<integrator::MDIntegrator> integrator;
      std::string file_name;
      std::string author;
      std::string author_email;
      std::string creator;
      std::string creator_version;
      int iomode; // 0: serial, 1: N-to-1, 2: N-to-N; real 0 not there now
      boost::python::list data_to_store;  // python list: can pass either 'all' and all particle data structure info
      // is stored or ['pid', 'mass', 'position'] and only pid, mass and position of the particles will be stored: pid always stored
      // Default is: 'all'
      //info_table_ini datas;
      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
      bool append; //append to existing trajectory file or create a new one
      real length_factor;  // for example
      std::string length_unit; // length unit: {could be LJ, nm, A} it is just for user info
      bool sort_pids;
      int writers;
      uint16_t STORE_OPTIONS;
      h5md_file ilfile;
      h5md_particles_group part_group;
    };
  }
}

#endif

