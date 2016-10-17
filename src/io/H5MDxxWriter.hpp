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

// ESPP_CLASS
#ifndef _IO_H5MDXXWRITER_HPP
#define _IO_H5MDXXWRITER_HPP

// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "ParticleAccess.hpp"  // keep python.hpp on top
#include "System.hpp"
#include "integrator/MDIntegrator.hpp"
#include "storage/Storage.hpp"
#include "io/FileBackup.hpp"
#include "esutil/Error.hpp"
#include "Version.hpp"
#include "template_py_to_cpp_containers.hpp"
#include "h5xx/h5xx.hpp"
#include <string>
#include <set>

#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/args.hpp>

typedef boost::array<int, 2> array_type_int_size_2;
typedef boost::multi_array<int, 2> array_2d_t;
typedef boost::multi_array<int, 1> array_1d_t;
typedef boost::multi_array<double, 2> array_2d_dbl_t;
typedef boost::array<std::string, 11> array_type_string_size_11;

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



    class H5MDxxWriter : public ParticleAccess {

    public:

    	H5MDxxWriter(shared_ptr<System> system,
    	                shared_ptr<integrator::MDIntegrator> _integrator,
    	                std::string _file_name,
    					std::string _author,
    					std::string _author_email,
    					std::string _creator,
    					std::string _creator_version,
						int _iomode,
    	                boost::python::list _data_to_store,
    	                bool _unfolded,
    	                real _length_factor,
    	                std::string _length_unit,
    	                bool _append):
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
    	                          append(_append) {

    		  	  //setLengthUnit(_length_unit);
    	          //info_table_ini table_store;
    			  //printf("datas before setting: pid %d\n", datas.pid);
    	          set_init_table(&datas, data_to_store);
    	          //printf("datas after setting: pid %d\n", datas.pid);
    	          // we can create the file here, but we need to create with MPI in parallel...


    	          if (iomode == 1 || iomode == 0) {

    	  			if (system->comm->rank() == 0  && !append) {
    	  				FileBackup backup(file_name);
    	  			}
    	          } else if (iomode == 2) {

    	          	if (!append) {
    	  				int rank = system->comm->rank();
    	  				size_t filename_length = file_name.length();
    	  				std::string suffix = file_name.substr(filename_length-3, 3);
    	  				std::string base_filename = file_name.substr(0,filename_length-3);
    	  				std::string rankstring = static_cast<std::ostringstream*>( &(std::ostringstream() << rank) )->str();
    	  				std::string final_name = base_filename + "_" + rankstring + suffix;
    	  				FileBackup backup(final_name);
    	          	}

    	          }


    	      	int rank = system->comm->rank();
    	      	int mpi_ranks = system->comm->size();
    	      	//string rankstring = static_cast<ostringstream*>( &(ostringstream() << rank) )->str();

    	      	int ierr;
    	      	int myN = system->storage->getNRealParticles();  // could avoid this call by putting a counter in the Cells loop below...
    	      	int maxN;   // maximal number of particles one processor has
    	      	int totalN; // total number of particles all processors have

    	      	boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
    	      	boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace

    	      	int* array_nparticles = new int [mpi_ranks];   // to write contiguous in file

    	      	boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguous writing


			 // will create a file and the particles group
				MPI_Info info;
				MPI_Info_create(&info);
				const char* hint_stripe = "striping_unit";
				const char* stripe_value = "4194304";
				MPI_Info_set(info, (char*)hint_stripe, (char*)stripe_value); // 4MB stripe. cast to avoid spurious warnings
				//add my detection of fs...don't know where it went!!!!!
				 //h5xx::file fileh = file(filename, MPI_COMM_WORLD, info, h5xx::file::out); // this should append ....?
				hdf5_file = h5xx::file(file_name, MPI_COMM_WORLD, info, h5xx::file::trunc); // default

				h5xx::group h5md_group = h5xx::group(hdf5_file, "h5md");

				array_type_int_size_2 version_attr = {{1, 0}};
				h5xx::write_attribute(h5md_group, "version", version_attr);

				h5xx::group author_group = h5xx::group(h5md_group, "author");
				h5xx::write_attribute(author_group, "name", author);

				if (!author_email.empty()) {
				  h5xx::write_attribute(author_group, "author_email", author_email);
				}

				h5xx::group creator_group = h5xx::group(h5md_group, "creator");
				h5xx::write_attribute(creator_group, "name", creator);
				h5xx::write_attribute(creator_group, "version", creator_version);

				h5xx::group particles_group = h5xx::group(h5md_group, "particles");
				h5xx::group observables_group = h5xx::group(h5md_group, "observables");
				h5xx::group parameters_group = h5xx::group(h5md_group, "parameters");
				// create a group inside the "Particles group" --> this will hold all the subgroups pids, positions, etc.

				atoms_group = h5xx::group(particles_group, "atoms");
				// here ends the creation of the file and of the particles group

				// can create the particles related datasets here ???!!!

				std::vector<size_t> dataspace_2d_dims;
				dataspace_2d_dims[0] = totalN;
				dataspace_2d_dims[1] = 3;
				std::vector<size_t> dataspace_1d_dims;
				dataspace_1d_dims[0] = totalN;
				dataspace_1d_dims[0] = 1;

				std::vector<size_t> chunk_dims;
				chunk_dims[0] = 128;
				chunk_dims[1] = 128;


				h5xx::dataspace space_1d(dataspace_1d_dims);
				h5xx::dataspace space_2d(dataspace_2d_dims);

				array_type_string_size_11 datasets_names {{"positions", "pids", "types", "velocities", "forces", "masses", "charges", "states", "drifts", "lambdas", "lambdaDerivs"}};

				//boost::array<size_t, myN> pids;
				//boost::array<size_t, myN> types;
				array_2d_dbl_t positions(boost::extents[myN][3]);
				array_2d_dbl_t velocities(boost::extents[myN][3]);
				array_2d_dbl_t forces(boost::extents[myN][3]);
				//boost::array<double, myN> charges;
				//boost::array<double, myN> masses;
				//boost::array<int, myN> states;
				//boost::array<double, myN> drifts;
				//boost::array<double, myN> lambdas;
				//boost::array<double, myN> lambdaDerivs;
				std::vector<size_t> pids = std::vector<size_t>(myN);
				std::vector<size_t> types = std::vector<size_t>(myN);
				std::vector<double> charges = std::vector<double>(myN);
				std::vector<double> masses = std::vector<double>(myN);
				std::vector<int> states = std::vector<int>(myN);
				std::vector<double> drifts = std::vector<double>(myN);
				std::vector<double> lambdas = std::vector<double>(myN);
				std::vector<double> lambdaDerivs = std::vector<double>(myN);

				if (datas.position == 1 || datas.all == 1) 	 	h5xx::create_dataset(atoms_group, datasets_names[0], array_2d_dbl_t, space_2d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.pid == 1 || datas.all == 1)      	 	h5xx::create_dataset(atoms_group, datasets_names[1], std::vector<size_t>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.type == 1 || datas.all == 1)     	 	h5xx::create_dataset(atoms_group, datasets_names[2], std::vector<size_t>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.velocity == 1 || datas.all == 1) 	 	h5xx::create_dataset(atoms_group, datasets_names[3], array_2d_dbl_t, space_2d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.force == 1 || datas.all == 1)    	 	h5xx::create_dataset(atoms_group, datasets_names[4], array_2d_dbl_t, space_2d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.mass == 1 || datas.all == 1)     	 	h5xx::create_dataset(atoms_group, datasets_names[5], std::vector<double>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.charge == 1 || datas.all == 1)   	 	h5xx::create_dataset(atoms_group, datasets_names[6], std::vector<double>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.state == 1 || datas.all == 1)      	h5xx::create_dataset(atoms_group, datasets_names[7], std::vector<int>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.drift == 1 || datas.all == 1)      	h5xx::create_dataset(atoms_group, datasets_names[8], std::vector<double>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.lambda == 1 || datas.all == 1)     	h5xx::create_dataset(atoms_group, datasets_names[9], std::vector<double>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
				if (datas.lambdaDeriv == 1 || datas.all == 1)   h5xx::create_dataset(atoms_group, datasets_names[10], std::vector<double>, space_1d, h5xx::policy::storage::chunked(chunk_dims));






    	        }


      ~H5MDxxWriter() {}

      void perform_action(){
        write();
      }

      /*
       * essential helper routines
       *
       */
      void open();
      void sort_by_pid();
      void flush_file_stable_storage();
      void close_file();
      void create_box();
      void create_time_depend_datasets();
      void write_time_dependent_datasets();

      void create_time_independent_datasets();
      void write_time_indipendent_datasets();






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
      std::string getCreatorVersion(){return creator;}
      void setCreatorVersion(std::string v){creator_version = v;}

      inline int getIomode(){return iomode;}
      inline void setIomode(int v){iomode = v;}
      // set and get data_to_store


      bool getUnfolded(){return unfolded;}
      void setUnfolded(bool v){unfolded = v;}
      bool getAppend(){return append;}
      void setAppend(bool v){append = v;}

	  std::set<std::string> getSetfrompythonlist(boost::python::list data_to_store) {

		  return python_list_to_set<std::string>(data_to_store);

	  };

      void set_init_table(info_table_ini* table_init, boost::python::list initial_list){
    	  //printf("Are we called ?\n");
       std::set<std::string> ciao = getSetfrompythonlist(initial_list);

        if (ciao.find("all") != ciao.end())
        {
        	//printf("Should NOT be inside ALLLLLLL!!!!!\n");
            table_init->all = 1;
        } else {

            if (ciao.find("pid") != ciao.end()) {table_init->pid = 1;}
            if (ciao.find("type") != ciao.end()) {table_init->type = 1;}
            if (ciao.find("mass") != ciao.end()) {table_init->mass = 1;}
            if (ciao.find("charge") != ciao.end()) {table_init->charge = 1;}
            if (ciao.find("lambda") != ciao.end()) {table_init->lambda = 1;}
            if (ciao.find("drift") != ciao.end()) {table_init->drift = 1;}
            if (ciao.find("lambdaDeriv") != ciao.end()) {table_init->lambdaDeriv = 1;}
            if (ciao.find("state") != ciao.end()) {table_init->state = 1;}
            if (ciao.find("position") != ciao.end()) {table_init->position = 1;}
            if (ciao.find("velocity") != ciao.end()) {table_init->velocity = 1;}
            if (ciao.find("force") != ciao.end()) {table_init->force = 1;}

        }

        //return std::set<std::string>

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

      // integrator we need to know an integration step
      shared_ptr<integrator::MDIntegrator> integrator;

      std::string file_name;


      std::string author;
      std::string author_email;
      std::string creator;
      std::string creator_version;

      int iomode; // 0: serial, 1: N-to-1, 2: N-to-N; real 0 not there now
      boost::python::list data_to_store;  // python list: can pass either 'all' and all particle data structure info
      // is stored or ['pid', 'mass', 'position'] and only pid, mass and position of the particles will be stored
      // Default is: 'all'
      info_table_ini datas;

      bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
      bool append; //append to existing trajectory file or create a new one
      real length_factor;  // for example
      std::string length_unit; // length unit: {could be LJ, nm, A} it is just for user info
      h5xx::file hdf5_file;
      h5xx::group atoms_group;
    };
  }
}

#endif

