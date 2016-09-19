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


// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "H5MDFile.hpp"  // keep python.hpp on top
#include "storage/Storage.hpp"
#include "System.hpp"
#include "storage/DomainDecomposition.hpp"
#include "bc/BC.hpp"
#include "analysis/ConfigurationExt.hpp"
#include "analysis/ConfigurationsExt.hpp"
#include "iterator/CellListIterator.hpp"

#include "boost/mpi.hpp"
#include "boost/mpi/communicator.hpp"
#include "mpi.h"

#include <fstream>
#include <sstream>

#ifdef HDF5_LAYER
    #include "hdf5.h"
    #include "hdf5_hl.h"
#endif

#ifdef H5MD_LAYER
    #include "ch5md.h"
#endif

//#include "H5Cpp.h"

using namespace espressopp;
using namespace espressopp::analysis;
using namespace std;


namespace espressopp {
  namespace io {


   void H5MDFile::write_n_to_1(){

    	shared_ptr<System> system = getSystem();
    	char *ch_f_name = new char[file_name.length() + 1];
    	strcpy(ch_f_name, file_name.c_str());
    	int rank = system->comm->rank();
    	int mpi_ranks = system->comm->size();
    	string rankstring = static_cast<ostringstream*>( &(ostringstream() << rank) )->str();

    	int ierr;
		int myN = system->storage->getNRealParticles();  // could avoid this call by putting a counter in the Cells loop below...
		int maxN;   // maximal number of particles one processor has
		int totalN; // total number of particles all processors have

		boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
		boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace

		int* array_nparticles = new int [mpi_ranks];   // to write contiguos in file

		boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguos writing

		char dataSetName[256];
		int RANK = 2;


		typedef struct {
			size_t pid;
			size_t type;
			real mass;
			real charge;
			real lambda; // adress flag: 0 means no adress
			real drift;
			real lambdaDeriv;
			int state;
			double x[3];
			double v[3];
			double force[3];
		} particle_info;

		CellList realCells = system->storage->getRealCells();
		//MPI_Info info  = MPI_INFO_NULL;  // if on normal laptop
        MPI_Info info;
        MPI_Info_create(&info);
        const char* hint_stripe = "striping_unit";
        const char* stripe_value = "4194304";
        MPI_Info_set(info, (char*)hint_stripe, (char*)stripe_value); // 4MB stripe. cast to avoid spurious warnings

        // my detect fs and set appropriate MPI hints!




        // create type for array-like objects, like coordinates, vel and force
		hsize_t dimearr[1] = {3};
		hid_t loctype = H5Tarray_create1(H5T_NATIVE_DOUBLE, 1, dimearr, NULL);

		// create the HDF5 compound datatype from the struct
		hid_t particle_record = H5Tcreate (H5T_COMPOUND, sizeof(particle_info));
		H5Tinsert(particle_record, "pid", HOFFSET(particle_info,pid), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "type", HOFFSET(particle_info,type), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "mass", HOFFSET(particle_info,mass), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "charge", HOFFSET(particle_info,charge), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "lambda", HOFFSET(particle_info,lambda), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "drift", HOFFSET(particle_info,drift), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "lambdaDeriv", HOFFSET(particle_info,lambdaDeriv), H5T_NATIVE_DOUBLE);
		H5Tinsert(particle_record, "state", HOFFSET(particle_info,state), H5T_NATIVE_INT);
		H5Tinsert(particle_record, "x", HOFFSET(particle_info,x), loctype);
		H5Tinsert(particle_record, "v", HOFFSET(particle_info,v), loctype);
		H5Tinsert(particle_record, "force", HOFFSET(particle_info,force), loctype);
        // create the hdf5 data type to store according to user preferences
        //


		hid_t acc_template;
		hid_t file_id, dset_id;
		hid_t filespace, memspace;
		herr_t status;

		hsize_t count[RANK];
		hsize_t offset[RANK];

		hsize_t dimsf[RANK];


		acc_template = H5Pcreate(H5P_FILE_ACCESS);
		H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, info);
		file_id = H5Fcreate(ch_f_name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template); // change this to create the file with H5MD create file call
		//get version;
        //Version();
        //std::string espressopp_version = Version.info();


        // h5md_file h5md_create_file (ch_f_name, author_or_username_on_the_machine, author_email_default_to_null, prog_name, espressopp_version.c_str())
		assert(file_id > 0);

        // h5md_file file;
        // h5md_particles_group particles;
        // particles = h5md_create_particles_group(file, "atoms");
        // // Create a time-dependent dataset
        //     // There is no data yet in "pos"
        // dims[0] = myN;
        // dims[1] = 3;
        // particles.position = h5md_create_time_data(particles.group, "position", 2, dims, H5T_NATIVE_DOUBLE, NULL);

		H5Pclose(acc_template);

		particle_info* particles_u  = new particle_info [myN];

		int i = 0;
		assert( i == 0);

		  if( unfolded ){
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();
			  Int3D& img = cit->image();
			  Real3D L = system->bc->getBoxL();
			  particles_u[i].pid = cit->id();
			  particles_u[i].type = cit->type();
			  particles_u[i].mass = cit->mass();
			  particles_u[i].charge = cit->q();
			  particles_u[i].lambda = cit->lambda();
			  particles_u[i].drift = cit->drift();
			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
			  particles_u[i].state = cit->state();
			  particles_u[i].x[0] = pos[0] + img[0] * L[0];
			  particles_u[i].x[1] = pos[1] + img[1] * L[1];
			  particles_u[i].x[2] = pos[2] + img[2] * L[2];
			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];
			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];

			  i++;
			}
		  }
		  else{
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();
			  particles_u[i].pid = cit->id();
			  particles_u[i].type = cit->type();
			  particles_u[i].mass = cit->mass();
			  particles_u[i].charge = cit->q();
			  particles_u[i].lambda = cit->lambda();
			  particles_u[i].drift = cit->drift();
			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
			  particles_u[i].state = cit->state();
			  particles_u[i].x[0] = pos[0];
			  particles_u[i].x[1] = pos[1];
			  particles_u[i].x[2] = pos[2];
			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];
			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];

			  i++;
			}
		  }

        hsize_t di[RANK];
        di[0] = myN;
        di[1] = 6;

        hsize_t di2[RANK];
		di2[0] = myN;
		di2[1] = H5Tget_size(particle_record);


        dimsf[0] = totalN;
        dimsf[1] = 1;

		filespace = H5Screate_simple(RANK, dimsf, NULL); //create filespace: check h5md call

		sprintf(dataSetName, "Particles");
		dset_id = H5Dcreate2(file_id, dataSetName, particle_record, filespace,
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // check h5md dataset creation func call!

        H5Sclose(filespace);

        // logic to write a file contiguosly!
		int sumup = 0;
		count[0] = myN;
		count[1] = dimsf[1];
		if (rank == 0) {
		offset[0] = rank;
		} else {

			for(int L=0; L<rank; L++) {
				sumup += array_nparticles[L];
			}

			offset[0] = sumup;
		}

		offset[1] = 0;

		memspace = H5Screate_simple(RANK, count, NULL);


		filespace = H5Dget_space(dset_id);
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL); // check with H5MD



		hid_t plist_id = H5Pcreate(H5P_DATASET_XFER); // keep it, h5md doesn't do parallel
		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // important, include in h5md
		//H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // to write datasets independently, usually not suggested!

		status = H5Dwrite(dset_id, particle_record, memspace, filespace, plist_id, particles_u); // check h5md: for parallel we have to pass preference list defined plist_id

		assert( status != -1);

		H5Dclose(dset_id);
		H5Sclose(filespace);
		H5Tclose(particle_record);
		H5Sclose(memspace);

		H5Pclose(plist_id);
		H5Fclose(file_id);
        // h5md_close_file(h5md_file file)
		//MPI_Info_free(&info);

		delete [] ch_f_name;
		delete [] array_nparticles;
		delete [] particles_u;
	}




    void H5MDFile::write_n_to_n(){

    	shared_ptr<System> system = getSystem();
    	int rank = system->comm->rank();

    	size_t filename_length = file_name.length();
    	string suffix = file_name.substr(filename_length-3, 3);
    	string base_filename = file_name.substr(0,filename_length-3);
    	string rankstring = static_cast<ostringstream*>( &(ostringstream() << rank) )->str();
    	std::string final_name = base_filename + "_" + rankstring + suffix;

    	char *ch_f_name = new char[final_name.length() + 1];
    	strcpy(ch_f_name, final_name.c_str());

    	hid_t acc_template;
		hid_t dataspace, memspace, dset, file_id;
		herr_t status;
		int ierr;

		int myN = system->storage->getNRealParticles();  // could avoid this call by putting a counter in the Cells loop below...
		int maxN;   // maximal number of particles one processor has
		int totalN; // total number of particles all processors have

		char dataSetName[256];
		int RANK = 2;

		CellList realCells = system->storage->getRealCells();

		typedef struct {
			size_t pid;
			size_t type;
			real mass;
			real charge;
			real lambda; // adress flag: 0 means no adress
			real drift;
			real lambdaDeriv;
			int state;
			double x[3];
			double v[3];
			double force[3];
		} particle_info;

		  hsize_t     dimearr[1] = {3};
		  hid_t loctype = H5Tarray_create1(H5T_NATIVE_DOUBLE, 1, dimearr, NULL);

		  hid_t particle_record = H5Tcreate (H5T_COMPOUND, sizeof(particle_info));
		  H5Tinsert(particle_record, "pid", HOFFSET(particle_info,pid), H5T_NATIVE_INT);
		  H5Tinsert(particle_record, "type", HOFFSET(particle_info,type), H5T_NATIVE_INT);
		  H5Tinsert(particle_record, "mass", HOFFSET(particle_info,mass), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "charge", HOFFSET(particle_info,charge), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "lambda", HOFFSET(particle_info,lambda), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "drift", HOFFSET(particle_info,drift), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "lambdaDeriv", HOFFSET(particle_info,lambdaDeriv), H5T_NATIVE_DOUBLE);
		  H5Tinsert(particle_record, "state", HOFFSET(particle_info,state), H5T_NATIVE_INT);

		  H5Tinsert(particle_record, "x", HOFFSET(particle_info,x), loctype);
		  H5Tinsert(particle_record, "v", HOFFSET(particle_info,v), loctype);
		  H5Tinsert(particle_record, "force", HOFFSET(particle_info,force), loctype);

		file_id = H5Fcreate(ch_f_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		assert(file_id > 0);

		particle_info* particles_u  = new particle_info [myN];

		int i = 0;
		assert( i == 0);

		  if( unfolded ){
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();
			  Int3D& img = cit->image();
			  Real3D L = system->bc->getBoxL();

			  particles_u[i].pid = cit->id();
			  particles_u[i].type = cit->type();
			  particles_u[i].mass = cit->mass();
			  particles_u[i].charge = cit->q();
			  particles_u[i].lambda = cit->lambda();
			  particles_u[i].drift = cit->drift();
			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
			  particles_u[i].state = cit->state();
			  particles_u[i].x[0] = pos[0] + img[0] * L[0];
			  particles_u[i].x[1] = pos[1] + img[1] * L[1];
			  particles_u[i].x[2] = pos[2] + img[2] * L[2];

			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];

			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];



			  i++;
			}
		  }
		  else{
			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

			  Real3D& pos = cit->position();
			  Real3D& vel = cit->velocity();
			  Real3D& force = cit->force();

			  particles_u[i].pid = cit->id();
			  particles_u[i].type = cit->type();
			  particles_u[i].mass = cit->mass();
			  particles_u[i].charge = cit->q();
			  particles_u[i].lambda = cit->lambda();
			  particles_u[i].drift = cit->drift();
			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
			  particles_u[i].state = cit->state();
			  particles_u[i].x[0] = pos[0];
			  particles_u[i].x[1] = pos[1];
			  particles_u[i].x[2] = pos[2];

			  particles_u[i].v[0] = vel[0];
			  particles_u[i].v[1] = vel[1];
			  particles_u[i].v[2] = vel[2];

			  particles_u[i].force[0] = force[0];
			  particles_u[i].force[1] = force[1];
			  particles_u[i].force[2] = force[2];


			  i++;
			}
		  }

        hsize_t di[RANK];
        di[0] = myN;
        di[1] = 6;

        hsize_t di2[RANK];
		di2[0] = myN;
		di2[1] = H5Tget_size(particle_record);

		hsize_t test[1];
		test[0] = myN;

		dataspace = H5Screate_simple(1, test, NULL);

		sprintf(dataSetName, "Particles");
		dset = H5Dcreate2(file_id, dataSetName, particle_record, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, particle_record, H5S_ALL, H5S_ALL, H5P_DEFAULT, particles_u);
		assert( status != -1);

		hid_t dataspace_id, attribute_id;
		hsize_t     dims;

		dims = 1;
		dataspace_id = H5Screate_simple(1, &dims, NULL);
		double timestep = integrator->getTimeStep();
		long long step = integrator->getStep();
		double attr_data_timestep = timestep;
		long long attr_data_step = step;

        //disabling momentarily
        //attribute_id = H5Acreate(dset, "timestep", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		//status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &attr_data_timestep);
		//attribute_id = H5Acreate(dset, "step", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		//status = H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data_step);
		//H5Aclose(attribute_id);
		H5Sclose(dataspace);
		H5Dclose(dset);

		H5Tclose(particle_record);


		H5Fclose(file_id);

		delete [] ch_f_name;
		delete [] particles_u;
      }



    void H5MDFile::write(){

    	int iomodus = getIomode();
    	if (iomodus == 1 || iomodus == 0) {
    		write_n_to_1();
		} else if (iomodus == 2) {
			write_n_to_n();
		}
    }


    // Python wrapping
    void H5MDFile::registerPython() {

      using namespace espressopp::python;

      class_<H5MDFile, bases<ParticleAccess>, boost::noncopyable >
      ("io_H5MDFile", init< shared_ptr< System >,
                           shared_ptr< integrator::MDIntegrator >,
                           std::string,
						   int,
                           boost::python::list,
                           bool,
                           real,
                           std::string ,
                           bool>())
        .add_property("filename", &H5MDFile::getFilename,
                                  &H5MDFile::setFilename)
		.add_property("iomode", &H5MDFile::getIomode,
										&H5MDFile::setIomode)
        .add_property("unfolded", &H5MDFile::getUnfolded,
                                  &H5MDFile::setUnfolded)
        .add_property("length_factor", &H5MDFile::getLengthFactor,
                                       &H5MDFile::setLengthFactor)
        .add_property("length_unit", &H5MDFile::getLengthUnit,
                                     &H5MDFile::setLengthUnit)
        .add_property("append", &H5MDFile::getAppend,
                                  &H5MDFile::setAppend)
        .def("write", &H5MDFile::write)
      ;
    }
  }
}

