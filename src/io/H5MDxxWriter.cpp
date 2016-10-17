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
#include "H5MDxxWriter.hpp"  // keep python.hpp on top
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
    #include "ch5md.h"
//#endif

using namespace espressopp;
using namespace std;
using namespace boost::python;


namespace espressopp {
  namespace io {


  void H5MDxxWriter::open() {}
  void H5MDxxWriter::sort_by_pid() {}
  void H5MDxxWriter::close_file() {

	  // move here closing of all the elements and the file
	  //H5Gclose(atoms.group);
	  //h5md_close_file(the_File);

  }
  void H5MDxxWriter::flush_file_stable_storage() {

	  // move here flush to stable storage to avoid corruptions
	  //H5Fflush()


  }


   void H5MDxxWriter::write_n_to_1(){

	   // move the creating file part in the constructor!!!!
	   // Separate the close...i.e. do not close the file here! It's a write method, not a close!

	shared_ptr<System> system = getSystem();
	//long int the_seed = system->rng->seed_;

	//char *ch_f_name = new char[file_name.length() + 1];
	//strcpy(ch_f_name, file_name.c_str());
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

	h5md_file the_File = h5md_create_file (file_name.c_str(), author.c_str(), author_email.c_str(), creator.c_str(), creator_version.c_str(), 1);

	const char *boundary[] = {"periodic", "periodic", "periodic"};

	Real3D L = system->bc->getBoxL();
	double box_edges[3];
	box_edges[0] = L[0];
	box_edges[1] = L[1];
	box_edges[2] = L[2];
	long long step = integrator->getStep();
	real time_ = step * integrator->getTimeStep();
	std::string stepstring = static_cast<std::ostringstream*>(&(std::ostringstream() << step))->str();
	std::string fin = "step_" + stepstring;

	h5md_particles_group atoms = h5md_create_particles_group(the_File, "Atoms");



	//hid_t atoms = H5Gcreate(atomsi.group, fin.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//h5md_particles_group atoms = h5md_create_particles_group()
	//h5md_create_parameters_group(the_File, "parameters");
	//system->getSkin();
	//system->storage->getInt3DCellGrid();
	//int number_interactions = system->getNumberOfInteractions();

	//real maximum_cutoff = system->maxCutoff;

	h5md_create_box(&atoms, 3, boundary, false, box_edges, NULL);

	//double* coordin = new double [myN];


	int RANK = 2;
	int dims[RANK];
	dims[0] = totalN;
	dims[1] = 3;

//	std::vector<size_t> dataspace_2d_dims;
//	dataspace_2d_dims[0] = totalN;
//	dataspace_2d_dims[1] = 3;
//	std::vector<size_t> dataspace_1d_dims;
//	dataspace_1d_dims[0] = totalN;
//	dataspace_1d_dims[0] = 1;
//
//	std::vector<size_t> chunk_dims;
//	chunk_dims[0] = 128;
//	chunk_dims[1] = 128;


//	h5xx::dataspace space_1d(dataspace_1d_dims);
//	h5xx::dataspace space_2d(dataspace_2d_dims);
//
//	array_type_string_size_8 datasets_names {{"positions", "pids", "types", "velocities", "forces", "masses", "charges", "states", "drifts", "lambdas", "lambdaDerivs"}};
//
//	boost::array<size_t, myN> pids;
//	boost::array<size_t, myN> types;
	array_2d_dbl_t positions(boost::extents[myN][3]);
	array_2d_dbl_t velocities(boost::extents[myN][3]);
	array_2d_dbl_t forces(boost::extents[myN][3]);
//	boost::array<double, myN> charges;
//	boost::array<double, myN> masses;
//	boost::array<int, myN> states;
//	boost::array<double, myN> drifts;
//	boost::array<double, myN> lambdas;
//	boost::array<double, myN> lambdaDerivs;





	std::vector<size_t> pids = std::vector<size_t>(myN);
	std::vector<size_t> types = std::vector<size_t>(myN);
	std::vector<double> charges = std::vector<double>(myN);
	std::vector<double> masses = std::vector<double>(myN);
	std::vector<int> states = std::vector<int>(myN);
	std::vector<double> drifts = std::vector<double>(myN);
	std::vector<double> lambdas = std::vector<double>(myN);
	std::vector<double> lambdaDerivs = std::vector<double>(myN);

//
//	if (datas.position == 1 || datas.all == 1) 	 	h5xx::create_dataset(ciao, datasets_names[0], array_2d_dbl_t, space_2d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.pid == 1 || datas.all == 1)      	 	h5xx::create_dataset(ciao, datasets_names[1], boost::array<size_t, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.type == 1 || datas.all == 1)     	 	h5xx::create_dataset(ciao, datasets_names[2], boost::array<size_t, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.velocity == 1 || datas.all == 1) 	 	h5xx::create_dataset(ciao, datasets_names[3], array_2d_dbl_t, space_2d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.force == 1 || datas.all == 1)    	 	h5xx::create_dataset(ciao, datasets_names[4], array_2d_dbl_t, space_2d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.mass == 1 || datas.all == 1)     	 	h5xx::create_dataset(ciao, datasets_names[5], boost::array<double, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.charge == 1 || datas.all == 1)   	 	h5xx::create_dataset(ciao, datasets_names[6], boost::array<double, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.state == 1 || datas.all == 1)      	h5xx::create_dataset(ciao, datasets_names[7], boost::array<int, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.drift == 1 || datas.all == 1)      	h5xx::create_dataset(ciao, datasets_names[8], boost::array<double, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.lambda == 1 || datas.all == 1)     	h5xx::create_dataset(ciao, datasets_names[9], boost::array<double, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//	if (datas.lambdaDeriv == 1 || datas.all == 1)   h5xx::create_dataset(ciao, datasets_names[10], boost::array<double, myN>, space_1d, h5xx::policy::storage::chunked(chunk_dims));
//

	CellList realCells = system->storage->getRealCells();

	int i = 0;
	assert( i == 0);

	  if( unfolded ){
		for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

		  if (datas.pid == 1 || datas.all == 1) {
			  pids[i] = cit->id();

		  }
		  if (datas.position == 1 || datas.all == 1) {
		  Real3D& pos = cit->position();
		  Int3D& img = cit->image();
		  positions[i][0] = pos[0] + img[0] * L[0];
		  positions[i][1] = pos[1] + img[1] * L[1];
		  positions[i][2] = pos[2] + img[2] * L[2];
		  }
		  if (datas.velocity == 1 || datas.all == 1) {
			  Real3D& vel = cit->velocity();
			  velocities[i][0] = vel[0];
			  velocities[i][1] = vel[1];
			  velocities[i][2] = vel[2];

		  }

		  if (datas.type == 1 || datas.all == 1) {

			  types[i] = cit->type();
		  }
		  if (datas.force == 1 || datas.all == 1) {
			  Real3D& force = cit->force();
			  forces[i][0] = force[0];
			  forces[i][1] = force[1];
			  forces[i][2] = force[2];
		  }
		  if (datas.mass == 1 || datas.all == 1) {

			  masses[i] = cit->mass();
		  }
		  if (datas.state == 1) {

			  states[i] = cit->state();
		  }

		  if (datas.drift == 1 || datas.all == 1) {

			  drifts[i] = cit->drift();
		  }

		  if (datas.charge == 1 || datas.all == 1) {

			  charges[i] = cit->q();
		  }

		  if (datas.lambda == 1 || datas.all == 1) {

			  lambdas[i] = cit->lambda();
		  }

		  if (datas.lambdaDeriv == 1 || datas.all == 1) {

			  lambdaDerivs[i] = cit->lambdaDeriv();
		  }

		  i++;
		}
	  }
	  else{
		for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

		  if (datas.pid == 1 || datas.all == 1) {
			  pids[i] = cit->id();

		  }
		  if (datas.position == 1 || datas.all == 1) {
		  Real3D& pos = cit->position();
		  positions[i][0] = pos[0];
		  positions[i][1] = pos[1];
		  positions[i][2] = pos[2];
		  }
		  if (datas.velocity == 1 || datas.all == 1) {
			  Real3D& vel = cit->velocity();
			  velocities[i][0] = vel[0];
			  velocities[i][1] = vel[1];
			  velocities[i][2] = vel[2];

		  }

		  if (datas.type == 1 || datas.all == 1) {

			  types[i] = cit->type();
		  }
		  if (datas.force == 1 || datas.all == 1) {
			  Real3D& force = cit->force();
			  forces[i][0] = force[0];
			  forces[i][1] = force[1];
			  forces[i][2] = force[2];
		  }
		  if (datas.mass == 1 || datas.all == 1) {

			  masses[i] = cit->mass();
		  }
		  if (datas.state == 1 || datas.all == 1) {

			  states[i] = cit->state();
		  }

		  if (datas.drift == 1 || datas.all == 1) {

			  drifts[i] = cit->drift();
		  }

		  if (datas.charge == 1 || datas.all == 1) {

			  charges[i] = cit->q();
		  }

		  if (datas.lambda == 1 || datas.all == 1) {

			  lambdas[i] = cit->lambda();
		  }

		  if (datas.lambdaDeriv == 1 || datas.all == 1) {

			  lambdaDerivs[i] = cit->lambdaDeriv();
		  }

		  i++;
		}
	  }

//

//	  // Create a time-dependent dataset
//	  	// There is no data yet in "pos"
//	  	int RANK = 2;
////	  	//hsize_t dims[RANK];
//	  	int dims[RANK];
//
	  	hsize_t count[RANK];
		hsize_t offset[RANK];
		hsize_t dimsf[RANK];
		//        dimsf[0] = totalN;
		//        dimsf[1] = 1;




	  	// As done in #xx
	  	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	  	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

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


		// --- define slice and write dataset
		std::vector<size_t> slice_offset;
		slice_offset.push_back(0);  // stack matrix blocks
		slice_offset.push_back(sumup);
		std::vector<size_t> slice_count;
		slice_count.push_back(1);
		slice_count.push_back(myN);
		h5xx::slice data_slice(slice_offset, slice_count);
		//
		if (datas.position == 1  || datas.all == 1)     h5xx::write_dataset(atoms_group, "positions", positions, data_slice);
		if (datas.pid == 1  || datas.all == 1) h5xx::write_dataset(atoms_group, "pids", pids, data_slice);


//	if (datas.position == 1  || datas.all == 1)   h5md_append(atoms.position, positions, step, time_, plist_id, offset[0], myN);
//	if (datas.pid == 1  || datas.all == 1)        h5md_append(atoms.id, pids, step, time_, plist_id, offset[0], myN);
//	if (datas.type == 1 || datas.all == 1)     	  h5md_append(atoms.species, types, step, time_, plist_id, offset[0], myN);
//	if (datas.velocity == 1 || datas.all == 1) 	  h5md_append(atoms.velocity, velocities, step, time_, plist_id, offset[0], myN);
//	if (datas.mass == 1 || datas.all == 1)     	  h5md_append(atoms.mass, masses, step, time_, plist_id, offset[0], myN);
//	if (datas.charge == 1 || datas.all == 1)   	  h5md_append(atoms.charge, charges, step, time_, plist_id, offset[0], myN);
//	if (datas.state == 1 || datas.all == 1)    	  h5md_append(atoms.state, states, step, time_, plist_id, offset[0], myN);
//	if (datas.drift == 1 || datas.all == 1)       h5md_append(atoms.drift, drifts, step, time_, plist_id, offset[0], myN);
//	if (datas.force == 1 || datas.all == 1)       h5md_append(atoms.force, forces, step, time_, plist_id, offset[0], myN);
//	if (datas.lambda == 1 || datas.all == 1)      h5md_append(atoms.lambda, lambdas, step, time_, plist_id, offset[0], myN);
//	if (datas.lambdaDeriv == 1 || datas.all == 1) h5md_append(atoms.lambdaDeriv, lambdaDerivs, step, time_, plist_id, offset[0], myN);

  	H5Pclose(plist_id);

//	if (datas.position == 1 || datas.all == 1) 	  h5md_close_element(atoms.position);
//	if (datas.pid == 1 || datas.all == 1) 		  h5md_close_element(atoms.id);
//	if (datas.type == 1 || datas.all == 1) 		  h5md_close_element(atoms.species);
//	if (datas.velocity == 1 || datas.all == 1) 	  h5md_close_element(atoms.velocity);
//	if (datas.force == 1 || datas.all == 1) 	  h5md_close_element(atoms.force);
//	if (datas.mass == 1 || datas.all == 1) 		  h5md_close_element(atoms.mass);
//	if (datas.charge == 1 || datas.all == 1) 	  h5md_close_element(atoms.charge);
//	if (datas.state == 1 || datas.all == 1) 	  h5md_close_element(atoms.state);
//	if (datas.drift == 1 || datas.all == 1) 	  h5md_close_element(atoms.drift);
//	if (datas.lambda == 1 || datas.all == 1)      h5md_close_element(atoms.lambda);
//	if (datas.lambdaDeriv == 1 || datas.all == 1) h5md_close_element(atoms.lambdaDeriv);

//	if (!append) {
//		H5Gclose(atoms.group);
//		h5md_close_file(the_File);
//	}

	delete [] array_nparticles;

	}




    void H5MDxxWriter::write_n_to_n(){

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



    void H5MDxxWriter::write(){

    	int iomodus = getIomode();
    	if (iomodus == 1 || iomodus == 0) {
    		write_n_to_1();
		} else if (iomodus == 2) {
			write_n_to_n();
		}
    }



    void create_box() {



    }
          void create_time_depend_datasets(h5xx::group group_name, std::string dataset_name, int chunk_dimension ) {
        	  // create chunked dataset
        	  //h5xx::create_dataset(group);



          }
          void write_time_dependent_datasets() {



          }

          void create_time_independent_datasets() {



          }
          void write_time_indipendent_datasets() {



          }

    // Python wrapping
    void H5MDxxWriter::registerPython() {

      using namespace espressopp::python;

      class_<H5MDxxWriter, bases<ParticleAccess>, boost::noncopyable >
      ("io_H5MDxxWriter", init< shared_ptr< System >,
                           shared_ptr< integrator::MDIntegrator >,
                           std::string,
						   std::string,
					       std::string,
						   std::string,
						   std::string,
						   int,
                           boost::python::list,
                           bool,
                           real,
                           std::string ,
                           bool>())
        .add_property("filename", &H5MDxxWriter::getFilename,
                                  &H5MDxxWriter::setFilename)
	    .add_property("author", &H5MDxxWriter::getFilename,
								  &H5MDxxWriter::setFilename)
	    .add_property("author_email", &H5MDxxWriter::getFilename,
								  &H5MDxxWriter::setFilename)
	    .add_property("creator", &H5MDxxWriter::getFilename,
								  &H5MDxxWriter::setFilename)
		.add_property("creator_version", &H5MDxxWriter::getFilename,
								  &H5MDxxWriter::setFilename)
	    .add_property("iomode", &H5MDxxWriter::getIomode,
								&H5MDxxWriter::setIomode)
		.add_property("data_to_store", &H5MDxxWriter::getDataToStore,
										&H5MDxxWriter::set_init_table)

        .add_property("unfolded", &H5MDxxWriter::getUnfolded,
                                  &H5MDxxWriter::setUnfolded)
        .add_property("length_factor", &H5MDxxWriter::getLengthFactor,
                                       &H5MDxxWriter::setLengthFactor)
        .add_property("length_unit", &H5MDxxWriter::getLengthUnit,
                                     &H5MDxxWriter::setLengthUnit)
        .add_property("append", &H5MDxxWriter::getAppend,
                                  &H5MDxxWriter::setAppend)
        .def("write", &H5MDxxWriter::write)
      ;
    }
  }
}

