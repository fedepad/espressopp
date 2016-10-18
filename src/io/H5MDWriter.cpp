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
#include "H5MDWriter.hpp"  // keep python.hpp on top
//#include "storage/Storage.hpp"
//#include "System.hpp"
#include "storage/DomainDecomposition.hpp"
//#include "bc/BC.hpp"

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


  void H5MDWriter::open() {}
  void H5MDWriter::sort_by_pid() {}
  void H5MDWriter::close() {
	 // check if we need to sort the pids...
	  if (sort_pids) sort_by_pid();


	if (datas.position == 1 || datas.all == 1) 	  h5md_close_element(part_group.position);
	if (datas.image == 1 || datas.all == 1) 	  h5md_close_element(part_group.image);
	if (datas.pid == 1 || datas.all == 1) 		  h5md_close_element(part_group.id);
	if (datas.type == 1 || datas.all == 1) 		  h5md_close_element(part_group.species);
	if (datas.velocity == 1 || datas.all == 1) 	  h5md_close_element(part_group.velocity);
	if (datas.force == 1 || datas.all == 1) 	  h5md_close_element(part_group.force);
	if (datas.mass == 1 || datas.all == 1) 		  h5md_close_element(part_group.mass);
	if (datas.charge == 1 || datas.all == 1) 	  h5md_close_element(part_group.charge);
	if (datas.state == 1 || datas.all == 1) 	  h5md_close_element(part_group.state);
	if (datas.drift == 1 || datas.all == 1) 	  h5md_close_element(part_group.drift);
	if (datas.lambda == 1 || datas.all == 1)      h5md_close_element(part_group.lambda);
	if (datas.lambdaDeriv == 1 || datas.all == 1) h5md_close_element(part_group.lambdaDeriv);
	h5md_close_element(part_group.box_edges);
	H5Gclose(part_group.group);
    //h5md_close_element(part_group.);
	h5md_close_file(ilfile);

  }

  void H5MDWriter::flush_file_stable_storage() {

	  // move here flush to stable storage to avoid corruptions
	  //H5Fflush()


  }


   void H5MDWriter::write_n_to_1(){

	   // move the creating file part in the constructor!!!!
	   // Separate the close...i.e. do not close the file here! It's a write method, not a close!

	shared_ptr<System> system = getSystem();
//	//long int the_seed = system->rng->seed_;
//
//	char *ch_f_name = new char[file_name.length() + 1];
//	strcpy(ch_f_name, file_name.c_str());
	int rank = system->comm->rank();
	int mpi_ranks = system->comm->size();
//	string rankstring = static_cast<ostringstream*>( &(ostringstream() << rank) )->str();
//
//	int ierr;
	int myN = system->storage->getNRealParticles();
	int maxN;   // maximal number of particles one processor has
	int totalN; // total number of particles all processors have
//
	boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
	boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace
//
	int* array_nparticles = new int [mpi_ranks];   // to write contiguous in file
//
	boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguous writing
//
//	h5md_file the_File = h5md_create_file (file_name.c_str(), author.c_str(), author_email.c_str(), creator.c_str(), creator_version.c_str(), 1);
//
//	const char *boundary[] = {"periodic", "periodic", "periodic"};
//
	Real3D L = system->bc->getBoxL();
//	double box_edges[3];
//	box_edges[0] = L[0];
//	box_edges[1] = L[1];
//	box_edges[2] = L[2];
	long long step = integrator->getStep();
	real time_ = step * integrator->getTimeStep();
//	std::string stepstring = static_cast<std::ostringstream*>(&(std::ostringstream() << step))->str();
//	std::string fin = "step_" + stepstring;
//
//	h5md_particles_group atoms = h5md_create_particles_group(the_File, "Atoms");
//
//
//
//	//hid_t atoms = H5Gcreate(atomsi.group, fin.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//	//h5md_particles_group atoms = h5md_create_particles_group()
//	//h5md_create_parameters_group(the_File, "parameters");
//	//system->getSkin();
//	//system->storage->getInt3DCellGrid();
//	//int number_interactions = system->getNumberOfInteractions();
//
//	//real maximum_cutoff = system->maxCutoff;
//
//	h5md_create_box(&atoms, 3, boundary, false, box_edges, NULL);
//
//	//double* coordin = new double [myN];
//
//
////		if (datas.pid == 1  || datas.all == 1)        size_t pids[myN];
////		if (datas.position == 1 || datas.all == 1)    double coordina[myN][3];
////		if (datas.velocity == 1 || datas.all == 1)    double velocities[myN][3];
////		if (datas.type == 1 || datas.all == 1)        size_t types[myN];
////		if (datas.force == 1 || datas.all == 1)       double forces[myN][3];
////		if (datas.charge == 1 || datas.all == 1)      double charges[myN];
////		if (datas.mass == 1 || datas.all == 1)        double masses[myN];
////		if (datas.state == 1 || datas.all == 1)       int states[myN];
////		if (datas.drift == 1 || datas.all == 1)       double drifts[myN];
////		if (datas.lambda == 1 || datas.all == 1)      double lambdas[myN];
////		if (datas.lambdaDeriv == 1 || datas.all == 1) double lambdaDerivs[myN];

	int RANK = 2;
//	int dims[RANK];
//	dims[0] = totalN;
//	dims[1] = 3;
//
//
//
//
//	if (datas.position == 1 || datas.all == 1) 	  atoms.position = h5md_create_time_data(atoms.group, "positions", RANK, dims, H5T_NATIVE_DOUBLE, NULL);
//	if (datas.pid == 1 || datas.all == 1)      	  atoms.id = h5md_create_time_data(atoms.group, "pids", 1, &totalN, H5T_NATIVE_INT, NULL);
//	if (datas.type == 1 || datas.all == 1)     	  atoms.species = h5md_create_time_data(atoms.group, "types", 1, &totalN, H5T_NATIVE_INT, NULL);
//	if (datas.velocity == 1 || datas.all == 1) 	  atoms.velocity = h5md_create_time_data(atoms.group, "velocities", RANK, dims, H5T_NATIVE_DOUBLE, NULL);
//	if (datas.force == 1 || datas.all == 1)    	  atoms.force = h5md_create_time_data(atoms.group, "forces", RANK, dims, H5T_NATIVE_DOUBLE, NULL);
//	if (datas.mass == 1 || datas.all == 1)     	  atoms.mass = h5md_create_time_data(atoms.group, "masses", 1, &totalN, H5T_NATIVE_DOUBLE, NULL);
//	if (datas.charge == 1 || datas.all == 1)   	  atoms.charge = h5md_create_time_data(atoms.group, "charges", 1, &totalN, H5T_NATIVE_DOUBLE, NULL);
//	if (datas.state == 1 || datas.all == 1)       atoms.state = h5md_create_time_data(atoms.group, "states", 1, &totalN, H5T_NATIVE_INT, NULL);
//	if (datas.drift == 1 || datas.all == 1)       atoms.drift = h5md_create_time_data(atoms.group, "drifts", 1, &totalN, H5T_NATIVE_DOUBLE, NULL);
//	if (datas.lambda == 1 || datas.all == 1)      atoms.lambda = h5md_create_time_data(atoms.group, "lambdas", 1, &totalN, H5T_NATIVE_DOUBLE, NULL);
//	if (datas.lambdaDeriv == 1 || datas.all == 1) atoms.lambdaDeriv = h5md_create_time_data(atoms.group, "lambdaDerivs", 1, &totalN, H5T_NATIVE_DOUBLE, NULL);

//	size_t pids[myN];
//	double coordina[myN][3];
//	double velocities[myN][3];
//	double forces[myN][3];
//	double charges[myN];
//	double masses[myN];
//	int states[myN];
//	double drifts[myN];
//	double lambdas[myN];
//	double lambdaDerivs[myN];
//	size_t types[myN];


	size_t* pids = new size_t [myN];
	size_t* types = new size_t [myN];
	double* coordina = new double [myN*3];
	int* images = new int [myN*3];
	double* velocities = new double [myN*3];
	double* forces = new double [myN*3];
	double* charges = new double [myN];
	double* masses = new double [myN];
	double* drifts = new double [myN];
	double* lambdas = new double [myN];
	double* lambdaDerivs = new double [myN];
	int* states = new int [myN];




	CellList realCells = system->storage->getRealCells();

	int i = 0;
	assert( i == 0);

	  if( unfolded ){
		for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

		  //if (datas.pid == 1 || datas.all == 1) { // always pids
			  pids[i] = cit->id();

		  //}
		  if (datas.position == 1 || datas.all == 1) {
			  Real3D& pos = cit->position();
			  Int3D& img = cit->image();
	//		  coordina[i][0] = pos[0] + img[0] * L[0];
	//		  coordina[i][1] = pos[1] + img[1] * L[1];
	//		  coordina[i][2] = pos[2] + img[2] * L[2];
			  coordina[i*3]   = pos[0] + img[0] * L[0];
			  coordina[i*3+1] = pos[1] + img[1] * L[1];
			  coordina[i*3+2] = pos[2] + img[2] * L[2];
		  }
		  //if (datas.image == 1 || datas.all == 1) {
		  		  Int3D& img = cit->image();
		  		  images[i*3]   = img[0];
		  		  images[i*3+1] = img[1];
		  		  images[i*3+2] = img[2];
		  //		  }
		  if (datas.velocity == 1 || datas.all == 1) {
			  Real3D& vel = cit->velocity();
//			  velocities[i][0] = vel[0];
//			  velocities[i][1] = vel[1];
//			  velocities[i][2] = vel[2];

			  velocities[i*3] = vel[0];
			  velocities[i*3+1] = vel[1];
			  velocities[i*3+2] = vel[2];
		  }

		  if (datas.type == 1 || datas.all == 1) {

			  types[i] = cit->type();
		  }
		  if (datas.force == 1 || datas.all == 1) {
			  Real3D& force = cit->force();
//			  forces[i][0] = force[0];
//			  forces[i][1] = force[1];
//			  forces[i][2] = force[2];
			  forces[i*3] = force[0];
			  forces[i*3+1] = force[1];
			  forces[i*3+2] = force[2];
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

		  //if (datas.pid == 1 || datas.all == 1) { //always pids
			  pids[i] = cit->id();

		  //}
		  if (datas.position == 1 || datas.all == 1) {
		  Real3D& pos = cit->position();
//		  coordina[i][0] = pos[0];
//		  coordina[i][1] = pos[1];
//		  coordina[i][2] = pos[2];
		  coordina[i*3] = pos[0];
		  coordina[i*3+1] = pos[1];
		  coordina[i*3+2] = pos[2];
		  }
		  if (datas.velocity == 1 || datas.all == 1) {
			  Real3D& vel = cit->velocity();
//			  velocities[i][0] = vel[0];
//			  velocities[i][1] = vel[1];
//			  velocities[i][2] = vel[2];
				velocities[i*3] = vel[0];
				velocities[i*3+1] = vel[1];
				velocities[i*3+2] = vel[2];

		  }

		  if (datas.type == 1 || datas.all == 1) {

			  types[i] = cit->type();
		  }
		  if (datas.force == 1 || datas.all == 1) {
			  Real3D& force = cit->force();
//			  forces[i][0] = force[0];
//			  forces[i][1] = force[1];
//			  forces[i][2] = force[2];
			  forces[i*3] = force[0];
			  forces[i*3+1] = force[1];
			  forces[i*3+2] = force[2];
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

        // logic to write a file contiguously!
		int sumup = 0;
		//count[0] = myN;
		//count[1] = dimsf[1];
		if (rank == 0) {
			offset[1] = 0;
		} else {

			for(int L=0; L<rank; L++) {
				sumup += array_nparticles[L];
			}

			offset[1] = sumup;
		}

		offset[0] = 0;

	if (datas.position == 1  || datas.all == 1)   h5md_append(part_group.position, coordina, step, time_, plist_id, offset[1], myN);
	if (unfolded && (datas.image == 1  || datas.all == 1))   h5md_append(part_group.position, images, step, time_, plist_id, offset[1], myN);
	if (datas.pid == 1  || datas.all == 1)        h5md_append(part_group.id, pids, step, time_, plist_id, offset[1], myN);
	if (datas.type == 1 || datas.all == 1)     	  h5md_append(part_group.species, types, step, time_, plist_id, offset[1], myN);
	if (datas.velocity == 1 || datas.all == 1) 	  h5md_append(part_group.velocity, velocities, step, time_, plist_id, offset[1], myN);
	if (datas.mass == 1 || datas.all == 1)     	  h5md_append(part_group.mass, masses, step, time_, plist_id, offset[1], myN);
	if (datas.charge == 1 || datas.all == 1)   	  h5md_append(part_group.charge, charges, step, time_, plist_id, offset[1], myN);
	if (datas.state == 1 || datas.all == 1)    	  h5md_append(part_group.state, states, step, time_, plist_id, offset[1], myN);
	if (datas.drift == 1 || datas.all == 1)       h5md_append(part_group.drift, drifts, step, time_, plist_id, offset[1], myN);
	if (datas.force == 1 || datas.all == 1)       h5md_append(part_group.force, forces, step, time_, plist_id, offset[1], myN);
	if (datas.lambda == 1 || datas.all == 1)      h5md_append(part_group.lambda, lambdas, step, time_, plist_id, offset[1], myN);
	if (datas.lambdaDeriv == 1 || datas.all == 1) h5md_append(part_group.lambdaDeriv, lambdaDerivs, step, time_, plist_id, offset[1], myN);

  	H5Pclose(plist_id);

//	if (datas.position == 1 || datas.all == 1) 	  h5md_close_element(part_group.position);
//	if (datas.pid == 1 || datas.all == 1) 		  h5md_close_element(part_group.id);
//	if (datas.type == 1 || datas.all == 1) 		  h5md_close_element(part_group.species);
//	if (datas.velocity == 1 || datas.all == 1) 	  h5md_close_element(part_group.velocity);
//	if (datas.force == 1 || datas.all == 1) 	  h5md_close_element(part_group.force);
//	if (datas.mass == 1 || datas.all == 1) 		  h5md_close_element(part_group.mass);
//	if (datas.charge == 1 || datas.all == 1) 	  h5md_close_element(part_group.charge);
//	if (datas.state == 1 || datas.all == 1) 	  h5md_close_element(part_group.state);
//	if (datas.drift == 1 || datas.all == 1) 	  h5md_close_element(part_group.drift);
//	if (datas.lambda == 1 || datas.all == 1)      h5md_close_element(part_group.lambda);
//	if (datas.lambdaDeriv == 1 || datas.all == 1) h5md_close_element(part_group.lambdaDeriv);

//	if (!append) {
//		H5Gclose(atoms.group);
//		h5md_close_file(the_File);
//	}

	delete [] array_nparticles;
	delete [] pids;
	delete [] types;
	delete [] coordina;
	delete [] images;
	delete [] velocities;
	delete [] forces;
	delete [] charges;
	delete [] masses;
	delete [] states;
	delete [] lambdas;
	delete [] lambdaDerivs;
	delete [] drifts;



//    	shared_ptr<System> system = getSystem();
//    	char *ch_f_name = new char[file_name.length() + 1];
//    	strcpy(ch_f_name, file_name.c_str());
//    	int rank = system->comm->rank();
//    	int mpi_ranks = system->comm->size();
//    	string rankstring = static_cast<ostringstream*>( &(ostringstream() << rank) )->str();
//
//    	int ierr;
//		int myN = system->storage->getNRealParticles();  // could avoid this call by putting a counter in the Cells loop below...
//		int maxN;   // maximal number of particles one processor has
//		int totalN; // total number of particles all processors have
//
//		boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
//		boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace
//
//		int* array_nparticles = new int [mpi_ranks];   // to write contiguous in file
//
//		boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguous writing
//
//		char dataSetName[256];
//		int RANK = 2;
//
//
//		typedef struct {
//			size_t pid;
//			size_t type;
//			real mass;
//			real charge;
//			real lambda; // adress flag: 0 means no adress
//			real drift;
//			real lambdaDeriv;
//			int state;
//			double x[3];
//			double v[3];
//			double force[3];
//		} particle_info;
//
//		CellList realCells = system->storage->getRealCells();
//		//MPI_Info info  = MPI_INFO_NULL;  // if on normal laptop
//        MPI_Info info;
//        MPI_Info_create(&info);
//        const char* hint_stripe = "striping_unit";
//        const char* stripe_value = "4194304";
//        MPI_Info_set(info, (char*)hint_stripe, (char*)stripe_value); // 4MB stripe. cast to avoid spurious warnings
//
//        // my detect fs and set appropriate MPI hints! Can move this logic into ch5md
//
//
//
//
//        // create type for array-like objects, like coordinates, vel and force
//		hsize_t dimearr[1] = {3};
//		hid_t loctype = H5Tarray_create1(H5T_NATIVE_DOUBLE, 1, dimearr, NULL);
//
//		// create the HDF5 compound datatype from the struct
//		hid_t particle_record = H5Tcreate (H5T_COMPOUND, sizeof(particle_info));
//		H5Tinsert(particle_record, "pid", HOFFSET(particle_info,pid), H5T_NATIVE_INT);
//		H5Tinsert(particle_record, "type", HOFFSET(particle_info,type), H5T_NATIVE_INT);
//		H5Tinsert(particle_record, "mass", HOFFSET(particle_info,mass), H5T_NATIVE_DOUBLE);
//		H5Tinsert(particle_record, "charge", HOFFSET(particle_info,charge), H5T_NATIVE_DOUBLE);
//		H5Tinsert(particle_record, "lambda", HOFFSET(particle_info,lambda), H5T_NATIVE_DOUBLE);
//		H5Tinsert(particle_record, "drift", HOFFSET(particle_info,drift), H5T_NATIVE_DOUBLE);
//		H5Tinsert(particle_record, "lambdaDeriv", HOFFSET(particle_info,lambdaDeriv), H5T_NATIVE_DOUBLE);
//		H5Tinsert(particle_record, "state", HOFFSET(particle_info,state), H5T_NATIVE_INT);
//		H5Tinsert(particle_record, "x", HOFFSET(particle_info,x), loctype);
//		H5Tinsert(particle_record, "v", HOFFSET(particle_info,v), loctype);
//		H5Tinsert(particle_record, "force", HOFFSET(particle_info,force), loctype);
//        // create the hdf5 data type to store according to user preferences
//        //
//
//
//		hid_t acc_template;
//		hid_t file_id, dset_id;
//		hid_t filespace, memspace;
//		herr_t status;
//
//		hsize_t count[RANK];
//		hsize_t offset[RANK];
//
//		hsize_t dimsf[RANK];
//
//
//		acc_template = H5Pcreate(H5P_FILE_ACCESS);
//		H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, info);
//		file_id = H5Fcreate(ch_f_name, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template); // change this to create the file with H5MD create file call
//		//h5md_file h5md_create_file (const char *filename, const char *author, const char *author_email, const char *creator, const char *creator_version, int parallel);
//		//h5md_file h5md_create_file (file_name.c_str(), author.c_str(), NULL, creator.c_str(), creator_version.c_str(), 1);
//
//		//get version;
//        //Version();
//        //std::string espressopp_version = Version.info();
//
//
//        // h5md_file h5md_create_file (ch_f_name, author_or_username_on_the_machine, author_email_default_to_null, prog_name, espressopp_version.c_str())
//		assert(file_id > 0);
//
//        // h5md_file file;
//        // h5md_particles_group particles;
//        // particles = h5md_create_particles_group(file, "atoms");
//        // // Create a time-dependent dataset
//        //     // There is no data yet in "pos"
//        // dims[0] = myN;
//        // dims[1] = 3;
//        // particles.position = h5md_create_time_data(particles.group, "position", 2, dims, H5T_NATIVE_DOUBLE, NULL);
//
//		//H5Pclose(acc_template);
//
//		particle_info* particles_u  = new particle_info [myN];
//
//		int i = 0;
//		assert( i == 0);
//
//		  if( unfolded ){
//			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
//			  Real3D& pos = cit->position();
//			  Real3D& vel = cit->velocity();
//			  Real3D& force = cit->force();
//			  Int3D& img = cit->image();
//			  Real3D L = system->bc->getBoxL();
//			  particles_u[i].pid = cit->id();
//			  particles_u[i].type = cit->type();
//			  particles_u[i].mass = cit->mass();
//			  particles_u[i].charge = cit->q();
//			  particles_u[i].lambda = cit->lambda();
//			  particles_u[i].drift = cit->drift();
//			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
//			  particles_u[i].state = cit->state();
//			  particles_u[i].x[0] = pos[0] + img[0] * L[0];
//			  particles_u[i].x[1] = pos[1] + img[1] * L[1];
//			  particles_u[i].x[2] = pos[2] + img[2] * L[2];
//			  particles_u[i].v[0] = vel[0];
//			  particles_u[i].v[1] = vel[1];
//			  particles_u[i].v[2] = vel[2];
//			  particles_u[i].force[0] = force[0];
//			  particles_u[i].force[1] = force[1];
//			  particles_u[i].force[2] = force[2];
//
//			  i++;
//			}
//		  }
//		  else{
//			for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {
//			  Real3D& pos = cit->position();
//			  Real3D& vel = cit->velocity();
//			  Real3D& force = cit->force();
//			  particles_u[i].pid = cit->id();
//			  particles_u[i].type = cit->type();
//			  particles_u[i].mass = cit->mass();
//			  particles_u[i].charge = cit->q();
//			  particles_u[i].lambda = cit->lambda();
//			  particles_u[i].drift = cit->drift();
//			  particles_u[i].lambdaDeriv = cit->lambdaDeriv();
//			  particles_u[i].state = cit->state();
//			  particles_u[i].x[0] = pos[0];
//			  particles_u[i].x[1] = pos[1];
//			  particles_u[i].x[2] = pos[2];
//			  particles_u[i].v[0] = vel[0];
//			  particles_u[i].v[1] = vel[1];
//			  particles_u[i].v[2] = vel[2];
//			  particles_u[i].force[0] = force[0];
//			  particles_u[i].force[1] = force[1];
//			  particles_u[i].force[2] = force[2];
//
//			  i++;
//			}
//		  }
//
//        hsize_t di[RANK];
//        di[0] = myN;
//        di[1] = 6;
//
//        hsize_t di2[RANK];
//		di2[0] = myN;
//		di2[1] = H5Tget_size(particle_record);
//
//
//        dimsf[0] = totalN;
//        dimsf[1] = 1;
//
//		filespace = H5Screate_simple(RANK, dimsf, NULL); //create filespace: check h5md call
//
//		sprintf(dataSetName, "Particles");
//		dset_id = H5Dcreate2(file_id, dataSetName, particle_record, filespace,
//                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); // check h5md dataset creation func call!
//
//        H5Sclose(filespace);
//
//        // logic to write a file contiguosly!
//		int sumup = 0;
//		count[0] = myN;
//		count[1] = dimsf[1];
//		if (rank == 0) {
//		offset[0] = rank;
//		} else {
//
//			for(int L=0; L<rank; L++) {
//				sumup += array_nparticles[L];
//			}
//
//			offset[0] = sumup;
//		}
//
//		offset[1] = 0;
//
//		memspace = H5Screate_simple(RANK, count, NULL);
//
//
//		filespace = H5Dget_space(dset_id);
//		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL); // check with H5MD
//
//
//
//		hid_t plist_id = H5Pcreate(H5P_DATASET_XFER); // keep it, h5md doesn't do parallel
//		H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); // important, include in h5md
//		//H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);  // to write datasets independently, usually not suggested!
//
//		status = H5Dwrite(dset_id, particle_record, memspace, filespace, plist_id, particles_u); // check h5md: for parallel we have to pass preference list defined plist_id
//
//		assert( status != -1);
//
//		H5Dclose(dset_id);
//		H5Sclose(filespace);
//		H5Tclose(particle_record);
//		H5Sclose(memspace);
//
//		H5Pclose(plist_id);
//		H5Fclose(file_id);
//        // h5md_close_file(h5md_file file)
//		//MPI_Info_free(&info);
//
//		delete [] ch_f_name;
//		delete [] array_nparticles;
//		delete [] particles_u;
	}




    void H5MDWriter::write_n_to_n(){

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



    void H5MDWriter::write(){

    	int iomodus = getIomode();
    	if (iomodus == 1 || iomodus == 0) {
    		write_n_to_1();
		} else if (iomodus == 2) {
			write_n_to_n();
		}
    }


    // Python wrapping
    void H5MDWriter::registerPython() {

      using namespace espressopp::python;

      class_<H5MDWriter, bases<ParticleAccess>, boost::noncopyable >
      ("io_H5MDWriter", init< shared_ptr< System >,
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
						   bool,
                           bool>())
        .add_property("filename", &H5MDWriter::getFilename,
                                  &H5MDWriter::setFilename)
	    .add_property("author", &H5MDWriter::getFilename,
								  &H5MDWriter::setFilename)
	    .add_property("author_email", &H5MDWriter::getFilename,
								  &H5MDWriter::setFilename)
	    .add_property("creator", &H5MDWriter::getFilename,
								  &H5MDWriter::setFilename)
		.add_property("creator_version", &H5MDWriter::getFilename,
								  &H5MDWriter::setFilename)
	    .add_property("iomode", &H5MDWriter::getIomode,
								&H5MDWriter::setIomode)
		.add_property("data_to_store", &H5MDWriter::getDataToStore,
										&H5MDWriter::set_init_table)

        .add_property("unfolded", &H5MDWriter::getUnfolded,
                                  &H5MDWriter::setUnfolded)
        .add_property("length_factor", &H5MDWriter::getLengthFactor,
                                       &H5MDWriter::setLengthFactor)
        .add_property("length_unit", &H5MDWriter::getLengthUnit,
                                     &H5MDWriter::setLengthUnit)
		 .add_property("sort_pids", &H5MDWriter::getSortPids,
									  &H5MDWriter::setSortPids)
        .add_property("append", &H5MDWriter::getAppend,
                                  &H5MDWriter::setAppend)
        .def("dump", &H5MDWriter::write)
		.def("close", &H5MDWriter::close)
      ;
    }
  }
}

