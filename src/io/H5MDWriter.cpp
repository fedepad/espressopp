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
#include "storage/DomainDecomposition.hpp"
#include "iterator/CellListIterator.hpp"

#include <fstream>
#include <sstream>

using namespace espressopp;
using namespace std;
using namespace boost::python;


namespace espressopp {
  namespace io {


  void H5MDWriter::sort_by_pid() {}
  void H5MDWriter::close() {
	 // check if we need to sort the pids...
  if (sort_pids) sort_by_pid();

  h5md_close_element(part_group.id);
	if (STORE_OPTIONS & (IMAGE | POSITION | ALL)) 	  h5md_close_element(part_group.position);
	if (STORE_OPTIONS & (IMAGE | ALL)) 	  h5md_close_element(part_group.image);
	//if (STORE_OPTIONS & (PID | ALL)) 		  h5md_close_element(part_group.id);
	if (STORE_OPTIONS & (TYPE | ALL)) 		  h5md_close_element(part_group.species);
	if (STORE_OPTIONS & (VELOCITY | ALL)) 	  h5md_close_element(part_group.velocity);
	if (STORE_OPTIONS & (FORCE | ALL)) 	  h5md_close_element(part_group.force);
	if (STORE_OPTIONS & (MASS | ALL)) 		  h5md_close_element(part_group.mass);
	if (STORE_OPTIONS & (CHARGE | ALL)) 	  h5md_close_element(part_group.charge);
	if (STORE_OPTIONS & (STATE | ALL)) 	  h5md_close_element(part_group.state);
	if (STORE_OPTIONS & (DRIFT | ALL)) 	  h5md_close_element(part_group.drift);
	if (STORE_OPTIONS & (LAMBDA | ALL))      h5md_close_element(part_group.lambda);
	if (STORE_OPTIONS & (LAMBDADERIV | ALL)) h5md_close_element(part_group.lambdaDeriv);
	h5md_close_element(part_group.box_edges);
	H5Gclose(part_group.group);
	h5md_close_file(ilfile);

  }

  void H5MDWriter::flush_file_stable_storage() {

	  // move here flush to stable storage to avoid corruptions
	  //H5Fflush()


  }


   void H5MDWriter::write_n_to_1(){

	shared_ptr<System> system = getSystem();

	int rank = system->comm->rank();
	int mpi_ranks = system->comm->size();
	int myN = system->storage->getNRealParticles();
	//int maxN;   // maximal number of particles one processor has
	//int totalN; // total number of particles all processors have

	//boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
	//boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace: flaw if totalN changes, should move the dataspace creation here, in the write routine

	int* array_nparticles = new int [mpi_ranks];   // to write contiguous in file

	boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguous writing: every mpi rank knows how many particles the other ranks have

	Real3D L = system->bc->getBoxL();
	long long step = integrator->getStep();
	real time_ = step * integrator->getTimeStep();

	int RANK = 2;

	int* pids = new int [myN];
	int* types = new int [myN];
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
	//assert( i == 0);
	if (unfolded) {
		for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

			Real3D& pos = cit->position();
			Int3D& img = cit->image();
			Real3D& vel = cit->velocity();
			Real3D& force = cit->force();
			pids[i] = cit->id();
			types[i] = cit->type();
			masses[i] = cit->mass();
			states[i] = cit->state();
			drifts[i] = cit->drift();
			charges[i] = cit->q();
			lambdas[i] = cit->lambda();
			velocities[i*3]   = vel[0];
			velocities[i*3+1] = vel[1];
			velocities[i*3+2] = vel[2];
			forces[i*3]   = force[0];
			forces[i*3+1] = force[1];
			forces[i*3+2] = force[2];
			coordina[i*3]   = pos[0] + img[0] * L[0];
			coordina[i*3+1] = pos[1] + img[1] * L[1];
			coordina[i*3+2] = pos[2] + img[2] * L[2];
			images[i*3]   = img[0];
			images[i*3+1] = img[1];
			images[i*3+2] = img[2];

			++i;
	  }
	  //i++;
	  //++i;
	} else {
		for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

			Real3D& pos = cit->position();
			Int3D& img = cit->image();
			Real3D& vel = cit->velocity();
			Real3D& force = cit->force();
			pids[i] = cit->id();
			types[i] = cit->type();
			masses[i] = cit->mass();
			states[i] = cit->state();
			drifts[i] = cit->drift();
			charges[i] = cit->q();
			lambdas[i] = cit->lambda();
			velocities[i*3]   = vel[0];
			velocities[i*3+1] = vel[1];
			velocities[i*3+2] = vel[2];
			forces[i*3]   = force[0];
			forces[i*3+1] = force[1];
			forces[i*3+2] = force[2];
			coordina[i*3]   = pos[0];
			coordina[i*3+1] = pos[1];
			coordina[i*3+2] = pos[2];
			images[i*3]   = img[0];
			images[i*3+1] = img[1];
			images[i*3+2] = img[2];

			++i;
		}
		//++i;
	}


	  	hsize_t count[RANK]; // not used, directly pass myN
		hsize_t offset[RANK];

	  	// As done in #xxPR HDF5File
	  	//hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	  	//H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // logic to write a file contiguously!
		int sumup = 0;
		if (rank == 0) {
			offset[1] = 0;
		} else {

			for(int L=0; L<rank; L++) {
				sumup += array_nparticles[L];
			}

			offset[1] = sumup;
		}

		offset[0] = 0;

	h5md_append(part_group.id, pids, step, time_, offset[1], myN, 1); // always write pids
	if (STORE_OPTIONS & (IMAGE | POSITION | ALL))   	 h5md_append(part_group.position, coordina, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (IMAGE | ALL))   				 h5md_append(part_group.image, images, step, time_, offset[1], myN, 1);
	//if (STORE_OPTIONS & (PID | ALL))        			 h5md_append(part_group.id, pids, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (TYPE | ALL))     	  			 h5md_append(part_group.species, types, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (VELOCITY | ALL)) 	  			 h5md_append(part_group.velocity, velocities, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (MASS | ALL))     	  			 h5md_append(part_group.mass, masses, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (CHARGE | ALL))   	 			 h5md_append(part_group.charge, charges, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (STATE | ALL))    	  			 h5md_append(part_group.state, states, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (DRIFT | ALL))       			 h5md_append(part_group.drift, drifts, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (FORCE | ALL))       			 h5md_append(part_group.force, forces, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (LAMBDA | ALL))      			 h5md_append(part_group.lambda, lambdas, step, time_, offset[1], myN, 1);
	if (STORE_OPTIONS & (LAMBDADERIV | ALL)) 			 h5md_append(part_group.lambdaDeriv, lambdaDerivs, step, time_, offset[1], myN, 1);

  	//H5Pclose(plist_id);

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

	}




    void H5MDWriter::write_n_to_n(){

    	shared_ptr<System> system = getSystem();

    	int rank = system->comm->rank();
    	int mpi_ranks = system->comm->size();
    	int myN = system->storage->getNRealParticles();
    	//int maxN;   // maximal number of particles one processor has
    	//int totalN; // total number of particles all processors have

    	//boost::mpi::all_reduce(*system->comm, myN, maxN, boost::mpi::maximum<int>());
    	//boost::mpi::all_reduce(*system->comm, myN, totalN, std::plus<int>());  // to create the dataspace: flaw if totalN changes, should move the dataspace creation here, in the write routine

    	int* array_nparticles = new int [mpi_ranks];   // to write contiguous in file

    	boost::mpi::all_gather(*system->comm, myN, array_nparticles);  // needed for contiguous writing: every mpi rank knows how many particles the other ranks have

    	Real3D L = system->bc->getBoxL();
    	long long step = integrator->getStep();
    	real time_ = step * integrator->getTimeStep();

    	int RANK = 2;

    	int* pids = new int [myN];
    	int* types = new int [myN];
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
    	//assert( i == 0);

    	  if( unfolded ){
    		for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

    			  pids[i] = cit->id();

    		  if (STORE_OPTIONS & (POSITION | IMAGE | ALL)) {
    			  Real3D& pos = cit->position();
    			  Int3D& img = cit->image();

    			  coordina[i*3]   = pos[0] + img[0] * L[0];
    			  coordina[i*3+1] = pos[1] + img[1] * L[1];
    			  coordina[i*3+2] = pos[2] + img[2] * L[2];
    		  }
    		  if (STORE_OPTIONS & (IMAGE | ALL)) {
			  Int3D& img = cit->image();

			  images[i*3]   = img[0];
			  images[i*3+1] = img[1];
			  images[i*3+2] = img[2];
		  }
    		  if (STORE_OPTIONS & (VELOCITY | ALL)) {
    			  Real3D& vel = cit->velocity();

    			  velocities[i*3]   = vel[0];
    			  velocities[i*3+1] = vel[1];
    			  velocities[i*3+2] = vel[2];
    		  }

    		  if (STORE_OPTIONS & (TYPE | ALL)) {

    			  types[i] = cit->type();
    		  }
    		  if (STORE_OPTIONS & (FORCE | ALL)) {
    			  Real3D& force = cit->force();

    			  forces[i*3]   = force[0];
    			  forces[i*3+1] = force[1];
    			  forces[i*3+2] = force[2];
    		  }
    		  if (STORE_OPTIONS & (MASS | ALL)) {

    			  masses[i] = cit->mass();
    		  }
    		  if (STORE_OPTIONS & (STATE | ALL)) {

    			  states[i] = cit->state();
    		  }

    		  if (STORE_OPTIONS & (DRIFT | ALL)) {

    			  drifts[i] = cit->drift();
    		  }

    		  if (STORE_OPTIONS & (CHARGE | ALL)) {

    			  charges[i] = cit->q();
    		  }

    		  if (STORE_OPTIONS & (LAMBDA | ALL)) {

    			  lambdas[i] = cit->lambda();
    		  }

    		  if (STORE_OPTIONS & (LAMBDADERIV | ALL)) {

    			  lambdaDerivs[i] = cit->lambdaDeriv();
    		  }

    		  i++;
    		}
    	  }
    	  else{
    		for(iterator::CellListIterator cit(realCells); !cit.isDone(); ++cit) {

    			  pids[i] = cit->id();

    		  if (STORE_OPTIONS & (POSITION | IMAGE | ALL)) {
    		  	Real3D& pos = cit->position();

    		  	coordina[i*3]   = pos[0];
    		  	coordina[i*3+1] = pos[1];
    		  	coordina[i*3+2] = pos[2];
    		  }
    		  if (STORE_OPTIONS & (IMAGE | ALL)) {
			Int3D& img = cit->image();

			images[i*3]   = img[0];
			images[i*3+1] = img[1];
			images[i*3+2] = img[2];
		  }
    		  if (STORE_OPTIONS & (VELOCITY | ALL)) {
    			Real3D& vel = cit->velocity();

    			velocities[i*3]   = vel[0];
    			velocities[i*3+1] = vel[1];
    			velocities[i*3+2] = vel[2];

    		  }

    		  if (STORE_OPTIONS & (TYPE | ALL)) {

    			  types[i] = cit->type();
    		  }
    		  if (STORE_OPTIONS & (FORCE | ALL)) {
    			  Real3D& force = cit->force();

    			  forces[i*3]   = force[0];
    			  forces[i*3+1] = force[1];
    			  forces[i*3+2] = force[2];
    		  }
    		  if (STORE_OPTIONS & (MASS | ALL)) {

    			  masses[i] = cit->mass();
    		  }
    		  if (STORE_OPTIONS & (STATE | ALL)) {

    			  states[i] = cit->state();
    		  }

    		  if (STORE_OPTIONS & (DRIFT | ALL)) {

    			  drifts[i] = cit->drift();
    		  }

    		  if (STORE_OPTIONS & (CHARGE | ALL)) {

    			  charges[i] = cit->q();
    		  }

    		  if (STORE_OPTIONS & (LAMBDA | ALL)) {

    			  lambdas[i] = cit->lambda();
    		  }

    		  if (STORE_OPTIONS & (LAMBDADERIV | ALL)) {

    			  lambdaDerivs[i] = cit->lambdaDeriv();
    		  }

    		  i++;
    		}
    	  }

		hsize_t count[RANK]; // not used, directly pass myN
		hsize_t offset[RANK];

		// logic to write a file contiguously!
		int sumup = 0;
		if (rank == 0) {
			offset[1] = 0;
		} else {

			for(int L=0; L<rank; L++) {
				sumup += array_nparticles[L];
			}

		offset[1] = sumup;
		}

	if (STORE_OPTIONS & (POSITION | ALL))   			 h5md_append(part_group.position, coordina, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (IMAGE | ALL))   				 h5md_append(part_group.image, images, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (PID | ALL))        			 h5md_append(part_group.id, pids, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (TYPE | ALL))     	  			 h5md_append(part_group.species, types, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (VELOCITY | ALL)) 	  			 h5md_append(part_group.velocity, velocities, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (MASS | ALL))     	  			 h5md_append(part_group.mass, masses, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (CHARGE | ALL))   	 			 h5md_append(part_group.charge, charges, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (STATE | ALL))    	  			 h5md_append(part_group.state, states, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (DRIFT | ALL))       			 h5md_append(part_group.drift, drifts, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (FORCE | ALL))       			 h5md_append(part_group.force, forces, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (LAMBDA | ALL))      			 h5md_append(part_group.lambda, lambdas, step, time_, offset[1], myN, 0);
	if (STORE_OPTIONS & (LAMBDADERIV | ALL)) 			 h5md_append(part_group.lambdaDeriv, lambdaDerivs, step, time_, offset[1], myN, 0);

      	//H5Pclose(plist_id);

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
      }



    void H5MDWriter::write() {

    	if (iomode == 1 || iomode == 0) {
    		write_n_to_1();
    	} else if (iomode == 2) {
		write_n_to_n();
		}
    }


    // Python wrapping
    void H5MDWriter::registerPython() {

      using namespace espressopp::python;

      class_<H5MDWriter, bases<ParticleAccess>, boost::noncopyable >
      ("io_H5MDWriter", init< const shared_ptr< System >&,
                           const shared_ptr< integrator::MDIntegrator >&,
                           const std::string&,
			   const std::string&,
			   const std::string&,
			   const std::string&,
			   const std::string&,
			   const int&,
                           const boost::python::list&,
                           const bool&,
                           const real&,
                           const std::string&,
			   const bool&,
                           const bool&,
			   const int&>())
        .add_property("filename", &H5MDWriter::getFilename,
                                  &H5MDWriter::setFilename)
	.add_property("author", &H5MDWriter::getAuthor,
				  &H5MDWriter::setAuthor)
	.add_property("author_email", &H5MDWriter::getAuthorEmail,
				  &H5MDWriter::setAuthorEmail)
	.add_property("creator", &H5MDWriter::getCreator,
				  &H5MDWriter::setCreator)
	.add_property("creator_version", &H5MDWriter::getCreatorVersion,
					  &H5MDWriter::setCreatorVersion)
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
	.add_property("writers", &H5MDWriter::getAggregators,
				&H5MDWriter::setAggregators)
        .def("dump", &H5MDWriter::write)
	.def("close", &H5MDWriter::close)
      ;
    }
  }
}

