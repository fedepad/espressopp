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
#ifndef _IO_H5MD_HPP
#define _IO_H5MD_HPP

// http://stackoverflow.com/questions/10056393/g-with-python-h-how-to-compile
// http://realmike.org/blog/2012/07/08/embedding-python-tutorial-part-1/
#include "ParticleAccess.hpp"  // keep python.hpp on top
#include "integrator/MDIntegrator.hpp"
#include "io/FileBackup.hpp"
#include "esutil/Error.hpp"
#include "template_py_to_cpp_containers.hpp"
#include <string>
#include <set>

//#include <boost/python/def.hpp>
//#include <boost/python/module.hpp>
//#include <boost/python/args.hpp>


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


    //void read_H5MD();
    //void write_H5MD();





    class H5MD : public ParticleAccess {

    public:

//      H5MD(shared_ptr<System> system,
//              shared_ptr<integrator::MDIntegrator> _integrator,
//              std::string _file_name,
//			  int _iomode,
//              boost::python::list _data_to_store,
//              bool _unfolded,
//              real _length_factor,
//              std::string _length_unit,
//              bool _append) :
//                        ParticleAccess(system),
//                        integrator(_integrator),
//                        file_name( _file_name ),
//						iomode(_iomode),
//                        data_to_store(_data_to_store),
//                        unfolded(_unfolded),
//                        length_factor(_length_factor),
//                        append(_append){
//        setLengthUnit(_length_unit);
//        info_table_ini table_store;
//        set_init_table(table_store, data_to_store);
//
//
//        if (iomode == 1 || iomode == 0) {
//
//			if (system->comm->rank() == 0  && !append) {
//				FileBackup backup(file_name);
//			}
//        } else if (iomode == 2) {
//
//        	if (!append) {
//				int rank = system->comm->rank();
//				size_t filename_length = file_name.length();
//				std::string suffix = file_name.substr(filename_length-3, 3);
//				std::string base_filename = file_name.substr(0,filename_length-3);
//				std::string rankstring = static_cast<std::ostringstream*>( &(std::ostringstream() << rank) )->str();
//				std::string final_name = base_filename + "_" + rankstring + suffix;
//				FileBackup backup(final_name);
//        	}
//
//        }
//      }

      H5MD() {}
      ~H5MD() {}


      /*
       * To make it work with ExtAnalyze (Inherit from ParticleAccess)
       *
       */
      void perform_action(){
        write(shared_ptr<System> ,
                shared_ptr<integrator::MDIntegrator> ,
                std::string ,
    			int ,
                boost::python::list ,
                bool ,
                real ,
                std::string );
      }

      /*
       * Specialized write routines for most I/O common schemes
       *
       */

      void write_n_to_1(shared_ptr<System> system,
              shared_ptr<integrator::MDIntegrator> _integrator,
              std::string _file_name,

              boost::python::list _data_to_store,
              bool _unfolded,
              real _length_factor,
              std::string _length_unit);
      void write_n_to_n(shared_ptr<System> system,
              shared_ptr<integrator::MDIntegrator> _integrator,
              std::string _file_name,

              boost::python::list _data_to_store,
              bool _unfolded,
              real _length_factor,
              std::string _length_unit);

      /*
       * Specialized read routines for most I/O common schemes
       *
       */
      void read_n_to_n(std::string filename);
      void read_n_to_1(std::string filename);

      /*
       * For read or writes
       *
       */
      void distribute_data_to_procs();

      /*
       * I/O routines to be exported to Python
       *
       */
      void write(shared_ptr<System> system,
    		                shared_ptr<integrator::MDIntegrator> _integrator,
    		                std::string _file_name,
							int _iomode,
    		                boost::python::list _data_to_store,
    		                bool _unfolded,
    		                real _length_factor,
    		                std::string _length_unit,
    		                bool _append);
      void read(std::string filename);


      /*
       * detect filesystem type where the given file/files live on
       *  (to set appropriate MPI hints)
       *
       */
       void detect_fs(std::string filename);

      /*
       * printing routines (after file read...)
       *
       */
      void print_number_particles();
      void print_number_timesteps();
      void print_pids_position(bool sorted);
      void print_pids_velocity(bool sorted);
      void print_pids(bool sorted);
      //void sort_by_pid();



//      std::string getFilename(){return file_name;}
//      void setFilename(std::string v){file_name = v;}
//      inline int getIomode(){return iomode;}
//      inline void setIomode(int v){iomode = v;}
//      bool getUnfolded(){return unfolded;}
//      void setUnfolded(bool v){unfolded = v;}
//      bool getAppend(){return append;}
//      void setAppend(bool v){append = v;}

       std::set<std::string> getSetfrompythonlist(boost::python::list data_to_store) {

        return python_list_to_set<std::string>(data_to_store);

       };

       void set_init_table(info_table_ini table_init, boost::python::list initial_list){

        std::set<std::string> ciao = getSetfrompythonlist(initial_list);


        if (ciao.find("all") != ciao.end())
        {
            table_init.all = 1;
        } else {

            if (ciao.find("pid") != ciao.end()) {table_init.pid = 1;}
            if (ciao.find("type") != ciao.end()) {table_init.type = 1;}
            if (ciao.find("mass") != ciao.end()) {table_init.mass = 1;}
            if (ciao.find("charge") != ciao.end()) {table_init.charge = 1;}
            if (ciao.find("lambda") != ciao.end()) {table_init.lambda = 1;}
            if (ciao.find("drift") != ciao.end()) {table_init.drift = 1;}
            if (ciao.find("lambdaDeriv") != ciao.end()) {table_init.lambdaDeriv = 1;}
            if (ciao.find("state") != ciao.end()) {table_init.state = 1;}
            if (ciao.find("position") != ciao.end()) {table_init.position = 1;}
            if (ciao.find("velocity") != ciao.end()) {table_init.velocity = 1;}
            if (ciao.find("force") != ciao.end()) {table_init.force = 1;}

        }


       };

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
      //shared_ptr<integrator::MDIntegrator> integrator;

      //std::string file_name;
      //int iomode; // 0: serial, 1: N-to-1, 2: N-to-N; real 0 not there now
      //boost::python::list data_to_store;  // python list: can pass either 'all' and all particle data structure info
      //// is stored or ['pid', 'mass', 'position'] and only pid, mass and position of the particles will be stored
      //info_table_ini datas;

      //bool unfolded;  // one can choose folded or unfolded coordinates, by default it is folded
      //bool append; //append to existing trajectory file or create a new one
      //real length_factor;  // for example
      //std::string length_unit; // length unit: {could be LJ, nm, A} it is just for user info
    };
  }
}

#endif

