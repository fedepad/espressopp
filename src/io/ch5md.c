/* This file is part of ch5md
 * Copyright (C) 2013-2015 Pierre de Buyl
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD license. See the LICENSE file for details.
 */

//#include "hdf5.h"
#include "ch5md.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <sys/vfs.h>


#if OP_SYSTEM == OSX
    #define USING_OSX
        #undef USING_LINUX
        #undef USING_WINDOWS
#elif OP_SYSTEM == LINUX
        #define USING_LINUX
        #undef USING_OSX
        #undef USING_WINDOWS
#elif OP_SYSTEM == WINDOWS
        #define USING_WINDOWS
        #undef USING_LINUX
        #undef USING_OSX
#endif

//for OSX
#ifdef USING_OSX
	#include <sys/param.h>
	#include <sys/mount.h>
// end OSX
#endif

#ifdef USING_LINUX
	#include <sys/param.h>
	#include <sys/mount.h>
// end OSX
#endif

//#ifdef _PARALLEL_V
    #include "mpi.h"
//#endif


#define MIN_CHUNK 10
#define MAX_CHUNK 256
#define MAX_RANK 5


h5md_file h5md_create_file (const char *filename, const char *author, const char *author_email, const char *creator, const char *creator_version, int parallel, const char* aggregators)  // parallel default, otherwise fall back to serial lib...-> parallel = 0
{

  h5md_file file;
  hid_t g, g1;
  hid_t a, s, t;
  hsize_t dims[1];
  herr_t status;
  hid_t acc_template;
  file.version[0] = 1;
  file.version[1] = 0;

  if (parallel == 0 ) {

    file.id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file.id > 0);

  } else if (parallel == 1) {
	//printf("In parallel mode!\n");
    acc_template = H5Pcreate(H5P_FILE_ACCESS);

    MPI_Info info;
    MPI_Info_create(&info);

    const char* hint_stripe = "striping_unit";
    const char* stripe_value = "4194304";
    /* Disables ROMIO's data-sieving */
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");

     /* Enable ROMIO's collective buffering */
    MPI_Info_set(info, "romio_cb_read", "enable");
    MPI_Info_set(info, "romio_cb_write", "enable");


    // collective buffering tuning
    MPI_Info_set(info, "cb_nodes", (char*)aggregators);  // number of aggregators: how many MPI ranks perform the write...
    //MPI_Info_set(info, "cb_buffer_size", "disable"); //


    MPI_Info_set(info, (char*)hint_stripe, (char*)stripe_value); // 4MB stripe. cast to avoid spurious warnings. OK for Lustre


    //add my detection of fs...don't know where it went!!!!!
    H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, info);
    file.id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, acc_template);
    assert(file.id > 0);
    struct statfs filesystem_info;
	statfs(filename, &filesystem_info);
	if (filesystem_info.f_type == 0x0BD00BD0) {
		printf("Lustre baby!\n");
		// set approriate hints for Lustre
	} else if (filesystem_info.f_type == 0x47504653) {
		printf("GPFS baby!\n");
		// set approriate hints for Lustre
	} else if (filesystem_info.f_type == 0xEF53) {

		printf("Ext2,ext3 or ext4!!\n");
	} else if (filesystem_info.f_type == 0x4244) {
		printf("HFS, probably on a mac!!\n");

	} else if (filesystem_info.f_type == 0x6969) {
		printf("NFS!!\n");

	} else if (filesystem_info.f_type == 0xFF534D42) {
		printf("CIFS!!\n");

	} else if (filesystem_info.f_type == 0x58465342) {
		printf("XFS!!\n");

	}

    //should free the MPI_Info object
    //MPI_Info_free(&info);

  }

  g = H5Gcreate(file.id, "h5md", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  dims[0] = 2;
  s = H5Screate_simple(1, dims, NULL);
  a = H5Acreate(g, "version", H5T_NATIVE_INT, s, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(a, H5T_NATIVE_INT, file.version);
  status = H5Aclose(a);
  status = H5Sclose(s);

  s = H5Screate(H5S_SCALAR);

  g1 = H5Gcreate(g, "author", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  t = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(t, strlen(author));
  a = H5Acreate(g1, "name", t, s, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(a, t, author);
  status = H5Aclose(a);
  status = H5Tclose(t);
  if (NULL!=author_email) {
    t = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(t, strlen(author_email));
    a = H5Acreate(g1, "author_email", t, s, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(a, t, author_email);
    status = H5Aclose(a);
    status = H5Tclose(t);
  }
  status = H5Gclose(g1);

  g1 = H5Gcreate(g, "creator", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  t = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(t, strlen(creator));
  a = H5Acreate(g1, "name", t, s, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(a, t, creator);
  status = H5Aclose(a);
  status = H5Tclose(t);
  t = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(t, strlen(creator_version));
  a = H5Acreate(g1, "version", t, s, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(a, t, creator_version);
  status = H5Aclose(a);
  status = H5Tclose(t);
  status = H5Gclose(g1);

  status = H5Gclose(g); // to make it work in parallel, we have to close it. Without closing it seems to work in serial case...????!!!

  file.particles = H5Gcreate(file.id, "particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  file.observables = H5Gcreate(file.id, "observables", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  file.parameters = H5Gcreate(file.id, "parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (parallel == 1) H5Pclose(acc_template);

  return file;
}

int h5md_close_file(h5md_file file) {
  H5Gclose(file.particles);
  H5Gclose(file.observables);
  H5Gclose(file.parameters);

  //Any opened objects before file close?
//  printf("Inside this loop before\n");
//  int norphans = H5Fget_obj_count(file.id, H5F_OBJ_ALL);
//  if (norphans > 1) { /* expect 1 for the file we have not closed */
//	 int i;
//	 H5O_info_t info;
//	 char name[64];
//	 hid_t * objects = calloc(norphans, sizeof(hid_t));
//	 H5Fget_obj_ids(file.id, H5F_OBJ_ALL, -1, objects);
//	 for (i=0; i<norphans; i++) {
//		 H5Oget_info(objects[i], &info);
//		 H5Iget_name(objects[i], name, 64);
//		 printf("%d of %zd things still open: %d with name %s of type %d\n", i, norphans, objects[i], name, info.type);
//		 //if (info.type == H5O_TYPE_GROUP) H5Gclose(objects[i]);
//	 }
//}

  //hid_t grid;
  //H5Gopen2(grid, "/h5md", H5P_DEFAULT);
  //should close h5md group
 // H5Gclose(grid);
  //H5Gclose(33554432);

  //H5Gclose(33554436);


  H5Fclose(file.id);


//  int norphans = H5Fget_obj_count(file.id, H5F_OBJ_ALL);
//  if (norphans > 1) { /* expect 1 for the file we have not closed */
//	 printf("Inside this loop after\n");
//	 int i;
//	 H5O_info_t info[norphans];
//	 char name[64];
//	 hid_t * objects = calloc(norphans, sizeof(hid_t));
//	 H5Fget_obj_ids(file.id, H5F_OBJ_ALL, -1, objects);
//	 for (i=0; i<norphans; i++) {
//		 H5Oget_info(objects[i], &info[i]);
//		 H5Iget_name(objects[i], name, 64);
//		 printf("%d of %zd things still open: %d with name %s of type %d\n", i, norphans, objects[i], name, info[i].type);
//		 printf("objects id: %d\n", objects[i]);
//		 if (info[i].type == H5O_TYPE_GROUP) {
//
//			 printf("What the hell is going on!!!!!\n");
//			 H5Gopen2(objects[i], name, H5P_DEFAULT);
//			 H5Gclose(objects[i]);
//		 }
//		 //if (info.type == H5O_TYPE_GROUP) H5Gclose(objects[i]);
//	 }
//}









//  int norphans = H5Fget_obj_count(file.id, H5F_OBJ_ALL);
//    if (norphans > 1) { /* expect 1 for the file we have not closed */
//  	 int i;
//  	 H5O_info_t info;
//  	 char name[64];
//  	 hid_t * objects = calloc(norphans, sizeof(hid_t));
//  	 H5Fget_obj_ids(file.id, H5F_OBJ_ALL, -1, objects);
//  	 for (i=0; i<norphans; i++) {
//  		 H5Oget_info(objects[i], &info);
//  		 H5Iget_name(objects[i], name, 64);
//  		 printf("%d of %zd things still open: %d with name %s of type %d\n", i, norphans, objects[i], name, info.type);
//  	 }
//  }

  return 0;
}

hid_t h5md_open_file (const char *filename)
{
  hid_t file;
  hid_t g;
  hid_t a, s;
  hsize_t dims[1];
  int version[2];
  int version_ok;

  file = H5Fopen(filename,  H5F_ACC_RDWR, H5P_DEFAULT); // change it to access it in parallel

  g = H5Gopen(file, "h5md", H5P_DEFAULT);

  version_ok = false;
  dims[0] = 2;
  a = H5Aopen(g, "version", H5P_DEFAULT);
  s = H5Aget_space(a);
  if (!(H5Sis_simple(s)>0)) {
    printf("H5MD version is not a simple dataspace");
    H5Sclose(s);
    H5Aclose(a);
    H5Gclose(g);
    H5Fclose(file);
  } else {
    if (H5Sget_simple_extent_ndims(s)==1) {
      H5Sget_simple_extent_dims(s, dims, NULL);
      if (dims[0]==2) {
	H5Aread(a, H5T_NATIVE_INT, version);
	if ( (version[0]==1) && (version[1]==0) ) {version_ok = true;}
      }
    }
  }

  H5Aclose(a);
  H5Sclose(s);
  H5Gclose(g);

  if (!version_ok) H5Fclose(file);

  return file;

}

h5md_particles_group h5md_create_particles_group(h5md_file file, const char *name)
{
  h5md_particles_group group;

  group.group = H5Gcreate(file.particles, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  return group;
}

hid_t h5md_open_particles_group(hid_t particles, const char *name)
{
  hid_t group;

  group = H5Gopen(particles, name, H5P_DEFAULT);

  return group;
}

h5md_element h5md_create_time_data(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, h5md_element *link)
{

  h5md_element td;

  hid_t spc, plist;
  hsize_t dims[MAX_RANK], max_dims[MAX_RANK], chunks[MAX_RANK];
  herr_t status;

  int i;

  dims[0] = 0 ;
  max_dims[0] = H5S_UNLIMITED ;
  for (i=0; i<rank; i++) {
    dims[i+1] = int_dims[i];
    max_dims[i+1] = int_dims[i];
  }
  chunks[0] = 1 ;
  if (MAX_CHUNK<int_dims[0]/4) {
    chunks[1] = MAX_CHUNK;
  } else {
    chunks[1] = int_dims[0];
  }
  for (i=1; i<rank; i++) {
    chunks[i+1]=int_dims[i];
  }

  td.group = H5Gcreate(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (NULL==link) {
    spc = H5Screate_simple( 1 , dims, max_dims);
    plist = H5Pcreate(H5P_DATASET_CREATE);
    status = H5Pset_chunk(plist, 1, chunks);
    td.time = H5Dcreate(td.group, "time", H5T_NATIVE_DOUBLE, spc, H5P_DEFAULT, plist, H5P_DEFAULT);
    td.step = H5Dcreate(td.group, "step", H5T_NATIVE_INT, spc, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(spc);
    td.link=NULL;
  } else {
    td.link = link;
    status = H5Lcreate_hard(link->group, "step", td.group, "step", H5P_DEFAULT, H5P_DEFAULT);
    status = H5Lcreate_hard(link->group, "time", td.group, "time", H5P_DEFAULT, H5P_DEFAULT);
  }

//  chunks[0] = 128;
//  chunks[1] = 128;
  spc = H5Screate_simple( rank+1 , dims, max_dims) ;
  plist = H5Pcreate(H5P_DATASET_CREATE);
  status = H5Pset_chunk(plist, rank+1, chunks);
//  int mpi_rankk;
//  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rankk);
//  printf("chunks[0]: %d, chunks[1]: %d; from MPI rank: %d\n", chunks[0], chunks[1], mpi_rankk);
//  const int fill_value_pids = -1;
//  if (strcmp(name, "pids") == 0) {
//	  H5Pset_fill_value(plist, H5T_NATIVE_INT, &fill_value_pids);
//  }
  td.value = H5Dcreate(td.group, "value", datatype, spc, H5P_DEFAULT, plist, H5P_DEFAULT);
  H5Pclose(plist);
  status = H5Sclose(spc);

  td.datatype = datatype;
  td.is_time = true;

  return td;

}

int h5md_close_element(h5md_element e)
{
  herr_t status;

  if (!e.is_time) return 0;

  if (e.link==NULL) {
    status = H5Dclose(e.step);
    status = H5Dclose(e.time);
  }
  status = H5Dclose(e.value);
  status = H5Gclose(e.group);

  return 0;

}

h5md_element h5md_create_fixed_data_simple(hid_t loc, const char *name, int rank, int int_dims[], hid_t datatype, void *data)
{

  h5md_element fd;

  hid_t spc;
  hsize_t dims[H5S_MAX_RANK];
  herr_t status;
  int i;

  for (i=0; i<rank; i++) {
    dims[i] = int_dims[i];
  }

  spc = H5Screate_simple(rank, dims, NULL);
  fd.value = H5Dcreate(loc, name, datatype, spc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Sclose(spc);
  status = H5Dwrite(fd.value, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status = H5Dclose(fd.value);

  fd.is_time = false;

  return fd;

}

h5md_element h5md_create_fixed_data_scalar(hid_t loc, const char *name, hid_t datatype, void *data)
{

  h5md_element fd;

  hid_t spc;
  herr_t status;

  spc = H5Screate(H5S_SCALAR);
  fd.value = H5Dcreate(loc, name, datatype, spc, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Sclose(spc);
  status = H5Dwrite(fd.value, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  status = H5Dclose(fd.value);

  fd.is_time = false;

  return fd;

}

int h5md_extend_by_one(hid_t dset, hsize_t *dims) {

  hid_t file_space;
  int rank;
  hsize_t maxdims[H5S_MAX_RANK];

  // Get dataset information
  file_space = H5Dget_space(dset);
  rank = H5Sget_simple_extent_ndims(file_space);
  if (rank > H5S_MAX_RANK) {
    return CH5MD_RANK_ERROR;
  }
  H5Sget_simple_extent_dims(file_space, dims, maxdims);
  H5Sclose(file_space);

  // Extend dimensions by one
  dims[0] = dims[0]+1;
  H5Dset_extent(dset, dims);

  return 0;

}

int h5md_append(h5md_element e, void *data, int step, double time, hid_t plist_id, int offset, int data_size, int shared_file) {

  hid_t mem_space, file_space;
  int i, rank;
  hsize_t dims[H5S_MAX_RANK];
  hsize_t start[H5S_MAX_RANK], count[H5S_MAX_RANK];

  // If not a time-dependent H5MD element, do nothing
  if (!e.is_time) return 0;

  if (NULL==e.link) {
    h5md_extend_by_one(e.step, dims);

    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;

    // Select the hyperslab
    file_space = H5Dget_space(e.step);
    rank = H5Sget_simple_extent_ndims(file_space);
    mem_space = H5Screate_simple(rank-1, dims+1, NULL);
    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;
    for (i=1 ; i<rank ; i++) {
      start[i] = 0;
      count[i] = dims[i];
    }
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    if (shared_file == 1) {
    	//H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    	H5Dwrite(e.step, H5T_NATIVE_INT, mem_space, file_space, plist_id, (void *)&step);
    } else {
    	H5Dwrite(e.step, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, (void *)&step);
    }
    H5Sclose(file_space);
    H5Sclose(mem_space);

    h5md_extend_by_one(e.time, dims);

    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;

    // Select the hyperslab
    file_space = H5Dget_space(e.time);
    rank = H5Sget_simple_extent_ndims(file_space);
    mem_space = H5Screate_simple(rank-1, dims+1, NULL);
    // Define hyperslab selection
    start[0] = dims[0]-1;
    count[0] = 1;
    for (i=1 ; i<rank ; i++) {
      start[i] = 0;
      count[i] = dims[i];
    }
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    if (shared_file == 1) {
    	//H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
    	H5Dwrite(e.time, H5T_NATIVE_DOUBLE, mem_space, file_space, plist_id, (void *)&time);
    } else {
    	H5Dwrite(e.time, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, (void *)&time);
    }
    H5Sclose(file_space);
    H5Sclose(mem_space);
  }

  h5md_extend_by_one(e.value, dims);
  dims[1] = data_size;

  // Select the hyperslab
  file_space = H5Dget_space(e.value);
  rank = H5Sget_simple_extent_ndims(file_space);
  mem_space = H5Screate_simple(rank-1, dims+1, NULL);

  // Define hyperslab selection
  start[0] = dims[0]-1;
  count[0] = 1;
  for (i=1 ; i<rank ; i++) {
    start[i] = 0;
    count[i] = dims[i];
  }

  // reuse logic in #
  if (shared_file == 1) {
	  start[1] = offset;
	  count[1] = data_size;

	  //start[0] = offset;

	//  int mpi_rankk;
	//  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rankk);
	//  printf("start[0]: %d, start[1]: %d, count[0]: %d, count[1]: %d; from MPI rank: %d\n", start[0], start[1], count[0], count[1], mpi_rankk);
	  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
	  H5Dwrite(e.value, e.datatype, mem_space, file_space, plist_id, data);
  } else {
	  H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
	  H5Dwrite(e.value, e.datatype, mem_space, file_space, H5P_DEFAULT, data);

  }
  H5Sclose(file_space);
  H5Sclose(mem_space);

  return 0;
}

int h5md_create_box(h5md_particles_group *group, int dim, const char *boundary[], bool is_time, double value[], h5md_element *link)
{
  hid_t spc, att, t;
  hsize_t dims[1];
  int int_dims[1];
  herr_t status;
  int i;
  size_t boundary_length, tmp;
  //char *tmp_boundary;

  // Create box group
  group->box = H5Gcreate(group->group, "box", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create dimension attribute
  spc = H5Screate(H5S_SCALAR);
  att = H5Acreate(group->box, "dimension", H5T_NATIVE_INT, spc, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att, H5T_NATIVE_INT, &dim);
  status = H5Aclose(att);
  status = H5Sclose(spc);

  // Compute the size of the string type for boundary
  dims[0] = dim;
  boundary_length=0;
  for (i=0; i<dim; i++) {
    tmp = strlen(boundary[i])+1;
    if (tmp>boundary_length) {
      boundary_length=tmp;
    }
  }
  char *tmp_boundary = malloc(dim*sizeof(char)*boundary_length);
  for (i=0; i<dim; i++) {
    strcpy((tmp_boundary+i*boundary_length), boundary[i]);
  }
  // Create boundary attribute
  t = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(t, boundary_length);
  spc = H5Screate_simple(1, dims, NULL);
  att = H5Acreate(group->box, "boundary", t, spc, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(att, t, tmp_boundary);
  free(tmp_boundary);
  status = H5Aclose(att);
  status = H5Sclose(spc);
  status = H5Tclose(t);

  // Create edges
  // Check if the box is time-dependent or not
  int_dims[0]=dim;
  if (is_time) {
    group->box_edges = h5md_create_time_data(group->box, "edges", 1, int_dims, H5T_NATIVE_DOUBLE, link);
  } else {
    if (NULL!=value) {
      group->box_edges = h5md_create_fixed_data_simple(group->box, "edges", 1, int_dims, H5T_NATIVE_DOUBLE, value);
    }
  }

  status = H5Gclose(group->box);

  return 0;
}

int h5md_write_string_attribute(hid_t loc, const char *obj_name,
    const char *att_name, const char *value)
{
  hid_t obj;
  hid_t s, t, a;
  herr_t status;

#ifdef _PARALLEL_V
  obj = H5Oopen(loc, obj_name, H5P_DEFAULT);  // change the acces propertu list for parallel IO mode?
#else
  obj = H5Oopen(loc, obj_name, H5P_DEFAULT);
#endif

  t = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(t, strlen(value));
  s = H5Screate(H5S_SCALAR);
  a = H5Acreate(obj, att_name, t, s, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(a, t, value);
  status = H5Aclose(a);
  status = H5Sclose(s);
  status = H5Tclose(t);

  status = H5Oclose(obj);

  return 0;
}

//int h5md_write_double_attribute(hid_t loc, const char *obj_name,
//    const char *att_name, const char *value)
//{
//  hid_t obj;
//  hid_t s, t, a;
//  herr_t status;
//
//#ifdef _PARALLEL_V
//  obj = H5Oopen(loc, obj_name, H5P_DEFAULT);  // change the acces propertu list for parallel IO mode?
//#else
//  obj = H5Oopen(loc, obj_name, H5P_DEFAULT);
//#endif
//
//  t = H5Tcopy(H5T_C_S1);
//  status = H5Tset_size(t, strlen(value));
//  s = H5Screate(H5S_SCALAR);
//  a = H5Acreate(obj, att_name, t, s, H5P_DEFAULT, H5P_DEFAULT);
//  status = H5Awrite(a, t, value);
//  status = H5Aclose(a);
//  status = H5Sclose(s);
//  status = H5Tclose(t);
//
//  status = H5Oclose(obj);
//
//  return 0;
//}
//
//int h5md_write_string_long(hid_t loc, const char *obj_name,
//    const char *att_name, const char *value)
//{
//  hid_t obj;
//  hid_t s, t, a;
//  herr_t status;
//
//#ifdef _PARALLEL_V
//  obj = H5Oopen(loc, obj_name, H5P_DEFAULT);  // change the acces propertu list for parallel IO mode?
//#else
//  obj = H5Oopen(loc, obj_name, H5P_DEFAULT);
//#endif
//
//  t = H5Tcopy(H5T_C_S1);
//  status = H5Tset_size(t, strlen(value));
//  s = H5Screate(H5S_SCALAR);
//  a = H5Acreate(obj, att_name, t, s, H5P_DEFAULT, H5P_DEFAULT);
//  status = H5Awrite(a, t, value);
//  status = H5Aclose(a);
//  status = H5Sclose(s);
//  status = H5Tclose(t);
//
//  status = H5Oclose(obj);
//
//  return 0;
//}
