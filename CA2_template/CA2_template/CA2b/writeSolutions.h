#ifndef WRITESOLUTIONS_H_
#define WRITESOLUTIONS_H_
#include <hdf5.h>
#include <deal.II/lac/vector.h>
using namespace dealii;

void writeSolutionsToFileCA2(const Vector<double>& U, std::string solutionTag){
  hid_t file_id, dataset_id, dataspace_id;
  herr_t      status;
  hsize_t dimens_1d;
  std::cout << "writing solution for Coding Assignment 2 to file : " << solutionTag.c_str() << ".h5" << std::endl;

 //create HDF5 file
  std::string solutionFileName(solutionTag); solutionFileName +=  ".h5";
  file_id = H5Fcreate (solutionFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  
  //write solution vector to data set
  dimens_1d = U.size();
  std::vector<double> solutionVector(U.size());
  for (unsigned int i=0; i<U.size(); i++){
    solutionVector[i]=U(i);
  } 
  dataspace_id = H5Screate_simple(1, &dimens_1d, NULL);
  dataset_id = H5Dcreate(file_id, "U", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &solutionVector[0]);  
  status = H5Sclose(dataspace_id);  
  status = H5Dclose(dataset_id);

  //close HDF5 file
  status = H5Fclose(file_id);

}

#endif
