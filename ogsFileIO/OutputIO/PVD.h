
#pragma once

#include <string>
#include <vector>

typedef struct
{
	double timestep;
	std::string vtk_file;
} VTK_Info;

class PVDFileIO
{
public:
	std::vector<VTK_Info> vec_dataset;
	std::string pvd_file_name;
	std::string pvd_vtk_file_name_base;
    std::string pvd_vtk_file_path_base;

public:
    PVDFileIO(void){}
	virtual ~PVDFileIO(void){}

protected:
	bool WriteHeaderOfPVD(std::fstream &fin);
	bool WriteEndOfPVD(std::fstream &fin);
	bool WriteDatasetOfPVD(std::fstream &fin, double timestep, const std::string &vtkfile);

public:
	bool InitializePVD(const std::string &file_base_namepcs_type_name);
	bool UpdatePVD(const std::string &pvdfile, const std::vector<VTK_Info> &vec_vtk);

};
