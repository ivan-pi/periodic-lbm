#include <fstream>

extern "C" {

void c_output_vtk_polydata(
    const char* filename, 
    int n, 
    const double* p, 
    const double* rho, 
    const double* ux, 
    const double* uy) {

  std::ofstream vtkfile(filename);

  vtkfile << "# vtk DataFile Version 2.0\n";
  vtkfile << "meshless lb data\n";
  vtkfile << "ASCII\n";
  vtkfile << "DATASET POLYDATA\n";

  vtkfile << "POINTS " << n << " float\n";
  
  for (int i = 0; i < n; i++) {
    const double *p_ = &p[2*i];
    vtkfile << p_[0] << " " << p_[1] << " 0\n";
  }

  vtkfile << "POINT_DATA " << n << '\n';

  vtkfile << "VECTORS Velocity float\n";

  for (int i = 0; i < n; i++) {
    vtkfile << ux[0] << " " << uy[1] << " 0\n";
  }

  vtkfile << "SCALARS Density float 1\n";
  vtkfile << "LOOKUP_TABLE default\n";

  for (int i = 0; i < n; i++) {
    vtkfile << rho[i] << '\n';
  }

  vtkfile.close();
}

} // extern "C"
