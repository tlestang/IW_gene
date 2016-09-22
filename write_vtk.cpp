#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

void write_fluid_vtk(int time, int Dx, int Dy, double **T, double ***u, const char* folderName)
// Writes the normalized velocities u and v. (u/Uref and v/Uref)
{

  /// Create filename
  stringstream fileName;

  fileName << "fluid_t" << time << ".vtk";

  string output_filename = string(folderName) + "/vtk_fluid/" + fileName.str();

  ofstream output_file;
  double rho0 = 1.0; double u0 = 0.01;

  /// Open file
  output_file.open(output_filename.c_str());

  /// Write VTK header
  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "fluid_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET RECTILINEAR_GRID\n";
  output_file << "DIMENSIONS " << Dx << " " << Dy << " 1" << "\n";
  output_file << "X_COORDINATES " << Dx << " float\n";
  for(int i = 0; i < Dx; ++i)
    output_file << i << " ";
  output_file << "\n";
  output_file << "Y_COORDINATES " << Dy  << " float\n";
  for(int j = 0; j < Dy ; ++j)
    output_file << j  << " ";
  output_file << "\n";
  output_file << "Z_COORDINATES " << 1 << " float\n";
  output_file << 0 << "\n";
  output_file << "POINT_DATA " << (Dx) * (Dy) << "\n";

  /// Write density difference
  output_file << "SCALARS temperature float 1\n";
  output_file << "LOOKUP_TABLE default\n";
  for(int Y =0; Y < Dy ; ++Y)
    for(int X = 0; X < Dx; ++X)
      output_file <<  T[X][Y] << "\n";

  /// Write velocity
  output_file << "VECTORS velocity_vector float\n";
  for(int Y = 0; Y < Dy ; ++Y)
    for(int X = 0; X < Dx; ++X)
      output_file << u[X][Y][0] << " " << u[X][Y][1] << " 0\n";

  /// Close file
  output_file.close();
}
