/*------------------------------------------------------ 
This file is part of the LAMMPS-LBM-PSM project.

LAMMPS-LBM-PSM is an open-source project distributed
under the GNU General Public License.

See the README and License file in the top-level 
LAMMPS-LBM-PSM directory for more details.

Tim Najuch
------------------------------------------------------*/

#ifdef FIX_CLASS

FixStyle(lbm-psm-vti,WriteVTI)

#else

#ifndef WRITEVTI_H
#define WRITEVTI_H

#include <algorithm>
#include <fstream> 
#include <functional>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
#include <sstream>
#include <string>
#include <filesystem>

#include "comm.h"
#include "fix.h"
#include "modify.h"

#include "fix_LBM_PSM.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace std;

namespace fs = std::filesystem;

class WriteVTI : public Fix{
  private:
    int nx, ny, nz;
    int decomposition[3];
    std::string fileName;
    bool binary_;
    bool useDouble_;

  public:
    WriteVTI(class LAMMPS *, int, char **);
    ~WriteVTI();
    int setmask();
    void init();
    void pre_force(int);

    void write_vti_wrapper(std::string name_, int timestep, bool binary, bool useDouble,
                         vector<double> &x_,
                         vector<double> &y_,
                         vector<double> &z_,
                         vector<double> &B_, double B0_,
                         vector<double> &rho_, double rho0_,
                         vector<double> &u_, double u0_);


    template <typename T> void execute_write_vti(std::string name_, int timestep, bool binary,
                         vector<double> &x_,
                         vector<double> &y_,
                         vector<double> &z_,
                         vector<double> &B_, double B0_,
                         vector<double> &rho_, double rho0_,
                         vector<double> &u_, double u0_);

    class fix_LBM_PSM *fixLBMPSM;

};


template <typename T>
void WriteVTI::execute_write_vti(std::string name_, int timestep, bool binary,
                                 vector<double> &x_, vector<double> &y_,
                                 vector<double> &z_, vector<double> &B_, double B0_,
                                 vector<double> &rho_, double rho0_, vector<double> &u_, double u0_)
{
    int rank = comm->me;
    int size = comm->nprocs;
    int env = 1; 

    std::string typeStr = (sizeof(T) == 8) ? "Float64" : "Float32";
    size_t scalarSize   = sizeof(T);

    // Create folder for output 
    std::string folder = "flow_field_output_data";
    if (rank == 0 && timestep == 0) {
        if (fs::exists(folder)) fs::remove_all(folder);
        fs::create_directory(folder);
    }
    MPI_Barrier(world);

    // Indices and dimensions of lattice grid
    std::vector<int> pC = fixLBMPSM->dynamics->get_procCoordinates();
    int ip = pC[0], jp = pC[1], kp = pC[2];
    int nxL = fixLBMPSM->dynamics->get_nxLocal(ip);
    int nyL = fixLBMPSM->dynamics->get_nyLocal(jp);
    int nzL = (domain->dimension == 3) ? fixLBMPSM->dynamics->get_nzLocal(kp) : 1;

    // VTI PointData requires adjacent pieces to share exactly one boundary node.
    // Non-leftmost processors include their left ghost cell as the first output point;
    // that ghost cell holds the same value as the adjacent processor's last real cell,
    // so the shared node is consistent across pieces.
    int iStart = (ip > 0) ? env - 1 : env;
    int jStart = (jp > 0) ? env - 1 : env;
    int kStart = (domain->dimension == 3) ? ((kp > 0) ? env - 1 : env) : 0;

    int nx_i = nxL - env - iStart;
    int ny_i = nyL - env - jStart;
    int nz_i = (domain->dimension == 3) ? (nzL - env - kStart) : 1;

    // offset = sum of real-cell counts of all preceding processors, minus 1 for the shared node
    int offsetX = 0; for (int i = 0; i < ip; i++) offsetX += (fixLBMPSM->dynamics->get_nxLocal(i) - 2 * env);
    if (ip > 0) offsetX--;
    int offsetY = 0; for (int j = 0; j < jp; j++) offsetY += (fixLBMPSM->dynamics->get_nyLocal(j) - 2 * env);
    if (jp > 0) offsetY--;
    int offsetZ = 0;
    if (domain->dimension == 3) {
        for (int k = 0; k < kp; k++) offsetZ += (fixLBMPSM->dynamics->get_nzLocal(k) - 2 * env);
        if (kp > 0) offsetZ--;
    }

    int nxTot = 0; for (int i = 0; i < decomposition[0]; i++) nxTot += (fixLBMPSM->dynamics->get_nxLocal(i) - 2 * env);
    int nyTot = 0; for (int j = 0; j < decomposition[1]; j++) nyTot += (fixLBMPSM->dynamics->get_nyLocal(j) - 2 * env);
    int nzTot = 1;
    if (domain->dimension == 3) {
        nzTot = 0; for (int k = 0; k < decomposition[2]; k++) nzTot += (fixLBMPSM->dynamics->get_nzLocal(k) - 2 * env);
    }

    // Set the output path
    std::ostringstream oss;
    oss << name_ << "." << std::setw(8) << std::setfill('0') << timestep;
    std::string baseFileName = oss.str();
    std::string baseTimePath = folder + "/" + baseFileName;

    // Clip the solid fraction at one and get also density + velocity for output (_o)
    int nP = nx_i * ny_i * nz_i;
    std::vector<T> B_o(nP), rho_o(nP), u_o(nP * 3);

    int id = 0;
    int zEnd = (domain->dimension == 3) ? nzL - env : 1;

    for (int k = kStart; k < zEnd; k++) {
        for (int j = jStart; j < nyL - env; j++) {
            for (int i = iStart; i < nxL - env; i++) {
                int idx = i * nyL * nzL + j * nzL + k;
                auto pD = fixLBMPSM->dynamics->getParticleDataOnLatticeNode(idx);
                double Bv = (pD.solidFraction[0] + pD.solidFraction[1]);
                
                // Casting for T (float or double) // if double not needed (e.g. for detailed postprocessing), can save space
                B_o[id]     = static_cast<T>((Bv > 1.0 ? 1.0 : Bv) * B0_);
                rho_o[id]   = static_cast<T>(rho_[idx] * rho0_);
                u_o[3*id+0] = static_cast<T>(u_[3*idx+0] * u0_);
                u_o[3*id+1] = static_cast<T>(u_[3*idx+1] * u0_);
                u_o[3*id+2] = (domain->dimension == 3) ? static_cast<T>(u_[3*idx+2] * u0_) : static_cast<T>(0.0);
                id++;
            }
        }
    }

    // Physical lattice spacing (uniform grid: dx = dy = dz).
    // Global VTI origin = physical position of global node (0,0,0) = proc 0's first output point.
    // Each processor derives this consistently without MPI:
    //   origin = coordinate_of_local_first_output_node - offset * dx
    // origin_global being negative is handled correctly since x_/y_/z_ store full physical coords.
    double dx = fixLBMPSM->unitConversion->get_dx();
    int idx_first = iStart * nyL * nzL + jStart * nzL + kStart;
    double originX = x_[idx_first] - offsetX * dx;
    double originY = y_[idx_first] - offsetY * dx;
    double originZ = (domain->dimension == 3) ? z_[idx_first] - offsetZ * dx : 0.0;

    // Write the local vti files. One file per core
    std::string vfname = baseTimePath + "_rank" + std::to_string(rank) + ".vti";
    int xE = offsetX + nx_i - 1;
    int yE = offsetY + ny_i - 1;
    int zE_ext = (domain->dimension == 3) ? (offsetZ + nz_i - 1) : 0;

    if (binary) {
        std::ofstream vti(vfname, std::ios::out | std::ios::binary);
        uint32_t b_sz = (uint32_t)(nP * scalarSize);
        uint32_t u_sz = (uint32_t)(nP * 3 * scalarSize);

        vti << "<?xml version=\"1.0\"?>\n<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
        vti << " <ImageData WholeExtent=\"0 " << nxTot-1 << " 0 " << nyTot-1 << " 0 " << nzTot-1 << "\""
            << " Origin=\"" << originX << " " << originY << " " << originZ << "\""
            << " Spacing=\"" << dx << " " << dx << " " << dx << "\">\n";
        vti << "  <Piece Extent=\"" << offsetX << " " << xE << " " << offsetY << " " << yE << " " << offsetZ << " " << zE_ext << "\">\n";
        vti << "   <PointData>\n";
        vti << "    <DataArray type=\"" << typeStr << "\" Name=\"SolidFraction\" format=\"appended\" offset=\"0\"/>\n";
        vti << "    <DataArray type=\"" << typeStr << "\" Name=\"Density\" format=\"appended\" offset=\"" << (4 + b_sz) << "\"/>\n";
        vti << "    <DataArray type=\"" << typeStr << "\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << (2*(4 + b_sz)) << "\"/>\n";
        vti << "   </PointData>\n  </Piece>\n </ImageData>\n <AppendedData encoding=\"raw\">_";
        vti.flush();

        vti.write((char*)&b_sz, 4); vti.write((char*)B_o.data(), b_sz);
        vti.write((char*)&b_sz, 4); vti.write((char*)rho_o.data(), b_sz);
        vti.write((char*)&u_sz, 4); vti.write((char*)u_o.data(), u_sz);
        vti << "\n </AppendedData>\n</VTKFile>";
    } else {
        std::ofstream vti(vfname);
        vti << "<?xml version=\"1.0\"?>\n<VTKFile type=\"ImageData\" version=\"0.1\">\n";
        vti << " <ImageData WholeExtent=\"0 " << nxTot-1 << " 0 " << nyTot-1 << " 0 " << nzTot-1 << "\""
            << " Origin=\"" << originX << " " << originY << " " << originZ << "\""
            << " Spacing=\"" << dx << " " << dx << " " << dx << "\">\n";
        vti << "  <Piece Extent=\"" << offsetX << " " << xE << " " << offsetY << " " << yE << " " << offsetZ << " " << zE_ext << "\">\n";
        vti << "   <PointData>\n";
        
        auto ascii = [&](std::string nm, int comp, const std::vector<T>& data) {
            vti << "    <DataArray type=\"" << typeStr << "\" Name=\"" << nm << "\" NumberOfComponents=\"" << comp << "\" format=\"ascii\">\n";
            for(size_t i=0; i<data.size(); i++) vti << data[i] << (i%10==9 ? "\n" : " ");
            vti << "\n    </DataArray>\n";
        };
        ascii("SolidFraction", 1, B_o);
        ascii("Density", 1, rho_o);
        ascii("Velocity", 3, u_o);
        vti << "   </PointData>\n  </Piece>\n </ImageData>\n</VTKFile>";
    }

    // Write the master pvti file on core 0
    int lEx[6] = {offsetX, xE, offsetY, yE, offsetZ, zE_ext};
    std::vector<int> allEx;
    if (rank == 0) allEx.resize(6 * size);
    MPI_Gather(lEx, 6, MPI_INT, allEx.data(), 6, MPI_INT, 0, world);

    if (rank == 0) {
        std::ofstream pvti(baseTimePath + ".pvti");
        pvti << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PImageData\" version=\"0.1\">\n";
        pvti << " <PImageData WholeExtent=\"0 " << nxTot-1 << " 0 " << nyTot-1 << " 0 " << nzTot-1 << "\""
             << " Origin=\"" << originX << " " << originY << " " << originZ << "\""
             << " Spacing=\"" << dx << " " << dx << " " << dx << "\">\n";
        pvti << "  <PPointData>\n";
        pvti << "   <PDataArray type=\"" << typeStr << "\" Name=\"SolidFraction\"/>\n";
        pvti << "   <PDataArray type=\"" << typeStr << "\" Name=\"Density\"/>\n";
        pvti << "   <PDataArray type=\"" << typeStr << "\" Name=\"Velocity\" NumberOfComponents=\"3\"/>\n";
        pvti << "  </PPointData>\n";
        for (int i = 0; i < size; i++) {
            int* e = &allEx[i*6];
            pvti << "  <Piece Extent=\"" << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << " " << e[4] << " " << e[5] 
                 << "\" Source=\"" << baseFileName << "_rank" << i << ".vti\"/>\n";
        }
        pvti << " </PImageData>\n</VTKFile>";
    }
}

#endif
#endif
