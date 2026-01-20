#include "System.h"
#include "Storage.h"
#include <cmath>
#include <limits>

Storage::Storage(std::weak_ptr<System> system_) {
    system = system_;

    std::ofstream statesOutput("Temp.sta");

    if (auto SYSTEM = system.lock()) {
        statesOutput << "node_count " << SYSTEM->coordInfoVecs.nodeLocX.size() << '\n';
        statesOutput << "edge_count " << SYSTEM->coordInfoVecs.num_edges << '\n';
        statesOutput << "elem_count " << SYSTEM->coordInfoVecs.num_triangles << '\n';
    }

    statesOutput.close();
}

// Helper function to sanitize values for VTK output
// Replaces NaN and Inf with safe default values
static double sanitize_value(double val, double default_val = 0.0) {
    if (std::isnan(val) || std::isinf(val)) {
        return default_val;
    }
    return val;
}

void Storage::print_VTK_File(void) {

    std::shared_ptr<System> SYSTEM = system.lock();
    if (!SYSTEM) return;

    // ----------------- filename -----------------
    iteration += 1;
    int digits = (iteration > 0) ? static_cast<int>(std::ceil(std::log10(static_cast<double>(iteration + 1)))) : 1;
    std::string format = ".vtk";
    std::string Number;

    if (digits <= 1) {
        Number = "0000" + std::to_string(iteration - 1);
    } else if (digits == 2) {
        Number = "000" + std::to_string(iteration - 1);
    } else if (digits == 3) {
        Number = "00" + std::to_string(iteration - 1);
    } else {
        Number = "0" + std::to_string(iteration - 1);
    }

    std::string initial = "Animation_realistic/01_08_26/_strain_test_case_2_";
    std::string Filename = initial + Number + format;

    std::ofstream ofs(Filename.c_str());
    if (!ofs) {
        std::cerr << "Storage: Failed to open " << Filename << std::endl;
        return;
    }

    // ----------------- points -----------------
    const int numParticles = SYSTEM->generalParams.maxNodeCount;

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Point representing Sub_cellular elem model\n";
    ofs << "ASCII\n\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";

    ofs << "POINTS " << numParticles + 1 << " float\n";
    
    for (int i = 0; i < numParticles; ++i) {
        double xPos = sanitize_value(SYSTEM->coordInfoVecs.nodeLocX[i]);
        double yPos = sanitize_value(SYSTEM->coordInfoVecs.nodeLocY[i]);
        double zPos = sanitize_value(SYSTEM->coordInfoVecs.nodeLocZ[i]);
        ofs << std::setprecision(5) << std::fixed
            << xPos << " " << yPos << " " << zPos << "\n";
    }
    
    ofs << std::setprecision(5) << std::fixed
        << sanitize_value(SYSTEM->ljInfoVecs.LJ_PosX) << " "
        << sanitize_value(SYSTEM->ljInfoVecs.LJ_PosY) << " "
        << sanitize_value(SYSTEM->ljInfoVecs.LJ_PosZ) << "\n";

    // ----------------- cells -----------------
    const int numEdgesTotal = SYSTEM->coordInfoVecs.num_edges;

    auto edge_valid = [&](int e) {
        if (e < 0 || e >= numEdgesTotal) return false;
        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[e];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[e];
        return (idA != INT_MAX && idB != INT_MAX &&
                idA != -INT_MAX && idB != -INT_MAX &&
                idA >= 0 && idB >= 0 &&
                idA < numParticles && idB < numParticles);
    };

    int activeEdges = 0;
    for (int e = 0; e < numEdgesTotal; ++e) {
        if (edge_valid(e)) ++activeEdges;
    }

    const int numCells = activeEdges + 1;
    const int numNumsInCells = 3 * activeEdges + 2;

    ofs << "CELLS " << numCells << " " << numNumsInCells << "\n";

    for (int edge = 0; edge < numEdgesTotal; ++edge) {
        if (!edge_valid(edge)) continue;
        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];
        ofs << "2 " << idA << " " << idB << "\n";
    }
    ofs << "1 " << numParticles << "\n";

    ofs << "CELL_TYPES " << numCells << "\n";
    for (int i = 0; i < activeEdges; ++i) {
        ofs << "3\n";
    }
    ofs << "1\n";

    // ----------------- cell data -----------------
    ofs << "CELL_DATA " << numCells << "\n";
    ofs << "SCALARS Strain double 1\n";
    ofs << "LOOKUP_TABLE default\n";

    for (int edge = 0; edge < numEdgesTotal; ++edge) {
        if (!edge_valid(edge)) continue;

        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];

        double L0 = 1.0;
        if (edge < (int)SYSTEM->linearSpringInfoVecs.edge_rest_length.size()) {
            L0 = SYSTEM->linearSpringInfoVecs.edge_rest_length[edge];
        }
        if (L0 <= 0 || std::isnan(L0) || std::isinf(L0)) L0 = 1.0;

        double xL = sanitize_value(SYSTEM->coordInfoVecs.nodeLocX[idA]);
        double yL = sanitize_value(SYSTEM->coordInfoVecs.nodeLocY[idA]);
        double zL = sanitize_value(SYSTEM->coordInfoVecs.nodeLocZ[idA]);
        double xR = sanitize_value(SYSTEM->coordInfoVecs.nodeLocX[idB]);
        double yR = sanitize_value(SYSTEM->coordInfoVecs.nodeLocY[idB]);
        double zR = sanitize_value(SYSTEM->coordInfoVecs.nodeLocZ[idB]);

        double dx = xL - xR;
        double dy = yL - yR;
        double dz = zL - zR;
        double L1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (L1 <= 0 || std::isnan(L1) || std::isinf(L1)) L1 = L0;

        double strain = sanitize_value((L1 - L0) / L0);
        ofs << std::fixed << std::setprecision(6) << strain << "\n";
    }
    ofs << "0.0\n";

    ofs << "SCALARS SpringTension double 1\n";
    ofs << "LOOKUP_TABLE default\n";

    double k = SYSTEM->linearSpringInfoVecs.spring_constant;
    if (k <= 0 || std::isnan(k) || std::isinf(k)) k = 1.0;

    for (int edge = 0; edge < numEdgesTotal; ++edge) {
        if (!edge_valid(edge)) continue;

        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];

        double L0 = 1.0;
        if (edge < (int)SYSTEM->linearSpringInfoVecs.edge_rest_length.size()) {
            L0 = SYSTEM->linearSpringInfoVecs.edge_rest_length[edge];
        }
        if (L0 <= 0 || std::isnan(L0) || std::isinf(L0)) L0 = 1.0;

        double xL = sanitize_value(SYSTEM->coordInfoVecs.nodeLocX[idA]);
        double yL = sanitize_value(SYSTEM->coordInfoVecs.nodeLocY[idA]);
        double zL = sanitize_value(SYSTEM->coordInfoVecs.nodeLocZ[idA]);
        double xR = sanitize_value(SYSTEM->coordInfoVecs.nodeLocX[idB]);
        double yR = sanitize_value(SYSTEM->coordInfoVecs.nodeLocY[idB]);
        double zR = sanitize_value(SYSTEM->coordInfoVecs.nodeLocZ[idB]);

        double dx = xL - xR;
        double dy = yL - yR;
        double dz = zL - zR;
        double L1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (L1 <= 0 || std::isnan(L1) || std::isinf(L1)) L1 = L0;

        double tension = sanitize_value(std::fabs(k * (L1 - L0)));
        ofs << std::fixed << std::setprecision(6) << tension << "\n";
    }
    ofs << "0.0\n";

    ofs.close();
    std::cout << "Storage: Wrote " << Filename << " with " << activeEdges << " edges" << std::endl;
}

void Storage::storeVariables(void) {
    std::shared_ptr<System> SYSTEM = system.lock();
    if (!SYSTEM) return;

    iteration2 += 1;
    int digits = (iteration2 > 0) ? static_cast<int>(std::ceil(std::log10(static_cast<double>(iteration2 + 1)))) : 1;
    std::string format = ".sta";
    std::string Number;

    if (digits <= 1) {
        Number = "0000" + std::to_string(iteration2 - 1);
    } else if (digits == 2) {
        Number = "000" + std::to_string(iteration2 - 1);
    } else if (digits == 3) {
        Number = "00" + std::to_string(iteration2 - 1);
    } else {
        Number = "0" + std::to_string(iteration2 - 1);
    }

    std::string initial = "Variables_realistic/YB_cellwall4_newinitialmesh6_";
    std::string Filename = initial + Number + format;

    std::ofstream ofs(Filename.c_str());
    if (!ofs) return;

    unsigned numParticles = SYSTEM->generalParams.maxNodeCountLJ;
    ofs << std::setprecision(5) << std::fixed
        << "number of LJ particles " << numParticles << "\n";

    for (int i = 0; i < SYSTEM->generalParams.maxNodeCount; ++i) {
        double x = sanitize_value(SYSTEM->coordInfoVecs.nodeLocX[i]);
        double y = sanitize_value(SYSTEM->coordInfoVecs.nodeLocY[i]);
        double z = sanitize_value(SYSTEM->coordInfoVecs.nodeLocZ[i]);
        ofs << "<node>" << x << " " << y << " " << z << "</node>\n";
    }

    for (int i = 0; i < SYSTEM->coordInfoVecs.num_triangles; ++i) {
        ofs << "<elem> " << SYSTEM->coordInfoVecs.triangles2Nodes_1[i]+1 << " " 
            << SYSTEM->coordInfoVecs.triangles2Nodes_2[i]+1 << " " 
            << SYSTEM->coordInfoVecs.triangles2Nodes_3[i]+1 << " </elem>\n";
    }

    for (int i = 0; i < SYSTEM->coordInfoVecs.num_triangles; ++i) {
        ofs << "<elem2edge> " << SYSTEM->coordInfoVecs.triangles2Edges_1[i]+1 << " " 
            << SYSTEM->coordInfoVecs.triangles2Edges_2[i]+1 << " " 
            << SYSTEM->coordInfoVecs.triangles2Edges_3[i]+1 << " </elem2edge>\n";
    }

    for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; ++i) {
        ofs << "<edgeinfo> " << SYSTEM->coordInfoVecs.edges2Nodes_1[i]+1 << " " 
            << SYSTEM->coordInfoVecs.edges2Nodes_2[i]+1 << " </edgeinfo>\n";
    }

    for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; ++i) {
        ofs << "<edge2elem> " << SYSTEM->coordInfoVecs.edges2Triangles_1[i]+1 << " " 
            << SYSTEM->coordInfoVecs.edges2Triangles_2[i]+1 << " </edge2elem>\n";
    }

    ofs.close();
}
