#include "System.h"
#include "Storage.h"
#include <cmath>        // for std::ceil, std::log10

Storage::Storage(std::weak_ptr<System> system_) {
    system = system_;

    std::ofstream statesOutput("Temp.sta");

    if (auto SYSTEM = system.lock()) {
        statesOutput << "node_count " << SYSTEM->coordInfoVecs.nodeLocX.size() << '\n';
        statesOutput << "edge_count " << SYSTEM->coordInfoVecs.num_edges      << '\n';
        statesOutput << "elem_count " << SYSTEM->coordInfoVecs.num_triangles  << '\n';
    }

    statesOutput.close();
}

void Storage::print_VTK_File(void) {

    std::shared_ptr<System> SYSTEM = system.lock();
    if (!SYSTEM) return;

    // ----------------- filename -----------------
    iteration += 1;
    int digits = static_cast<int>(std::ceil(std::log10(std::max(iteration, 1))));
    std::string format = ".vtk";
    std::string Number;

    if (digits <= 1) {
        Number = "0000" + std::to_string(iteration - 1);
    } else if (digits == 2) {
        Number = "000" + std::to_string(iteration - 1);
    } else if (digits == 3) {
        Number = "00" + std::to_string(iteration - 1);
    } else { // digits >= 4
        Number = "0" + std::to_string(iteration - 1);
    }

    // choose your base path/name here
    std::string initial  = "Animation_realistic/12_22_25/_clean_code_test_";
    std::string Filename = initial + Number + format;

    std::ofstream ofs(Filename.c_str());
    if (!ofs) return;

    // ----------------- points -----------------
    const int numParticles = SYSTEM->generalParams.maxNodeCount;

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Point representing Sub_cellular elem model\n";
    ofs << "ASCII\n\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";

    ofs << "POINTS " << numParticles + 1 << " float\n";
    for (int i = 0; i < numParticles; ++i) {
        double xPos = SYSTEM->coordInfoVecs.nodeLocX[i];
        double yPos = SYSTEM->coordInfoVecs.nodeLocY[i];
        double zPos = SYSTEM->coordInfoVecs.nodeLocZ[i];
        ofs << std::setprecision(5) << std::fixed
            << xPos << " " << yPos << " " << zPos << "\n";
    }
    ofs << std::setprecision(5) << std::fixed
        << SYSTEM->ljInfoVecs.LJ_PosX << " "
        << SYSTEM->ljInfoVecs.LJ_PosY << " "
        << SYSTEM->ljInfoVecs.LJ_PosZ << "\n";

    // ----------------- cells (edges + 1 LJ cell) -----------------
    const int numEdgesTotal = SYSTEM->coordInfoVecs.num_edges;

    auto edge_valid = [&](int e) {
        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[e];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[e];
        return (idA != INT_MAX && idB != INT_MAX &&
                idA != -INT_MAX && idB != -INT_MAX &&
                idA >= 0 && idB >= 0);
    };

    // first pass: count valid edges
    int activeEdges = 0;
    for (int e = 0; e < numEdgesTotal; ++e) {
        if (edge_valid(e)) ++activeEdges;
    }

    const int numCells       = activeEdges + 1;           // +1 for LJ particle
    const int numNumsInCells = 3 * activeEdges + 2;       // 3 per edge cell, 2 for LJ cell

    ofs << "CELLS " << numCells << " " << numNumsInCells << "\n";

    // second pass: write only valid edges
    for (int edge = 0; edge < numEdgesTotal; ++edge) {
        if (!edge_valid(edge)) continue;
        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];
        ofs << "2 " << idA << " " << idB << "\n";  // line cell with 2 points
    }

    // final cell: LJ point
    ofs << "1 " << SYSTEM->generalParams.maxNodeCount << "\n";

    // cell types
    ofs << "CELL_TYPES " << numCells << "\n";
    for (int i = 0; i < activeEdges; ++i) {
        ofs << "3\n";   // VTK_LINE
    }
    ofs << "1\n";       // VTK_VERTEX for LJ particle

//    // ----------------- cell data: strain per edge + LJ dummy -----------------
//    ofs << "CELL_DATA " << numCells << "\n";
//    ofs << "SCALARS Strain double\n";
//    ofs << "LOOKUP_TABLE default\n";
//
//    for (int edge = 0; edge < numEdgesTotal; ++edge) {
//        if (!edge_valid(edge)) continue;
//
//        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
//        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];
//
//        double L0 = SYSTEM->linearSpringInfoVecs.edge_rest_length[edge];
//        double xL = SYSTEM->coordInfoVecs.nodeLocX[idA];
//        double yL = SYSTEM->coordInfoVecs.nodeLocY[idA];
//        double zL = SYSTEM->coordInfoVecs.nodeLocZ[idA];
//        double xR = SYSTEM->coordInfoVecs.nodeLocX[idB];
//        double yR = SYSTEM->coordInfoVecs.nodeLocY[idB];
//        double zR = SYSTEM->coordInfoVecs.nodeLocZ[idB];
//
//        double L1 = std::sqrt((xL - xR)*(xL - xR) +
//                              (yL - yR)*(yL - yR) +
//                              (zL - zR)*(zL - zR));
//        double strain = (L1 - L0) / L0;
//        ofs << std::fixed << strain << "\n";
//    }
//
//    // dummy value for LJ cell
//    ofs << std::fixed << 0.1 << "\n";
//
//    ofs.close();
    
    // ----------------- cell data: strain per edge + LJ dummy -----------------
    ofs << "CELL_DATA " << numCells << "\n";
    
    // 1) Strain
    ofs << "SCALARS Strain double\n";
    ofs << "LOOKUP_TABLE default\n";
    
    std::vector<double> strain_values;
    strain_values.reserve(activeEdges);
    
    for (int edge = 0; edge < numEdgesTotal; ++edge) {
        if (!edge_valid(edge)) continue;
    
        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];
    
        double L0 = SYSTEM->linearSpringInfoVecs.edge_rest_length[edge]; // This is incorrect. Must look at what the original rest length was for each spring. 
        double xL = SYSTEM->coordInfoVecs.nodeLocX[idA];
        double yL = SYSTEM->coordInfoVecs.nodeLocY[idA];
        double zL = SYSTEM->coordInfoVecs.nodeLocZ[idA];
        double xR = SYSTEM->coordInfoVecs.nodeLocX[idB];
        double yR = SYSTEM->coordInfoVecs.nodeLocY[idB];
        double zR = SYSTEM->coordInfoVecs.nodeLocZ[idB];
    
        double L1 = std::sqrt((xL - xR)*(xL - xR) +
                              (yL - yR)*(yL - yR) +
                              (zL - zR)*(zL - zR));
        double strain = (L1 - L0) / L0;
    
        strain_values.push_back(strain);
        ofs << std::fixed << strain << "\n";
    }
    
    // dummy for LJ cell
    ofs << std::fixed << 0.1 << "\n";
    
    // 2) Spring tension proxy = |strain|
    ofs << "SCALARS SpringTension double\n";
    ofs << "LOOKUP_TABLE default\n";
    for (double edge = 0; edge < numEdgesTotal; ++edge){
        
        int idA = SYSTEM->coordInfoVecs.edges2Nodes_1[edge];
        int idB = SYSTEM->coordInfoVecs.edges2Nodes_2[edge];
    
        double L0 = SYSTEM->linearSpringInfoVecs.edge_rest_length[edge]; // This is incorrect. Must look at what the original rest length was for each spring. 
        double xL = SYSTEM->coordInfoVecs.nodeLocX[idA];
        double yL = SYSTEM->coordInfoVecs.nodeLocY[idA];
        double zL = SYSTEM->coordInfoVecs.nodeLocZ[idA];
        double xR = SYSTEM->coordInfoVecs.nodeLocX[idB];
        double yR = SYSTEM->coordInfoVecs.nodeLocY[idB];
        double zR = SYSTEM->coordInfoVecs.nodeLocZ[idB];
    
        double L1 = std::sqrt((xL - xR)*(xL - xR) +
                              (yL - yR)*(yL - yR) +
                              (zL - zR)*(zL - zR));
        
        //double tension = std::fabs(s);   // or keep sign if you like
        
        double k = SYSTEM->linearSpringInfoVecs.spring_constant;  // at a later point we will have modified the spring constants to be according to the individual springs. 
        double tension = std::fabs(k * (L1 - L0));

        ofs << std::fixed << tension << "\n";
    }
    // dummy for LJ cell
    ofs << std::fixed << 0.0 << "\n";
    
    ofs.close();

}

void Storage::storeVariables(void) {
    std::shared_ptr<System> SYSTEM = system.lock();
    if (!SYSTEM) return;

    iteration2 += 1;
    int digits = static_cast<int>(std::ceil(std::log10(std::max(iteration2, 1))));
    std::string format = ".sta";
    std::string Number;

    if (digits <= 1) {
        Number = "0000" + std::to_string(iteration2 - 1);
    } else if (digits == 2) {
        Number = "000" + std::to_string(iteration2 - 1);
    } else if (digits == 3) {
        Number = "00" + std::to_string(iteration2 - 1);
    } else { // digits >= 4
        Number = "0" + std::to_string(iteration2 - 1);
    }

    std::string initial  = "Variables_realistic/YB_cellwall4_newinitialmesh6_";
    std::string Filename = initial + Number + format;

    std::ofstream ofs(Filename.c_str());
    if (!ofs) return;

    unsigned numParticles = SYSTEM->generalParams.maxNodeCountLJ;
    ofs << std::setprecision(5) << std::fixed
        << "number of LJ particles " << numParticles << "\n";

    // nodes
    for (int i = 0; i < SYSTEM->generalParams.maxNodeCount; ++i) {
        double x = SYSTEM->coordInfoVecs.nodeLocX[i];
        double y = SYSTEM->coordInfoVecs.nodeLocY[i];
        double z = SYSTEM->coordInfoVecs.nodeLocZ[i];
        ofs << "<node>" << x << " " << y << " " << z << "</node>\n";
    }

    // triangles -> nodes
    for (int i = 0; i < SYSTEM->coordInfoVecs.num_triangles; ++i) {
        int t2n_1 = SYSTEM->coordInfoVecs.triangles2Nodes_1[i];
        int t2n_2 = SYSTEM->coordInfoVecs.triangles2Nodes_2[i];
        int t2n_3 = SYSTEM->coordInfoVecs.triangles2Nodes_3[i];
        ofs << "<elem> " << t2n_1+1 << " " << t2n_2+1 << " " << t2n_3+1 << " </elem>\n";
    }

    // triangles -> edges
    for (int i = 0; i < SYSTEM->coordInfoVecs.num_triangles; ++i) {
        int e1 = SYSTEM->coordInfoVecs.triangles2Edges_1[i];
        int e2 = SYSTEM->coordInfoVecs.triangles2Edges_2[i];
        int e3 = SYSTEM->coordInfoVecs.triangles2Edges_3[i];
        ofs << "<elem2edge> " << e1+1 << " " << e2+1 << " " << e3+1 << " </elem2edge>\n";
    }

    // edges -> nodes
    for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; ++i) {
        int n1 = SYSTEM->coordInfoVecs.edges2Nodes_1[i];
        int n2 = SYSTEM->coordInfoVecs.edges2Nodes_2[i];
        ofs << "<edgeinfo> " << n1+1 << " " << n2+1 << " </edgeinfo>\n";
    }

    // edges -> triangles
    for (int i = 0; i < SYSTEM->coordInfoVecs.num_edges; ++i) {
        int t1 = SYSTEM->coordInfoVecs.edges2Triangles_1[i];
        int t2 = SYSTEM->coordInfoVecs.edges2Triangles_2[i];
        ofs << "<edge2elem> " << t1+1 << " " << t2+1 << " </edge2elem>\n";
    }

    ofs.close();
}
