/*
Types of nodes to include in double layer model data structure:

 - nodes - apical and basal
 - edgeinfo - apical and basal
 - elements
 - elem2edge
 - edge2elem
 
 The new data structure has been created. 
 
 I think the next step here is printing out the edge lengths with whether they are apical basal or otherwise and checking to see if they are truly at rest. Esp the vertical ones. Those seem to have the most tension at the start. Need to get this equilibrium structure down. 

*/ 
#include "XMLParser.h"
#include "System.h"
#include <iostream>
#include <cstdio>

bool XMLParser::parseFile(const std::string& filename, SystemBuilder& builder) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    
    std::cout<<"data struct file name = "<< filename.c_str()<<std::endl;
    if (!result) {
        std::cerr << "XMLParser: Error parsing file: " << result.description() << std::endl;
        return false;
    }
    
    pugi::xml_node root = doc.child("data");
    if (!root) {
        std::cerr << "XMLParser: No <data> root found." << std::endl;
        return false;
    }
    
        // ---------- Count nodes per group BEFORE building ----------
    pugi::xml_node apicalNodes = root.child("apicalNodes");
    pugi::xml_node bodyNodes   = root.child("bodyNodes");
    pugi::xml_node basalNodes  = root.child("basalNodes");
    
    auto countChildren = [](pugi::xml_node parent, const char* tag) {
        int c = 0;
        if (parent) for (pugi::xml_node n = parent.child(tag); n; n = n.next_sibling(tag)) ++c;
        return c;
    };
    
    const int A = countChildren(apicalNodes, "node"); // apical count per layer
    const int B = countChildren(bodyNodes,   "node"); // all body nodes (across all body sublayers)
    const int C = countChildren(basalNodes,  "node"); // basal count per layer
    
    if (A == 0) {
        std::cerr << "XMLParser: No <apicalNodes> or empty." << std::endl;
        return false;
    }
    if (C != 0 && C != A) {
        std::cerr << "XMLParser: Basal nodes (" << C << ") must equal apical nodes (" << A << ")." << std::endl;
        return false;
    }
    
    const int totalNodes = A + B + C;
    
    // Your formula: body layer count from apical count
    int N_body = 0;
    if (B > 0) {
        if (B % A != 0) {
            std::cerr << "XMLParser: bodyNodes (" << B << ") not a multiple of apicalNodes (" << A << ")." << std::endl;
            return false;
        }
        N_body = B / A;  // == (totalNodes - 2*A)/A as you specified
    }
    const int N_layers = N_body + 2; // Basal..Body..Apical
    
    builder.defaultLayers = N_layers;
    
    std::cout << "[Layers] A=" << A << " B=" << B << " C=" << C
              << " -> N_body=" << N_body << " total layers=" << N_layers << std::endl;
    
    
    // ---------------------------// currently all my constants are taken as the default ones and not the ones from the data structure. Should'nt be the case. I need to add a condition where they are taking in the data structure constants also. 
    // Process <settings> node
    pugi::xml_node settings = root.child("settings");
    if (settings) {
        if (auto p = settings.child("Tau")) {
            builder.defaultTau = p.text().as_double();
            std::cout << "Setting tau: " << builder.defaultTau << std::endl;
        }
        if (auto p = settings.child("KBT")) {
            builder.defaultKBT = p.text().as_double();
            std::cout << "Setting kbt: " << builder.defaultKBT << std::endl;
        }
        if (auto p = settings.child("Linear_Const")) {
            builder.defaultLinear_Const = p.text().as_double();
            std::cout << "Setting linear const: " << builder.defaultLinear_Const << std::endl;
        }
        if (auto p = settings.child("dt")) {
            builder.defaultdt = p.text().as_double();
            std::cout << "Setting dt: "<< builder.defaultdt << std::endl;
        }
        if (auto p = settings.child("Tf")) {
            builder.defaultTf = p.text().as_double();
            std::cout <<"Setting Tf:" << builder.defaultTf << std::endl;
        }
        if (auto p = settings.child("lambda_iso_in_DV_center")){
            builder.defaultlambda_iso_in_DV_center = p.text().as_double();
            std::cout << "Setting lambda_iso_in_DV_center:" << builder.defaultlambda_iso_in_DV_center << std::endl;
        }
        if (auto p = settings.child("lambda_aniso_in_DV_center")){
            builder.defaultlambda_aniso_in_DV_center = p.text().as_double();
            std::cout << "Setting lambda_aniso_in_DV_center:" << builder.defaultlambda_aniso_in_DV_center << std::endl;
        }
        if (auto p = settings.child("lambda_iso_in_DV_edge")){
            builder.defaultlambda_iso_in_DV_edge = p.text().as_double();
            std::cout << "Setting lambda_iso_in_DV_edge:" << builder.defaultlambda_iso_in_DV_edge << std::endl;
        }
        if (auto p = settings.child("lambda_aniso_in_DV_edge")){
            builder.defaultlambda_aniso_in_DV_edge = p.text().as_double();
            std::cout << "Setting lambda_aniso_in_DV_edge:" << builder.defaultlambda_aniso_in_DV_edge << std::endl;
        }
        if (auto p = settings.child("lambda_iso_out_DV_center")){
            builder.defaultlambda_iso_out_DV_center = p.text().as_double();
            std::cout << "Setting lambda_iso_out_DV_center:" << builder.defaultlambda_iso_out_DV_center << std::endl;
        }
        if (auto p = settings.child("lambda_aniso_out_DV_center")){
            builder.defaultlambda_aniso_out_DV_center = p.text().as_double();
            std::cout << "Setting lambda_aniso_out_DV_center:" << builder.defaultlambda_aniso_out_DV_center << std::endl;
        }
        if (auto p = settings.child("lambda_iso_out_DV_edge")){
            builder.defaultlambda_iso_out_DV_edge = p.text().as_double();
            std::cout << "Setting lambda_iso_out_DV_edge:" << builder.defaultlambda_iso_out_DV_edge << std::endl;
        }
        if (auto p = settings.child("lambda_ansio_out_DV_edge")){
            builder.defaultlambda_aniso_out_DV_edge = p.text().as_double();
            std::cout << "Setting lambda_aniso_out_DV_edge:" << builder.defaultlambda_aniso_out_DV_edge << std::endl;
        }
        if (auto p = settings.child("tol")){
            builder.defaulttol = p.text().as_double();
            std::cout << "Setting tol:" << builder.defaulttol << std::endl;
        }
        if (auto p = settings.child("Rmin")){
            builder.defaultRmin = p.text().as_double();
            std::cout << "Setting Rmin:" << builder.defaultRmin << std::endl;
        }
        if (auto p = settings.child("Area_Const")) {
            builder.defaultArea_Const = p.text().as_double();
            std::cout << "Setting area const: " << builder.defaultArea_Const << std::endl;
        }
        if (auto p = settings.child("Bend_Const")) {
            builder.defaultBending_Const = p.text().as_double();
            std::cout << "Setting bending const: " << builder.defaultBending_Const << std::endl;
        }
        if (auto p = settings.child("LJ_Eps")) {// not needed
            builder.defaultLJ_Eps = p.text().as_double();
            std::cout << "Setting lj eps: " << builder.defaultLJ_Eps << std::endl;
        }
        if (auto p = settings.child("LJ_Rmin")) {// not needed
            builder.defaultLJ_Rmin = p.text().as_double();
            std::cout << "Setting lj rmin: " << builder.defaultLJ_Rmin << std::endl;
        }
        if (auto p = settings.child("LJ_Rmax")) {// not needed
            builder.defaultLJ_Rmax = p.text().as_double();
            std::cout << "Setting lj rmax: " << builder.defaultLJ_Rmax << std::endl;
        }
        if (auto p = settings.child("LJ_Const")) {// not needed
            builder.defaultLJ_Const = p.text().as_double();
            std::cout << "Setting lj const: " << builder.defaultLJ_Const << std::endl;
        }
        if (auto p = settings.child("LJ_X")) {// not needed
            builder.defaultLJ_X = p.text().as_double();
            std::cout << "Setting lj x: " << builder.defaultLJ_X << std::endl;
        }
        if (auto p = settings.child("LJ_Y")) {// not needed
            builder.defaultLJ_Y = p.text().as_double();
            std::cout << "Setting lj y: " << builder.defaultLJ_Y << std::endl;
        }
        if (auto p = settings.child("LJ_Z")) {// not needed
            builder.defaultLJ_Z = p.text().as_double();
            std::cout << "Setting lj z: " << builder.defaultLJ_Z << std::endl;
        }
        //Additional settings:
        // Initial Bending angle:
        // Iso strain 
        // Aniso strain
        // 
        }
    
    // ---------------------------
    // Process <nodes> section
        
    
    // Helper: map a 1-based global node id to its layer index under the new convention.
    // Ordering is assumed Apical [1..A], then N_body blocks of A nodes, then Basal [last A].
    auto layerOfNode1Based = [A, N_body, totalNodes](unsigned gid) -> int {
        // Safety
        if (gid == 0 || (int)gid > totalNodes) return -1;
    
        const int apicalStart = 1;
        const int apicalEnd   = A;
        const int bodyStart   = A + 1;
        const int bodyEnd     = A + N_body * A;
        const int basalStart  = bodyEnd + 1;
        const int basalEnd    = bodyEnd + A;
    
        if ((int)gid >= apicalStart && (int)gid <= apicalEnd) {
            return N_body + 1;         // Apical
        } else if ((int)gid >= basalStart && (int)gid <= basalEnd) {
            return 0;                  // Basal
        } else if ((int)gid >= bodyStart && (int)gid <= bodyEnd) {
            // Which body sublayer? 1..N_body
            int offset = (int)gid - bodyStart;       // 0-based into body block
            return (offset / A) + 1;                 // Body layer index (1..N_body)
        }
        return -1;
    };
    
    // Small helper for coordinates parsing
    auto parseXYZ = [](const char* text, double& x, double& y, double& z) -> bool {
        return 3 == std::sscanf(text, "%lf %lf %lf", &x, &y, &z);
    };
    // ---------- Build nodes with correct layer flags ----------
    if (apicalNodes) {
        for (pugi::xml_node node = apicalNodes.child("node"); node; node = node.next_sibling("node")) {
            double x, y, z;
            if (!parseXYZ(node.text().as_string(), x, y, z)) {
                std::cerr << "XMLParser: parse apical node error" << std::endl; return false;
            }
            builder.addNode(x, y, z);
            builder.addLayerFlag_node(N_body + 1); // Apical layer = N_body+1
        }
        std::cout << "Apical nodes parsed successfully." << std::endl;
    }
    
    if (bodyNodes) {
        int idxWithinBody = 0; // how many body nodes seen
        for (pugi::xml_node node = bodyNodes.child("node"); node; node = node.next_sibling("node"), ++idxWithinBody) {
            double x, y, z;
            if (!parseXYZ(node.text().as_string(), x, y, z)) {
                std::cerr << "XMLParser: parse body node error" << std::endl; return false;
            }
            // Body sublayer index: 1..N_body
            int bodyLayer = (idxWithinBody / A) + 1;
            builder.addNode(x, y, z);
            builder.addLayerFlag_node(bodyLayer);
        }
        std::cout << "Body nodes parsed successfully." << std::endl;
    } else if (B > 0) {
        std::cerr << "XMLParser: Declared body nodes implied by counts, but <bodyNodes> missing." << std::endl;
        return false;
    }
    
    if (basalNodes) {
        for (pugi::xml_node node = basalNodes.child("node"); node; node = node.next_sibling("node")) {
            double x, y, z;
            if (!parseXYZ(node.text().as_string(), x, y, z)) {
                std::cerr << "XMLParser: parse basal node error" << std::endl; return false;
            }
            builder.addNode(x, y, z);
            builder.addLayerFlag_node(0); // Basal
        }
        std::cout << "Basal nodes parsed successfully." << std::endl;
    }

     else {
        // Single-layer model. Made layerflag 1 because the strain calculations all take place within the upperhem. 
        pugi::xml_node nodes = root.child("nodes");
        if (nodes) {
            for (pugi::xml_node node = nodes.child("node"); node; node = node.next_sibling("node")) {
                double x, y, z;
                int layerflag = 1;
                const char* text = node.text().as_string();
                if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
                    std::cerr << "XMLParser: parse node error" << std::endl;
                    return false;
                }
                builder.addNode(x, y, z);
                builder.addLayerFlag_node(layerflag);
            }
            std::cout << "Nodes parsed successfully." << std::endl;
        }
    }
    
    
    // ---------------------------
    // Process <edgeinfos> section
    // Everything below here will have apical, body, basal and vertical parts. 
    auto addEdgeWithLayer = [&](unsigned from1, unsigned to1, double* rest) {
        // from1, to1 are 1-based in the XML; builder uses 0-based
        const int L1 = layerOfNode1Based(from1);
        const int L2 = layerOfNode1Based(to1);
        int edgeLayer = (L1 == L2 ? L1 : -1); // vertical if crossing layers
    
        if (rest) builder.addEdge(from1 - 1, to1 - 1, *rest);
        else      builder.addEdge(from1 - 1, to1 - 1);
        builder.addLayerFlag_edge(edgeLayer);
    };
    
    // ---------- Apical edgeinfos ----------
    if (pugi::xml_node apical_edgeinfos = root.child("apical_edgeinfos")) {
        for (pugi::xml_node edge = apical_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to; double rl;
            int n = std::sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &rl);
            if (n == 3) addEdgeWithLayer(from, to, &rl);
            else if (n == 2) addEdgeWithLayer(from, to, nullptr);
            else { std::cerr << "XMLParser: parse apical edge error\n"; return false; }
        }
        std::cout << "Apical Edges parsed successfully." << std::endl;
    }
    
    // ---------- Body edgeinfos (may include edges for multiple body sublayers) ----------
    if (pugi::xml_node body_edgeinfos = root.child("body_edgeinfos")) {
        for (pugi::xml_node edge = body_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to; double rl;
            int n = std::sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &rl);
            if (n == 3) addEdgeWithLayer(from, to, &rl);
            else if (n == 2) addEdgeWithLayer(from, to, nullptr);
            else { std::cerr << "XMLParser: parse body edge error\n"; return false; }
        }
        std::cout << "Body Edges parsed successfully." << std::endl;
    }
    
    // ---------- Basal edgeinfos ----------
    if (pugi::xml_node basal_edgeinfos = root.child("basal_edgeinfos")) {
        for (pugi::xml_node edge = basal_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to; double rl;
            int n = std::sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &rl);
            if (n == 3) addEdgeWithLayer(from, to, &rl);
            else if (n == 2) addEdgeWithLayer(from, to, nullptr);
            else { std::cerr << "XMLParser: parse basal edge error\n"; return false; }
        }
        std::cout << "Basal Edges parsed successfully." << std::endl;
    }
    
    // ---------- Vertical edgeinfos (explicit vertical list gets layer -1) ----------
    if (pugi::xml_node vertical_edgeinfos = root.child("vertical_edgeinfos")) {
        for (pugi::xml_node edge = vertical_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to; double rl;
            int n = std::sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &rl);
            int vertical = -1;
            if (n == 3) { builder.addEdge(from - 1, to - 1, rl); builder.addLayerFlag_edge(vertical); }
            else if (n == 2) { builder.addEdge(from - 1, to - 1); builder.addLayerFlag_edge(vertical); }
            else { std::cerr << "XMLParser: parse vertical edge error\n"; return false; }
        }
        std::cout << "Vertical Edges parsed successfully." << std::endl;
    }

    else{
        pugi::xml_node edgeinfos = root.child("edgeinfos");
        if (edgeinfos) {
            for (pugi::xml_node edge = edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
                unsigned int from, to;
                double restLength;
                int layerflag = 1;
                int count = sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &restLength);
                if (count == 3) {
                    builder.addEdge(from - 1, to - 1, restLength);
                    builder.addLayerFlag_edge(layerflag);
                } else if (count == 2) {
                    builder.addEdge(from - 1, to - 1);
                    builder.addLayerFlag_edge(layerflag);
                } else {
                    std::cerr << "XMLParser: parse edge error" << std::endl;
                    return false;
                }
            }
            std::cout << "Edges parsed successfully." << std::endl;
        }
    }
    
    // ---------------------------
    // Process <elems> section (elements/triangles)
    pugi::xml_node elems = root.child("elems");
    if (elems) {
        for (pugi::xml_node elem = elems.child("elem"); elem; elem = elem.next_sibling("elem")) {
            unsigned int a, b, c;
            if (3 != sscanf(elem.text().as_string(), "%u %u %u", &a, &b, &c)) {
                std::cerr << "XMLParser: parse elem error" << std::endl;
                return false;
            }
            if(a-1 != b-1 + 7 || a-1 != c-1 + 7 || b-1 != c-1 + 7){
                builder.addElement(a - 1, b - 1, c - 1); // This needs to check layerflag as well. If it is a vertical edge it cannot have the triangle added.
            }// FIX THIS 
             
        }
        std::cout << "Elements parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Process <elem2edges> mapping, if present.
    pugi::xml_node elem2edges = root.child("elem2edges");
    if (elem2edges) {
        for (pugi::xml_node mapping = elem2edges.child("elem2edge"); mapping; mapping = mapping.next_sibling("elem2edge")) {
            unsigned int e1, e2, e3;
            if (3 != sscanf(mapping.text().as_string(), "%u %u %u", &e1, &e2, &e3)) {
                std::cerr << "XMLParser: parse elem2edge error" << std::endl;
                return false;
            }
            builder.addElement2Edge(e1 - 1, e2 - 1, e3 - 1); // This one can have vertical edges included in this. 
        }
        std::cout << "Element-to-edge mappings parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Process <edge2elems> mapping, if present.
    pugi::xml_node edge2elems = root.child("edge2elems");
    if (edge2elems) {
        for (pugi::xml_node mapping = edge2elems.child("edge2elem"); mapping; mapping = mapping.next_sibling("edge2elem")) {
            unsigned int elemID, edgeID;
            if (2 != sscanf(mapping.text().as_string(), "%u %u", &elemID, &edgeID)) {
                std::cerr << "XMLParser: parse edge2elem error" << std::endl;
                return false;
            }
            builder.addEdge2Elem(elemID - 1, edgeID - 1); // This one can have vertical edges included in it. 
        }
        std::cout << "Edge-to-element mappings parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Optionally process other sections (capsidnodes, fixed nodes, nndatas, etc.)
    pugi::xml_node capsidnodes = root.child("capsidnodes");
    if (capsidnodes) {
        for (pugi::xml_node node = capsidnodes.child("node"); node; node = node.next_sibling("node")) {
            double x, y, z;
            if (3 != sscanf(node.text().as_string(), "%lf %lf %lf", &x, &y, &z)) {
                std::cerr << "XMLParser: parse capsid node error" << std::endl;
                return false;
            }
            builder.addCapsidNode(x, y, z);
        }
        std::cout << "Capsid nodes parsed successfully." << std::endl;
    }
    
    pugi::xml_node fixnodes = root.child("fixed");
    if (fixnodes) {
        for (pugi::xml_node node = fixnodes.child("node"); node; node = node.next_sibling("node")) {
            int id = node.text().as_int();
            builder.fixNodes(id - 1);
        }
        std::cout << "Fixed nodes parsed successfully." << std::endl;
    }
    
    pugi::xml_node nndatas = root.child("nndatas");
    if (nndatas) {
        for (pugi::xml_node nndata = nndatas.child("nndata"); nndata; nndata = nndata.next_sibling("nndata")) {
            double d1, d2, d3, d4, d5, d6, d7, d8, d9;
            if (9 != sscanf(nndata.text().as_string(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                            &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9)) {
                std::cerr << "XMLParser: parse nndata error" << std::endl;
                return false;
            }
            builder.addNndata(d1, d2, d3, d4, d5, d6, d7, d8, d9);
        }
        std::cout << "nndata parsed successfully." << std::endl;
    }
    
    return true;
}
