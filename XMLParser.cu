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
    // This will be divided into apical and basal inputs.
    
    // Check if this is a double-layer model (using <apicalNodes> and <basalNodes>) or single layer (<nodes>)
    pugi::xml_node apicalNodes = root.child("apicalNodes");
    if (apicalNodes) {
        // Double-layered model
        for (pugi::xml_node node = apicalNodes.child("node"); node; node = node.next_sibling("node")) {
            double x, y, z;
            int layerflag = 1;
            const char* text = node.text().as_string();
            if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
                std::cerr << "XMLParser: parse apical node error" << std::endl;
                return false;
            }
            // Assume SystemBuilder has an addApicalNode method.
            builder.addNode(x, y, z);
            builder.addLayerFlag_node(layerflag);
        }
        std::cout << "Apical nodes parsed successfully." << std::endl;
        
        pugi::xml_node bodyNodes = root.child("bodyNodes");
        if (bodyNodes) {
            for (pugi::xml_node node = bodyNodes.child("node"); node; node = node.next_sibling("node")) {
                double x, y, z;
                int layerflag = 0;
                const char* text = node.text().as_string();
                if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
                    std::cerr << "XMLParser: parse body node error" << std::endl;
                    return false;
                }
                // Assume SystemBuilder has an addBodyNode method.
                builder.addNode(x, y, z);
                builder.addLayerFlag_node(layerflag);
            }
            std::cout << "Body nodes parsed successfully." << std::endl;
        } else {
            std::cerr << "XMLParser: No <bodyNodes> found in double-layered model." << std::endl;
            return false;
        }
        
        pugi::xml_node basalNodes = root.child("basalNodes");
        if (basalNodes) {
            for (pugi::xml_node node = basalNodes.child("node"); node; node = node.next_sibling("node")) {
                double x, y, z;
                int layerflag = -1;
                const char* text = node.text().as_string();
                if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
                    std::cerr << "XMLParser: parse basal node error" << std::endl;
                    return false;
                }
                // Assume SystemBuilder has an addBasalNode method.
                builder.addNode(x, y, z);
                builder.addLayerFlag_node(layerflag);
            }
            std::cout << "Basal nodes parsed successfully." << std::endl;
        } else {
            std::cerr << "XMLParser: No <basalNodes> found in double-layered model." << std::endl;
            return false;
        }
    } else {
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
    // Everything below here will have apical, basal and vertical parts. 
    
    pugi::xml_node apical_edgeinfos = root.child("apical_edgeinfos");
    if (apical_edgeinfos) {
        for (pugi::xml_node edge = apical_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
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
                std::cerr << "XMLParser: parse apical edge error" << std::endl;
                return false;
            }
        }
        std::cout << "Apical Edges parsed successfully." << std::endl;
    }

    pugi::xml_node body_edgeinfos = root.child("body_edgeinfos");
    if (body_edgeinfos) {
        for (pugi::xml_node edge = body_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to;
            double restLength;
            int layerflag = 0;
            int count = sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &restLength);
            if (count == 3) {
                builder.addEdge(from - 1, to - 1, restLength);
                builder.addLayerFlag_edge(layerflag);
            } else if (count == 2) {
                builder.addEdge(from - 1, to - 1);
                builder.addLayerFlag_edge(layerflag);
            } else {
                std::cerr << "XMLParser: parse body edge error" << std::endl;
                return false;
            }
        }
        std::cout << "Body Edges parsed successfully." << std::endl;
    }
    
    pugi::xml_node basal_edgeinfos = root.child("basal_edgeinfos");
    if (basal_edgeinfos) {
        for (pugi::xml_node edge = basal_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to;
            double restLength;
            int layerflag = -1;
            int count = sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &restLength);
            if (count == 3) {
                builder.addEdge(from - 1, to - 1, restLength);
                builder.addLayerFlag_edge(layerflag);
            } else if (count == 2) {
                builder.addEdge(from - 1, to - 1);
                builder.addLayerFlag_edge(layerflag);
            } else {
                std::cerr << "XMLParser: parse basal edge error" << std::endl;
                return false;
            }
        }
        std::cout << "Basal Edges parsed successfully." << std::endl;
    }
    
    pugi::xml_node vertical_edgeinfos = root.child("vertical_edgeinfos");
    if (vertical_edgeinfos) {
        for (pugi::xml_node edge = vertical_edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to;
            double restLength;
            int layerflag = 2;
            int count = sscanf(edge.text().as_string(),"%u %u %lf",&from,&to,&restLength);
            if (count == 3) {
                builder.addEdge(from - 1, to - 1, restLength);
                builder.addLayerFlag_edge(layerflag);
            } else if (count == 2) {
               //builder.addEdge(from -1, to -1);
              // std::cout<< "from = "<< from -1<< " to = "<<to -1 <<std::endl; 
                builder.addEdge(from - 1, to - 1);
                builder.addLayerFlag_edge(layerflag);
            } else {
                std::cerr << "XMLParser: parse vertical edge error" << std::endl;
                return false;
            }
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
            builder.addElement(a - 1, b - 1, c - 1);
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
            builder.addElement2Edge(e1 - 1, e2 - 1, e3 - 1);
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
            builder.addEdge2Elem(elemID - 1, edgeID - 1);
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
