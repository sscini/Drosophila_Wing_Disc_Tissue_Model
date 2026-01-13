#ifndef XMLPARSER_H_
#define XMLPARSER_H_

#include <string>
#include "SystemBuilder.h"
#include "SystemStructures.h"
#include "System.h"
#include "pugixml/include/pugixml.hpp"

// XMLParser encapsulates the logic for reading an XML file (such as Data_structure_circle_double_layer.xml)
// and populating the SystemBuilder with settings, nodes, edges (with optional rest lengths),
// elements, and other data.
class XMLParser {
public:
    // Parses the XML file and configures the builder.
    // Returns true if parsing succeeded, false otherwise.
    static bool parseFile(const std::string& filename, SystemBuilder& builder);
};

#endif  // XMLPARSER_H_
