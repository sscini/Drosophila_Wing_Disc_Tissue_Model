#include <iomanip>
#include <string>
#include <memory>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <inttypes.h>
#include <cstddef>

#include "System.h"
#include "SystemBuilder.h"
#include "Storage.h"
#include "pugixml/include/pugixml.hpp"
#include "XMLParser.h"   // new header for XML parsing

// Updated createSystem: simply calls XMLParser::parseFile to configure the builder,
// then creates and returns the System.
std::shared_ptr<System> createSystem(const char* schemeFile, std::shared_ptr<SystemBuilder> builder) {
    // Use XMLParser to populate the builder with settings, nodes, edges, etc.
    if (!XMLParser::parseFile(schemeFile, *builder)) {
        std::cout << "Error parsing XML file: " << schemeFile << std::endl;
        return nullptr;
    }
    // Create the system based on builder data.
    std::shared_ptr<System> virus_system = builder->createSystem();
    return virus_system;
}

std::string generateOutputFileName(std::string inputFileName)
{
    time_t now;
    const int MAX_DATE = 64;
    char theDate[MAX_DATE];
    theDate[0] = '\0';

    now = time(nullptr);

    if (now != -1) {
        strftime(theDate, MAX_DATE, "_%Y.%m.%d_%H-%M-%S", gmtime(&now));
        return inputFileName + theDate;
    }
    return "";
}

void run(int argc, char** argv) {
    time_t t0, t1;
    t0 = time(0);

    double forceStep = 0.0;
    double timestep = 0.001;
    int solve_time = 10000;
    bool time_found = true;

    for (int i = -1; i < argc - 1; i++) {
        std::string arg = argv[i];
        int pos = arg.find('=');
        std::string key = arg.substr(0, pos);
        std::string val = arg.substr(pos + 1);

        std::cout << "argc: " << argc << std::endl;
        std::cout << "arg: " << arg << std::endl;
        std::cout << "pos: " << pos << std::endl;
        std::cout << "key: " << key << std::endl;
        std::cout << "val: " << val << std::endl;

        if (key == "-dt") {
            time_found = true;
            timestep = std::atof(val.c_str());
            std::cout << "setting timestep: " << timestep << std::endl;
            continue;
        }
        if (key == "-solve_time") {
            solve_time = std::atof(val.c_str());
            std::cout << "setting solve time: " << solve_time << std::endl;
            continue;
        }
    }

    auto builder = std::make_shared<SystemBuilder>(timestep, solve_time);

    // Use the last argument as the XML scheme file name.
    std::shared_ptr<System> system = createSystem(argv[argc - 1], builder);
    if (!system) {
        std::cerr << "Failed to create system." << std::endl;
        return;
    }

    auto storage = std::make_shared<Storage>(system);
    system->assignStorage(storage);
    system->set_weak_builder(builder);

    std::cout << "solving system in main" << std::endl;
    system->solveSystem();

    t1 = time(0);  // current time at the end of solving the system.
    int total, hours, min, sec;
    total = difftime(t1, t0);
    hours = total / 3600;
    min = (total % 3600) / 60;
    sec = (total % 3600) % 60;
    std::cout << "Total time hh: " << hours << " mm:" << min << " ss:" << sec << "\n";
}

int main(int argc, char** argv)
{
    std::cout << argc << std::endl;
    run(argc - 2, argv + 2);
    return 0;
}
