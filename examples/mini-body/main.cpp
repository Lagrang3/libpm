#include "integration.hpp"
#include "parser.hpp"
#include "state.hpp"
#include "error.hpp"
#include <cmath>
#include <iostream>
#include <pm.hpp>

#include <boost/program_options.hpp>

state mystate;

auto comma_separated_to_vector(std::string text)
{
    std::vector<double> ret;

    while (text.size())
    {
        auto pos = std::find(text.begin(), text.end(), ',') - text.begin();
        auto first = text.substr(0, pos);
        ret.push_back(std::stod(first));
        text = text.substr(pos);
        if (text[0] == ',')
            text = text.substr(1);
    }

    return ret;
}

namespace po = boost::program_options;

int main(int argc, char** argv)
{
    
    try{
        // read command line options
        po::options_description usage("Allowed options");
        usage.add_options()
            ("help","produce help message")
            ("paramfile",po::value<std::string>(),"parameter file");
        
        po::variables_map vm;
        po::store(po::parse_command_line(argc,argv,usage),vm);
        po::notify(vm);
        
        if(vm.count("help") or vm.count("paramfile")==0)
        {
            std::cerr << usage << '\n';
            return 1;
        }
        
        
        // read paramfile
        ParamfileParser par_options( vm["paramfile"].as<std::string>() );

        mystate.BoxSize = std::stod(par_options["BoxSize"]);
        mystate.HubbleParam = std::stod(par_options["HubbleParam"]);
        mystate.Omega0 = std::stod(par_options["Omega0"]);
        mystate.zvect = comma_separated_to_vector(par_options["Redshifts"]);

        // read snapshot
        std::cout << "Reading snapshot file\n";
        mystate.read(par_options["Input"]);
        
        // run from z_begin to z_end
        std::cout << "Initializing grid\n";
        PM::grid<3, float, PM::CIC_filter, PM::CIC_filter> mygrid(
            std::stoi(par_options["PMGRID"]));

        // update force at the begining, for the KDK evolution
        std::cout << "Periodic boundary conditions\n";
        mystate.enforce_PBC();
        std::cout << "First force update\n";
        mygrid(mystate.F, mystate.X);

        const double afinal = 1 / (mystate.zvect.back() + 1);
        std::cout << "Main loop\n";
        while (mystate.a < afinal)
        {
            std::cout<< "    a: " << mystate.a 
                << ", final a: " << afinal << std::endl;
            double dt = 0.001 / mystate.get_H();

            leapfrogKDK_evolve(mygrid, mystate, dt);
            mystate.advance_a(dt);
            
            mystate.enforce_PBC();
        }
        std::cout << "done Main loop\n";

        // output snapshot
        std::cout << "Writing output snapshot\n";
        mystate.write(par_options["Output"]);
        
        return 0;
    }
    catch(std::exception& e)
    {
        std::cerr << "Exception! " << e.what() << '\n';
    }
    catch(...)
    {
        std::cerr << "Unknown exception!\n";
        
    }
    return 1;
}
