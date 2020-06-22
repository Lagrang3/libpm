#include "integration.hpp"
#include "parser.hpp"
#include "state.hpp"
#include <cmath>
#include <iostream>
#include <pm.hpp>

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

int main(int argc, char** argv)
{
    // read command line options
    CmdParser cmd_options(argc, argv);

    // read paramfile
    std::string paramfile = cmd_options["paramfile"];
    ParamfileParser par_options(paramfile);

    mystate.BoxSize = std::stod(par_options["BoxSize"]);
    mystate.HubbleParam = std::stod(par_options["HubbleParam"]);
    mystate.Omega0 = std::stod(par_options["Omega0"]);
    mystate.zvect = comma_separated_to_vector(par_options["Redshifts"]);

    // read snapshot
    mystate.read(par_options["Input"]);

    // run from z_begin to z_end
    PM::grid<3, float, PM::Gaussian_filter, PM::Gaussian_filter> mygrid(
        std::stoi(par_options["PMGRID"]));

    // update force at the begining, for the KDK evolution
    mystate.enforce_PBC();
    mygrid(mystate.F, mystate.X);

    const double afinal = 1 / (mystate.zvect.back() + 1);
    while (mystate.a < afinal)
    {
        double dt = 0.001 / mystate.get_H();

        leapfrogKDK_evolve(mygrid, mystate, dt);
        mystate.advance_a(dt);
    }

    // output snapshot
    mystate.write(par_options["Output"]);

    return 0;
}
