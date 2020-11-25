#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <pm.hpp>
#include <pm/field.hpp>
#include <pm/filters.hpp>

#include <cmath>
#include <fstream>
#include <string>
#include <unistd.h>

namespace mpi = boost::mpi;

auto Triangle_shape = [](double z) { return 3 * std::min(z, 1 - z); };
auto Sine_shape = [](double z) {
    static const double pi = acos(-1.0);
    return sin(2 * pi * z) + 3 * sin(2 * pi * z * 2);
};

template <class filter_t, class callable>
void write_pts(mpi::communicator com,
               int N_grid,
               int N_eval,
               callable F,
               std::string fname,
               std::string shape)
{
    if (com.rank() == 0)
    {
        std::ofstream fd(shape + "_exact.dat");
        for (int i = 0; i < N_eval; ++i)
        {
            double x = i * 1.0 / N_eval;
            fd << x << " " << F(x) << '\n';
        }
    }
    if (com.rank() == 0)
    {
        std::ofstream fd(shape + "_points.dat");
        for (int i = 0; i < N_grid; ++i)
        {
            double x = i * 1.0 / N_grid;
            fd << x << " " << F(x) << '\n';
        }
    }

    PM::Field<double, filter_t, filter_t> phi(com, N_grid, {1, 1});

    for (int i = 0; i < phi.extents()[0]; ++i)
        for (int j = 0; j < phi.extents()[1]; ++j)
            for (int k = 0; k < phi.extents()[2]; ++k)
                phi(i, j, k) = F(k * 1.0 / N_grid);
    phi.update_ghosts();

    if (com.rank() == 0)
    {
        std::ofstream fd(shape + "_" + fname + ".dat");
        int i = 0, j = 0;
        for (int k = 0; k < N_eval; ++k)
        {
            double x = i, y = j, z = k * (N_grid * 1.0 / N_eval);
            fd << z / N_grid << " " << phi.interpolate(x, y, z) << '\n';
        }
    }
}

int main(int narg, char** args)
{
    chdir(args[1]);
    mpi::environment env;
    mpi::communicator world;

    const int N_grid = 20, N_eval = 2000;
    write_pts<PM::filters::CIC>(world, N_grid, N_eval, Triangle_shape, "cic",
                                "triangle");
    write_pts<PM::filters::NGP>(world, N_grid, N_eval, Triangle_shape, "ngp",
                                "triangle");
    write_pts<PM::filters::TSC>(world, N_grid, N_eval, Triangle_shape, "tsc",
                                "triangle");
    write_pts<PM::filters::PCS>(world, N_grid, N_eval, Triangle_shape, "pcs",
                                "triangle");
    write_pts<PM::filters::Gaussian>(world, N_grid, N_eval, Triangle_shape,
                                     "gaussian", "triangle");

    write_pts<PM::filters::CIC>(world, N_grid, N_eval, Sine_shape, "cic",
                                "sine");
    write_pts<PM::filters::NGP>(world, N_grid, N_eval, Sine_shape, "ngp",
                                "sine");
    write_pts<PM::filters::TSC>(world, N_grid, N_eval, Sine_shape, "tsc",
                                "sine");
    write_pts<PM::filters::PCS>(world, N_grid, N_eval, Sine_shape, "pcs",
                                "sine");
    write_pts<PM::filters::Gaussian>(world, N_grid, N_eval, Sine_shape,
                                     "gaussian", "sine");

    return 0;
}
