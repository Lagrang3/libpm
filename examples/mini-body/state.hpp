#pragma once

#include <cmath>
#include <gadget.hpp>
#include <vector>

class state
{
   public:
    /* inmutable parameter */
    double HubbleParam, BoxSize, Gm, Omega0, OmegaLambda;
    std::vector<int64_t> ids;
    std::vector<double> zvect;

    /* dynamical parameters */
    double a, H;
    std::vector<float> X, P, F;

    void set_NumPart(const int N)
    {
        X.resize(3 * N);
        P.resize(3 * N);
        F.resize(3 * N);
        ids.resize(N);
    }
    auto get_NumPart() const { return ids.size(); }

    /* write to snapshot */
    void write(const std::string filename) const
    {
        gadget::osnapshot<1> snap(filename);
        auto& header = snap.get_raw_header();
        header.BoxSize = BoxSize;
        header.num_files = 1;
        header.npartTotal[1] = header.npart[1] = get_NumPart();
        header.redshift = 1.0 / a - 1;
        header.Omega0 = Omega0;
        header.OmegaLambda = OmegaLambda;
        // header.mass[1] = 1.0; OJO

        snap.write_header();

        {
            auto X_tmp{X};
            for (auto& x : X_tmp)
            {
                x *= BoxSize;
            }
            snap.write_block("POS ", X_tmp.begin(), X_tmp.end());
        }
        {
            const double LH_a2_inv = BoxSize * HubbleParam / (a * a);
            auto P_tmp{P};
            for (auto& v : P_tmp)
            {
                v *= LH_a2_inv;
            }
            snap.write_block("VEL ", P_tmp.begin(), P_tmp.end());
        }
        snap.write_block("ID  ", ids.begin(), ids.end());
    }

    /* read from snapshot */
    void read(const std::string filename)
    {
        gadget::isnapshot<1> snap(filename);

        set_NumPart(snap.npart(1));

        const double pi = acos(-1);

        // deduced quantities
        OmegaLambda = 1 - Omega0;
        Gm = 3 * Omega0 / (8 * pi * get_NumPart());
        a = 1 / (zvect[0] + 1);

        snap.scan_block("POS ", X.begin());
        snap.scan_block("VEL ", P.begin());
        snap.scan_block("ID  ", ids.begin());

        const double BoxSize_inv = 1 / BoxSize;
        for (auto& x : X)
        {
            x *= BoxSize_inv;
        }

        const double a2_LH_inv = a * a / BoxSize / HubbleParam;
        for (auto& v : P)
        {
            v *= a2_LH_inv;
        }
    }

    double get_H() const { return sqrt(OmegaLambda + Omega0 / (a * a * a)); }

    void advance_a(const double dt)
    {
        H = sqrt(OmegaLambda + Omega0 / (a * a * a));
        a += dt * a * H;
    }

    /* enforce Periodic-Boundary-Condition */
    void enforce_PBC()
    {
        for (auto& x : X)
            x -= std::floor(x);
    }
};
