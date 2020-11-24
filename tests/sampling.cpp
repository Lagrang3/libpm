#define BOOST_TEST_MODULE Mass Sampling
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <pm/field.hpp>
#include <pm/filters.hpp>
#include <vector>

using std::to_string;

template <class>
struct Decorator;

template <class R, class... Args>
struct Decorator<R(Args...)>
{
    std::function<R(Args...)> f_;

    Decorator(std::function<R(Args...)> f) : f_{f} {}

    R operator()(Args... args)
    {
        try
        {
            f_(args...);
        }
        catch (std::exception& e)
        {
            std::stringstream mes;
            mes << "Exception thrown: " << e.what();
            BOOST_TEST_MESSAGE(mes.str());
            throw;
        }
    }
};

template <class R, class... Args>
Decorator<R(Args...)> makeDecorator(R (*f)(Args...))
{
    return Decorator<R(Args...)>(std::function<R(Args...)>(f));
}

namespace ut = boost::unit_test;
namespace mpi = boost::mpi;

struct fixture
{
    fixture() {}

    ~fixture() {}

    static mpi::environment env;
    static mpi::communicator world;
};

template <class T, class filter_t>
struct sample_checker
{
    filter_t W;
    const int N;
    std::vector<T> data;

    T& operator()(int i, int j, int k)
    {
        i = PM::utilities::modulo(i, N);
        j = PM::utilities::modulo(j, N);
        k = PM::utilities::modulo(k, N);
        return data[k + N * (j + N * i)];
    }
    const T& operator()(int i, int j, int k) const
    {
        i = PM::utilities::modulo(i, N);
        j = PM::utilities::modulo(j, N);
        k = PM::utilities::modulo(k, N);
        return data[k + N * (j + N * i)];
    }

    sample_checker(int N) : N{N}, data(N * N * N) {}

    void clear() { std::fill(data.begin(), data.end(), T{}); }

    T sum() const
    {
        T s{};
        for (auto x : data)
            s += x;
        return s;
    }
    void sample(std::array<double, 3> p) { sample(p[0], p[1], p[2]); }
    void sample(double x, double y, double z)
    {
        x = PM::utilities::modulo(x, N);
        y = PM::utilities::modulo(y, N);
        z = PM::utilities::modulo(z, N);

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k)
                {
                    T w{};
                    for (int fx = -1; fx <= 1; ++fx)
                        for (int fy = -1; fy <= 1; ++fy)
                            for (int fz = -1; fz <= 1; ++fz)
                                w += W(x - i + fx * N) * W(y - j + fy * N) *
                                     W(z - k + fz * N);

                    (*this)(i, j, k) += w;
                }
    }
};

template <class mpiFun, class serFun>
void add_and_check(mpiFun& F,
                   serFun& F2,
                   std::array<double, 3> p,
                   double tolerance)
{
    const int N = F.global_extents();

    F2.sample(p);

    if (fixture::world.rank() == 0)
        F.sample({p});
    else
        F.sample({});

    auto ser = F.serialize();
    boost::math::fpc::close_at_tolerance<double> close(tolerance);

    if (fixture::world.rank() == 0)
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k)
                {
                    const double a = ser[k + N * (j + i * N)], b = F2(i, j, k);
                    if (not close(a, b))
                        throw std::runtime_error(
                            "densities diverge at (" + to_string(i) + "," +
                            to_string(j) + "," + to_string(k) + ") " + "a{" +
                            to_string(a) + "} b{" + to_string(b) + "}");
                }
}

mpi::environment fixture::env;
mpi::communicator fixture::world;

BOOST_TEST_GLOBAL_FIXTURE(fixture);

typedef boost::mpl::
    list<PM::filters::NGP, PM::filters::CIC, PM::filters::TSC, PM::filters::PCS>
        test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(Sampling, filter_t, test_types)
{
    const int N_grid = 10;
    const double tolerance = 0.0001;  // percentage
    double sum = 0;
    PM::Field<double, filter_t, filter_t> F(fixture::world, N_grid, {2, 2});
    sample_checker<double, filter_t> F2(N_grid);
    F.clear();
    F2.clear();

    BOOST_CHECK_NO_THROW(makeDecorator(
        add_and_check<decltype(F), decltype(F2)>)(F, F2, {0, 0, 0}, tolerance));
    sum += 1.0;
    BOOST_CHECK_CLOSE(F.sum(), sum, tolerance);
    BOOST_CHECK_CLOSE(F2.sum(), sum, tolerance);

    BOOST_CHECK_NO_THROW(makeDecorator(
        add_and_check<decltype(F), decltype(F2)>)(F, F2, {1, 1, 1}, tolerance));
    sum += 1.0;
    BOOST_CHECK_CLOSE(F.sum(), sum, tolerance);
    BOOST_CHECK_CLOSE(F2.sum(), sum, tolerance);

    BOOST_CHECK_NO_THROW(makeDecorator(
        add_and_check<decltype(F), decltype(F2)>)(F, F2, {0, 0, 0}, tolerance));
    sum += 1.0;
    BOOST_CHECK_CLOSE(F.sum(), sum, tolerance);
    BOOST_CHECK_CLOSE(F2.sum(), sum, tolerance);

    BOOST_CHECK_NO_THROW(
        makeDecorator(add_and_check<decltype(F), decltype(F2)>)(
            F, F2, {10, 10, 10}, tolerance));
    sum += 1.0;
    BOOST_CHECK_CLOSE(F.sum(), sum, tolerance);
    BOOST_CHECK_CLOSE(F2.sum(), sum, tolerance);

    for (int i = 0; i < N_grid; ++i)
        for (int j = 0; j < N_grid; ++j)
            for (int k = 0; k < N_grid; ++k)
            {
                BOOST_CHECK_NO_THROW(
                    makeDecorator(add_and_check<decltype(F), decltype(F2)>)(
                        F, F2, {i * 1.0, j * 1.0, k * 1.0}, tolerance));
                sum += 1.0;
                BOOST_CHECK_CLOSE(F.sum(), sum, tolerance);
                BOOST_CHECK_CLOSE(F2.sum(), sum, tolerance);
            }
}
