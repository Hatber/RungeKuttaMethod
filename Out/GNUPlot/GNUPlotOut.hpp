#ifndef GNUPLOTOUT_HPP
#define GNUPLOTOUT_HPP

#include "../AOut.hpp"
#include "lib/gnuplot_i.hpp"

class GNUPlotOut : public f_AOut {
public:
    GNUPlotOut();
    GNUPlotOut(std::size_t coordCount);

    void out(const std::string& name);

private:
    std::vector< Gnuplot > _plots;
};

#endif // GNUPLOTOUT_HPP
