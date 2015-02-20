#include "GNUPlotOut.hpp"

#include <iostream>

GNUPlotOut::GNUPlotOut() : f_AOut(0) { }
GNUPlotOut::GNUPlotOut(std::size_t coordCount) : f_AOut(coordCount) { }

void GNUPlotOut::out(const std::string& name) {
    Gnuplot plot;
    plot.set_title(name.c_str());
    plot.set_style("lines");
    plot.cmd("set terminal png size 1024,768");
    plot.cmd("set output \"" + name + ".png\"");
    plot.plot_xy(_y, _x);

    //std::cin.get();
}
