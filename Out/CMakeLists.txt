add_library(thermal_out SHARED
  AOut.cpp
  GNUPlot/GNUPlotOut.cpp
)

install(TARGETS thermal_out LIBRARY DESTINATION lib)
