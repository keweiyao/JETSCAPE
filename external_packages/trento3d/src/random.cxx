// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "random.h"


namespace trento3d { namespace random {

// Seed random number generator from hardware device.
Engine engine{std::random_device{}()};

}}  // namespace trento3d::random
