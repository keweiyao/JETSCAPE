// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#pragma once

// forward declarations

namespace boost {

namespace filesystem {
class path;
}

namespace math {}

namespace program_options {
class variables_map;
}

}  // namespace boost

namespace fs = boost::filesystem;
namespace po = boost::program_options;
namespace math = boost::math;

using VarMap = po::variables_map;

namespace trento3d {
class Event;
class Nucleus;
}
