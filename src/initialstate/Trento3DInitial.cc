/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <algorithm>
#include <functional>
#include <string>
#include "JetScapeLogger.h"

#include "Trento3DInitial.h"

namespace Jetscape {

// Register the module with the base class
RegisterJetScapeModule<Trento3DInitial> Trento3DInitial::reg("Trento3DInitial");

namespace {
/// @brief Tokenize a string.  The tokens will be separated by each non-quoted
///        space or equal character.  Empty tokens are removed.
///
/// @param input The string to tokenize.
///
/// @return Vector of tokens.
std::vector<std::string> tokenize(const std::string &input) {
  typedef boost::escaped_list_separator<char> separator_type;
  separator_type separator("\\",    // The escape characters.
                           "= ",    // The separator characters.
                           "\"\'"); // The quote characters.

  // Tokenize the intput.
  boost::tokenizer<separator_type> tokens(input, separator);

  // Copy non-empty tokens from the tokenizer into the result.
  std::vector<std::string> result;
  copy_if(tokens.begin(), tokens.end(), std::back_inserter(result),
          !boost::bind(&std::string::empty, _1));
  return result;
}
} // end namespace

// See header for explanation.
Trento3DInitial::Trento3DInitial() : InitialState() { SetId("Trento3D"); }

Trento3DInitial::~Trento3DInitial() = default;

void Trento3DInitial::InitTask() {
  JSINFO << " Initialzie TRENTo-3D initial condition ";

  // TRENTO-3D OPTION DESK
  using namespace trento3d;
  using OptDesc = po::options_description;
  using VecStr = std::vector<std::string>;
  OptDesc main_opts{};
  main_opts.add_options()
    ("projectile",
     po::value<VecStr>()
      ->required()
      ->notifier( // use a lambda to verify there are exactly two projectiles
        [](const VecStr &projectiles) {
          if (projectiles.size() != 2)
            throw po::required_option{"projectile"};
          }),
     "projectile symbols")
    ("number-events", po::value<int>()->default_value(1),
     "number of events");

  // Make all main arguments positional.
  po::positional_options_description positional_opts{};
  positional_opts
    .add("projectile", 2)
    .add("number-events", 1);

  using VecPath = std::vector<fs::path>;
  OptDesc general_opts{"general options"};
  general_opts.add_options()
    ("help,h", "show this help message and exit")
    ("version", "print version information and exit")
    ("bibtex", "print bibtex entry and exit")
    // ("default-config", "print a config file with default settings and exit")
    ("config-file,c", po::value<VecPath>()->value_name("FILE"),
     "configuration file\n(can be passed multiple times)");

  OptDesc output_opts{"output options"};
  output_opts.add_options()
    ("quiet,q", po::bool_switch(),
     "do not print event properties to stdout")
    ("output,o", po::value<fs::path>()->value_name("PATH"),
     "HDF5 file or directory for text files")
    ("no-header", po::bool_switch(),
     "do not write headers to text files")
    ("ncoll", po::bool_switch(),
     "calculate binary collisions");

  OptDesc phys_opts{"physical options"};
  phys_opts.add_options()
    // Random seed
    ("random-seed",
     po::value<int64_t>()->value_name("INT")->default_value(-1, "auto"),
     "random seed")
    // Nulcear configuration parameters
    ("form-width",
     po::value<double>()->value_name("FLOAT")->default_value(.5, "0.5"),
     "form-factor width of inelastic collisions [fm] (0.35, 1.0)")
    ("nucleon-width,w",
     po::value<double>()->value_name("FLOAT")->default_value(.5, "0.5"),
     "Gaussian nucleon width [fm] (0.35, 1.4)")
    ("constit-width,v",
     po::value<double>()->value_name("FLOAT")->default_value(.5, "same"),
     "Gaussian constituent width [fm]")
    ("constit-number,m",
     po::value<double>()->value_name("FLOAT")->default_value(1., "1"),
     "Number of constituents in the nucleon")
    ("nucleon-min-dist,d",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum nucleon-nucleon distance [fm]")
    ("b-min",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum impact parameter [fm]")
    ("b-max",
     po::value<double>()->value_name("FLOAT")->default_value(-1., "auto"),
     "maximum impact parameter [fm]")
    ("npart-min",
     po::value<int>()->value_name("INT")->default_value(0, "0"),
     "minimum Npart cut")
    ("npart-max",
     po::value<int>()->value_name("INT")->default_value(
     std::numeric_limits<int>::max(), "INT_MAX"), "maximum Npart cut")
    ("mult-min",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "minimum midrapidity dET/detas cut")
    ("mult-max",
     po::value<double>()->value_name("FLOAT")->default_value(
     std::numeric_limits<double>::max(), "DOUBLE_MAX"), "maxmimum midrapidity dET/detas cut")
    ("etas-shift",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "eta at grid center, use if Ebeam1 != Ebeam2")
    ("mult-etas-low",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "eta lower bound for calculating multiplicity")
    ("mult-etas-high",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0"),
     "eta upper bound for calculating multiplicity")
    ("sqrts,s",
     po::value<double>()->value_name("FLOAT")->default_value(2760., "2760."),
     "CoM energy of collision [GeV], determines cross-section, ybeam, etc")
    ("fluctuation,k",
     po::value<double>()->value_name("FLOAT")->default_value(0.3, "0.3"),
     "fluctuation of energy sharing of central and forward fireball (0.1, 0.6)")
    ("kT-min,t",
     po::value<double>()->value_name("FLOAT")->default_value(.4, ".4"),
     "eta_max ~ ln(sqrts/kT-min) [GeV], (0.2, 1.0)")
    ("reduced-thickness,p",
     po::value<double>()->value_name("FLOAT")->default_value(0., "0."),
     "reduced thickness parameter, !!currently set to 0 internally!!")
    ("shape-alpha",
     po::value<double>()->value_name("FLOAT")->default_value(4., "4."),
     "fragmentation region a: [-ln(x)]^a x^b, (3, 5)")
    ("shape-beta",
     po::value<double>()->value_name("FLOAT")->default_value(.5, ".5"),
     "fragmentation region b: [-ln(x)]^a x^b, (0.2, 1.0)")
    ("mid-power",
     po::value<double>()->value_name("FLOAT")->default_value(0.4, "0.4"),
     "power of ETmid ~ #1*sqrts^#2, (0.3, 0.48)")
    ("mid-norm",
     po::value<double>()->value_name("FLOAT")->default_value(.3, ".3"),
     "norm of ETmid ~ #1*sqrts^#2, (0.2, 0.45)")
    ("flatness",
     po::value<double>()->value_name("FLOAT")->default_value(1.5, "1.5"),
     "flatness parameter of the central plateau, (1.0, 1.5)");

  OptDesc grid_opts{"grid options"};
  grid_opts.add_options()
    ("grid-max",
     po::value<double>()->value_name("FLOAT")->default_value(10., "10.0"),
     "xy max [fm]\n(grid extends from -max to +max)")
    ("grid-step",
     po::value<double>()->value_name("FLOAT")->default_value(0.2, "0.2"),
     "step size [fm]")
    ("nsteps-etas",
     po::value<int>()->value_name("INT")->default_value(1, "1"),
     "number of grid steps in eta_s");

  // Make a meta-group containing all the option groups except the main
  // positional options (don't want the auto-generated usage info for those).
  OptDesc usage_opts{};
  usage_opts.add(general_opts)
      .add(output_opts)
      .add(phys_opts)
      .add(grid_opts);

  // Now a meta-group containing _all_ options.
  OptDesc all_opts{};
  all_opts.add(usage_opts).add(main_opts);

  // Will be used several times.
  const std::string usage_str{
      "usage: trento-3 [options] projectile projectile [number-events = 1]\n"};

  // NOW LETS FILL IN THE OPTION DESK
  auto phy_opts = GetXMLElement({"IS", "Trento3D", "PhysicsInputs"});
  auto cut_opts = GetXMLElement({"IS", "Trento3D", "CutInputs"});
  auto trans_opts = GetXMLElement({"IS", "Trento3D", "TransInputs"});
  auto longi_opts = GetXMLElement({"IS", "Trento3D", "LongiInputs"});

  double xymax = GetXMax(), dxy = GetXStep();
  double etamax = GetZMax(), deta = GetZStep();

  auto random_seed = (*GetMt19937Generator())();
  //TEMPORARY FOR TESTING
  //auto random_seed = 1;
  //TEMPORARY
  JSINFO << "Random seed used for Trento-3D " << random_seed;

  std::string proj(phy_opts->Attribute("projectile"));
  std::string targ(phy_opts->Attribute("target"));
  double sqrts = std::atof(phy_opts->Attribute("sqrts"));
  double mid_norm = std::atof(phy_opts->Attribute("mid-norm"));
  double mid_power = std::atof(phy_opts->Attribute("mid-power"));

  int cen_low = std::atoi(cut_opts->Attribute("centrality-low"));
  int cen_high = std::atoi(cut_opts->Attribute("centrality-high"));
  double etas_shift = std::atof(cut_opts->Attribute("etas-shift"));
  double mult_etas_low = std::atof(cut_opts->Attribute("mult-etas-low"));
  double mult_etas_high = std::atof(cut_opts->Attribute("mult-etas-high"));

  double p = std::atof(trans_opts->Attribute("reduced-thickness"));
  double u = std::atof(trans_opts->Attribute("form-width"));
  double w = std::atof(trans_opts->Attribute("nucleon-width"));
  double d = std::atof(trans_opts->Attribute("nucleon-min-dist"));
  int    m = std::atoi(trans_opts->Attribute("constit-number"));
  double v = std::atof(trans_opts->Attribute("constit-width"));

  double k        = std::atof(longi_opts->Attribute("fluctuation"));
  double kT_min   = std::atof(longi_opts->Attribute("kT-min"));
  double alpha    = std::atof(longi_opts->Attribute("shape-alpha"));
  double beta     = std::atof(longi_opts->Attribute("shape-beta"));
  double flatness = std::atof(longi_opts->Attribute("flatness"));

  std::string options1 =  // options affecting centrality table
      " --sqrts "             + std::to_string(sqrts) +
      " --reduced-thickness " + std::to_string(p) + //TODO//need backward compat work on trento-3 side//XXX//
      " --form-width "        + std::to_string(u) +
      " --nucleon-width "     + std::to_string(w) +
      " --nucleon-min-dist "  + std::to_string(d) +
      " --constit-number "    + std::to_string(m) +
      " --constit-width "     + std::to_string(v) +
      " --fluctuation "       + std::to_string(k) + //TODO//need backward compat work on trento-3 side//XXX//
      " --kT-min "            + std::to_string(kT_min) +
      " --shape-alpha "       + std::to_string(alpha) +
      " --shape-beta "        + std::to_string(beta) +
      " --flatness "          + std::to_string(flatness) +
      " --mid-norm "          + std::to_string(mid_norm) +
      " --mid-power "         + std::to_string(mid_power) +
      " --nsteps-etas "       + std::to_string(static_cast<int>(std::round( 2. * GetZMax() / GetZStep() ))) +//TODO//is this right for now? should revamp JETSCAPE config to do this better//XXX//
      " --etas-shift "        + std::to_string(etas_shift) +
      " --mult-etas-low "     + std::to_string(mult_etas_low) +
      " --mult-etas-high "    + std::to_string(mult_etas_high) +
      " --quiet ";
  std::string options2 =
      " --random-seed "       + std::to_string(random_seed) +
      " --ncoll "             + // calcualte # of binary collision
      " --grid-max "          + std::to_string(xymax) +
      " --grid-step "         + std::to_string(dxy);

  // Handle centrality table, not normzlized, default grid, 2D (fast) !!!
  std::string cmd_basic = proj + " " + targ + " 10000 --random-seed 1 " + options1;  // random seed shouldn't matter, but specify a fixed seed anyway for reproducibility
  VarMap var_map_basic{};
  po::store(po::command_line_parser(tokenize(cmd_basic))
                .options(all_opts)
                .positional(positional_opts)
                .run(),
            var_map_basic);

  std::string options_cut = "";
  if (cen_low == 0 && cen_high == 100) {
    JSINFO << "TRENTo-3D Minimum Biased Mode Generates 0-100(%) of nuclear "
              "inelastic cross-section";
  } else {
    auto Ecut = GenCenTab(proj, targ, var_map_basic, cen_low, cen_high);
    double Ehigh = Ecut.first;
    double Elow = Ecut.second;

    JSINFO << "The total energy density cut for centrality = [" << cen_low
           << ", " << cen_high << "] (%) is:";
    JSINFO << Elow << "<dE/deta(eta=0)<" << Ehigh;
    options_cut = " --mult-max " + std::to_string(Ehigh) + " --mult-min " +
                  std::to_string(Elow);
    // Set trento configuration
  }
  std::string cmd =
      proj + " " + targ + " 1 " + options1 + options2 + options_cut;
  JSINFO << cmd;
  VarMap var_map{};
  po::store(po::command_line_parser(tokenize(cmd))
                .options(all_opts)
                .positional(positional_opts)
                .run(),
            var_map);
  Trento3DGen_ = std::make_shared<trento3d::Collider>(var_map);
  SetRanges(xymax, xymax, etamax);
  SetSteps(dxy, dxy, deta);
  JSINFO << "TRENTo-3D set";
}

bool compare_E(trento3d::records r1, trento3d::records r2) {
  return r1.mult > r2.mult;
}

std::pair<double, double> Trento3DInitial::GenCenTab(std::string proj,
                                                     std::string targ,
                                                     VarMap var_map, int cL,
                                                     int cH) {
  // Terminate for nonsense
  if (cL < 0 || cL > 100 || cH < 0 || cH > 100 || cH < cL) {
    JSWARN << "Wrong centrality cuts! To be terminated.";
    exit(-1);
  }
  // These are all the parameters that could change the shape of centrality tables
  // Normalization prefactor parameter is factorized
  // They form a table header
  trento3d::Collider another_collider(var_map);
  auto sqrts          = var_map["sqrts"].as<double>();
  auto form_width     = var_map["form-width"].as<double>();
  auto nucleon_width  = var_map["nucleon-width"].as<double>();
  auto constit_width  = var_map["constit-width"].as<double>();
  auto constit_number = var_map["constit-number"].as<double>();
  auto dmin           = var_map["nucleon-min-dist"].as<double>();
  auto kT_min         = var_map["kT-min"].as<double>();
  auto shape_alpha    = var_map["shape-alpha"].as<double>();
  auto shape_beta     = var_map["shape-beta"].as<double>();
  auto mid_norm       = var_map["mid-norm"].as<double>();
  auto mid_power      = var_map["mid-power"].as<double>();
  auto fluctuation    = var_map["fluctuation"].as<double>();
  auto flatness       = var_map["flatness"].as<double>();
  auto nsteps_etas    = var_map["nsteps-etas"].as<int>();  // nsteps_etas doesn't affect the continuous function dET/detas, but does affect at which etas it's evaluated, and thus the reported multiplicity
  auto etas_shift     = var_map["etas-shift"].as<double>();
  auto mult_etas_low  = var_map["mult-etas-low"].as<double>();
  auto mult_etas_high = var_map["mult-etas-high"].as<double>();
  char buffer[512] = {};  // initializer ensures there's a terminating NUL
  std::snprintf(buffer, sizeof(buffer) - 1,
                "%s_%s_%.1f_u%.4f_w%.4f_v%.4f_m%.4f_d%.4f_kT%.4f_a%.4f_b%.4f_N%.4f_pmid%.4f_k%.4f_f%.4f_netas%d_etassh%.4f_multetas%.4f_%.4f",
                 proj.c_str(),
                    targ.c_str(),
                       sqrts,
                             form_width,
                                   nucleon_width,
                                         constit_width,
                                               constit_number,
                                                     dmin,
                                                            kT_min,
                                                                  shape_alpha,
                                                                        shape_beta,
                                                                              mid_norm,
                                                                                       mid_power,
                                                                                             fluctuation,
                                                                                                   flatness,
                                                                                                             nsteps_etas,
                                                                                                                      etas_shift,
                                                                                                                                   mult_etas_low,
                                                                                                                                        mult_etas_high);
  std::string header(buffer);
  JSINFO << "TRENTo-3D centrality table header: " << header;
  // Create headering string hash tage for these parameter combination
  // Use this tag as a unique table filename for this specific parameter set
  std::hash<std::string> hash_function;
  size_t header_hash = hash_function(header);
  JSINFO << "Hash tag for this header: " << header_hash;
  // create dir incase it does not exist
  boost::filesystem::path datadir{"trento3d_data"};
  boost::filesystem::create_directory(datadir);
  char filename[512];
  std::sprintf(filename, "%zu", header_hash);
  // Step1: check it a table exist
  auto filepath = (datadir / filename).string();
  std::ifstream infile(filepath);
  double Etab[101];
  double buff1, buff2;
  std::string line;
  if (infile.good()) {
    JSINFO << "The required centrality table exists. Load the table.";
    int i = 0;
    while (std::getline(infile, line)) {
      if (line[0] != '#') {
        std::istringstream iss(line);
        iss >> buff1 >> buff2 >> Etab[i];
        i++;
      }
    }
    infile.close();
  } else {
    JSINFO << "TRENTo-3D is generating new centrality table for this new "
              "parameter set";
    JSINFO << "It may take 10(s) to 1(min).";

    another_collider.run_events();
    // Get all records and sort according to totoal energy
    auto event_records = another_collider.all_records();
    std::sort(event_records.begin(), event_records.end(), compare_E);
    // write centrality table
    int nstep = std::ceil(event_records.size() / 100);
    std::ofstream fout(filepath);
    fout << "# projectile target sqrts "
           << "form_width nucleon_width constit_width constit_number dmin kT_min shape_alpha shape_beta mid_norm mid_power fluctuation flatness "
           << "nsteps_etas etas_shift mult_etas_low mult_etas_high" << std::endl
         << "# " << proj << " " << targ << " " << sqrts
          << " " << form_width << " " << nucleon_width << " " << constit_width << " " << constit_number << " " << dmin
          << " " << kT_min << " " << shape_alpha << " " << shape_beta << " " << mid_norm << " " << mid_power << " " << fluctuation << " " << flatness
          << " " << nsteps_etas << " " << etas_shift << " " << mult_etas_low << " " << mult_etas_high
          << std::endl
         << "# cen_L  cen_H  un-normalized total density" << std::endl;
    Etab[0] = 1e10;
    for (int i = 1; i < 100; i += 1) {
      auto ee = event_records[i * nstep];
      fout << i - 1 << " " << i << " " << ee.mult << std::endl;
      Etab[i] = ee.mult;
    }
    auto ee = event_records.back();
    fout << 99 << " " << 100 << " " << ee.mult << std::endl;
    Etab[100] = ee.mult;
    fout.close();
  }
  JSINFO << "#########" << Etab[cL] << " " << Etab[cH];
  return std::make_pair(Etab[cL], Etab[cH]);
}

std::map<int,double>QUICKHACK(const std::map<int,std::vector<double>>a){//TODO//this is a quick hack for a 2D (N x N x 1) test run; remove it//XXX//
 std::map<int,double>b;
 for(const auto&iter:a)b[iter.first]=iter.second[0];
 return b;}
void Trento3DInitial::Exec() {
  JSINFO << " Exec TRENTo-3D initial condition ";
  Trento3DGen_->run_events();

  JSINFO << " TRENTo-3D event info: ";
  auto tmp_event = Trento3DGen_->expose_event();
  info_.impact_parameter = Trento3DGen_->all_records().back().b;
  info_.num_participant = tmp_event.npart();
  info_.num_binary_collisions = tmp_event.ncoll();
  info_.total_entropy = tmp_event.multiplicity();
  info_.ecc = QUICKHACK(tmp_event.ecc_mag());
  info_.psi = QUICKHACK(tmp_event.ecc_ang());
  info_.xmid =
      -GetXMax() + tmp_event.mass_center_index()[0].first * tmp_event.dxy();//TODO//this [0] is a quick hack for a 2D (N x N x 1) test run//XXX//
  info_.ymid =
      -GetYMax() + tmp_event.mass_center_index()[0].second * tmp_event.dxy();//TODO//this [0] is a quick hack for a 2D (N x N x 1) test run//XXX//
  JSINFO << "b\tnpart\tncoll\tET\t(x-com, y-com) (fm)";
  JSINFO << info_.impact_parameter << "\t" << info_.num_participant << "\t"
         << info_.num_binary_collisions << "\t" << info_.total_entropy << "\t"
         << "(" << info_.xmid << ", " << info_.ymid << ")";

  JSINFO << " Load TRENTo-3D density and ncoll density to JETSCAPE memory ";
  auto density_field = tmp_event.Density();
  auto ncoll_field = tmp_event.TAB2D_grid();
  JSINFO << density_field.num_elements() << " density elements";
  for (int i = 0; i < density_field.num_elements(); i++) {
    entropy_density_distribution_.push_back(density_field.data()[i]);  // Event::Density().data() is in z-major order (zyx = {000, 001, ..., 010, 011, ... 100, 101, ..., 110, 111, ...})
  }
  JSINFO << ncoll_field.num_elements() << " ncoll elements";
  for (int i = 0; i < ncoll_field.num_elements(); i++) {
    num_of_binary_collisions_.push_back(ncoll_field.data()[i]);
  }
  JSINFO << " TRENTo-3D event generated and loaded ";
}

void Trento3DInitial::Clear() {
  VERBOSE(2) << " : Finish creating initial condition ";
  entropy_density_distribution_.clear();
  num_of_binary_collisions_.clear();
}

} // end namespace Jetscape
