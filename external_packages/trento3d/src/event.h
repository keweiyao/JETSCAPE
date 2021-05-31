// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#pragma once

#include <functional>
#include <map>

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
using boost::math::interpolators::cardinal_cubic_b_spline;

#include "fwd_decl.h"
//#include "random_field.h"

namespace trento3d {

class NucleonCommon;
class NucleonData;

/// \rst
/// The primary computation class, responsible for constructing nuclear
/// thickness functions and calculating event observables.  Designed to be
/// created once and used many times by repeatedly calling ``compute()``.
/// Stores its observables internally and provides inspector methods.
///
/// Example::
///
///   Event event{var_map};
///   for (int n = 0; n < nevents; ++n) {
///     event.compute(nucleusA, nucleusB, nucleon_profile);
///     do_something(
///       event.npart(),
///       event.multiplicity(),
///       event.eccentricity(),
///       event.reduced_thickness_grid()
///     );
///   }
///
/// \endrst
class Event {
 public:
  /// Instantiate from the configuration.
  explicit Event(const VarMap& var_map);

  /// \rst
  /// Compute thickness functions and event observables for a pair of
  /// ``Nucleus`` objects and a ``NucleonProfile``.  The nuclei must have
  /// already sampled nucleon positions and participants before passing to this
  /// function.
  /// \endrst
  void compute(const Nucleus& nucleusA, const Nucleus& nucleusB,
               const NucleonCommon& nucleon_common);

  // Alias for a two-dimensional thickness grid.
  using Grid = boost::multi_array<double, 2>;

  // Alias for a three-dimensional thickness grid.
  using Grid3D = boost::multi_array<double, 3>;


  /// Number of nucleon participants.
  const int& npart() const
  { return npart_; }
  const int& npartA() const
  { return npartA_; }
  const int& npartB() const
  { return npartB_; }

  /// WK: Number of binary collision.
  const int& ncoll() const
  { return ncoll_; }

  /// Used by Collider to store the number of binary collisions it computes inside this Event,
  /// so that it can be exposed to external code.
  void set_ncoll(int ncoll)
  { ncoll_ = ncoll; }

  const std::vector<std::pair<double, double>> mass_center_index() const {
      std::vector<std::pair<double, double>> vec;
      for (int iz = 0; iz < nsteps_etas_; ++iz)
          vec.push_back(std::make_pair(ixcm_[iz], iycm_[iz]));
      return vec;
  }

  /// \rst
  /// Multiplicity---or more specifically, total integrated reduced thickness.  May be interpreted
  /// as `dS/d\eta` or `dE/d\eta` at midrapidity.
  /// \endrst
  const double& multiplicity() const
  { return dET_detas_[std::floor(nsteps_etas_/2)]; }


  const std::vector<double> & dET_detas() const
  { return dET_detas_; }

  /// \rst
  /// Eccentricity harmonics `\varepsilon_n` for *n* = 2--5.
  /// Returns a map of `(n : \varepsilon_n)` pairs, so e.g.::
  ///
  ///   double e2 = event.eccentricity().at(2);
  ///
  /// \endrst
  const std::map<int, std::vector<double> > & ecc_mag() const
  { return ecc_mag_; }

  const std::map<int, std::vector<double> > & ecc_ang() const
  { return ecc_ang_; }

  /// Density grid
  const Grid3D & Density() const
  { return Density_; }

  /// returns grid steps
  const double& dxy() const
  { return dxy_; }

  /// WK: The TAB grid for hard process vertex sampling
  const Grid& TAB2D_grid() const
  { return TAB2D_; }

  /// WK: clear and increase TAB
  void clear_TAB2D(void);
  void accumulate_TAB2D(NucleonData& A, NucleonData& B, NucleonCommon& profile);

  double central_profile(double eta) const{
     if (std::abs(eta)>eta_max_) return 0.;
     double u = eta*eta/2./eta_max_;
     return std::exp(-std::pow(u, flatness_))
           *std::pow(1.-std::pow(eta/eta_max_, 4), 4);
  }
  double center_of_mass_eta(double ta, double tb) const{
     return .5*std::log((ta*Pplus_+tb*Pminus_)
                       /(tb*Pplus_+ta*Pminus_));
  }
  double frag_profile(double x) const{
     return std::pow(-std::log(x), shape_alpha_)
          * std::pow(x, shape_beta_+1) * std::exp(-2*kT_min_/sqrts_/x);
  }

 private:
  /// Compute a nuclear thickness function (TA or TB) onto a grid for a given
  /// nucleus and nucleon profile.  This destroys any data previously contained
  /// by the grid.
  void compute_nuclear_thickness(
      const Nucleus& nucleus, const NucleonCommon& nucleon_common, 
                   Grid& TX, Grid& FX);

  /// Compute the reduced thickness function (TR) after computing TA and TB.
  void compute_density();

  /// Compute observables that require a second pass 
  /// over the reduced thickness grid.
  void compute_observables();

  /// Normalization factor.
  double norm_trento_, Nfrag_, xloss_, Pplus_, Pminus_;

  //const double scaling_p0_;
  const double mid_power_, mid_norm_, flatness_; 
  const double shape_alpha_, shape_beta_;
  const double kT_min_;
  const double sqrts_, nucleon_pabs_;
  const double eta_max_, eta_grid_max_;

  /// Number of grid steps
  const int nsteps_etas_;  
  const double detas_;
  const double dxy_;
  const int nsteps_;

  /// Grid maximum (half width).
  const double xymax_;
  
  /// Nuclear thickness grids TA, TB;
  Grid TA_, TB_, FA_, FB_, TAB2D_;

  // 3D grid for matter deposition density
  Grid3D Density_;

  /// Center of mass coordinates in "units" of grid index (not fm).
  /// as an array of rapidity
  std::vector<double> ixcm_, iycm_;

  /// Number of participants.
  int npart_, npartA_, npartB_;

  /// WK: Number of binary collisions.
  int ncoll_;

  /// Total matter production as an array of rapidity
  std::vector<double> dET_detas_;

  /// Eccentricity harmonics.
  std::map<int, std::vector<double> > ecc_mag_;
  /// participant plane angles
  std::map<int, std::vector<double> > ecc_ang_;
};

}  // namespace trento3d
