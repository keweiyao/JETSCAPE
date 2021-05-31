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

#ifndef TRENTO3DINITIAL_H
#define TRENTO3DINITIAL_H

#include <tuple>
#include <memory>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

#include "fwd_decl.h"
#include "JetScapeModuleBase.h"
#include "trento3d/src/collider.h"
#include "InitialState.h"
#include "JetScapeLogger.h"

using OptDesc = po::options_description;
using VarMap = po::variables_map;
using namespace trento3d;

namespace Jetscape {

typedef struct {
  double impact_parameter;
  double num_participant;
  double num_binary_collisions;
  double total_entropy;
  std::map<int, double> ecc; // order, eccentricity
  std::map<int, double> psi; // order, participant_plane (complex eccentricity argument)
  double xmid, ymid;
} EventInfo;

////////////////////////// Trento-3D Initial Condition Wrapper //////////////////////
class Trento3DInitial : public InitialState {
public:
  // Initialize from XML configuration
  Trento3DInitial();
  ~Trento3DInitial();

  //void Init();
  void Exec();
  void Clear();
  void InitTask();

  struct RangeFailure : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };
  EventInfo info_;

private:
  std::shared_ptr<trento3d::Collider> Trento3DGen_;
  std::pair<double, double> GenCenTab(std::string proj, std::string targ,
                                      VarMap var_map, int cL, int cH);

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<Trento3DInitial> reg;
};

} // end namespace Jetscape

#endif // TRENTO3DINITIAL_H
