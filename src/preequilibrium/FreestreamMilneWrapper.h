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

#ifndef FREESTREAMMILNEWRAPPER_H
#define FREESTREAMMILNEWRAPPER_H

#include <string>
#include "PreequilibriumDynamics.h"
#include "FreestreamMilne.h"

using namespace Jetscape;

class FreestreamMilneWrapper : public PreequilibriumDynamics {
private:
  //int mode; //!< records running mode
  FREESTREAMMILNE *fsmilne_ptr;

  std::string input_file_;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<FreestreamMilneWrapper> reg;

public:
  FreestreamMilneWrapper();
  ~FreestreamMilneWrapper();

  void InitializePreequilibrium(PreEquilibriumParameterFile parameter_list);
  void EvolvePreequilibrium();
};

#endif // FREESTREAMMILNEWRAPPER_H
