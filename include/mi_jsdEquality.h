#pragma once
#include "mi_equality.h"
#include "mi_impactModel.h"
#include "mi_osd.h"
#include "mi_utils.h"

namespace mc_impact
{

class mi_jsdEquality : public mi_equality
/** \brief Specify the joint space impulse equation: [M, 0, -J_i^T, ...] [\delta q_dot, I_1,\ldots, I_m ]^T
 */
{
public:
  mi_jsdEquality(const std::shared_ptr<mi_osd> & osdPtr,
                 const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels,
                 const std::map<std::string, endEffector> & endEffectors);
  ~mi_jsdEquality() {}

  inline std::string nameEq() const override
  {
    return "JointSpaceDynamicsImpulseEqualityConstraint";
  }

  /*!
    \param contactEe We use the end-effectors in contact to specify the jsd impulse equation
    */
  void update() override;

protected:
  void reset_() override;
  const std::map<std::string, std::shared_ptr<mi_impactModel>> & impactModels_;
  const std::map<std::string, endEffector> & endEffectors_;
  int nameToIndex_(const std::string & eeName);
};

} // namespace mc_impact
