#pragma once
#include "mi_equality.h"
#include "mi_impactModel.h"
#include "mi_osd.h"

namespace mc_impact
{

class mi_invOsdEquality : public mi_equality
{
public:
  mi_invOsdEquality(const std::shared_ptr<mi_osd> & osdPtr, const int & numEe);
  ~mi_invOsdEquality() {}

  inline std::string nameEq() const override
  {
    return "InverseOperationSpaceDynamicsImpulseEqualityConstraint";
  }
  void update() override;

protected:
  void reset_() override;
  const int numEe_;
};
} // namespace mc_impact
