# pragma once 
# include "mi_osd.h" 
# include "mi_impactModel.h"
# include "mi_equality.h" 

class mi_invOsdEquality: public mi_equality
{
  public: 
  mi_invOsdEquality(
		const std::shared_ptr<mi_osd> & osdPtr
	       );
  ~mi_invOsdEquality(){}
  
  inline std::string nameEq() const override
  {
    return "InverseOperationSpaceDynamicsImpulseEqualityConstraint";
  }
  void update() override;
  protected: 
  void reset_() override;
};

