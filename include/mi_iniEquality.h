# pragma once 
# include "mi_osd.h" 
# include "mi_equality.h"
# include "mi_impactModel.h"


class mi_iniEquality: public mi_equality
/** \brief initial condition given the predicted impact 
 * J_k \delta q_dot = kth end-effector velocity jump.
 */
{
  public: 
  mi_iniEquality(
		const std::shared_ptr<mi_osd> & osdPtr,
		const mi_impactModel * impactPtr 
	       );
  ~mi_iniEquality(){}
  
  inline std::string nameEq() const override 
  {
    return "InitialConditionEqualityConstraint";
  }
  void update() override;
  protected: 
  void reset_() override;
   
  const mi_impactModel* impactPtr_;
   };
