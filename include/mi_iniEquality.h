# pragma once 
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
		const std::shared_ptr<mi_impactModel> & impactPtr,
		const int & numEe
	       );
  ~mi_iniEquality(){}
  
  inline std::string nameEq() const override 
  {
    return "InitialConditionEqualityConstraint";
  }
  void update() override;
  inline const std::shared_ptr<mi_impactModel> & getImpactModel()
  {
    return impactPtr_; 
  }
  protected: 
  void reset_() override;

  const std::shared_ptr<mi_impactModel> & impactPtr_;

  const int numEe_;
};
