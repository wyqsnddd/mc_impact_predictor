# pragma once 
# include "mi_osd.h" 
# include "mi_impactModel.h"
# include "mi_equality.h" 

class mi_osdEquality: public mi_equality
/** \brief Specify the operation space impuplse equation: [\sum lambda_ki * J_i, 0,..., -I, ..] [\delta_q_dot, I_1,\ldots, I_m]^T.
 * for all the m end-effectors.
 */
{
  public: 
  mi_osdEquality(
		const std::shared_ptr<mi_osd> & osdPtr
	       );
  ~mi_osdEquality(){}
  
  inline std::string nameEq() const override
  {
    return "OperationSpaceDynamicsImpulseEqualityConstraint";
  }
  void update() override;
  protected: 
  void reset_() override;
};

