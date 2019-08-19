# pragma once 
# include "mi_osd.h" 
# include "mi_equality.h"

class mi_jsdEquality: public mi_equality
/** \brief Specify the joint space impulse equation: [M, 0, -J_i^T, ...] [\delta q_dot, I_1,\ldots, I_m ]^T
 */
{
  public: 
  mi_jsdEquality(
		const std::shared_ptr<mi_osd> & osdPtr,
		const std::map<std::string, Eigen::Vector3d> & impactNameAndNormals
		//const std::vector<std::string> & impactBodies
	       );
  ~mi_jsdEquality(){}
  
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
  std::vector<std::string> impactBodyNames_;
};

