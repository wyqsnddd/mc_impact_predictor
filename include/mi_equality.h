#pragma once
#include "mi_osd.h"

namespace mc_impact
{

class mi_equality
/** \brief base class for the equality constraints
 */
{
public:
  mi_equality(const std::shared_ptr<mi_osd> & osdPtr) : osdPtr_(osdPtr) {}
  ~mi_equality() {}

  virtual inline std::string nameEq() const
  {
    return "baseEqualityConstraint";
  }
  inline const Eigen::MatrixXd & AEq() const
  {
    return A_;
  }
  inline const Eigen::VectorXd & bEq() const
  {
    return b_;
  }

  inline int nrEq() const
  {
    return static_cast<int>(b_.size());
  }

  /*!
    \param contactEe We use the end-effectors in contact to specify the jsd impulse equation
    */
  virtual void update() = 0;
  inline void printInfo()
  {
    std::cout << nameEq() << " has Aeq size: (" << A_.rows() << ", " << A_.cols() << "), beq size: (" << b_.rows()
              << ", " << b_.cols() << ").  " << std::endl;
  }

protected:
  inline const std::shared_ptr<mi_osd> & getOsd_() const
  {
    return osdPtr_;
  }
  virtual void reset_() = 0;
  const std::shared_ptr<mi_osd> & osdPtr_;
  Eigen::MatrixXd A_;
  Eigen::VectorXd b_;
};
} // namespace mc_impact
