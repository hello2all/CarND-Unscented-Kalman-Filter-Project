#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);
  // Polar measurements to Cartesian
  Eigen::VectorXd Polar2Cartesian(const Eigen::VectorXd &z);
  // Calc NIS
  float NIS(const Eigen::VectorXd &z, const Eigen::VectorXd &z_1, const Eigen::MatrixXd &Si);
  double Normalize(double angle);

};

#endif /* TOOLS_H_ */
