#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
          || estimations.size() == 0){
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

VectorXd Tools::Polar2Cartesian(const Eigen::VectorXd &z){
  float rho = z(0);
  float phi = z(1);
  float rho_dot = z(2);

  float px = rho * cos(phi);
  float py = rho * sin(phi);
  float v = rho_dot;
  float phi_dot = 0;

  VectorXd x_ = VectorXd(5);
  x_ << px, py, v, phi, phi_dot;

  return x_;
}

float Tools::NIS(const Eigen::VectorXd &z, const Eigen::VectorXd &z_1, const Eigen::MatrixXd &Si){
  Eigen::VectorXd diff = z_1 - z;
  return diff.transpose() * Si * diff;
}
