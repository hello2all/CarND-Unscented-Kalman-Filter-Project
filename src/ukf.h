#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* previous time stamp
  double previous_timestamp_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* predicted measurement sigma points matrix
  MatrixXd Zsig_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  // Lidar measurement dimension
  int n_z_laser_;

  // Radar measurement dimension
  int n_z_radar_;

  // Measurement dimension;
  int n_z_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Augmented state vector
  VectorXd x_aug_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  // Tools
  Tools tools;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
  * Init
  */
  void Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
    Eigen::MatrixXd &H_in, Eigen::MatrixXd &R_laser_in, Eigen::MatrixXd &R_radar_in, Eigen::MatrixXd &Q_in);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  float UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  float UpdateRadar(MeasurementPackage meas_package);

  Eigen::MatrixXd AugmentedSigmaPoints(Eigen::VectorXd x, Eigen::MatrixXd P);

  Eigen::MatrixXd SigmaPointPrediction(MatrixXd Xsig_in, double delta_t);
  void PredictMeanAndCovariance(Eigen::VectorXd* x_out, Eigen::MatrixXd* P_out);
  Eigen::MatrixXd PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out);
  Eigen::MatrixXd PredictLidarMeasurement(VectorXd* z_out, MatrixXd* S_out);
  void LidarUpdateState(Eigen::VectorXd z, Eigen::VectorXd z_pred, Eigen::MatrixXd S, Eigen::MatrixXd Zsig);
  void RadarUpdateState(Eigen::VectorXd z, Eigen::VectorXd z_pred, Eigen::MatrixXd S, Eigen::MatrixXd Zsig);

};

#endif /* UKF_H */
