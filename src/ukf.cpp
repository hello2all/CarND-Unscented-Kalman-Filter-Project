#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // initialization state
  is_initialized_ = false;

  // previous time stamp
  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);

  H_ = MatrixXd(2, 5);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

  // Lidar measurement noise covariant matrix
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  // Radar measurement noise covariant matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_* std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  // Measurement dimension;
  n_z_ = n_z_laser_;

  // initial augmented state vector
  x_aug_ = VectorXd(n_aug_);

  // Predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  // the current NIS for radar
  NIS_radar_ = 0;

  // the current NIS for laser
  NIS_laser_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ = tools.Polar2Cartesian(meas_package.raw_measurements_);

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
    }

    // set px, py to very small float number in case they are 0
    if (x_(0) == 0){
      x_(0) = 0.00001;
    }
    if (x_(1) == 0){
      x_(1) = 0.00001;
    }

    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_;
    // cout << "init success" << endl;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

	//compute the time elapsed between the current and previous measurements
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;

	//predict
  // cout << "predict" << endl;
  Prediction(dt);

  // Update
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar updates
    // cout << "update radar" << endl;
    NIS_radar_ = UpdateRadar(meas_package);
  } else if(use_laser_){
    // Laser updates
    // cout << "update lidar" << endl;
    NIS_laser_ = UpdateLidar(meas_package);
  }

  // print the output
  // cout << "x_ = " << x_ << endl;
  // cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  MatrixXd Xsig_aug = AugmentedSigmaPoints(x_, P_);
  Xsig_pred_ = SigmaPointPrediction(Xsig_aug, delta_t);
  PredictMeanAndCovariance(&x_, &P_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
float UKF::UpdateLidar(MeasurementPackage meas_package) {

	VectorXd z_pred = H_ * x_;
	VectorXd y = meas_package.raw_measurements_ - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
  return tools.NIS(meas_package.raw_measurements_, z_pred, S.inverse());
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
float UKF::UpdateRadar(MeasurementPackage meas_package) {

  // normalize measurement in case measurement is not between -pi and pi
  meas_package.raw_measurements_(1) = tools.Normalize(meas_package.raw_measurements_(1));

  n_z_ = n_z_radar_;
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);
  MatrixXd Zsig = PredictRadarMeasurement(&z_pred, &S);
  RadarUpdateState(meas_package.raw_measurements_, z_pred, S, Zsig);
  return tools.NIS(meas_package.raw_measurements_, z_pred, S.inverse());

}
MatrixXd UKF::AugmentedSigmaPoints(VectorXd x, MatrixXd P) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //return augmented sigma points
  return Xsig_aug;

}

MatrixXd UKF::SigmaPointPrediction(MatrixXd Xsig_in, double delta_t) {

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_in(0,i);
    double p_y = Xsig_in(1,i);
    double v = Xsig_in(2,i);
    double yaw = Xsig_in(3,i);
    double yawd = Xsig_in(4,i);
    double nu_a = Xsig_in(5,i);
    double nu_yawdd = Xsig_in(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // angle normalization
    // yaw_p = tools.Normalize(yaw_p);
    // yawd_p = tools.Normalize(yawd_p);

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  //write result
  return Xsig_pred_;

}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3) = tools.Normalize(x_diff(3));

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  *x_out = x;
  *P_out = P;
}

MatrixXd UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = tools.Normalize(z_diff(1));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;

  //write result
  *z_out = z_pred;
  *S_out = S;

  return Zsig;
}


void UKF::RadarUpdateState(VectorXd z, VectorXd z_pred, MatrixXd S, MatrixXd Zsig) {

  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = tools.Normalize(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = tools.Normalize(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = tools.Normalize(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}