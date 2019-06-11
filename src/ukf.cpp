#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // size of state vector
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  //  size of augumented state vector
  n_aug_ = n_x_ + 2;

  // spreading parameter used for sigma-point matrix Xsig_aug
  // (see Prediction() below )
  lambda_ = 3 - n_aug_;

  // number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // matrix with predicted sigma points in columns
  Xsig_pred_= MatrixXd(n_x_, n_sig_);

  // weights used for updates
  weights_ = VectorXd(n_sig_);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);

  // dimension of measurement for lidar
  n_z_lidar = 2;

  // dimension of measurement for radar
  n_z_radar = 3;

  // measurement noise covariance matrix for lidar
  R_lidar = MatrixXd(n_z_lidar, n_z_lidar);
  R_lidar <<    std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;

  // measurement noise covariance matrix for radar
  R_radar = MatrixXd(n_z_radar, n_z_radar);
  R_radar <<    std_radr_ * std_radr_, 0, 0,
                0, std_radphi_ * std_radphi_, 0,
                0, 0, std_radrd_ * std_radrd_;

  // whether initialization is done or not
  is_initialized_ = false;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // if initialization is not done yet
  if(!is_initialized_){

    //  initialization of x_
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){

      // for radar
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float v = meas_package.raw_measurements_(2);
      // note that x_ contains p_x, p_y, v, yaw, yaw_dot
      x_ <<  rho * cos(phi), rho * sin(phi), v, 0, 0;

    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER){

      // for lidar
      float p_x = meas_package.raw_measurements_(0);
      float p_y = meas_package.raw_measurements_(1);
      x_ << p_x, p_y, 0, 0, 0;

    }

    // initialization of P_
    P_.fill(0.0);
    for(int i=0; i<n_x_; i++){
      P_(i, i) = 1.0;
    }

    // get time stamp
    time_us_ = meas_package.timestamp_;

    // initialization is done now
    is_initialized_ = true;

    return;

  }

  // compute time lapse (in second)
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  // predicts sigma points, the state vector, and state covariance matrix
  Prediction(dt);

  // update of measurement pack
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    // for radar
    UpdateRadar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    // for lidar
    UpdateLidar(meas_package);

  }


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // create augumented state vector
  // x_aug = (x_, longitudinal_acceleration noise, yaw acceleration noise)
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  // create augmented state covarnant matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) =  std_yawdd_ * std_yawdd_;

  // square root of P_aug
  MatrixXd L = P_aug.llt().matrixL();

  // create sigma-point matrix　
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;

  for(int i = 0; i< n_aug_; i++){

    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);

  }

  // predict sigma points
  for(int i=0; i < n_sig_; i++){

    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // predicted state values
    double px_p, py_p;

    if (fabs(yawd) > 0.001) { // avoid division by zero

        px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));

    }
    else {

        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);

    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    double delta_t2 = delta_t * delta_t;
    px_p += 0.5 * nu_a * delta_t2 * cos(yaw);
    py_p += 0.5 * nu_a * delta_t2 * sin(yaw);
    v_p += nu_a * delta_t;
    yaw_p += 0.5 * nu_yawdd * delta_t2;
    yawd_p += nu_yawdd * delta_t;

    //write the predicted sigma point into the columns of Xsig_pred_
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

  }

  // predict state mean
  x_.fill(0.0);
  for (int i=0; i < n_sig_; i++) {  //iterate over sigma points
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i=0; i < n_sig_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    Angle_Normalize(x_diff(3));

    P_ += weights_(i) * x_diff * x_diff.transpose();

  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //sigma-point matrix in measurement space
  MatrixXd Zsig = MatrixXd(n_z_lidar, n_sig_);

  //transform sigma points into measurement space
  for(int i=0; i < n_sig_; i++){

    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);

  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_lidar);
  z_pred.fill(0.0);
  for (int i=0; i < n_sig_; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_lidar, n_z_lidar);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();

  }
  //add measurement noise covariance matrix to S
  S += R_lidar;

  // matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_lidar);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++){

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    Angle_Normalize(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();

  }

  //Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // radar measurement
  VectorXd z = meas_package.raw_measurements_;

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += - K * S * K.transpose();

  // compute NIS
  NIS_lidar = z_diff.transpose() * S.inverse() * z_diff;
  // cout << NIS_lidar << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */

  //matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar, n_sig_);
  //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v  = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v_x = v * cos(yaw);
    double v_y = v * sin(yaw);

    double p_abs = sqrt(p_x * p_x + p_y * p_y);
    Zsig(0,i) = p_abs;  //rho
    Zsig(1,i) = atan2(p_y, p_x); //phi
    Zsig(2,i) = (p_x * v_x + p_y * v_y) / p_abs;  //rho_dot

  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar);
  z_pred.fill(0.0);
  for (int i=0; i < n_sig_; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar, n_z_radar);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    Angle_Normalize(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();

  }
  //add measurement noise covariance matrix to S
  S += R_radar;

  // matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++){

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    Angle_Normalize(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    Angle_Normalize(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();

  }

  //Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // radar measurement
  VectorXd z = meas_package.raw_measurements_;

  //residual
  VectorXd z_diff = z - z_pred;
  //angle normalization
  Angle_Normalize(z_diff(1));

  //update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += - K * S * K.transpose();

  // compute NIS
  NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
  // cout << NIS_radar << endl;

}

// for angle normalization into range [-pi, pi]
void UKF::Angle_Normalize(double angle){

  while(angle > M_PI) angle -= 2.0 * M_PI;
  while(angle < -M_PI) angle += 2.0 * M_PI;

}
