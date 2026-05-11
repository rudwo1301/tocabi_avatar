#pragma once
#include <Eigen/Dense>
#include <OsqpEigen/OsqpEigen.h>
#include <fstream>
#include "math_type_define.h"

class SgnPlanner
{
public:
    SgnPlanner(double dt);

    void updateMpcParam(const double& mpc_freq, const double& mpc_dt, const double& mpc_preview_window, const int& mpc_synchro_hz);
    void updateWalkingParam(const double& walking_tick_mpc, const double& com_start_tick_mpc, const int& current_step_num_mpc);
    void getInitWalkingParam(const double& w0, const double& t_start, const double& t_total, const int& total_step_num);
    void calculateMpcMatrix();
    void getPlannerStates(const Eigen::VectorXd& get_planner_state);
    void getFootStepSupportFrame(const Eigen::MatrixXd& get_foot_step_support_frame, const Eigen::MatrixXd& get_foot_step_support_frame_wo_offset);
    void calculateRef(const Eigen::MatrixXd& get_com_ref, const Eigen::MatrixXd& get_alpha_ref);
    void setInitialState(const Eigen::VectorXd& get_planner_state);
    void frameChange(Eigen::VectorXd* state_before_frame_change);

    Eigen::VectorXd runPlanner(const double& mpc_freq, const double& mpc_dt, const double& mpc_preview_window, const int& mpc_synchro_hz);

    int sgn_first_loop_ = 0;

    std::string sgn_file_path = "/home/econom2/data/sgn/";
    static constexpr int sgn_file_num_ = 30;

    bool sgn_data_save_ = true;
    std::ofstream e_sgn_data[sgn_file_num_];

private:
    double sgn_dt_;
    Eigen::Vector3d sgn_gravity_vec_ = Eigen::Vector3d(0.0, 0.0, -9.81);

    //mpc parameters
    double sgn_K1_, sgn_K2_, sgn_K4_, sgn_K5_;
    double sgn_Rw_, sgn_Rh_;

    double sgn_mpc_freq_;
    double sgn_mpc_dt_;
    double sgn_mpc_preview_window_;
    int sgn_mpc_synchro_hz_;
    int sgn_mpc_tick_;
    int sgn_matlab_tick_;
    int sgn_N_plan_mpc_;
    int sgn_state_num_;
    int sgn_input_num_w_;
    int sgn_input_num_h_;
    int sgn_input_num_;

    //walking parameters
    int sgn_current_step_num_mpc_;
    double sgn_w0_;
    double sgn_t_start_mpc_;
    double sgn_t_total_mpc_;
    int sgn_total_step_num_mpc_;
    double sgn_walking_tick_mpc_;

    //mpc variables
    Eigen::VectorXd sgn_planner_state_mpc_;
    Eigen::VectorXd sgn_planner_state_m1_mpc_;
    Eigen::VectorXd sgn_planner_state_p1_mpc_;
    Eigen::MatrixXd sgn_foot_step_support_frame_mpc_;
    Eigen::MatrixXd sgn_foot_step_support_frame_wo_offset_mpc_;
    Eigen::VectorXd sgn_com_ref_vec_;
    Eigen::VectorXd sgn_com_ref0_vec_;
    Eigen::MatrixXd sgn_com_ref_mat_;
    Eigen::MatrixXd sgn_com_ref_left_;
    Eigen::MatrixXd sgn_com_ref_right_;
    Eigen::VectorXd sgn_alpha_ref_vec_;
    Eigen::VectorXd sgn_lfoot_contact_schedule_;
    Eigen::VectorXd sgn_rfoot_contact_schedule_;
    Eigen::MatrixXd sgn_s_x_ref_, sgn_s_y_ref_, sgn_s_z_ref_;
    Eigen::MatrixXd sgn_contact_mask_ref_;
    Eigen::VectorXd sgn_X_bar_;
    Eigen::VectorXd sgn_U_bar_;

    //mpc matrix
    Eigen::SparseMatrix<double> sgn_D_total_;
    Eigen::SparseMatrix<double> sgn_SSz_;
    Eigen::SparseMatrix<double> sgn_SUw_;
    Eigen::SparseMatrix<double> sgn_HD_;
    Eigen::SparseMatrix<double> sgn_Hz_;
    Eigen::SparseMatrix<double> sgn_H_;

    Eigen::SparseMatrix<double> sgn_GX_;
    double sgn_GX_nonzeros_;
    Eigen::SparseMatrix<double> sgn_GU_;
    double sgn_GU_nonzeros_;
    std::vector<Eigen::Triplet<double>> sgn_GX_triplets_;
    std::vector<Eigen::Triplet<double>> sgn_GU_triplets_;

    Eigen::MatrixXd sgn_JX_;
    Eigen::MatrixXd sgn_JU_;

    Eigen::SparseMatrix<double> sgn_JXX_;
    double sgn_JXX_nonzeros_;
    Eigen::SparseMatrix<double> sgn_JUU_;
    double sgn_JUU_nonzeros_;

    Eigen::SparseMatrix<double> sgn_Hessian_;
    double sgn_Hessian_nonzeros_init_;
    double sgn_Hessian_nonzeros_tick_;
    double sgn_Hessian_nonzeros_total_;
    std::vector<Eigen::Triplet<double>> sgn_Hessian_triplets_init_;
    std::vector<Eigen::Triplet<double>> sgn_Hessian_triplets_tick_;
    std::vector<Eigen::Triplet<double>> sgn_Hessian_triplets_total_;

    Eigen::SparseMatrix<double> sgn_grad_X_;
    Eigen::SparseMatrix<double> sgn_grad_U_;
    double sgn_gradient_nonzeros_;
    Eigen::VectorXd sgn_gradient_;
    std::vector<Eigen::Triplet<double>> sgn_gradient_triplets_;

    Eigen::SparseMatrix<double> sgn_dZ_;
    Eigen::SparseMatrix<double> sgn_dU_;
};