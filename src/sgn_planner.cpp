#include "sgn_planner.h"
#include <iostream>

using namespace std;

// 생성자
SgnPlanner::SgnPlanner(double dt)
    : sgn_dt_(dt)
{
    if(sgn_dt_ > 1e10)
    {
        sgn_dt_ = 1/2000.0;
    }

    sgn_data_save_ = true;
    
    for(int i = 0; i < sgn_file_num_; i++)
    {
        e_sgn_data[i].open(sgn_file_path + "e_sgn_data" + to_string(i) + ".txt");
    }
}

void SgnPlanner::updateMpcParam(const double& mpc_freq, const double& mpc_dt, const double& mpc_preview_window, const int& mpc_synchro_hz)
{
    sgn_K1_ = 1e0;
    sgn_K2_ = 1e0;
    sgn_K4_ = 1e2;
    sgn_K5_ = 1e0;
    sgn_Rw_ = 1e-2;
    sgn_Rh_ = 1e-2;

    sgn_mpc_freq_ = mpc_freq;
    sgn_mpc_dt_ = mpc_dt;
    sgn_mpc_preview_window_ = mpc_preview_window;
    sgn_mpc_preview_window_ = 1.6;
    sgn_mpc_synchro_hz_ = mpc_synchro_hz;
    sgn_N_plan_mpc_ = sgn_mpc_preview_window_ * sgn_mpc_freq_;
    sgn_state_num_ = 3;
    sgn_input_num_w_ = 8;
    sgn_input_num_h_ = 1;
    sgn_input_num_ = sgn_input_num_w_ + sgn_input_num_h_;
}

void SgnPlanner::updateWalkingParam(const double& walking_tick_mpc, const double& com_start_tick_mpc, const int& current_step_num_mpc)
{
    sgn_matlab_tick_ = int(walking_tick_mpc - floor(sgn_mpc_dt_/sgn_dt_) + 1)/floor(sgn_mpc_dt_/sgn_dt_) + 1;
    cout << "sgn_matlab_tick_: " << sgn_matlab_tick_ << endl;

    sgn_walking_tick_mpc_ = walking_tick_mpc;
    sgn_mpc_tick_ = walking_tick_mpc - com_start_tick_mpc;
    sgn_current_step_num_mpc_ = current_step_num_mpc;
}

void SgnPlanner::getInitWalkingParam(const double& w0, const double& t_start, const double& t_total, const int& total_step_num)
{
    sgn_w0_ = w0;
    sgn_t_start_mpc_ = t_start;
    sgn_t_total_mpc_ = t_total;
    sgn_total_step_num_mpc_ = total_step_num;
}

void SgnPlanner::calculateMpcMatrix()
{
    sgn_D_total_.resize(sgn_state_num_*(sgn_N_plan_mpc_ - 1), sgn_state_num_*sgn_N_plan_mpc_);
    sgn_SSz_.resize(sgn_N_plan_mpc_, sgn_state_num_*sgn_N_plan_mpc_);

    for(int i = 0; i < sgn_N_plan_mpc_ - 1; i++)
    {           
        for(int j = 0; j < sgn_state_num_; j++)
        {
            sgn_D_total_.insert(sgn_state_num_ * i + j, sgn_state_num_ *  i + j) = -1.0;
            sgn_D_total_.insert(sgn_state_num_ * i + j, sgn_state_num_ * (i+1) + j) =  1.0;
        }

        sgn_SSz_.insert(i, sgn_state_num_*(i + 1) - 1 + 0) = 1.0;
    }

    sgn_SSz_.insert(sgn_N_plan_mpc_ - 1, sgn_state_num_*sgn_N_plan_mpc_ - 1) = 1.0;

    sgn_SUw_.resize(sgn_input_num_w_ * sgn_N_plan_mpc_, sgn_input_num_ * sgn_N_plan_mpc_);
    
    std::vector<Eigen::Triplet<double>> SUw_triplets;
    SUw_triplets.reserve(sgn_input_num_w_ * sgn_N_plan_mpc_);

    for(int k = 0; k < sgn_N_plan_mpc_; k++)
    {
        int row_offset = sgn_input_num_w_ * k;
        int col_offset = sgn_input_num_ * k;

        for(int i = 0; i < sgn_input_num_w_; i++)
        {
            SUw_triplets.emplace_back(row_offset + i, col_offset + i, 1.0);
        }
    }

    sgn_SUw_.setFromTriplets(SUw_triplets.begin(), SUw_triplets.end());

    sgn_HD_ = 2*sgn_K1_*(sgn_D_total_.transpose() * sgn_D_total_);
    sgn_Hz_ = 2*sgn_K2_*(sgn_SSz_.transpose() * sgn_SSz_);
    for(int i = 0; i < sgn_state_num_; i++)
    {
        //sgn_HD_.coeffRef(i,i) += 2*sgn_K1_;
    }

    sgn_H_  = sgn_HD_ + sgn_Hz_;
}

void SgnPlanner::getPlannerStates(const Eigen::VectorXd& get_planner_state)
{
    sgn_planner_state_mpc_.resize(sgn_state_num_);
    for (int i = 0; i < sgn_state_num_; i++)
    {
        sgn_planner_state_mpc_(i) = get_planner_state(3*i);
    }
}

void SgnPlanner::setInitialState(const Eigen::VectorXd& get_planner_state)
{
    if(sgn_first_loop_ == 0)
    {
        calculateMpcMatrix();
        getPlannerStates(get_planner_state);

        sgn_planner_state_m1_mpc_ = sgn_planner_state_mpc_;
        sgn_planner_state_p1_mpc_ = sgn_planner_state_mpc_;
        
        sgn_GX_.resize(sgn_state_num_ * sgn_N_plan_mpc_, sgn_state_num_ * sgn_N_plan_mpc_);
        sgn_GU_.resize(sgn_state_num_ * sgn_N_plan_mpc_, sgn_input_num_ * sgn_N_plan_mpc_);

        sgn_X_bar_.setZero(sgn_state_num_ * sgn_N_plan_mpc_);
        sgn_U_bar_.setZero(sgn_input_num_ * sgn_N_plan_mpc_);

        sgn_GX_nonzeros_ = sgn_state_num_ *  sgn_N_plan_mpc_
                         + sgn_state_num_ * (sgn_N_plan_mpc_ - 1)
                         + sgn_state_num_ * (sgn_N_plan_mpc_ - 2)
                         + sgn_N_plan_mpc_ - 1
                         + sgn_N_plan_mpc_ - 1;
        
        sgn_GX_.reserve(sgn_GX_nonzeros_);

        sgn_GU_nonzeros_ = (2 * sgn_input_num_w_ + sgn_state_num_ * sgn_input_num_h_) * sgn_N_plan_mpc_;

        sgn_GU_.reserve(sgn_GU_nonzeros_);

        sgn_JXX_nonzeros_ = sgn_state_num_ *  sgn_N_plan_mpc_
                          + sgn_state_num_ * (sgn_N_plan_mpc_ - 1)
                          + sgn_state_num_ * (sgn_N_plan_mpc_ - 1);

        sgn_JXX_.reserve(sgn_JXX_nonzeros_);

        sgn_JUU_nonzeros_ = (sgn_input_num_w_ * sgn_input_num_w_ + sgn_input_num_h_ * sgn_input_num_h_) * sgn_N_plan_mpc_;

        sgn_JUU_.reserve(sgn_JUU_nonzeros_);

        sgn_Hessian_nonzeros_init_ = sgn_JXX_nonzeros_ + sgn_JUU_nonzeros_;

        sgn_Hessian_nonzeros_tick_ = 2*sgn_GX_nonzeros_ + 2*sgn_GU_nonzeros_;

        sgn_Hessian_nonzeros_total_ = sgn_Hessian_nonzeros_init_ + sgn_Hessian_nonzeros_tick_;

        sgn_Hessian_triplets_init_.reserve(sgn_Hessian_nonzeros_init_);

        sgn_Hessian_triplets_tick_.reserve(sgn_Hessian_nonzeros_tick_);

        sgn_Hessian_.reserve(sgn_Hessian_nonzeros_total_);

        for(int i = 0; i < sgn_state_num_ * sgn_N_plan_mpc_; i++)
        {
            sgn_Hessian_triplets_init_.emplace_back(i, i, sgn_H_.coeff(i, i)); //sgn_JXX_ = sgn_H_;

            if (i < sgn_state_num_ * (sgn_N_plan_mpc_ - 1))
            {
                sgn_Hessian_triplets_init_.emplace_back(i + sgn_state_num_, i, sgn_H_.coeff(i + sgn_state_num_, i)); //sgn_JXX_ = sgn_H_;
                sgn_Hessian_triplets_init_.emplace_back(i, i + sgn_state_num_, sgn_H_.coeff(i, i + sgn_state_num_)); //sgn_JXX_ = sgn_H_;
            }
        }

        for(int i = sgn_state_num_ * sgn_N_plan_mpc_; i < (sgn_state_num_ + sgn_input_num_) * sgn_N_plan_mpc_; i++)
        {
            sgn_Hessian_triplets_init_.emplace_back(i, i, 2*sgn_Rw_); //sgn_JUU_ = 2*sgn_Rw_ * Identity + ...;
        }

        sgn_Hessian_.resize((sgn_state_num_ + sgn_input_num_ + sgn_state_num_) * sgn_N_plan_mpc_, (sgn_state_num_ + sgn_input_num_ + sgn_state_num_) * sgn_N_plan_mpc_);

        sgn_gradient_.setZero((sgn_state_num_ + sgn_input_num_ + sgn_state_num_) * sgn_N_plan_mpc_);

        sgn_first_loop_ = 1;
    }
}

void SgnPlanner::getFootStepSupportFrame(const Eigen::MatrixXd& get_foot_step_support_frame, const Eigen::MatrixXd& get_foot_step_support_frame_wo_offset)
{
    sgn_foot_step_support_frame_mpc_           = get_foot_step_support_frame;
    sgn_foot_step_support_frame_wo_offset_mpc_ = get_foot_step_support_frame_wo_offset;
}

void SgnPlanner::calculateRef(const Eigen::MatrixXd& get_com_ref, const Eigen::MatrixXd& get_alpha_ref)
{
    sgn_com_ref_vec_.setZero(sgn_state_num_ * sgn_N_plan_mpc_);

    sgn_com_ref_mat_.setZero(sgn_N_plan_mpc_, sgn_state_num_);

    sgn_alpha_ref_vec_.setZero(sgn_N_plan_mpc_);

    sgn_com_ref_left_ .setZero(sgn_N_plan_mpc_, sgn_state_num_);

    sgn_com_ref_right_.setZero(sgn_N_plan_mpc_, sgn_state_num_);

    sgn_com_ref0_vec_.setZero(sgn_state_num_);

    for(int i = 0; i < sgn_N_plan_mpc_; i++)
    {
        sgn_alpha_ref_vec_(i) = get_alpha_ref(sgn_mpc_tick_ + sgn_mpc_synchro_hz_*(i+1),0);

        for(int j = 0; j < sgn_state_num_; j++)
        {
            sgn_com_ref_vec_(sgn_state_num_*i + j) = get_com_ref(sgn_mpc_tick_ + sgn_mpc_synchro_hz_*(i+1),j);

            sgn_com_ref_mat_(i, j) = sgn_com_ref_vec_(sgn_state_num_*i + j);

            if(i == 0)
            {
                sgn_com_ref0_vec_(j) = get_com_ref(sgn_mpc_tick_ + sgn_mpc_synchro_hz_,j);

                sgn_com_ref_left_ .col(j).setConstant(max(sgn_foot_step_support_frame_wo_offset_mpc_(sgn_current_step_num_mpc_, j), 0.0));

                sgn_com_ref_right_.col(j).setConstant(min(sgn_foot_step_support_frame_wo_offset_mpc_(sgn_current_step_num_mpc_, j), 0.0));
            }
        }
    }

    sgn_lfoot_contact_schedule_ = sgn_alpha_ref_vec_.array().ceil();
    sgn_rfoot_contact_schedule_ = (Eigen::VectorXd::Constant(sgn_alpha_ref_vec_.size(), 1.0) - sgn_alpha_ref_vec_).array().ceil();

    sgn_contact_mask_ref_.resize(sgn_N_plan_mpc_, sgn_input_num_w_);
    sgn_contact_mask_ref_.col(0) = sgn_lfoot_contact_schedule_;
    sgn_contact_mask_ref_.col(1) = sgn_lfoot_contact_schedule_;
    sgn_contact_mask_ref_.col(2) = sgn_lfoot_contact_schedule_;
    sgn_contact_mask_ref_.col(3) = sgn_lfoot_contact_schedule_;
    sgn_contact_mask_ref_.col(4) = sgn_rfoot_contact_schedule_;
    sgn_contact_mask_ref_.col(5) = sgn_rfoot_contact_schedule_;
    sgn_contact_mask_ref_.col(6) = sgn_rfoot_contact_schedule_;
    sgn_contact_mask_ref_.col(7) = sgn_rfoot_contact_schedule_;

    double foot_width = 0.08, foot_front = 0.17, foot_back = 0.11;

    Eigen::VectorXd s_l_x_ref_front_left, s_l_x_ref_front_right, s_l_x_ref_back_left, s_l_x_ref_back_right;
    Eigen::VectorXd s_l_y_ref_front_left, s_l_y_ref_front_right, s_l_y_ref_back_left, s_l_y_ref_back_right;
    Eigen::VectorXd s_l_z_ref_front_left, s_l_z_ref_front_right, s_l_z_ref_back_left, s_l_z_ref_back_right;

    const int rows_size = sgn_com_ref_left_.rows();

    s_l_x_ref_front_left  = sgn_com_ref_left_.col(0)  + Eigen::VectorXd::Constant(rows_size, foot_front);
    s_l_x_ref_back_left   = sgn_com_ref_left_.col(0)  - Eigen::VectorXd::Constant(rows_size, foot_back);
    s_l_x_ref_back_right  = sgn_com_ref_left_.col(0)  - Eigen::VectorXd::Constant(rows_size, foot_back);
    s_l_x_ref_front_right = sgn_com_ref_left_.col(0)  + Eigen::VectorXd::Constant(rows_size, foot_front);

    s_l_y_ref_front_left  = sgn_com_ref_left_.col(1)  + Eigen::VectorXd::Constant(rows_size, foot_width);
    s_l_y_ref_back_left   = sgn_com_ref_left_.col(1)  + Eigen::VectorXd::Constant(rows_size, foot_width);
    s_l_y_ref_back_right  = sgn_com_ref_left_.col(1)  - Eigen::VectorXd::Constant(rows_size, foot_width);
    s_l_y_ref_front_right = sgn_com_ref_left_.col(1)  - Eigen::VectorXd::Constant(rows_size, foot_width);

    s_l_z_ref_front_left  = sgn_com_ref_left_.col(2);
    s_l_z_ref_back_left   = sgn_com_ref_left_.col(2);
    s_l_z_ref_back_right  = sgn_com_ref_left_.col(2);
    s_l_z_ref_front_right = sgn_com_ref_left_.col(2);

    Eigen::VectorXd s_r_x_ref_front_left, s_r_x_ref_front_right, s_r_x_ref_back_left, s_r_x_ref_back_right;
    Eigen::VectorXd s_r_y_ref_front_left, s_r_y_ref_front_right, s_r_y_ref_back_left, s_r_y_ref_back_right;
    Eigen::VectorXd s_r_z_ref_front_left, s_r_z_ref_front_right, s_r_z_ref_back_left, s_r_z_ref_back_right;

    s_r_x_ref_front_left  = sgn_com_ref_right_.col(0) + Eigen::VectorXd::Constant(rows_size, foot_front);
    s_r_x_ref_back_left   = sgn_com_ref_right_.col(0) - Eigen::VectorXd::Constant(rows_size, foot_back);
    s_r_x_ref_back_right  = sgn_com_ref_right_.col(0) - Eigen::VectorXd::Constant(rows_size, foot_back);
    s_r_x_ref_front_right = sgn_com_ref_right_.col(0) + Eigen::VectorXd::Constant(rows_size, foot_front);

    s_r_y_ref_front_left  = sgn_com_ref_right_.col(1) + Eigen::VectorXd::Constant(rows_size, foot_width);
    s_r_y_ref_back_left   = sgn_com_ref_right_.col(1) + Eigen::VectorXd::Constant(rows_size, foot_width);
    s_r_y_ref_back_right  = sgn_com_ref_right_.col(1) - Eigen::VectorXd::Constant(rows_size, foot_width);
    s_r_y_ref_front_right = sgn_com_ref_right_.col(1) - Eigen::VectorXd::Constant(rows_size, foot_width);

    s_r_z_ref_front_left  = sgn_com_ref_right_.col(2);
    s_r_z_ref_back_left   = sgn_com_ref_right_.col(2);
    s_r_z_ref_back_right  = sgn_com_ref_right_.col(2);
    s_r_z_ref_front_right = sgn_com_ref_right_.col(2);
    
    sgn_s_x_ref_.resize(sgn_N_plan_mpc_, sgn_input_num_w_);
    sgn_s_y_ref_.resize(sgn_N_plan_mpc_, sgn_input_num_w_);
    sgn_s_z_ref_.resize(sgn_N_plan_mpc_, sgn_input_num_w_);

    sgn_s_x_ref_.col(0) = sgn_lfoot_contact_schedule_.array() * s_l_x_ref_front_left.array();
    sgn_s_x_ref_.col(1) = sgn_lfoot_contact_schedule_.array() * s_l_x_ref_back_left.array();
    sgn_s_x_ref_.col(2) = sgn_lfoot_contact_schedule_.array() * s_l_x_ref_back_right.array();
    sgn_s_x_ref_.col(3) = sgn_lfoot_contact_schedule_.array() * s_l_x_ref_front_right.array();
    sgn_s_x_ref_.col(4) = sgn_rfoot_contact_schedule_.array() * s_r_x_ref_front_left.array();
    sgn_s_x_ref_.col(5) = sgn_rfoot_contact_schedule_.array() * s_r_x_ref_back_left.array();
    sgn_s_x_ref_.col(6) = sgn_rfoot_contact_schedule_.array() * s_r_x_ref_back_right.array();
    sgn_s_x_ref_.col(7) = sgn_rfoot_contact_schedule_.array() * s_r_x_ref_front_right.array();

    sgn_s_y_ref_.col(0) = sgn_lfoot_contact_schedule_.array() * s_l_y_ref_front_left.array();
    sgn_s_y_ref_.col(1) = sgn_lfoot_contact_schedule_.array() * s_l_y_ref_back_left.array();
    sgn_s_y_ref_.col(2) = sgn_lfoot_contact_schedule_.array() * s_l_y_ref_back_right.array();
    sgn_s_y_ref_.col(3) = sgn_lfoot_contact_schedule_.array() * s_l_y_ref_front_right.array();
    sgn_s_y_ref_.col(4) = sgn_rfoot_contact_schedule_.array() * s_r_y_ref_front_left.array();
    sgn_s_y_ref_.col(5) = sgn_rfoot_contact_schedule_.array() * s_r_y_ref_back_left.array();
    sgn_s_y_ref_.col(6) = sgn_rfoot_contact_schedule_.array() * s_r_y_ref_back_right.array();
    sgn_s_y_ref_.col(7) = sgn_rfoot_contact_schedule_.array() * s_r_y_ref_front_right.array();

    sgn_s_z_ref_.col(0) = sgn_lfoot_contact_schedule_.array() * s_l_z_ref_front_left.array();
    sgn_s_z_ref_.col(1) = sgn_lfoot_contact_schedule_.array() * s_l_z_ref_back_left.array();
    sgn_s_z_ref_.col(2) = sgn_lfoot_contact_schedule_.array() * s_l_z_ref_back_right.array();
    sgn_s_z_ref_.col(3) = sgn_lfoot_contact_schedule_.array() * s_l_z_ref_front_right.array();
    sgn_s_z_ref_.col(4) = sgn_rfoot_contact_schedule_.array() * s_r_z_ref_front_left.array();
    sgn_s_z_ref_.col(5) = sgn_rfoot_contact_schedule_.array() * s_r_z_ref_back_left.array();
    sgn_s_z_ref_.col(6) = sgn_rfoot_contact_schedule_.array() * s_r_z_ref_back_right.array();
    sgn_s_z_ref_.col(7) = sgn_rfoot_contact_schedule_.array() * s_r_z_ref_front_right.array();

    if(sgn_data_save_)
    {
        e_sgn_data[0]  << sgn_com_ref_mat_.  col(1).transpose()   << endl;
        e_sgn_data[1]  << sgn_com_ref_left_. col(1).transpose()   << endl;
        e_sgn_data[2]  << sgn_com_ref_right_.col(1).transpose()   << endl;
        e_sgn_data[3]  << sgn_alpha_ref_vec_.transpose()          << endl;
        e_sgn_data[4]  << sgn_lfoot_contact_schedule_.transpose() << endl;
        e_sgn_data[5]  << sgn_rfoot_contact_schedule_.transpose() << endl;
        e_sgn_data[6]  << sgn_s_y_ref_.col(0).transpose() << endl;
        e_sgn_data[7]  << sgn_s_y_ref_.col(1).transpose() << endl;
        e_sgn_data[8]  << sgn_s_y_ref_.col(2).transpose() << endl;
        e_sgn_data[9]  << sgn_s_y_ref_.col(3).transpose() << endl;
        e_sgn_data[10] << sgn_s_y_ref_.col(4).transpose() << endl;
        e_sgn_data[11] << sgn_s_y_ref_.col(5).transpose() << endl;
        e_sgn_data[12] << sgn_s_y_ref_.col(6).transpose() << endl;
        e_sgn_data[13] << sgn_s_y_ref_.col(7).transpose() << endl;
        e_sgn_data[14] << sgn_w0_ << "," << 9.81/(sgn_w0_*sgn_w0_) << endl;
        e_sgn_data[15] << sgn_com_ref_mat_.  col(2).transpose() << endl;
    }
}

Eigen::VectorXd SgnPlanner::runPlanner(const double& mpc_freq, const double& mpc_dt, const double& mpc_preview_window, const int& mpc_synchro_hz)
{   
    sgn_X_bar_.setZero();
    
    sgn_GX_triplets_.clear();
    sgn_GX_triplets_.reserve(sgn_GX_nonzeros_);

    sgn_GU_triplets_.clear();
    sgn_GU_triplets_.reserve(sgn_GU_nonzeros_);

    sgn_Hessian_triplets_tick_.clear();
    sgn_Hessian_triplets_tick_.reserve(sgn_Hessian_nonzeros_tick_);

    sgn_Hessian_triplets_total_.clear();
    sgn_Hessian_triplets_total_.reserve(sgn_Hessian_nonzeros_total_);

    Eigen::MatrixXd s_ref_calc;
    s_ref_calc.resize(sgn_state_num_, sgn_input_num_w_);

    sgn_contact_mask_ref_.setZero(sgn_N_plan_mpc_, sgn_input_num_w_*sgn_N_plan_mpc_);

    Eigen::MatrixXd sgn_contact_mask_calc;
    sgn_contact_mask_calc.resize(sgn_N_plan_mpc_, sgn_input_num_w_);
    sgn_contact_mask_calc.setZero();

    int FIP_or_VHIP = 1; // 0: FIP, 1: VHIP
    
    for (int i = 0; i < sgn_N_plan_mpc_; i++)
    {
        s_ref_calc.row(0) = sgn_s_x_ref_.row(i);
        s_ref_calc.row(1) = sgn_s_y_ref_.row(i);
        s_ref_calc.row(2) = sgn_s_z_ref_.row(i);

        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 0) = sgn_lfoot_contact_schedule_(i);
        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 1) = sgn_lfoot_contact_schedule_(i);
        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 2) = sgn_lfoot_contact_schedule_(i);
        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 3) = sgn_lfoot_contact_schedule_(i);
        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 4) = sgn_rfoot_contact_schedule_(i);
        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 5) = sgn_rfoot_contact_schedule_(i);
        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 6) = sgn_rfoot_contact_schedule_(i);
        sgn_contact_mask_ref_(i, sgn_input_num_w_ * i + 7) = sgn_rfoot_contact_schedule_(i);

        if(i == 0)
        {
            if(FIP_or_VHIP == 0)
            {
                sgn_X_bar_.segment(sgn_state_num_ * i, sgn_state_num_) = 2*sgn_planner_state_mpc_
                                                                       - sgn_planner_state_m1_mpc_
                                                                       + sgn_w0_*sgn_w0_*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_planner_state_mpc_ - s_ref_calc*sgn_U_bar_.segment(sgn_input_num_ * i, sgn_input_num_w_))
                                                                       + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
            }
            else
            {
                sgn_X_bar_.segment(sgn_state_num_ * i, sgn_state_num_) = 2*sgn_planner_state_mpc_
                                                                       - sgn_planner_state_m1_mpc_
                                                                       + ((GRAVITY + sgn_U_bar_(sgn_input_num_*(i+1)-1))/sgn_planner_state_mpc_(sgn_state_num_-1))*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_planner_state_mpc_ - s_ref_calc*sgn_U_bar_.segment(sgn_input_num_ * i, sgn_input_num_w_))
                                                                       + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
            }
        }
        else if (i == 1)
        {
            if(FIP_or_VHIP == 0)
            {
                sgn_X_bar_.segment(sgn_state_num_ * i, sgn_state_num_) = 2*sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_)
                                                                       - sgn_planner_state_mpc_
                                                                       + sgn_w0_*sgn_w0_*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_) - s_ref_calc*sgn_U_bar_.segment(sgn_input_num_ * i, sgn_input_num_w_))
                                                                       + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
            }
            else
            {
                sgn_X_bar_.segment(sgn_state_num_ * i, sgn_state_num_) = 2*sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_)
                                                                       - sgn_planner_state_mpc_
                                                                       + ((GRAVITY + sgn_U_bar_(sgn_input_num_*(i+1)-1))/sgn_X_bar_(sgn_state_num_ * i - 1))*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_) - s_ref_calc*sgn_U_bar_.segment(sgn_input_num_ * i, sgn_input_num_w_))
                                                                       + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
            }
        }
        else
        {
            if(FIP_or_VHIP == 0)
            {
                sgn_X_bar_.segment(sgn_state_num_ * i, sgn_state_num_) = 2*sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_)
                                                                       - sgn_X_bar_.segment(sgn_state_num_ * (i-2), sgn_state_num_)
                                                                       + sgn_w0_*sgn_w0_*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_) - s_ref_calc*sgn_U_bar_.segment(sgn_input_num_ * i, sgn_input_num_w_))
                                                                       + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
            }
            else
            {
                sgn_X_bar_.segment(sgn_state_num_ * i, sgn_state_num_) = 2*sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_)
                                                                       - sgn_X_bar_.segment(sgn_state_num_ * (i-2), sgn_state_num_)
                                                                       + ((GRAVITY + sgn_U_bar_(sgn_input_num_*(i+1)-1))/sgn_X_bar_(sgn_state_num_ * i - 1))*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_X_bar_.segment(sgn_state_num_ * (i-1), sgn_state_num_) - s_ref_calc*sgn_U_bar_.segment(sgn_input_num_ * i, sgn_input_num_w_))
                                                                       + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
            }
        }

        double com_x   = sgn_X_bar_(sgn_state_num_*i + 0);
        double com_y   = sgn_X_bar_(sgn_state_num_*i + 1);
        double com_z   = sgn_X_bar_(sgn_state_num_*i + 2);

        double input_zx = s_ref_calc.row(0)*sgn_U_bar_.segment(sgn_input_num_*i, sgn_input_num_w_);
        double input_zy = s_ref_calc.row(1)*sgn_U_bar_.segment(sgn_input_num_*i, sgn_input_num_w_);
        double input_zz = s_ref_calc.row(2)*sgn_U_bar_.segment(sgn_input_num_*i, sgn_input_num_w_);

        double input_hdd = sgn_U_bar_(sgn_input_num_*i + sgn_input_num_ - 1);

        int GX_row_index = 0;
        int GX_col_index = 0;
        int GX_to_Hessian_row_offset = (sgn_state_num_ + sgn_input_num_) * sgn_N_plan_mpc_;
        int GX_to_Hessian_col_offset = 0;
        double GX_value = 0.0;

        int GU_row_index = 0;
        int GU_col_index = 0;
        int GU_to_Hessian_row_offset = (sgn_state_num_ + sgn_input_num_) * sgn_N_plan_mpc_;
        int GU_to_Hessian_col_offset =  sgn_state_num_                   * sgn_N_plan_mpc_;
        double GU_value = 0.0;

        //GX tick
        for (int j = 0; j < sgn_state_num_; j++)
        {
            GX_row_index = sgn_state_num_ * i + j;
            GX_col_index = sgn_state_num_ * i + j;

            GX_value = 1.0;

            sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index, GX_value);
            //sgn_Hessian GX  part
            sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index, GX_value);
            //sgn_Hessian GX' part
            sgn_Hessian_triplets_tick_.emplace_back(GX_col_index, GX_to_Hessian_row_offset + GX_row_index, GX_value);

            if (i < sgn_N_plan_mpc_ - 1)
            {
                //FIP
                if(FIP_or_VHIP == 0)
                {
                    GX_row_index = sgn_state_num_ * (i + 1) + j;
                    GX_col_index = sgn_state_num_ *  i      + j;
                    GX_value = - (2 + sgn_w0_ * sgn_w0_ * sgn_mpc_dt_ * sgn_mpc_dt_);

                    sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index, GX_value);
                    //sgn_Hessian GX  part
                    sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index, GX_value);
                    //sgn_Hessian GX' part
                    sgn_Hessian_triplets_tick_.emplace_back(GX_col_index, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                }

                //VHIP
                if(FIP_or_VHIP == 1)
                {
                    GX_row_index = sgn_state_num_ * (i + 1) + j;
                    GX_col_index = sgn_state_num_ *  i;
                    if(j == 0)
                    {
                        //x axis
                        GX_value = - 2 - ((GRAVITY + input_hdd)/com_z) * sgn_mpc_dt_ * sgn_mpc_dt_;
                        sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 0, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 0, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 0, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                        //y axis
                        //GX_value = 0.0;
                        //sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 1, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 1, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 1, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                        //z axis
                        GX_value = (com_x - input_zx) * ((GRAVITY + input_hdd)/(com_z*com_z)) * sgn_mpc_dt_ * sgn_mpc_dt_;
                        sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 2, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 2, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 2, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                    }
                    else if(j == 1)
                    {   
                        //x axis
                        //GX_value = 0.0;
                        //sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 0, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 0, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 0, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                        //y axis
                        GX_value = - 2 - ((GRAVITY + input_hdd)/com_z) * sgn_mpc_dt_ * sgn_mpc_dt_;
                        sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 1, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 1, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 1, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                        //z axis
                        GX_value = (com_y - input_zy) * ((GRAVITY + input_hdd)/(com_z*com_z)) * sgn_mpc_dt_ * sgn_mpc_dt_;
                        sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 2, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 2, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 2, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                    }
                    else if(j == 2)
                    {
                        //x axis
                        //GX_value = 0.0;
                        //sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 0, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 0, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 0, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                        //y axis
                        //GX_value = 0.0;
                        //sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 1, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 1, GX_value);
                        //sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 1, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                        //z axis
                        GX_value = - 2 - (input_zz)*((GRAVITY + input_hdd)/(com_z*com_z)) * sgn_mpc_dt_ * sgn_mpc_dt_;
                        sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index + 2, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index + 2, GX_value);
                        sgn_Hessian_triplets_tick_.emplace_back(GX_col_index + 2, GX_to_Hessian_row_offset + GX_row_index, GX_value);
                    }
                }
            }

            if (i < sgn_N_plan_mpc_ - 2)
            {
                GX_row_index = sgn_state_num_ * (i + 2) + j;
                GX_col_index = sgn_state_num_ *  i      + j;

                GX_value = 1.0;

                sgn_GX_triplets_.emplace_back(GX_row_index, GX_col_index, GX_value);
                //sgn_Hessian GX  part
                sgn_Hessian_triplets_tick_.emplace_back(GX_to_Hessian_row_offset + GX_row_index, GX_col_index, GX_value);
                //sgn_Hessian GX' part
                sgn_Hessian_triplets_tick_.emplace_back(GX_col_index, GX_to_Hessian_row_offset + GX_row_index, GX_value);
            }
        }

        //GU tick
        for (int j = 0; j < sgn_input_num_w_; j++)
        {
            GU_row_index = sgn_state_num_ * i;
            GU_col_index = sgn_input_num_ * i + j;

            //FIP
            if(FIP_or_VHIP == 0)
            {
                GU_value = sgn_w0_ * sgn_w0_ * sgn_mpc_dt_ * sgn_mpc_dt_ * sgn_s_x_ref_(i, j);
                sgn_GU_triplets_.emplace_back(GU_row_index + 0, GU_col_index, GU_value);
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 0, GU_to_Hessian_col_offset + GU_col_index, GU_value); //sgn_Hessian GU  part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 0, GU_value); //sgn_Hessian GU' part

                GU_value = sgn_w0_ * sgn_w0_ * sgn_mpc_dt_ * sgn_mpc_dt_ * sgn_s_y_ref_(i, j);
                sgn_GU_triplets_.emplace_back(GU_row_index + 1, GU_col_index, GU_value);
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 1, GU_to_Hessian_col_offset + GU_col_index, GU_value); //sgn_Hessian GU  part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 1, GU_value); //sgn_Hessian GU' part

                GU_value = sgn_w0_ * sgn_w0_ * sgn_mpc_dt_ * sgn_mpc_dt_ * sgn_s_z_ref_(i, j);
                sgn_GU_triplets_.emplace_back(GU_row_index + 2, GU_col_index, GU_value);
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 2, GU_to_Hessian_col_offset + GU_col_index, GU_value); //sgn_Hessian GU  part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 2, GU_value); //sgn_Hessian GU' part
            }

            //VHIP
            if(FIP_or_VHIP == 1)
            {
                GU_value = sgn_s_x_ref_(i,j)*((GRAVITY + input_hdd)/com_z)*sgn_mpc_dt_ * sgn_mpc_dt_;
                sgn_GU_triplets_.emplace_back(GU_row_index + 0, GU_col_index, GU_value);
                //sgn_Hessian GU  part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 0, GU_to_Hessian_col_offset + GU_col_index, GU_value);
                //sgn_Hessian GU' part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 0, GU_value);

                GU_value = sgn_s_y_ref_(i,j)*((GRAVITY + input_hdd)/com_z)*sgn_mpc_dt_ * sgn_mpc_dt_;
                sgn_GU_triplets_.emplace_back(GU_row_index + 1, GU_col_index, GU_value);
                //sgn_Hessian GU  part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 1, GU_to_Hessian_col_offset + GU_col_index, GU_value);
                //sgn_Hessian GU' part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 1, GU_value);

                GU_value = sgn_s_z_ref_(i,j)*((GRAVITY + input_hdd)/com_z)*sgn_mpc_dt_ * sgn_mpc_dt_;
                sgn_GU_triplets_.emplace_back(GU_row_index + 2, GU_col_index, GU_value);
                //sgn_Hessian GU  part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 2, GU_to_Hessian_col_offset + GU_col_index, GU_value);
                //sgn_Hessian GU' part
                sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 2, GU_value);
            }
        }

        //VHIP hdd part
        if(sgn_input_num_h_ != 0 && FIP_or_VHIP == 1)
        {
            GU_row_index = sgn_state_num_ *  i;
            GU_col_index = sgn_input_num_ * (i + 1) - 1;

            GU_value = - (com_x - input_zx)*sgn_mpc_dt_ * sgn_mpc_dt_ / com_z;
            sgn_GU_triplets_.emplace_back(GU_row_index + 0, GU_col_index, GU_value);
            //sgn_Hessian GU  part
            sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 0, GU_to_Hessian_col_offset + GU_col_index, GU_value);
            //sgn_Hessian GU' part
            sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 0, GU_value);

            GU_value = - (com_y - input_zy)*sgn_mpc_dt_ * sgn_mpc_dt_ / com_z;
            sgn_GU_triplets_.emplace_back(GU_row_index + 1, GU_col_index, GU_value);
            //sgn_Hessian GU  part
            sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 1, GU_to_Hessian_col_offset + GU_col_index, GU_value); //sgn_Hessian GU  part
            //sgn_Hessian GU' part
            sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 1, GU_value); //sgn_Hessian GU' part

            GU_value = - (com_z - input_zz)*sgn_mpc_dt_ * sgn_mpc_dt_ / com_z;
            sgn_GU_triplets_.emplace_back(GU_row_index + 2, GU_col_index, GU_value);
            //sgn_Hessian GU  part
            sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_row_offset + GU_row_index + 2, GU_to_Hessian_col_offset + GU_col_index, GU_value); //sgn_Hessian GU  part
            //sgn_Hessian GU' part
            sgn_Hessian_triplets_tick_.emplace_back(GU_to_Hessian_col_offset + GU_col_index, GU_to_Hessian_row_offset + GU_row_index + 2, GU_value); //sgn_Hessian GU' part
        }
    }
    
    Eigen::VectorXd sgn_com_ref_vel_vec; sgn_com_ref_vel_vec.setZero(sgn_state_num_ * (sgn_N_plan_mpc_ - 1));
    sgn_com_ref_vel_vec = sgn_D_total_ * sgn_com_ref_vec_;

    Eigen::VectorXd sgn_com_v0_ref, sgn_com_v0;
    sgn_com_v0_ref = sgn_com_ref_vec_.segment(0, sgn_state_num_) - sgn_com_ref0_vec_;
    sgn_com_v0     = sgn_X_bar_.segment(0, sgn_state_num_) - sgn_planner_state_mpc_;

    sgn_JX_ = sgn_H_ * (sgn_X_bar_ - sgn_com_ref_vec_);

    //soft constraint
    Eigen::VectorXd JU_soft_constraint;
    JU_soft_constraint.setZero(sgn_input_num_ * sgn_N_plan_mpc_);
    Eigen::VectorXd Gamma = sgn_U_bar_;
    double epsilon = 1e-1;
    for(int i = 0; i < Gamma.size(); i++)
    {
        double g = Gamma(i);

        if(g >= epsilon)
        {
            JU_soft_constraint(i) = 0.0;
        }
        else if(g <= -epsilon)
        {
            JU_soft_constraint(i) = - (1.0/(2.0*epsilon))*g*g
                                    + g
                                    - (epsilon/2.0);
                                    
            double val = sgn_K5_ * (-(1.0/epsilon)*g + 1.0);

            sgn_Hessian_triplets_tick_.emplace_back(sgn_state_num_ * sgn_N_plan_mpc_ + i, sgn_state_num_ * sgn_N_plan_mpc_ + i, val); //sgn_Hessian soft constraint part
        }
        else
        {
            JU_soft_constraint(i) = 2.0*g;

            sgn_Hessian_triplets_tick_.emplace_back(sgn_state_num_ * sgn_N_plan_mpc_ + i, sgn_state_num_ * sgn_N_plan_mpc_ + i, 2.0*sgn_K5_); //sgn_Hessian soft constraint part
        }
    }
    
    Eigen::MatrixXd sgn_Hessian_contact_mask_ref_SUw_calc;
    sgn_Hessian_contact_mask_ref_SUw_calc = sgn_contact_mask_ref_ * sgn_SUw_;

    Eigen::MatrixXd sgn_Hessian_contact_mask_ref_SUw_calcTcalc;
    sgn_Hessian_contact_mask_ref_SUw_calcTcalc = sgn_Hessian_contact_mask_ref_SUw_calc.transpose() * sgn_Hessian_contact_mask_ref_SUw_calc;

    sgn_JU_ = 2*sgn_Rw_ * sgn_U_bar_
            + sgn_K4_ * sgn_Hessian_contact_mask_ref_SUw_calcTcalc * sgn_U_bar_
            - sgn_K4_ * sgn_Hessian_contact_mask_ref_SUw_calc.transpose() * Eigen::VectorXd::Ones(sgn_N_plan_mpc_)
            + sgn_K5_ * JU_soft_constraint;
    
    //Hessian
    for(int i = sgn_state_num_ * sgn_N_plan_mpc_; i < (sgn_state_num_ + sgn_input_num_) * sgn_N_plan_mpc_; i++)
    {
        for(int j = sgn_state_num_ * sgn_N_plan_mpc_; j < (sgn_state_num_ + sgn_input_num_) * sgn_N_plan_mpc_; j++)
        {
            if(sgn_Hessian_contact_mask_ref_SUw_calcTcalc(i - sgn_state_num_ * sgn_N_plan_mpc_, j - sgn_state_num_ * sgn_N_plan_mpc_) != 0.0)
            {
                sgn_Hessian_triplets_tick_.emplace_back(i, j, sgn_K4_ * sgn_Hessian_contact_mask_ref_SUw_calcTcalc.coeff(i - sgn_state_num_ * sgn_N_plan_mpc_, j - sgn_state_num_ * sgn_N_plan_mpc_)); //sgn_Hessian JUU part contact mask
            }
        }
    }
    
    sgn_GX_.setFromTriplets(sgn_GX_triplets_.begin(), sgn_GX_triplets_.end());
    sgn_GU_.setFromTriplets(sgn_GU_triplets_.begin(), sgn_GU_triplets_.end());

    sgn_Hessian_triplets_total_.insert(sgn_Hessian_triplets_total_.end(), sgn_Hessian_triplets_init_.begin(), sgn_Hessian_triplets_init_.end());
    sgn_Hessian_triplets_total_.insert(sgn_Hessian_triplets_total_.end(), sgn_Hessian_triplets_tick_.begin(), sgn_Hessian_triplets_tick_.end());

    sgn_Hessian_.setFromTriplets(sgn_Hessian_triplets_total_.begin(), sgn_Hessian_triplets_total_.end());
    
    //gradient
    Eigen::VectorXd sgn_G0;
    sgn_G0.setZero(sgn_state_num_ * sgn_N_plan_mpc_);
    sgn_gradient_ << sgn_JX_, sgn_JU_, sgn_G0;

    if(sgn_matlab_tick_ == 4.4 * sgn_mpc_freq_) // 182
    {
        if(sgn_data_save_)
        {
            //e_sgn_data[sgn_file_num_ - 2] << Eigen::MatrixXd(sgn_Hessian_) << endl;
            //e_sgn_data[sgn_file_num_ - 3] << Eigen::MatrixXd(sgn_HD_) << endl;
            //e_sgn_data[sgn_file_num_ - 4] << Eigen::MatrixXd(sgn_Hz_) << endl;
            //e_sgn_data[sgn_file_num_ - 5] << Eigen::MatrixXd(sgn_H_)  << endl;
            //e_sgn_data[sgn_file_num_ - 6] << Eigen::MatrixXd(sgn_D_total_) << endl;
            //e_sgn_data[sgn_file_num_ - 7] << Eigen::MatrixXd(sgn_GX_) << endl;
            //e_sgn_data[sgn_file_num_ - 8] << Eigen::MatrixXd(sgn_GU_) << endl;
            //e_sgn_data[sgn_file_num_ - 9] << Eigen::MatrixXd(sgn_gradient_) << endl;
        }
    }

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(sgn_Hessian_);
    Eigen::VectorXd dZ = solver.solve(-sgn_gradient_);
    Eigen::VectorXd dX = dZ.segment(0, sgn_state_num_ * sgn_N_plan_mpc_);
    Eigen::VectorXd dU = dZ.segment(sgn_state_num_ * sgn_N_plan_mpc_, sgn_input_num_ * sgn_N_plan_mpc_);

    Eigen::VectorXd sgn_U_opt = sgn_U_bar_ + dU;

    Eigen::VectorXd sgn_u_apply = sgn_U_opt.segment(0, sgn_input_num_);
    
    s_ref_calc.row(0) = sgn_s_x_ref_.row(0);
    s_ref_calc.row(1) = sgn_s_y_ref_.row(0);
    s_ref_calc.row(2) = sgn_s_z_ref_.row(0);

    
    if(FIP_or_VHIP == 0)
    {
        sgn_planner_state_p1_mpc_ = 2*sgn_planner_state_mpc_
                                  - sgn_planner_state_m1_mpc_
                                  + sgn_w0_*sgn_w0_*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_planner_state_mpc_ - s_ref_calc*sgn_u_apply.segment(0, sgn_input_num_w_))
                                  + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
    }
    else
    {
        sgn_planner_state_p1_mpc_ = 2*sgn_planner_state_mpc_
                                  - sgn_planner_state_m1_mpc_
                                  + ((GRAVITY + sgn_u_apply(sgn_input_num_-1))/sgn_planner_state_mpc_(sgn_state_num_-1))*sgn_mpc_dt_*sgn_mpc_dt_*(sgn_planner_state_mpc_ - s_ref_calc*sgn_u_apply.segment(0, sgn_input_num_w_))
                                  + sgn_gravity_vec_*sgn_mpc_dt_*sgn_mpc_dt_;
    }
    
    sgn_planner_state_m1_mpc_ = sgn_planner_state_mpc_;
    sgn_planner_state_mpc_    = sgn_planner_state_p1_mpc_;

    sgn_U_bar_.setZero();
    sgn_U_bar_.segment(0, sgn_input_num_*(sgn_N_plan_mpc_ - 1)) = sgn_U_opt.segment(sgn_input_num_, sgn_input_num_*(sgn_N_plan_mpc_ - 1));

    Eigen::Vector3d calculated_grf;
    calculated_grf = s_ref_calc * sgn_u_apply.segment(0, sgn_input_num_w_);

    e_sgn_data[sgn_file_num_ - 1] << sgn_com_ref_mat_(0,0)        << "," << sgn_com_ref_mat_(0,1)        << "," << sgn_com_ref_mat_(0,2)        << ","
                                  << sgn_planner_state_mpc_(0)    << "," << sgn_planner_state_mpc_(1)    << "," << sgn_planner_state_mpc_(2)    << ","
                                  << sgn_planner_state_m1_mpc_(0) << "," << sgn_planner_state_m1_mpc_(1) << "," << sgn_planner_state_m1_mpc_(2) << ","
                                  << sgn_planner_state_p1_mpc_(0) << "," << sgn_planner_state_p1_mpc_(1) << "," << sgn_planner_state_p1_mpc_(2) << ","
                                  << calculated_grf(0)            << "," << calculated_grf(1)            << "," << calculated_grf(2)            << ","
                                  << endl;
    
    frameChange(&sgn_planner_state_mpc_);
    frameChange(&sgn_planner_state_m1_mpc_);
    
    return sgn_planner_state_mpc_;
}

void SgnPlanner::frameChange(Eigen::VectorXd* state_before_frame_change)
{
    if((sgn_walking_tick_mpc_ - (sgn_t_start_mpc_ + sgn_t_total_mpc_) >= - 1.0/(sgn_mpc_freq_*sgn_dt_)) && (sgn_current_step_num_mpc_ < sgn_total_step_num_mpc_ - 1))
    {
        Eigen::Vector3d var_after_step_change, var_before_step_change, frame_pos_diff;
        Eigen::Matrix3d frame_rot_diff;

        frame_rot_diff = DyrosMath::rotateWithZ(-sgn_foot_step_support_frame_wo_offset_mpc_(sgn_current_step_num_mpc_, 5));
        frame_pos_diff = sgn_foot_step_support_frame_wo_offset_mpc_.block(sgn_current_step_num_mpc_, 0, 1, 3).transpose();

        //com pos step change
        var_before_step_change = *state_before_frame_change;
        var_after_step_change = frame_rot_diff*(var_before_step_change - frame_pos_diff);

        *state_before_frame_change = var_after_step_change;
    }
}