#pragma once
#include <Eigen/Dense>
#include "math_type_define.h"
#include "wholebody_functions.h"
#include <fstream>

extern std::vector<std::ofstream> e_planner_data_txt_;

class fipmPlanner
{
    public:
        static const int num_ofstream = 10;

        fipmPlanner();

        //functions
        void ofstream_open();
        void ofstream_close();
        void IS_FIPM_CoM_Planner_MPC(double mpc_freq, double preview_window);
        void MPC_State_Step_Change();

        //variables
        //walking parameters
        int walking_tick_mpc_ = 0;
        int com_start_tick_mpc_ = 0;
        double t_total_const_mpc_ = 0;
        double t_total_mpc_ = 0;
        double t_start_mpc_ = 0;
        double w_ = 0;
        double b_ = 0;
        double hz_ = 2000.0;
        int total_step_num_mpc_ = 0;
        int current_step_num_mpc_ = 0;
        double zmp_x_max = 0.13;
        double zmp_x_min = 0.09;
        double zmp_y_max = 0.085;
        double zmp_y_min = 0.085;

        Eigen::VectorXd zmp_max_x_mpc_;
        Eigen::VectorXd zmp_min_x_mpc_;
        Eigen::VectorXd zmp_max_y_mpc_;
        Eigen::VectorXd zmp_min_y_mpc_;

        Eigen::MatrixXd ref_vrp_mpc_;

        Eigen::MatrixXd ref_zmp_wo_offset_mpc_;

        Eigen::MatrixXd foot_step_support_frame_mpc_;

        Eigen::MatrixXd foot_step_support_frame_offset_mpc_;

        //MPC parameters
        bool MPC_first_loop_ = 0;

        double thread3_hz_;

        double print_sec_;

        bool print_sec_bool_;

        Eigen::MatrixXd A_mpc_;
        Eigen::MatrixXd B_mpc_;
        Eigen::MatrixXd Ccp_mpc_;
        Eigen::MatrixXd Ccv_mpc_;
        Eigen::MatrixXd Cvp_mpc_;

        CQuadraticProgram QP_MPC_Planner_;

        Eigen::MatrixXd Pcps_plan_mpc_;
        Eigen::MatrixXd Pcvs_plan_mpc_;
        Eigen::MatrixXd Pvps_plan_mpc_;

        Eigen::MatrixXd Pcpu_plan_mpc_;
        Eigen::MatrixXd Pcvu_plan_mpc_;
        Eigen::MatrixXd Pvpu_plan_mpc_;

        Eigen::MatrixXd P_IS_step_mpc_;
        Eigen::MatrixXd b_IS_plan_mpc_;
        Eigen::MatrixXd b_IS_step_mpc_;
        Eigen::MatrixXd p_IS_step_mpc_;

        Eigen::MatrixXd SUx_plan_mpc_;
        Eigen::MatrixXd SUy_plan_mpc_;
        Eigen::MatrixXd SUz_plan_mpc_;

        Eigen::MatrixXd ssx_plan_mpc_;
        Eigen::MatrixXd ssy_plan_mpc_;
        Eigen::MatrixXd ssz_plan_mpc_;

        Eigen::MatrixXd Qmat_plan_mpc_;

        Eigen::MatrixXd Qxcalc_plan_mpc_;
        Eigen::MatrixXd Qycalc_plan_mpc_;
        Eigen::MatrixXd Qzcalc_plan_mpc_;

        Eigen::MatrixXd Qcalc_plan_mpc_;

        Eigen::MatrixXd gxcalc_plan_mpc_;
        Eigen::MatrixXd gycalc_plan_mpc_;
        Eigen::MatrixXd gzcalc_plan_mpc_;

        Eigen::MatrixXd gcalc_plan_mpc_;

        Eigen::VectorXd MPC_Planner_u_mpc_;
        Eigen::VectorXd MPC_Planner_SQP_du_mpc_;

        Eigen::VectorXd MPC_Planner_u_mpc_sep_;

        Eigen::VectorXd MPC_Planner_state_mpc_;

        Eigen::MatrixXd Planner_State_Prev_mpc_;

        Eigen::MatrixXd Pv_dot_ref_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_x_phi_N_plan_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_x_pi1_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_pi2_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_pi3_N_plan_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_x_ri1_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_ri2_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_ri3_N_plan_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_y_phi_N_plan_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_y_pi1_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_pi2_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_pi3_N_plan_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_y_ri1_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_ri2_N_plan_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_ri3_N_plan_mpc_;

        Eigen::MatrixXd const_SQP_phi_mpc_;
        Eigen::MatrixXd const_SQP_phi_mpc_calc_;
        Eigen::MatrixXd const_SQP_pi_mpc_;
        Eigen::MatrixXd const_SQP_ri_mpc_;
        Eigen::MatrixXd const_SQP_hi_mpc_;

        Eigen::MatrixXd SQP_deldel_Qcalc_plan_mpc_;
        Eigen::MatrixXd SQP_del_gcalc_plan_mpc_;

        Eigen::MatrixXd const_A_mpc_;
        Eigen::MatrixXd const_ub_mpc_;
        Eigen::MatrixXd const_lb_mpc_;

        Eigen::VectorXd MPC_Stabilizer_delf_mpc_x_;
        Eigen::VectorXd MPC_Stabilizer_delf_mpc_y_;
        
        const double equality_condition_eps_ = 1e-8;
};