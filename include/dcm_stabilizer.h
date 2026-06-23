#pragma once
#include <Eigen/Dense>
#include "math_type_define.h"
#include "wholebody_functions.h"
#include <fstream>

extern std::vector<std::ofstream> e_stabilizer_data_txt_;

class dcmStabilizer
{
    public:
        static const int num_ofstream = 10;

        dcmStabilizer();

        //functions
        void ofstream_open();
        void ofstream_close();
        void IS_FIPM_3D_DCM_Stabililzer_MPC(double mpc_freq, double mpc_preview_window);

        //variables
        //walking parameters
        int walking_tick_mpc_ = 0;
        int com_start_tick_mpc_ = 0;
        double t_total_const_mpc_ = 0;
        double t_total_mpc_ = 0;
        double t_start_mpc_ = 0;
        double w_ = 0;
        double b_ = 0;
        int total_step_num_mpc_ = 0;
        int current_step_num_mpc_ = 0;

        double zmp_x_max;
        double zmp_x_min;
        double zmp_y_max;
        double zmp_y_min;

        Eigen::VectorXd zmp_max_x_mpc_;
        Eigen::VectorXd zmp_min_x_mpc_;

        Eigen::VectorXd zmp_max_y_mpc_;
        Eigen::VectorXd zmp_min_y_mpc_;

        Eigen::MatrixXd ref_vrp_mpc_;

        Eigen::MatrixXd ref_zmp_wo_offset_mpc_;

        Eigen::MatrixXd foot_step_support_frame_mpc_;

        Eigen::MatrixXd foot_step_support_frame_offset_mpc_;

        Eigen::Vector3d com_measured_mpc_;
        Eigen::Vector3d com_dot_measured_mpc_;

        int step_time_adj_candidate_num_;

        //mpc parameters
        double hz_ = 2000.0;

        double thread3_hz_ = 0.0;

        bool MPC_first_loop_ = 0;

        double print_sec_;

        bool print_sec_bool_;
        
        Eigen::MatrixXd Ccp_mpc_;
        Eigen::MatrixXd Ccv_mpc_;
        Eigen::MatrixXd Cvp_mpc_;
        Eigen::MatrixXd Cdp_mpc_;

        Eigen::MatrixXd A_mpc_;
        Eigen::MatrixXd B_mpc_;

        Eigen::MatrixXd Pdps_stab_mpc_;
        Eigen::MatrixXd Pcps_stab_mpc_;
        Eigen::MatrixXd Pcvs_stab_mpc_;
        Eigen::MatrixXd Pvps_stab_mpc_;

        Eigen::MatrixXd Pdpu_stab_mpc_;
        Eigen::MatrixXd Pcpu_stab_mpc_;
        Eigen::MatrixXd Pcvu_stab_mpc_;
        Eigen::MatrixXd Pvpu_stab_mpc_;

        Eigen::MatrixXd P_IS_step_mpc_;
        Eigen::MatrixXd b_IS_plan_mpc_;
        Eigen::MatrixXd b_IS_step_mpc_;
        Eigen::MatrixXd p_IS_step_mpc_;

        CQuadraticProgram QP_MPC_Stabilizer_;

        Eigen::MatrixXd Qmat_stab_mpc_Q_;
        Eigen::MatrixXd Qmat_stab_mpc_R_;

        Eigen::MatrixXd Qcalc_stab_mpc_;

        Eigen::MatrixXd gcalc_stab_mpc_;
        Eigen::MatrixXd gxpcalc_stab_mpc_;
        Eigen::MatrixXd gypcalc_stab_mpc_;
        Eigen::MatrixXd gzpcalc_stab_mpc_;

        Eigen::VectorXd MPC_Stabilizer_u_mpc_;
        Eigen::VectorXd MPC_Stabilizer_SQP_du_mpc_;

        Eigen::VectorXd MPC_Stabilizer_alpha_mpc_;

        Eigen::VectorXd MPC_Stabilizer_aux_mpc_;
        Eigen::VectorXd MPC_Stabilizer_aux_mpc_x_;
        Eigen::VectorXd MPC_Stabilizer_aux_mpc_y_;

        Eigen::VectorXd MPC_Stabilizer_delf_mpc_;
        Eigen::VectorXd MPC_Stabilizer_delf_mpc_x_;
        Eigen::VectorXd MPC_Stabilizer_delf_mpc_y_;

        Eigen::MatrixXd ssx_stab_mpc_;
        Eigen::MatrixXd ssy_stab_mpc_;
        Eigen::MatrixXd ssz_stab_mpc_;

        Eigen::MatrixXd SUp_stab_mpc_;
        Eigen::MatrixXd SUpx_stab_mpc_;
        Eigen::MatrixXd SUpxp_stab_mpc_;
        Eigen::MatrixXd SUpy_stab_mpc_;
        Eigen::MatrixXd SUpyp_stab_mpc_;
        Eigen::MatrixXd SUpz_stab_mpc_;
        Eigen::MatrixXd SUpzp_stab_mpc_;

        Eigen::MatrixXd Sf1_stab_mpc_;

        Eigen::VectorXd MPC_Stabilizer_u_mpc_sep_;

        Eigen::VectorXd MPC_Stabilizer_state_mpc_;

        Eigen::MatrixXd Pv_dot_ref_mpc_;
        
        Eigen::MatrixXd IS_FIPM_SQP_x_phi1_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_phi2_N_stab_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_x_pi1_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_pi2_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_pi3_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_pi4_N_stab_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_x_ri1_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_ri2_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_x_ri3_N_stab_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_y_phi1_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_phi2_N_stab_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_y_pi1_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_pi2_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_pi3_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_pi4_N_stab_mpc_;

        Eigen::MatrixXd IS_FIPM_SQP_y_ri1_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_ri2_N_stab_mpc_;
        Eigen::MatrixXd IS_FIPM_SQP_y_ri3_N_stab_mpc_;

        Eigen::MatrixXd const_SQP_phi_mpc_;
        Eigen::MatrixXd const_SQP_phi_mpc_calc_;
        Eigen::MatrixXd const_SQP_pi_mpc_;
        Eigen::MatrixXd const_SQP_ri_mpc_;
        Eigen::MatrixXd const_SQP_hi_mpc_;

        Eigen::MatrixXd SQP_deldel_Qcalc_stab_mpc_;
        Eigen::MatrixXd SQP_del_gcalc_stab_mpc_;

        Eigen::MatrixXd Planner_State_Prev_mpc_;

        Eigen::MatrixXd const_A_mpc_;
        Eigen::MatrixXd const_ub_mpc_;
        Eigen::MatrixXd const_lb_mpc_;

        const double equality_condition_eps_ = 1e-8;
};