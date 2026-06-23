#include "fipm_planner.h"
#include <iostream>

using namespace std;

std::vector<std::ofstream> e_planner_data_txt_;

fipmPlanner::fipmPlanner()
{
    cout << "fipmPlanner constructor called" << endl;
}

void fipmPlanner::ofstream_open()
{
    if (!e_planner_data_txt_.empty() && e_planner_data_txt_[0].is_open()) {
        return;
    }

    e_planner_data_txt_.resize(num_ofstream);

    for (int i = 0; i < num_ofstream; ++i) {
        if (!e_planner_data_txt_[i].is_open()) {
            e_planner_data_txt_[i].open("/home/econom2-20/data/fipm_planner/e_data_" + std::to_string(i) + ".txt");
        }
    }
}

void fipmPlanner::ofstream_close()
{
    for (auto &file_stream : e_planner_data_txt_) {
        if (file_stream.is_open()) {
            file_stream.close();
        }
    }
}

void fipmPlanner::IS_FIPM_CoM_Planner_MPC(double mpc_freq, double preview_window)
{
    double wpvx, wpvy, wpvz;
    double wdvx, wdvy, wdvz;

    wpvx = 1e+2; wpvy = 1e+6; wpvz = 1e+2;
    wdvx = 1e+0; wdvy = 1e+0; wdvz = 1e+0;

    int mpc_tick = walking_tick_mpc_ - com_start_tick_mpc_;
    double mpc_synchro_hz = hz_/mpc_freq;
    const int N_plan_mpc = preview_window*mpc_freq;
    const int N_step = t_total_const_mpc_/mpc_synchro_hz;
    const int N_state = 3; //com position, com velocity, vrp position
    double dt_plan_mpc = 1/mpc_freq;

    double lambda_is_calc = exp(-w_*dt_plan_mpc);

    int input_num = 3 * N_plan_mpc;
    //              xyz vrp                    
    int const_num = 6 * N_plan_mpc  + 3;
    //              xyz vrp min,max   IS
    if(MPC_first_loop_ == 0)
    {
        cout << "Initialization of IS FIPM Planner MPC." << endl;
        b_ = 1/w_;
        ofstream_open();
        A_mpc_.resize(N_state,N_state);
        Eigen::MatrixXd A_mpc_cont; A_mpc_cont.resize(N_state, N_state); A_mpc_cont.setZero();
        A_mpc_cont <<     0, 1,      0,
                      w_*w_, 0, -w_*w_,
                          0, 0,      0;
        
        A_mpc_ = MatrixXd::Identity(N_state,N_state) 
               + A_mpc_cont*dt_plan_mpc 
               + A_mpc_cont*A_mpc_cont*dt_plan_mpc*dt_plan_mpc/(1*2) 
               + A_mpc_cont*A_mpc_cont*A_mpc_cont*dt_plan_mpc*dt_plan_mpc*dt_plan_mpc/(1*2*3)
               + A_mpc_cont*A_mpc_cont*A_mpc_cont*A_mpc_cont*dt_plan_mpc*dt_plan_mpc*dt_plan_mpc*dt_plan_mpc/(1*2*3*4);

        B_mpc_.resize(N_state,1);
        Eigen::MatrixXd B_mpc_cont; B_mpc_cont.resize(N_state,1); B_mpc_cont.setZero();
        B_mpc_cont << 0,
                      0,
                      1;
        
        B_mpc_ = B_mpc_cont*dt_plan_mpc 
               + A_mpc_cont*B_mpc_cont*dt_plan_mpc*dt_plan_mpc/(1*2) 
               + A_mpc_cont*A_mpc_cont*B_mpc_cont*dt_plan_mpc*dt_plan_mpc*dt_plan_mpc/(1*2*3)
               + A_mpc_cont*A_mpc_cont*A_mpc_cont*B_mpc_cont*dt_plan_mpc*dt_plan_mpc*dt_plan_mpc*dt_plan_mpc/(1*2*3*4);

        Ccp_mpc_.resize(1,N_state);
        Ccv_mpc_.resize(1,N_state);
        Cvp_mpc_.resize(1,N_state);

        Ccp_mpc_ << 1, 0, 0;
        Ccv_mpc_ << 0, 1, 0;
        Cvp_mpc_ << 0, 0, 1;

        Pcps_plan_mpc_.resize(N_plan_mpc,N_state);
        Pcvs_plan_mpc_.resize(N_plan_mpc,N_state);
        Pvps_plan_mpc_.resize(N_plan_mpc,N_state);

        Eigen::MatrixXd Ps_calc;
        Ps_calc.resize(N_state,N_state);
        Ps_calc = A_mpc_;
        
        Pcpu_plan_mpc_.setZero(N_plan_mpc, N_plan_mpc);
        Pcvu_plan_mpc_.setZero(N_plan_mpc, N_plan_mpc);
        Pvpu_plan_mpc_.setZero(N_plan_mpc, N_plan_mpc);
        
        P_IS_step_mpc_.setZero(N_step, N_step);

        b_IS_plan_mpc_.setZero(N_plan_mpc,1);
        b_IS_step_mpc_.setZero(N_step,1);

        p_IS_step_mpc_.setOnes(N_step,1);
        
        Eigen::MatrixXd Pu_calc, Pu_step_calc;
        Pu_calc.setZero(N_state,N_plan_mpc);
        Pu_step_calc.setZero(N_state,N_step);
        
        for(int i = 0; i < N_plan_mpc; i++)
        {
            Pcps_plan_mpc_.row(i) = Ccp_mpc_*Ps_calc;
            Pcvs_plan_mpc_.row(i) = Ccv_mpc_*Ps_calc;
            Pvps_plan_mpc_.row(i) = Cvp_mpc_*Ps_calc;
            Ps_calc = Ps_calc*A_mpc_;

            Pu_calc.col(i) = B_mpc_;
            Pcpu_plan_mpc_.row(i) = Ccp_mpc_*Pu_calc;
            Pcvu_plan_mpc_.row(i) = Ccv_mpc_*Pu_calc;
            Pvpu_plan_mpc_.row(i) = Cvp_mpc_*Pu_calc;
            b_IS_plan_mpc_(i,0) = pow(lambda_is_calc, i);
            
            if(i < N_step)
            {
                Pu_step_calc.col(i)   = B_mpc_;
                P_IS_step_mpc_.row(i) = Cvp_mpc_*Pu_step_calc;    
                b_IS_step_mpc_(i,0)   = pow(lambda_is_calc, i);
            }

            Pu_calc = A_mpc_*Pu_calc;
            Pu_step_calc = A_mpc_*Pu_step_calc;
        }

        SUx_plan_mpc_.setZero(N_plan_mpc, input_num); SUx_plan_mpc_ << MatrixXd::Identity(N_plan_mpc, N_plan_mpc), MatrixXd::Zero(N_plan_mpc, N_plan_mpc), MatrixXd::Zero(N_plan_mpc, N_plan_mpc);
        SUy_plan_mpc_.setZero(N_plan_mpc, input_num); SUy_plan_mpc_ << MatrixXd::Zero(N_plan_mpc, N_plan_mpc), MatrixXd::Identity(N_plan_mpc, N_plan_mpc), MatrixXd::Zero(N_plan_mpc, N_plan_mpc);
        SUz_plan_mpc_.setZero(N_plan_mpc, input_num); SUz_plan_mpc_ << MatrixXd::Zero(N_plan_mpc, N_plan_mpc), MatrixXd::Zero(N_plan_mpc, N_plan_mpc), MatrixXd::Identity(N_plan_mpc, N_plan_mpc);

        ssx_plan_mpc_.setZero(N_state, 3*N_state); ssx_plan_mpc_ << MatrixXd::Identity(N_state, N_state), MatrixXd::Zero(N_state, N_state), MatrixXd::Zero(N_state, N_state);
        ssy_plan_mpc_.setZero(N_state, 3*N_state); ssy_plan_mpc_ << MatrixXd::Zero(N_state, N_state), MatrixXd::Identity(N_state, N_state), MatrixXd::Zero(N_state, N_state);
        ssz_plan_mpc_.setZero(N_state, 3*N_state); ssz_plan_mpc_ << MatrixXd::Zero(N_state, N_state), MatrixXd::Zero(N_state, N_state), MatrixXd::Identity(N_state, N_state);

        Qmat_plan_mpc_.resize(N_plan_mpc, N_plan_mpc);
        Qmat_plan_mpc_.setIdentity();

        Qxcalc_plan_mpc_.setZero(N_plan_mpc, N_plan_mpc);
        Qxcalc_plan_mpc_ = Pvpu_plan_mpc_.transpose()*wpvx*Qmat_plan_mpc_*Pvpu_plan_mpc_ + wdvx*Qmat_plan_mpc_;

        Qycalc_plan_mpc_.setZero(N_plan_mpc, N_plan_mpc);
        Qycalc_plan_mpc_ = Pvpu_plan_mpc_.transpose()*wpvy*Qmat_plan_mpc_*Pvpu_plan_mpc_ + wdvy*Qmat_plan_mpc_;

        Qzcalc_plan_mpc_.setZero(N_plan_mpc, N_plan_mpc);
        Qzcalc_plan_mpc_ = Pvpu_plan_mpc_.transpose()*wpvz*Qmat_plan_mpc_*Pvpu_plan_mpc_ + wdvz*Qmat_plan_mpc_;

        Qcalc_plan_mpc_.setZero(input_num, input_num);
        Qcalc_plan_mpc_ = SUx_plan_mpc_.transpose()*Qxcalc_plan_mpc_*SUx_plan_mpc_ 
                        + SUy_plan_mpc_.transpose()*Qycalc_plan_mpc_*SUy_plan_mpc_
                        + SUz_plan_mpc_.transpose()*Qzcalc_plan_mpc_*SUz_plan_mpc_;
        
        gxcalc_plan_mpc_.setZero(input_num, N_plan_mpc);
        gxcalc_plan_mpc_ = SUx_plan_mpc_.transpose()*Pvpu_plan_mpc_.transpose()*wpvx*Qmat_plan_mpc_;

        gycalc_plan_mpc_.setZero(input_num, N_plan_mpc);
        gycalc_plan_mpc_ = SUy_plan_mpc_.transpose()*Pvpu_plan_mpc_.transpose()*wpvy*Qmat_plan_mpc_;

        gzcalc_plan_mpc_.setZero(input_num, N_plan_mpc);
        gzcalc_plan_mpc_ = SUz_plan_mpc_.transpose()*Pvpu_plan_mpc_.transpose()*wpvz*Qmat_plan_mpc_;

        QP_MPC_Planner_.InitializeProblemSize(input_num, const_num);

        MPC_Planner_u_mpc_.setZero(input_num);
        MPC_Planner_SQP_du_mpc_.setZero(input_num);

        MPC_Planner_u_mpc_sep_.setZero(3);

        Planner_State_Prev_mpc_.setZero(9, N_plan_mpc);

        Pv_dot_ref_mpc_.setZero(N_step, 3);
        
        zmp_max_x_mpc_.setZero(N_plan_mpc);
        zmp_min_x_mpc_.setZero(N_plan_mpc);

        zmp_max_y_mpc_.setZero(N_plan_mpc);
        zmp_min_y_mpc_.setZero(N_plan_mpc);

        IS_FIPM_SQP_x_phi_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 3*N_plan_mpc);
        IS_FIPM_SQP_x_pi1_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 3*N_state);
        IS_FIPM_SQP_x_pi2_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 1);
        IS_FIPM_SQP_x_pi3_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 1);
        IS_FIPM_SQP_x_ri1_N_plan_mpc_.setZero(3*N_state*N_plan_mpc,    3*N_state);
        IS_FIPM_SQP_x_ri2_N_plan_mpc_.setZero(1*N_plan_mpc,            3*N_state);
        IS_FIPM_SQP_x_ri3_N_plan_mpc_.setZero(1*N_plan_mpc,            3*N_state);
        IS_FIPM_SQP_y_phi_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 3*N_plan_mpc);
        IS_FIPM_SQP_y_pi1_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 3*N_state);
        IS_FIPM_SQP_y_pi2_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 1);
        IS_FIPM_SQP_y_pi3_N_plan_mpc_.setZero(3*N_plan_mpc*N_plan_mpc, 1);
        IS_FIPM_SQP_y_ri1_N_plan_mpc_.setZero(3*N_state*N_plan_mpc,    3*N_state);
        IS_FIPM_SQP_y_ri2_N_plan_mpc_.setZero(1*N_plan_mpc,            3*N_state);
        IS_FIPM_SQP_y_ri3_N_plan_mpc_.setZero(1*N_plan_mpc,            3*N_state);

        int calc_index  = 0;
        int calc_index2 = 0;
        Eigen::MatrixXd Si_mpc; Si_mpc.setZero(1, N_plan_mpc);

        for(int i = 0; i < N_plan_mpc; i++)
        {
            Si_mpc.setZero(1, N_plan_mpc);
            Si_mpc(0, i) = 1;

            IS_FIPM_SQP_x_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc) = (Si_mpc*Pvpu_plan_mpc_*SUx_plan_mpc_).transpose()*(Si_mpc*Pcpu_plan_mpc_*SUz_plan_mpc_)
                                                                                           - (Si_mpc*Pvpu_plan_mpc_*SUz_plan_mpc_).transpose()*(Si_mpc*Pcpu_plan_mpc_*SUx_plan_mpc_);

            IS_FIPM_SQP_y_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc) = (Si_mpc*Pvpu_plan_mpc_*SUy_plan_mpc_).transpose()*(Si_mpc*Pcpu_plan_mpc_*SUz_plan_mpc_)
                                                                                           - (Si_mpc*Pvpu_plan_mpc_*SUz_plan_mpc_).transpose()*(Si_mpc*Pcpu_plan_mpc_*SUy_plan_mpc_);

            IS_FIPM_SQP_x_pi1_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_state) = (Si_mpc*Pcpu_plan_mpc_*SUz_plan_mpc_).transpose()*(Si_mpc*Pvps_plan_mpc_*ssx_plan_mpc_)
                                                                                        + (Si_mpc*Pvpu_plan_mpc_*SUx_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssz_plan_mpc_)
                                                                                        - (Si_mpc*Pcpu_plan_mpc_*SUx_plan_mpc_).transpose()*(Si_mpc*Pvps_plan_mpc_*ssz_plan_mpc_)
                                                                                        - (Si_mpc*Pvpu_plan_mpc_*SUz_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssx_plan_mpc_);
            
            IS_FIPM_SQP_x_pi2_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1) = (GRAVITY*b_*b_*(Si_mpc*Pcpu_plan_mpc_*SUx_plan_mpc_)).transpose();

            IS_FIPM_SQP_x_pi3_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1) = (Si_mpc*(Pcpu_plan_mpc_ - Pvpu_plan_mpc_)*SUz_plan_mpc_).transpose();

            IS_FIPM_SQP_y_pi1_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_state) = (Si_mpc*Pcpu_plan_mpc_*SUz_plan_mpc_).transpose()*(Si_mpc*Pvps_plan_mpc_*ssy_plan_mpc_)
                                                                                       + (Si_mpc*Pvpu_plan_mpc_*SUy_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssz_plan_mpc_)
                                                                                       - (Si_mpc*Pcpu_plan_mpc_*SUy_plan_mpc_).transpose()*(Si_mpc*Pvps_plan_mpc_*ssz_plan_mpc_)
                                                                                       - (Si_mpc*Pvpu_plan_mpc_*SUz_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssy_plan_mpc_);

            IS_FIPM_SQP_y_pi2_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1) = (GRAVITY*b_*b_*(Si_mpc*Pcpu_plan_mpc_*SUy_plan_mpc_)).transpose();

            IS_FIPM_SQP_y_pi3_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1) = (Si_mpc*(Pcpu_plan_mpc_ - Pvpu_plan_mpc_)*SUz_plan_mpc_).transpose();

            calc_index += 3*N_plan_mpc;

            IS_FIPM_SQP_x_ri1_N_plan_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state) = (Si_mpc*Pvps_plan_mpc_*ssx_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssz_plan_mpc_)
                                                                                     - (Si_mpc*Pvps_plan_mpc_*ssz_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssx_plan_mpc_);

            IS_FIPM_SQP_x_ri2_N_plan_mpc_.block(i, 0, 1, 3*N_state)                  = GRAVITY*b_*b_*(Si_mpc*Pcps_plan_mpc_*ssx_plan_mpc_);
            
            IS_FIPM_SQP_x_ri3_N_plan_mpc_.block(i, 0, 1, 3*N_state)                  = (Si_mpc*(Pcps_plan_mpc_ - Pvps_plan_mpc_)*ssz_plan_mpc_);

            IS_FIPM_SQP_y_ri1_N_plan_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state) = (Si_mpc*Pvps_plan_mpc_*ssy_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssz_plan_mpc_)
                                                                                     - (Si_mpc*Pvps_plan_mpc_*ssz_plan_mpc_).transpose()*(Si_mpc*Pcps_plan_mpc_*ssy_plan_mpc_);

            IS_FIPM_SQP_y_ri2_N_plan_mpc_.block(i, 0, 1, 3*N_state)                  = GRAVITY*b_*b_*(Si_mpc*Pcps_plan_mpc_*ssy_plan_mpc_);

            IS_FIPM_SQP_y_ri3_N_plan_mpc_.block(i, 0, 1, 3*N_state)                  = (Si_mpc*(Pcps_plan_mpc_ - Pvps_plan_mpc_)*ssz_plan_mpc_);

            calc_index2 += 3*N_state;
        }

        cout << "Initialization of IS FIPM Planner MPC is completed." << endl;

        MPC_first_loop_ = 1;

        print_sec_ = 2.0;

        t_total_mpc_ = t_total_const_mpc_;
    }

    print_sec_bool_ = ((((walking_tick_mpc_ - int(mpc_synchro_hz) + 1)/int(mpc_synchro_hz))%int(print_sec_*thread3_hz_))==0);

    int matlab_tick = ((walking_tick_mpc_ - int(mpc_synchro_hz) + 1)/int(mpc_synchro_hz));

    Eigen::VectorXd Pv_x_ref(N_plan_mpc);
    Eigen::VectorXd Pv_y_ref(N_plan_mpc);
    Eigen::VectorXd Pv_z_ref(N_plan_mpc);

    for(int i = 0; i < N_plan_mpc; i++)
    {
        Pv_x_ref(i) = ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i+1),0);
        Pv_y_ref(i) = ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i+1),1);
        Pv_z_ref(i) = ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i+1),2);

        if(i < N_step)
        {
            Pv_dot_ref_mpc_(i,0) = (ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i + N_plan_mpc + 1),0) - ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i + N_plan_mpc + 0),0))/dt_plan_mpc;
            Pv_dot_ref_mpc_(i,1) = (ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i + N_plan_mpc + 1),1) - ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i + N_plan_mpc + 0),1))/dt_plan_mpc;
            Pv_dot_ref_mpc_(i,2) = (ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i + N_plan_mpc + 1),2) - ref_vrp_mpc_(mpc_tick + mpc_synchro_hz*(i + N_plan_mpc + 0),2))/dt_plan_mpc;
        }

        zmp_max_x_mpc_(i) = ref_zmp_wo_offset_mpc_(mpc_tick + mpc_synchro_hz*(i+1),0) + zmp_x_max;     
        zmp_min_x_mpc_(i) = ref_zmp_wo_offset_mpc_(mpc_tick + mpc_synchro_hz*(i+1),0) - zmp_x_min;

        zmp_max_y_mpc_(i) = ref_zmp_wo_offset_mpc_(mpc_tick + mpc_synchro_hz*(i+1),1) + zmp_y_max;
        zmp_min_y_mpc_(i) = ref_zmp_wo_offset_mpc_(mpc_tick + mpc_synchro_hz*(i+1),1) - zmp_y_min;
    }

    gcalc_plan_mpc_ = gxcalc_plan_mpc_*(Pvps_plan_mpc_*ssx_plan_mpc_*MPC_Planner_state_mpc_ - Pv_x_ref)
                    + gycalc_plan_mpc_*(Pvps_plan_mpc_*ssy_plan_mpc_*MPC_Planner_state_mpc_ - Pv_y_ref)
                    + gzcalc_plan_mpc_*(Pvps_plan_mpc_*ssz_plan_mpc_*MPC_Planner_state_mpc_ - Pv_z_ref);

    SQP_deldel_Qcalc_plan_mpc_ = Qcalc_plan_mpc_;

    int sqp_iter = 1;
    for(int s = 0; s < sqp_iter; s++)
    {
        if(print_sec_bool_)
            cout << "Planner SQP Iteration: " << s << endl;

        SQP_del_gcalc_plan_mpc_    = Qcalc_plan_mpc_*MPC_Planner_u_mpc_ + gcalc_plan_mpc_;

        QP_MPC_Planner_.EnableEqualityCondition(equality_condition_eps_);
        QP_MPC_Planner_.UpdateMinProblem(SQP_deldel_Qcalc_plan_mpc_, SQP_del_gcalc_plan_mpc_);
        QP_MPC_Planner_.DeleteSubjectToAx();
        QP_MPC_Planner_.DeleteSubjectToX();
    
        const_A_mpc_.setZero( const_num, input_num);
        const_ub_mpc_.setZero(const_num, 1);
        const_lb_mpc_.setZero(const_num, 1);

        int constraint_index = 0;
        int calc_index = 0;
        int calc_index2 = 0;
        ////VRP Constraint
        Eigen::MatrixXd Si_mpc; Si_mpc.setZero(1, N_plan_mpc);
        //for(int i = 0; i < N_plan_mpc; i++)
        for(int i = 0; i < 1; i++)
        {
            const_SQP_phi_mpc_.setZero(input_num, input_num);
            const_SQP_phi_mpc_calc_.setZero(input_num, input_num);
            const_SQP_pi_mpc_.setZero(input_num, 1);
            const_SQP_ri_mpc_.setZero(1, 1);

            Si_mpc.setZero(1, N_plan_mpc);
            Si_mpc(0, i) = 1;

            //X max
            const_SQP_phi_mpc_ = 0.5*(IS_FIPM_SQP_x_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc)
                                     +IS_FIPM_SQP_x_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc).transpose());

            const_SQP_pi_mpc_  = IS_FIPM_SQP_x_pi1_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_state)*MPC_Planner_state_mpc_

                               + IS_FIPM_SQP_x_pi2_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1)

                               - zmp_max_x_mpc_(i)*IS_FIPM_SQP_x_pi3_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1);

            const_SQP_ri_mpc_ = MPC_Planner_state_mpc_.transpose()*IS_FIPM_SQP_x_ri1_N_plan_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Planner_state_mpc_

                              + IS_FIPM_SQP_x_ri2_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                              - zmp_max_x_mpc_(i)*IS_FIPM_SQP_x_ri3_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                              - GRAVITY*b_*b_*zmp_max_x_mpc_(i)*MatrixXd::Identity(1,1);

            const_SQP_hi_mpc_ = MPC_Planner_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;

            const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - 1e+3*MatrixXd::Identity(1,1);
            constraint_index += 1;

            //X min
            const_SQP_phi_mpc_ = 0.5*(IS_FIPM_SQP_x_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc)
                                     +IS_FIPM_SQP_x_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc).transpose());

            const_SQP_pi_mpc_  = IS_FIPM_SQP_x_pi1_N_plan_mpc_.block (calc_index, 0, 3*N_plan_mpc, 3*N_state)*MPC_Planner_state_mpc_

                               + IS_FIPM_SQP_x_pi2_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1)

                               - zmp_min_x_mpc_(i)*IS_FIPM_SQP_x_pi3_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1);

            const_SQP_ri_mpc_  = MPC_Planner_state_mpc_.transpose()*IS_FIPM_SQP_x_ri1_N_plan_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Planner_state_mpc_

                               + IS_FIPM_SQP_x_ri2_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                               - zmp_min_x_mpc_(i)*IS_FIPM_SQP_x_ri3_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                               - GRAVITY*b_*b_*zmp_min_x_mpc_(i)*MatrixXd::Identity(1,1);

            const_SQP_hi_mpc_  = MPC_Planner_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;

            const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) = + 1e+3*MatrixXd::Identity(1,1);
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            constraint_index += 1;

            //Y direction
            //Y max
            const_SQP_phi_mpc_ = 0.5*(IS_FIPM_SQP_y_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc)
                                     +IS_FIPM_SQP_y_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc).transpose());

            const_SQP_pi_mpc_  = IS_FIPM_SQP_y_pi1_N_plan_mpc_.block (calc_index, 0, 3*N_plan_mpc, 3*N_state)*MPC_Planner_state_mpc_

                               + IS_FIPM_SQP_y_pi2_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1)

                               - zmp_max_y_mpc_(i)*IS_FIPM_SQP_y_pi3_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1);

            const_SQP_ri_mpc_  = MPC_Planner_state_mpc_.transpose()*IS_FIPM_SQP_y_ri1_N_plan_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Planner_state_mpc_

                               + IS_FIPM_SQP_y_ri2_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                               - zmp_max_y_mpc_(i)*IS_FIPM_SQP_y_ri3_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                               - GRAVITY*b_*b_*zmp_max_y_mpc_(i)*MatrixXd::Identity(1,1);

            const_SQP_hi_mpc_  = MPC_Planner_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;

            const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - 1e+3*MatrixXd::Identity(1,1);
            constraint_index += 1;

            //Y min
            const_SQP_phi_mpc_ = 0.5*(IS_FIPM_SQP_y_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc)
                                     +IS_FIPM_SQP_y_phi_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 3*N_plan_mpc).transpose());

            const_SQP_pi_mpc_  = IS_FIPM_SQP_y_pi1_N_plan_mpc_.block (calc_index, 0, 3*N_plan_mpc, 3*N_state)*MPC_Planner_state_mpc_

                               + IS_FIPM_SQP_y_pi2_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1)

                               - zmp_min_y_mpc_(i)*IS_FIPM_SQP_y_pi3_N_plan_mpc_.block(calc_index, 0, 3*N_plan_mpc, 1);

            const_SQP_ri_mpc_  = MPC_Planner_state_mpc_.transpose()*IS_FIPM_SQP_y_ri1_N_plan_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Planner_state_mpc_

                               + IS_FIPM_SQP_y_ri2_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                               - zmp_min_y_mpc_(i)*IS_FIPM_SQP_y_ri3_N_plan_mpc_.block(i, 0, 1, 3*N_state)*MPC_Planner_state_mpc_

                               - GRAVITY*b_*b_*zmp_min_y_mpc_(i)*MatrixXd::Identity(1,1);

            const_SQP_hi_mpc_  = MPC_Planner_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;

            const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) = + 1e+3*MatrixXd::Identity(1,1);
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            constraint_index += 1;

            calc_index  += 3*N_plan_mpc;
            calc_index2 += 3*N_state;
        }

        //IS EQ
        Eigen::MatrixXd Const_b_eq_;  Const_b_eq_.setZero(3,1);

        Const_b_eq_(0,0) = (w_/(1 - lambda_is_calc))*(MPC_Planner_state_mpc_(0) + MPC_Planner_state_mpc_(1)/w_ - MPC_Planner_state_mpc_(2))
                           -(pow(lambda_is_calc, N_plan_mpc)/(1 - pow(lambda_is_calc,N_step))*(b_IS_step_mpc_.transpose()*Pv_dot_ref_mpc_.col(0))(0,0));
        Const_b_eq_(1,0) = (w_/(1 - lambda_is_calc))*(MPC_Planner_state_mpc_(3) + MPC_Planner_state_mpc_(4)/w_ - MPC_Planner_state_mpc_(5))
                           -(pow(lambda_is_calc, N_plan_mpc)/(1 + pow(lambda_is_calc,N_step))*(b_IS_step_mpc_.transpose()*Pv_dot_ref_mpc_.col(1))(0,0));
        Const_b_eq_(2,0) = (w_/(1 - lambda_is_calc))*(MPC_Planner_state_mpc_(6) + MPC_Planner_state_mpc_(7)/w_ - MPC_Planner_state_mpc_(8));

        //X direction
        const_SQP_phi_mpc_.setZero(input_num, input_num);
        const_SQP_pi_mpc_ = (b_IS_plan_mpc_.transpose()*SUx_plan_mpc_).transpose();
        const_SQP_ri_mpc_ = - Const_b_eq_.block(0, 0, 1, 1);
        //const_SQP_hi_mpc_ = MPC_Planner_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;
        const_SQP_hi_mpc_ = const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;
    
        //const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_).transpose();
        const_A_mpc_.row(constraint_index) = (const_SQP_pi_mpc_).transpose();
        const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        constraint_index += 1;
    
        //Y direction
        const_SQP_phi_mpc_.setZero(input_num, input_num);
        const_SQP_pi_mpc_ = (b_IS_plan_mpc_.transpose()*SUy_plan_mpc_).transpose();
        const_SQP_ri_mpc_ = - Const_b_eq_.block(1, 0, 1, 1);
        //const_SQP_hi_mpc_ = MPC_Planner_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;
        const_SQP_hi_mpc_ = const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;

        //const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_).transpose();
        const_A_mpc_.row(constraint_index) = (const_SQP_pi_mpc_).transpose();
        const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        constraint_index += 1;
    
        //Z direction
        const_SQP_phi_mpc_.setZero(input_num, input_num);
        const_SQP_pi_mpc_ = (b_IS_plan_mpc_.transpose()*SUz_plan_mpc_).transpose();
        const_SQP_ri_mpc_ = - Const_b_eq_.block(2, 0, 1, 1);
        //const_SQP_hi_mpc_ = MPC_Planner_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;
        const_SQP_hi_mpc_ = const_SQP_pi_mpc_.transpose()*MPC_Planner_u_mpc_ + const_SQP_ri_mpc_;

        //const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Planner_u_mpc_ + const_SQP_pi_mpc_).transpose();
        const_A_mpc_.row(constraint_index) = (const_SQP_pi_mpc_).transpose();
        const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        constraint_index += 1;

        QP_MPC_Planner_.UpdateSubjectToAx(const_A_mpc_, const_lb_mpc_, const_ub_mpc_);

        if(QP_MPC_Planner_.SolveQPoases(100, MPC_Planner_SQP_du_mpc_))
        {
            if(print_sec_bool_)
                cout << "IS FIPM Planner MPC Solved" << endl;
        
            MPC_Planner_u_mpc_ = MPC_Planner_u_mpc_ + MPC_Planner_SQP_du_mpc_;

            if(s == sqp_iter - 1)
            {
                Planner_State_Prev_mpc_.row(0).transpose() = Pcps_plan_mpc_*MPC_Planner_state_mpc_.segment(0,3) + Pcpu_plan_mpc_*MPC_Planner_u_mpc_.segment(0*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(1).transpose() = Pcvs_plan_mpc_*MPC_Planner_state_mpc_.segment(0,3) + Pcvu_plan_mpc_*MPC_Planner_u_mpc_.segment(0*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(2).transpose() = Pvps_plan_mpc_*MPC_Planner_state_mpc_.segment(0,3) + Pvpu_plan_mpc_*MPC_Planner_u_mpc_.segment(0*N_plan_mpc, N_plan_mpc);

                Planner_State_Prev_mpc_.row(3).transpose() = Pcps_plan_mpc_*MPC_Planner_state_mpc_.segment(3,3) + Pcpu_plan_mpc_*MPC_Planner_u_mpc_.segment(1*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(4).transpose() = Pcvs_plan_mpc_*MPC_Planner_state_mpc_.segment(3,3) + Pcvu_plan_mpc_*MPC_Planner_u_mpc_.segment(1*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(5).transpose() = Pvps_plan_mpc_*MPC_Planner_state_mpc_.segment(3,3) + Pvpu_plan_mpc_*MPC_Planner_u_mpc_.segment(1*N_plan_mpc, N_plan_mpc);

                Planner_State_Prev_mpc_.row(6).transpose() = Pcps_plan_mpc_*MPC_Planner_state_mpc_.segment(6,3) + Pcpu_plan_mpc_*MPC_Planner_u_mpc_.segment(2*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(7).transpose() = Pcvs_plan_mpc_*MPC_Planner_state_mpc_.segment(6,3) + Pcvu_plan_mpc_*MPC_Planner_u_mpc_.segment(2*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(8).transpose() = Pvps_plan_mpc_*MPC_Planner_state_mpc_.segment(6,3) + Pvpu_plan_mpc_*MPC_Planner_u_mpc_.segment(2*N_plan_mpc, N_plan_mpc);

                MPC_Planner_state_mpc_.segment(0,3) = A_mpc_*MPC_Planner_state_mpc_.segment(0,3) + B_mpc_*MPC_Planner_u_mpc_(0*N_plan_mpc);
                MPC_Planner_state_mpc_.segment(3,3) = A_mpc_*MPC_Planner_state_mpc_.segment(3,3) + B_mpc_*MPC_Planner_u_mpc_(1*N_plan_mpc);
                MPC_Planner_state_mpc_.segment(6,3) = A_mpc_*MPC_Planner_state_mpc_.segment(6,3) + B_mpc_*MPC_Planner_u_mpc_(2*N_plan_mpc);

                MPC_Planner_u_mpc_sep_(0) = MPC_Planner_u_mpc_(0*N_plan_mpc);
                MPC_Planner_u_mpc_sep_(1) = MPC_Planner_u_mpc_(1*N_plan_mpc);
                MPC_Planner_u_mpc_sep_(2) = MPC_Planner_u_mpc_(2*N_plan_mpc);   
            }
        }
        else
        { 
            cout << "IS FIPM Planner MPC Not Solved" << endl;
            cout << matlab_tick << endl;

            MPC_Planner_u_mpc_ = MPC_Planner_u_mpc_ + MPC_Planner_SQP_du_mpc_;

            if(s == sqp_iter - 1)
            {
                Planner_State_Prev_mpc_.row(0).transpose() = Pcps_plan_mpc_*MPC_Planner_state_mpc_.segment(0,3) + Pcpu_plan_mpc_*MPC_Planner_u_mpc_.segment(0*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(1).transpose() = Pcvs_plan_mpc_*MPC_Planner_state_mpc_.segment(0,3) + Pcvu_plan_mpc_*MPC_Planner_u_mpc_.segment(0*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(2).transpose() = Pvps_plan_mpc_*MPC_Planner_state_mpc_.segment(0,3) + Pvpu_plan_mpc_*MPC_Planner_u_mpc_.segment(0*N_plan_mpc, N_plan_mpc);

                Planner_State_Prev_mpc_.row(3).transpose() = Pcps_plan_mpc_*MPC_Planner_state_mpc_.segment(3,3) + Pcpu_plan_mpc_*MPC_Planner_u_mpc_.segment(1*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(4).transpose() = Pcvs_plan_mpc_*MPC_Planner_state_mpc_.segment(3,3) + Pcvu_plan_mpc_*MPC_Planner_u_mpc_.segment(1*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(5).transpose() = Pvps_plan_mpc_*MPC_Planner_state_mpc_.segment(3,3) + Pvpu_plan_mpc_*MPC_Planner_u_mpc_.segment(1*N_plan_mpc, N_plan_mpc);

                Planner_State_Prev_mpc_.row(6).transpose() = Pcps_plan_mpc_*MPC_Planner_state_mpc_.segment(6,3) + Pcpu_plan_mpc_*MPC_Planner_u_mpc_.segment(2*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(7).transpose() = Pcvs_plan_mpc_*MPC_Planner_state_mpc_.segment(6,3) + Pcvu_plan_mpc_*MPC_Planner_u_mpc_.segment(2*N_plan_mpc, N_plan_mpc);
                Planner_State_Prev_mpc_.row(8).transpose() = Pvps_plan_mpc_*MPC_Planner_state_mpc_.segment(6,3) + Pvpu_plan_mpc_*MPC_Planner_u_mpc_.segment(2*N_plan_mpc, N_plan_mpc);

                MPC_Planner_state_mpc_.segment(0,3) = A_mpc_*MPC_Planner_state_mpc_.segment(0,3) + B_mpc_*MPC_Planner_u_mpc_(0*N_plan_mpc);
                MPC_Planner_state_mpc_.segment(3,3) = A_mpc_*MPC_Planner_state_mpc_.segment(3,3) + B_mpc_*MPC_Planner_u_mpc_(1*N_plan_mpc);
                MPC_Planner_state_mpc_.segment(6,3) = A_mpc_*MPC_Planner_state_mpc_.segment(6,3) + B_mpc_*MPC_Planner_u_mpc_(2*N_plan_mpc);

                MPC_Planner_u_mpc_sep_(0) = MPC_Planner_u_mpc_(0*N_plan_mpc);
                MPC_Planner_u_mpc_sep_(1) = MPC_Planner_u_mpc_(1*N_plan_mpc);
                MPC_Planner_u_mpc_sep_(2) = MPC_Planner_u_mpc_(2*N_plan_mpc);
            }
        }
    }

    
    if(e_planner_data_txt_[0].is_open())
    {
        e_planner_data_txt_[0] << walking_tick_mpc_         << "," << N_plan_mpc                << "," << 0                         << ","
                                << MPC_Planner_state_mpc_(0) << "," << MPC_Planner_state_mpc_(3) << "," << MPC_Planner_state_mpc_(6) << ","
                                << MPC_Planner_state_mpc_(1) << "," << MPC_Planner_state_mpc_(4) << "," << MPC_Planner_state_mpc_(7) << ","
                                << MPC_Planner_state_mpc_(2) << "," << MPC_Planner_state_mpc_(5) << "," << MPC_Planner_state_mpc_(8) << ","
                                << Pv_x_ref(0)               << "," << Pv_y_ref(0)               << "," << Pv_z_ref(0)               << ","
                                << zmp_max_x_mpc_(0)         << "," << zmp_max_y_mpc_(0)         << "," << 0                         << ","
                                << zmp_min_x_mpc_(0)         << "," << zmp_min_y_mpc_(0)         << "," << 0                         << ","
                                << endl;
    }
    else
    {
        cout << "Error: Unable to write to e_planner_data_txt_[0]. File is not open." << endl;
    }

    Eigen::VectorXd data_save_calc; data_save_calc.setZero(3*N_plan_mpc);
    data_save_calc << Pv_x_ref, Pv_y_ref, Pv_z_ref;
    e_planner_data_txt_[1] << data_save_calc.transpose() << endl;
    data_save_calc << Planner_State_Prev_mpc_.row(0).transpose(), Planner_State_Prev_mpc_.row(3).transpose(), Planner_State_Prev_mpc_.row(6).transpose();
    e_planner_data_txt_[2] << data_save_calc.transpose() << endl;
    data_save_calc << Planner_State_Prev_mpc_.row(2).transpose(), Planner_State_Prev_mpc_.row(5).transpose(), Planner_State_Prev_mpc_.row(8).transpose();
    e_planner_data_txt_[3] << data_save_calc.transpose() << endl;
    data_save_calc << zmp_max_x_mpc_, zmp_max_y_mpc_, 0*zmp_max_y_mpc_;
    e_planner_data_txt_[4] << data_save_calc.transpose() << endl;
}

void fipmPlanner::MPC_State_Step_Change()
{
    if((walking_tick_mpc_ - (t_start_mpc_ + t_total_mpc_) >= -hz_/thread3_hz_) && (current_step_num_mpc_ < total_step_num_mpc_ - 1))
    {
        Eigen::Vector3d var_after_step_change, var_before_step_change, frame_pos_diff;
        Eigen::Matrix3d frame_rot_diff;

        frame_rot_diff = DyrosMath::rotateWithZ(-foot_step_support_frame_mpc_(current_step_num_mpc_, 5));
        frame_pos_diff(0) = foot_step_support_frame_mpc_(current_step_num_mpc_,0);
        frame_pos_diff(1) = foot_step_support_frame_mpc_(current_step_num_mpc_,1);
        frame_pos_diff(2) = foot_step_support_frame_mpc_(current_step_num_mpc_,2);

        //com pos step change
        var_before_step_change(0) = MPC_Planner_state_mpc_(0);
        var_before_step_change(1) = MPC_Planner_state_mpc_(3);
        var_before_step_change(2) = MPC_Planner_state_mpc_(6);
        var_after_step_change = frame_rot_diff*(var_before_step_change - frame_pos_diff);

        MPC_Planner_state_mpc_(0) = var_after_step_change(0);
        MPC_Planner_state_mpc_(3) = var_after_step_change(1);
        MPC_Planner_state_mpc_(6) = var_after_step_change(2);

        //com vel step change
        var_before_step_change(0) = MPC_Planner_state_mpc_(1);
        var_before_step_change(1) = MPC_Planner_state_mpc_(4);
        var_before_step_change(2) = MPC_Planner_state_mpc_(7);
        var_after_step_change = frame_rot_diff*var_before_step_change;

        MPC_Planner_state_mpc_(1) = var_after_step_change(0);
        MPC_Planner_state_mpc_(4) = var_after_step_change(1);
        MPC_Planner_state_mpc_(7) = var_after_step_change(2);

        //vrp pos step change
        var_before_step_change(0) = MPC_Planner_state_mpc_(2);
        var_before_step_change(1) = MPC_Planner_state_mpc_(5);
        var_before_step_change(2) = MPC_Planner_state_mpc_(8);
        var_after_step_change = frame_rot_diff*(var_before_step_change - frame_pos_diff);

        MPC_Planner_state_mpc_(2) = var_after_step_change(0);
        MPC_Planner_state_mpc_(5) = var_after_step_change(1);
        MPC_Planner_state_mpc_(8) = var_after_step_change(2);
    }
}