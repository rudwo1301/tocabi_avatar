#include "dcm_stabilizer.h"
#include <iostream>

using namespace std;

std::vector<std::ofstream> e_stabilizer_data_txt_;

dcmStabilizer::dcmStabilizer()
{
    cout << "dcmStabilizer constructor called" << endl;
}

void dcmStabilizer::ofstream_open()
{
    if (!e_stabilizer_data_txt_.empty() && e_stabilizer_data_txt_[0].is_open()) {
        return;
    }

    e_stabilizer_data_txt_.resize(num_ofstream);

    for (int i = 0; i < num_ofstream; ++i) {
        if (!e_stabilizer_data_txt_[i].is_open()) {
            e_stabilizer_data_txt_[i].open("/home/econom2-20/data/dcm_stabilizer/e_data_" + std::to_string(i) + ".txt");
        }
    }
}

void dcmStabilizer::ofstream_close()
{
    for (auto &file_stream : e_stabilizer_data_txt_) {
        if (file_stream.is_open()) {
            file_stream.close();
        }
    }
}

void dcmStabilizer::IS_FIPM_3D_DCM_Stabilizer_MPC(double mpc_freq, double mpc_preview_window)
{
    double Q_dcm_x, Q_dcm_y, Q_dcm_z, R_dcm_x, R_dcm_y, R_dcm_z, R_dalp, R_df_x, R_df_y;

    Q_dcm_x = 1e-0; R_dcm_x = 1e-2; R_dalp = 1e+1; R_df_x = 1e-0;
    Q_dcm_y = 1e-0; R_dcm_y = 1e-2;                R_df_y = 1e+2; //need tuning 1e+2 - 1e+3
    Q_dcm_z = 9e-1; R_dcm_z = 1e-1; //for 0.9 step time

    int mpc_tick = walking_tick_mpc_ - com_start_tick_mpc_;
    double mpc_synchro_hz = hz_/mpc_freq;
    const int N_stab_mpc = mpc_preview_window*mpc_freq;
    const int N_step = t_total_const_mpc_/mpc_synchro_hz;
    const int N_state = 3;
    double dt_stab_mpc = 1/mpc_freq;

    double lambda_is_calc = exp(-w_*dt_stab_mpc);

    int input_num = 3 * N_stab_mpc      + 2            + 1 * step_time_adj_candidate_num_             + 2 * step_time_adj_candidate_num_;
    //              xyz VRP               delf           alpha                                          aux
    int const_num = 4 * N_stab_mpc + 3  + 4            + 1 * step_time_adj_candidate_num_ + 1         + 2 * step_time_adj_candidate_num_;
    //              VRP min,max      IS   delf min,max   alpha min,max                      alpha sum   alpha, delf, aux eq

    if(MPC_first_loop_ == 0)
    {
        cout << "Initialization of IS FIPM 3D DCM Stabilizer MPC." << endl;
        b_ = 1/w_;
        ofstream_open();

        Cdp_mpc_.resize(1, N_state);
        Cdp_mpc_ << 1, 1/w_, 0;

        Pdps_stab_mpc_.resize(N_stab_mpc, N_state);
        Pcps_stab_mpc_.resize(N_stab_mpc, N_state);
        Pcvs_stab_mpc_.resize(N_stab_mpc, N_state);
        Pvps_stab_mpc_.resize(N_stab_mpc, N_state);

        Eigen::MatrixXd Ps_calc;
        Ps_calc.resize(N_state,N_state);
        Ps_calc = A_mpc_;
        
        Pdpu_stab_mpc_.setZero(N_stab_mpc, N_stab_mpc);
        Pcpu_stab_mpc_.setZero(N_stab_mpc, N_stab_mpc);
        Pcvu_stab_mpc_.setZero(N_stab_mpc, N_stab_mpc);
        Pvpu_stab_mpc_.setZero(N_stab_mpc, N_stab_mpc);
        
        Eigen::MatrixXd Pu_calc;
        Pu_calc.setZero(N_state,N_stab_mpc);
        
        for(int i = 0; i < N_stab_mpc; i++)
        {
            Pdps_stab_mpc_.row(i) = Cdp_mpc_*Ps_calc;
            Pcps_stab_mpc_.row(i) = Ccp_mpc_*Ps_calc;
            Pcvs_stab_mpc_.row(i) = Ccv_mpc_*Ps_calc;
            Pvps_stab_mpc_.row(i) = Cvp_mpc_*Ps_calc;
            Ps_calc = Ps_calc*A_mpc_;

            Pu_calc.col(i) = B_mpc_;
            Pdpu_stab_mpc_.row(i) = Cdp_mpc_*Pu_calc;
            Pcpu_stab_mpc_.row(i) = Ccp_mpc_*Pu_calc;
            Pcvu_stab_mpc_.row(i) = Ccv_mpc_*Pu_calc;
            Pvpu_stab_mpc_.row(i) = Cvp_mpc_*Pu_calc;
            Pu_calc = A_mpc_*Pu_calc;
        }

        QP_MPC_Stabilizer_.InitializeProblemSize(input_num, const_num);

        Qmat_stab_mpc_Q_.resize(N_stab_mpc, N_stab_mpc); //DCM
        Qmat_stab_mpc_Q_.setIdentity();
        
        Qmat_stab_mpc_R_.resize(N_stab_mpc, N_stab_mpc); //VRP
        Qmat_stab_mpc_R_.setIdentity();

        Qcalc_stab_mpc_.setZero(input_num, input_num);

        gcalc_stab_mpc_.setZero(input_num, 1);
        gxpcalc_stab_mpc_.setZero(input_num, N_stab_mpc);
        gypcalc_stab_mpc_.setZero(input_num, N_stab_mpc);

        MPC_Stabilizer_u_mpc_.setZero(input_num); // VRP, delf, eps

        MPC_Stabilizer_SQP_du_mpc_.setZero(input_num);
        
        MPC_Stabilizer_alpha_mpc_.setZero(step_time_adj_candidate_num_);
        MPC_Stabilizer_alpha_mpc_(0) = 1;

        MPC_Stabilizer_aux_mpc_.setZero(2*step_time_adj_candidate_num_);
        MPC_Stabilizer_aux_mpc_x_.setZero(step_time_adj_candidate_num_);
        MPC_Stabilizer_aux_mpc_y_.setZero(step_time_adj_candidate_num_);

        MPC_Stabilizer_delf_mpc_.setZero(2);
        MPC_Stabilizer_delf_mpc_x_.setZero(1);
        MPC_Stabilizer_delf_mpc_y_.setZero(1);

        int input_index_calc = 0;

        SUp_stab_mpc_.setZero(3*N_stab_mpc, input_num);
        SUp_stab_mpc_.block  (0, input_index_calc, 3*N_stab_mpc, 3*N_stab_mpc) = MatrixXd::Identity(3*N_stab_mpc, 3*N_stab_mpc);

        SUpx_stab_mpc_.setZero(N_stab_mpc, 3*N_stab_mpc);
        SUpx_stab_mpc_.block  (0, 0*N_stab_mpc, N_stab_mpc, N_stab_mpc) = MatrixXd::Identity(N_stab_mpc, N_stab_mpc);
        SUpxp_stab_mpc_ = SUpx_stab_mpc_*SUp_stab_mpc_;
        SUpy_stab_mpc_.setZero(N_stab_mpc, 3*N_stab_mpc);
        SUpy_stab_mpc_.block  (0, 1*N_stab_mpc, N_stab_mpc, N_stab_mpc) = MatrixXd::Identity(N_stab_mpc, N_stab_mpc);
        SUpyp_stab_mpc_ = SUpy_stab_mpc_*SUp_stab_mpc_;
        SUpz_stab_mpc_.setZero(N_stab_mpc, 3*N_stab_mpc);
        SUpz_stab_mpc_.block  (0, 2*N_stab_mpc, N_stab_mpc, N_stab_mpc) = MatrixXd::Identity(N_stab_mpc, N_stab_mpc);
        SUpzp_stab_mpc_ = SUpz_stab_mpc_*SUp_stab_mpc_;

        input_index_calc += 3*N_stab_mpc;

        cout << "Selection Matrix VRP Complete" << endl;
        cout << "input_num: " << input_num << endl;
        cout << "input_index_calc: " << input_index_calc << endl << endl;

        SUf_stab_mpc_.setZero(2, input_num);
        SUf_stab_mpc_.block  (0, input_index_calc, 2, 2) = MatrixXd::Identity(2, 2);

        SUfx_stab_mpc_.setZero(1, 2);
        SUfx_stab_mpc_.block  (0, 0*1, 1, 1) = MatrixXd::Identity(1, 1);
        SUfy_stab_mpc_.setZero(1, 2);
        SUfy_stab_mpc_.block  (0, 1*1, 1, 1) = MatrixXd::Identity(1, 1);
        
        SUfxf_stab_mpc_ = SUfx_stab_mpc_*SUf_stab_mpc_;
        SUfyf_stab_mpc_ = SUfy_stab_mpc_*SUf_stab_mpc_;

        input_index_calc += 2*1;

        cout << "Selection Matrix delf Complete" << endl;
        cout << "input_num: "                    << input_num        << endl;
        cout << "input_index_calc: "             << input_index_calc << endl << endl;

        SUalpha_stab_mpc_.setZero(1 * step_time_adj_candidate_num_, input_num);
        SUalpha_stab_mpc_.block  (0, input_index_calc, 1 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_) = MatrixXd::Identity(1 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_);

        input_index_calc += 1 * step_time_adj_candidate_num_;

        cout << "Selection Matrix alpha Complete" << endl;
        cout << "input_num: "                     << input_num        << endl;
        cout << "input_index_calc: "              << input_index_calc << endl << endl;

        SUaux_stab_mpc_.setZero(2 * step_time_adj_candidate_num_, input_num);
        SUaux_stab_mpc_.block  (0, input_index_calc, 2 * step_time_adj_candidate_num_, 2 * step_time_adj_candidate_num_) = MatrixXd::Identity(2 * step_time_adj_candidate_num_, 2 * step_time_adj_candidate_num_);

        SUauxx_stab_mpc_.setZero(1 * step_time_adj_candidate_num_, 2 * step_time_adj_candidate_num_);
        SUauxx_stab_mpc_.block  (0, 0 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_) = MatrixXd::Identity(1 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_);
        SUauxy_stab_mpc_.setZero(1 * step_time_adj_candidate_num_, 2 * step_time_adj_candidate_num_);
        SUauxy_stab_mpc_.block  (0, 1 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_) = MatrixXd::Identity(1 * step_time_adj_candidate_num_, 1 * step_time_adj_candidate_num_);

        input_index_calc += 2 * step_time_adj_candidate_num_;

        cout << "Selection Matrix aux Complete" << endl;
        cout << "input_num: "                   << input_num        << endl;
        cout << "input_index_calc: "            << input_index_calc << endl << endl;

        Qcalc_stab_mpc_ = SUp_stab_mpc_.transpose()*(SUpx_stab_mpc_.transpose()*(Pdpu_stab_mpc_.transpose()*Q_dcm_x*MatrixXd::Identity(N_stab_mpc, N_stab_mpc)*Pdpu_stab_mpc_ + R_dcm_x*Qmat_stab_mpc_R_)*SUpx_stab_mpc_
                                                    +SUpy_stab_mpc_.transpose()*(Pdpu_stab_mpc_.transpose()*Q_dcm_y*MatrixXd::Identity(N_stab_mpc, N_stab_mpc)*Pdpu_stab_mpc_ + R_dcm_y*Qmat_stab_mpc_R_)*SUpy_stab_mpc_
                                                    +SUpz_stab_mpc_.transpose()*(Pdpu_stab_mpc_.transpose()*Q_dcm_z*MatrixXd::Identity(N_stab_mpc, N_stab_mpc)*Pdpu_stab_mpc_ + R_dcm_z*Qmat_stab_mpc_R_)*SUpz_stab_mpc_)*SUp_stab_mpc_
                                                    
                        + SUf_stab_mpc_.transpose()*(SUfx_stab_mpc_.transpose()*R_df_x*SUfx_stab_mpc_
                                                   + SUfy_stab_mpc_.transpose()*R_df_y*SUfy_stab_mpc_)*SUf_stab_mpc_;

        gxpcalc_stab_mpc_ = SUp_stab_mpc_.transpose()*SUpx_stab_mpc_.transpose()*Pdpu_stab_mpc_.transpose()*Q_dcm_x*MatrixXd::Identity(N_stab_mpc, N_stab_mpc);
        gypcalc_stab_mpc_ = SUp_stab_mpc_.transpose()*SUpy_stab_mpc_.transpose()*Pdpu_stab_mpc_.transpose()*Q_dcm_y*MatrixXd::Identity(N_stab_mpc, N_stab_mpc);
        gzpcalc_stab_mpc_ = SUp_stab_mpc_.transpose()*SUpz_stab_mpc_.transpose()*Pdpu_stab_mpc_.transpose()*Q_dcm_z*MatrixXd::Identity(N_stab_mpc, N_stab_mpc);

        gxdfcalc_stab_mpc_ = SUf_stab_mpc_.transpose()*SUfx_stab_mpc_.transpose()*R_df_x;
        gydfcalc_stab_mpc_ = SUf_stab_mpc_.transpose()*SUfy_stab_mpc_.transpose()*R_df_y;


        Sf1_stab_mpc_.setZero(N_stab_mpc, step_time_adj_candidate_num_);

        MPC_Stabilizer_u_mpc_sep_.setZero(3);

        IS_FIPM_SQP_x_phi1_N_stab_mpc_.setZero(input_num*N_stab_mpc, input_num);
        IS_FIPM_SQP_x_phi2_N_stab_mpc_.setZero(input_num*N_stab_mpc, N_stab_mpc);

        IS_FIPM_SQP_x_pi1_N_stab_mpc_.setZero (input_num*N_stab_mpc, 3*N_state);
        IS_FIPM_SQP_x_pi2_N_stab_mpc_.setZero(input_num*N_stab_mpc, 1);
        IS_FIPM_SQP_x_pi3_N_stab_mpc_.setZero(input_num*N_stab_mpc, 1);
        IS_FIPM_SQP_x_pi4_N_stab_mpc_.setZero(3*N_state*N_stab_mpc, N_stab_mpc);

        IS_FIPM_SQP_x_ri1_N_stab_mpc_.setZero (3*N_state*N_stab_mpc, 3*N_state);
        IS_FIPM_SQP_x_ri2_N_stab_mpc_.setZero(        1*N_stab_mpc, 3*N_state);
        IS_FIPM_SQP_x_ri3_N_stab_mpc_.setZero(        1*N_stab_mpc, 3*N_state);

        IS_FIPM_SQP_y_phi1_N_stab_mpc_.setZero(input_num*N_stab_mpc, input_num);
        IS_FIPM_SQP_y_phi2_N_stab_mpc_.setZero(input_num*N_stab_mpc, N_stab_mpc);

        IS_FIPM_SQP_y_pi1_N_stab_mpc_.setZero (input_num*N_stab_mpc, 3*N_state);
        IS_FIPM_SQP_y_pi2_N_stab_mpc_.setZero(input_num*N_stab_mpc, 1);
        IS_FIPM_SQP_y_pi3_N_stab_mpc_.setZero(input_num*N_stab_mpc, 1);
        IS_FIPM_SQP_y_pi4_N_stab_mpc_.setZero(3*N_state*N_stab_mpc, N_stab_mpc);

        IS_FIPM_SQP_y_ri1_N_stab_mpc_.setZero (3*N_state*N_stab_mpc, 3*N_state);
        IS_FIPM_SQP_y_ri2_N_stab_mpc_.setZero(        1*N_stab_mpc, 3*N_state);
        IS_FIPM_SQP_y_ri3_N_stab_mpc_.setZero(        1*N_stab_mpc, 3*N_state);
        
        int calc_index  = 0;
        int calc_index2 = 0;
        Eigen::MatrixXd Si_mpc; Si_mpc.setZero(1, N_stab_mpc);

        for(int i = 0; i < N_stab_mpc; i++)
        {
            Si_mpc.setZero(1, N_stab_mpc);
            Si_mpc(0, i) = 1;

            IS_FIPM_SQP_x_phi1_N_stab_mpc_.block(calc_index, 0, input_num, input_num)  = (Pvpu_stab_mpc_.row(i)*SUpx_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_)
                                                                                      - (Pvpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcpu_stab_mpc_.row(i)*SUpx_stab_mpc_*SUp_stab_mpc_);
                                                                                     
            IS_FIPM_SQP_x_phi2_N_stab_mpc_.block(calc_index, 0, input_num, N_stab_mpc)= ((Pvpu_stab_mpc_ - Pcpu_stab_mpc_).row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*Si_mpc;

            IS_FIPM_SQP_y_phi1_N_stab_mpc_.block(calc_index, 0, input_num, input_num)  = (Pvpu_stab_mpc_.row(i)*SUpy_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_)
                                                                                      - (Pvpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcpu_stab_mpc_.row(i)*SUpy_stab_mpc_*SUp_stab_mpc_);

            IS_FIPM_SQP_y_phi2_N_stab_mpc_.block(calc_index, 0, input_num, N_stab_mpc)= ((Pvpu_stab_mpc_ - Pcpu_stab_mpc_).row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*Si_mpc;

            IS_FIPM_SQP_x_pi1_N_stab_mpc_.block(calc_index,  0, input_num, 3*N_state)  = (Pcpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*(Pvps_stab_mpc_.row(i)*ssx_stab_mpc_)
                                                                                      + (Pvpu_stab_mpc_.row(i)*SUpx_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssz_stab_mpc_)
                                                                                      - (Pcpu_stab_mpc_.row(i)*SUpx_stab_mpc_*SUp_stab_mpc_).transpose()*(Pvps_stab_mpc_.row(i)*ssz_stab_mpc_)
                                                                                      - (Pvpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssx_stab_mpc_);

            IS_FIPM_SQP_x_pi2_N_stab_mpc_.block(calc_index, 0, input_num, 1)          = GRAVITY*b_*b_*(Si_mpc*Pcpu_stab_mpc_*SUpx_stab_mpc_*SUp_stab_mpc_).transpose();

            IS_FIPM_SQP_x_pi3_N_stab_mpc_.block(calc_index, 0, input_num, 1)          = ((Pcpu_stab_mpc_ - Pvpu_stab_mpc_).row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose();
            
            IS_FIPM_SQP_y_pi1_N_stab_mpc_.block(calc_index,  0, input_num, 3*N_state)  = (Pcpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*(Pvps_stab_mpc_.row(i)*ssy_stab_mpc_)
                                                                                      + (Pvpu_stab_mpc_.row(i)*SUpy_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssz_stab_mpc_)
                                                                                      - (Pcpu_stab_mpc_.row(i)*SUpy_stab_mpc_*SUp_stab_mpc_).transpose()*(Pvps_stab_mpc_.row(i)*ssz_stab_mpc_)
                                                                                      - (Pvpu_stab_mpc_.row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssy_stab_mpc_);

            IS_FIPM_SQP_y_pi2_N_stab_mpc_.block(calc_index, 0, input_num, 1)          = GRAVITY*b_*b_*(Pcpu_stab_mpc_.row(i)*SUpy_stab_mpc_*SUp_stab_mpc_).transpose();
            
            IS_FIPM_SQP_y_pi3_N_stab_mpc_.block(calc_index, 0, input_num, 1)          = ((Pcpu_stab_mpc_ - Pvpu_stab_mpc_).row(i)*SUpz_stab_mpc_*SUp_stab_mpc_).transpose();

            calc_index += input_num;

            IS_FIPM_SQP_x_ri1_N_stab_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)  = (Pvps_stab_mpc_.row(i)*ssx_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssz_stab_mpc_)
                                                                                      - (Pvps_stab_mpc_.row(i)*ssz_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssx_stab_mpc_);

            IS_FIPM_SQP_x_ri2_N_stab_mpc_.block(i, 0, 1, 3*N_state)                   = GRAVITY*b_*b_*(Pcps_stab_mpc_.row(i)*ssx_stab_mpc_);
            
            IS_FIPM_SQP_x_ri3_N_stab_mpc_.block(i, 0, 1, 3*N_state)                   = ((Pcps_stab_mpc_ - Pvps_stab_mpc_).row(i)*ssz_stab_mpc_);
            
            IS_FIPM_SQP_x_pi4_N_stab_mpc_.block(calc_index2, 0, 3*N_state, N_stab_mpc)= ((Pvps_stab_mpc_ - Pcps_stab_mpc_).row(i)*ssz_stab_mpc_).transpose()*Si_mpc;

            IS_FIPM_SQP_y_ri1_N_stab_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)  = (Pvps_stab_mpc_.row(i)*ssy_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssz_stab_mpc_)
                                                                                      - (Pvps_stab_mpc_.row(i)*ssz_stab_mpc_).transpose()*(Pcps_stab_mpc_.row(i)*ssy_stab_mpc_);

            IS_FIPM_SQP_y_ri2_N_stab_mpc_.block(i, 0, 1, 3*N_state)                   = GRAVITY*b_*b_*(Pcps_stab_mpc_.row(i)*ssy_stab_mpc_);

            IS_FIPM_SQP_y_ri3_N_stab_mpc_.block(i, 0, 1, 3*N_state)                   = ((Pcps_stab_mpc_ - Pvps_stab_mpc_).row(i)*ssz_stab_mpc_);

            IS_FIPM_SQP_y_pi4_N_stab_mpc_.block(calc_index2, 0, 3*N_state, N_stab_mpc)= ((Pvps_stab_mpc_ - Pcps_stab_mpc_).row(i)*ssz_stab_mpc_).transpose()*Si_mpc;

            calc_index2 += 3*N_state;
        }

        cout << "Initialization of IS 3D DCM MPC is completed" << endl;

        MPC_first_loop_ = 1;

        print_sec_ = 2.0;

        t_total_mpc_ = t_total_const_mpc_;
    }

    Sf1_stab_mpc_.setZero();

    print_sec_bool_ = ((((walking_tick_mpc_ - int(mpc_synchro_hz) + 1)/int(mpc_synchro_hz))%int(print_sec_*thread3_hz_))==0);

    int matlab_tick = ((walking_tick_mpc_ - int(mpc_synchro_hz) + 1)/int(mpc_synchro_hz));

    double step_x_norm = foot_step_support_frame_mpc_(current_step_num_mpc_, 0);
    double step_y_norm = foot_step_support_frame_mpc_(current_step_num_mpc_, 1);

    if(mpc_tick < hz_/mpc_freq)
    {
        MPC_Stabilizer_delf_mpc_x_(0) = step_x_norm;
        MPC_Stabilizer_delf_mpc_y_(0) = step_y_norm;
        MPC_Stabilizer_delf_mpc_ << MPC_Stabilizer_delf_mpc_x_, MPC_Stabilizer_delf_mpc_y_;

        MPC_Stabilizer_alpha_mpc_.setZero();
        MPC_Stabilizer_alpha_mpc_(0) = 1;

        MPC_Stabilizer_aux_mpc_x_ = MPC_Stabilizer_alpha_mpc_*MPC_Stabilizer_delf_mpc_x_;
        MPC_Stabilizer_aux_mpc_y_ = MPC_Stabilizer_alpha_mpc_*MPC_Stabilizer_delf_mpc_y_;
        MPC_Stabilizer_aux_mpc_ << MPC_Stabilizer_aux_mpc_x_, MPC_Stabilizer_aux_mpc_y_;
    }
    
    step_enable_bool_mpc_                 = (bool)((mpc_tick + com_start_tick_mpc_ > t_temp_- (t_dsp2_const_mpc_ + step_enable_fix_time_pre_*hz_))&&(mpc_tick                  < (1 - (bool)(current_step_num_mpc_))*t_temp_ + t_total_const_mpc_ - t_dsp2_const_mpc_ - step_enable_time_fwd_*hz_ - step_enable_fix_time_pre_*hz_));
    step_enable_bool_one_tick_before_mpc_ = (bool)((mpc_tick + com_start_tick_mpc_ > t_temp_- (t_dsp2_const_mpc_ + step_enable_fix_time_pre_*hz_))&&(mpc_tick + mpc_synchro_hz < (1 - (bool)(current_step_num_mpc_))*t_temp_ + t_total_const_mpc_ - t_dsp2_const_mpc_ - step_enable_time_fwd_*hz_ - step_enable_fix_time_pre_*hz_));

    Eigen::VectorXd zmp_max_x_mpc_step_calc; zmp_max_x_mpc_step_calc = zmp_max_x_mpc_;
    Eigen::VectorXd zmp_min_x_mpc_step_calc; zmp_min_x_mpc_step_calc = zmp_min_x_mpc_;

    Eigen::VectorXd zmp_max_y_mpc_step_calc; zmp_max_y_mpc_step_calc = zmp_max_y_mpc_;
    Eigen::VectorXd zmp_min_y_mpc_step_calc; zmp_min_y_mpc_step_calc = zmp_min_y_mpc_;

    for(int i = 0; i < N_stab_mpc; i++)
    {
        for(int j = 0; j < 1; j++)
        { 
            int dsp_length_calc           = int((t_dsp1_const_mpc_ + t_dsp2_const_mpc_)/mpc_synchro_hz + 0.5);
            int next_step_start_prev_tick = max(ceil(((j+1)*((1 - (bool)(current_step_num_mpc_))*t_temp_ + t_total_const_mpc_ - t_dsp2_const_mpc_) - mpc_tick)/mpc_synchro_hz), 0.0);
            bool  next_step_prev_bool     = (bool)(mpc_tick + mpc_synchro_hz*(i + 2) > ((1 - (bool)(current_step_num_mpc_))*t_temp_ + (j + 1)*t_total_mpc_ - t_dsp2_const_mpc_));
            bool nnext_step_prev_bool     = (bool)(mpc_tick + mpc_synchro_hz*(i + 2) > ((1 - (bool)(current_step_num_mpc_))*t_temp_ + (j + 2)*t_total_mpc_ - t_dsp2_const_mpc_));

            if(walking_tick_mpc_ < t_temp_)
            {
                next_step_prev_bool       = (bool)(mpc_tick + mpc_synchro_hz*(i + 1) > ((1 - (bool)(current_step_num_mpc_))*t_temp_ + (j + 1)*t_total_mpc_ - t_dsp2_const_mpc_));
            }

            if((walking_tick_mpc_ > t_temp_ - (t_dsp2_const_mpc_ + step_enable_fix_time_pre_*hz_)) && (next_step_start_prev_tick < N_stab_mpc))
            //if(0)
            {
                Sf1_stab_mpc_(i,0) = step_enable_bool_mpc_* next_step_prev_bool*(zmp_max_y_mpc_(i) - zmp_max_y_mpc_(int(t_dsp1_const_mpc_/mpc_synchro_hz)))/(MPC_Stabilizer_delf_mpc_y_(0));
                if(walking_tick_mpc_ < t_temp_)
                {
                    Sf1_stab_mpc_(i,0) = step_enable_bool_mpc_* next_step_prev_bool*(zmp_max_y_mpc_(i) - zmp_max_y_mpc_(int(t_dsp1_const_mpc_/mpc_synchro_hz) + step_time_adj_candidate_num_))/(MPC_Stabilizer_delf_mpc_y_(0));
                }

                zmp_max_x_mpc_step_calc(i)  = step_enable_bool_mpc_*(next_step_prev_bool* zmp_max_x_mpc_(max(0, next_step_start_prev_tick - 2)) + (1 - next_step_prev_bool)*zmp_max_x_mpc_(i)) + (1 - step_enable_bool_mpc_)*zmp_max_x_mpc_(i);
                zmp_min_x_mpc_step_calc(i)  = step_enable_bool_mpc_*(next_step_prev_bool* zmp_min_x_mpc_(max(0, next_step_start_prev_tick - 2)) + (1 - next_step_prev_bool)*zmp_min_x_mpc_(i)) + (1 - step_enable_bool_mpc_)*zmp_min_x_mpc_(i);

                zmp_max_y_mpc_step_calc(i)  = step_enable_bool_mpc_*(next_step_prev_bool* zmp_max_y_mpc_(max(0, next_step_start_prev_tick - 2)) + (1 - next_step_prev_bool)*zmp_max_y_mpc_(i)) + (1 - step_enable_bool_mpc_)*zmp_max_y_mpc_(i);
                zmp_min_y_mpc_step_calc(i)  = step_enable_bool_mpc_*(next_step_prev_bool* zmp_min_y_mpc_(max(0, next_step_start_prev_tick - 2)) + (1 - next_step_prev_bool)*zmp_min_y_mpc_(i)) + (1 - step_enable_bool_mpc_)*zmp_min_y_mpc_(i);
            }
        }
    }

    for(int i = 1; i < step_time_adj_candidate_num_; i++)
    {
        Sf1_stab_mpc_.col(i) << Sf1_stab_mpc_.col(i-1).segment(1, N_stab_mpc - 1), min((Sf1_stab_mpc_(N_stab_mpc - 1, i-1) + (Sf1_stab_mpc_(N_stab_mpc - 1, i-1) - Sf1_stab_mpc_(N_stab_mpc - 2, i-1))), 1.0);
    }

    Eigen::VectorXd dcm_refx, dcm_refy, dcm_refz;
    dcm_refx.setZero(N_stab_mpc);
    dcm_refy.setZero(N_stab_mpc);
    dcm_refz.setZero(N_stab_mpc);

    dcm_refx = Planner_State_Prev_mpc_.block(0, 0, 1, N_stab_mpc).transpose() + b_*Planner_State_Prev_mpc_.block(1, 0, 1, N_stab_mpc).transpose();
    dcm_refy = Planner_State_Prev_mpc_.block(3, 0, 1, N_stab_mpc).transpose() + b_*Planner_State_Prev_mpc_.block(4, 0, 1, N_stab_mpc).transpose();
    dcm_refz = Planner_State_Prev_mpc_.block(6, 0, 1, N_stab_mpc).transpose() + b_*Planner_State_Prev_mpc_.block(7, 0, 1, N_stab_mpc).transpose();

    MPC_Stabilizer_state_mpc_(0) = com_measured_mpc_(0);
    MPC_Stabilizer_state_mpc_(1) = com_dot_measured_mpc_(0);
    MPC_Stabilizer_state_mpc_(3) = com_measured_mpc_(1);
    MPC_Stabilizer_state_mpc_(4) = com_dot_measured_mpc_(1);
    MPC_Stabilizer_state_mpc_(6) = com_measured_mpc_(2);
    MPC_Stabilizer_state_mpc_(7) = com_dot_measured_mpc_(2);

    gcalc_stab_mpc_ = gxpcalc_stab_mpc_  *(Pdps_stab_mpc_*ssx_stab_mpc_*MPC_Stabilizer_state_mpc_ - dcm_refx)
                     +gypcalc_stab_mpc_  *(Pdps_stab_mpc_*ssy_stab_mpc_*MPC_Stabilizer_state_mpc_ - dcm_refy)
                     +gzpcalc_stab_mpc_  *(Pdps_stab_mpc_*ssz_stab_mpc_*MPC_Stabilizer_state_mpc_ - dcm_refz)
                     
                     +gxdfcalc_stab_mpc_ *(                                                       - MPC_Stabilizer_delf_mpc_x_)
                     +gydfcalc_stab_mpc_ *(                                                       - MPC_Stabilizer_delf_mpc_y_);;

    SQP_deldel_Qcalc_stab_mpc_ = Qcalc_stab_mpc_;

    int sqp_iter = 4;

    const_A_mpc_.setZero(const_num, input_num);
    const_lb_mpc_.setZero(const_num, 1);
    const_ub_mpc_.setZero(const_num, 1);

    Eigen::MatrixXd prev_state;    prev_state.setZero(9,N_stab_mpc);
    Eigen::MatrixXd prev_zmp;      prev_zmp.setZero(2, N_stab_mpc);

    std::vector<Eigen::MatrixXd> IS_FIPM_SQP_x_phi_max_vec; std::vector<Eigen::MatrixXd> IS_FIPM_SQP_x_phi_min_vec;
    std::vector<Eigen::MatrixXd> IS_FIPM_SQP_y_phi_max_vec; std::vector<Eigen::MatrixXd> IS_FIPM_SQP_y_phi_min_vec;
    std::vector<Eigen::MatrixXd> IS_FIPM_SQP_x_pi_max_vec;  std::vector<Eigen::MatrixXd> IS_FIPM_SQP_x_pi_min_vec;
    std::vector<Eigen::MatrixXd> IS_FIPM_SQP_y_pi_max_vec;  std::vector<Eigen::MatrixXd> IS_FIPM_SQP_y_pi_min_vec;
    std::vector<Eigen::MatrixXd> IS_FIPM_SQP_x_ri_max_vec;  std::vector<Eigen::MatrixXd> IS_FIPM_SQP_x_ri_min_vec;
    std::vector<Eigen::MatrixXd> IS_FIPM_SQP_y_ri_max_vec;  std::vector<Eigen::MatrixXd> IS_FIPM_SQP_y_ri_min_vec;

    for(int s = 0; s < sqp_iter; s++)
    {
        if(print_sec_bool_)
            cout << "Stabilizer SQP Iteration: " << s << endl;

        SQP_del_gcalc_stab_mpc_ = Qcalc_stab_mpc_*MPC_Stabilizer_u_mpc_ + gcalc_stab_mpc_;

        QP_MPC_Stabilizer_.EnableEqualityCondition(equality_condition_eps_);
        QP_MPC_Stabilizer_.UpdateMinProblem(SQP_deldel_Qcalc_stab_mpc_, SQP_del_gcalc_stab_mpc_);
        QP_MPC_Stabilizer_.DeleteSubjectToAx();
        QP_MPC_Stabilizer_.DeleteSubjectToX();

        int constraint_index = 0;

        int calc_index = 0;
        int calc_index2 = 0;

        for(int i = 0; i < N_stab_mpc; i++)
        {
            //X max
            if(s == 0)
            {
                const_SQP_phi_mpc_ = IS_FIPM_SQP_x_phi1_N_stab_mpc_.block(calc_index, 0, input_num, input_num)
                
                                   + IS_FIPM_SQP_x_phi2_N_stab_mpc_.block(calc_index, 0, input_num, N_stab_mpc)*Sf1_stab_mpc_*SUfx_stab_mpc_*SUf_stab_mpc_;

                const_SQP_phi_mpc_calc_ = 0.5*(const_SQP_phi_mpc_ + const_SQP_phi_mpc_.transpose());

                IS_FIPM_SQP_x_phi_max_vec.push_back(const_SQP_phi_mpc_calc_);

                const_SQP_pi_mpc_  = IS_FIPM_SQP_x_pi1_N_stab_mpc_.block (calc_index, 0, input_num, 3*N_state)*MPC_Stabilizer_state_mpc_

                                   + IS_FIPM_SQP_x_pi2_N_stab_mpc_.block(calc_index, 0, input_num, 1)

                                   - zmp_max_x_mpc_(i)*IS_FIPM_SQP_x_pi3_N_stab_mpc_.block(calc_index, 0, input_num, 1)
                                   
                                   + (MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_x_pi4_N_stab_mpc_.block(calc_index2, 0, 3*N_state, N_stab_mpc)*Sf1_stab_mpc_*SUfx_stab_mpc_*SUf_stab_mpc_).transpose()
                                   
                                   - (((GRAVITY*b_*b_)*MatrixXd::Identity(1,1)).transpose()*Sf1_stab_mpc_.row(i)*SUfx_stab_mpc_*SUf_stab_mpc_).transpose();
            
                IS_FIPM_SQP_x_pi_max_vec.push_back(const_SQP_pi_mpc_);

                const_SQP_ri_mpc_  = MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_x_ri1_N_stab_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Stabilizer_state_mpc_

                                   + IS_FIPM_SQP_x_ri2_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_

                                   - zmp_max_x_mpc_(i)*IS_FIPM_SQP_x_ri3_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_
                          
                                   - GRAVITY*b_*b_*zmp_max_x_mpc_(i)*MatrixXd::Identity(1,1);
            
                IS_FIPM_SQP_x_ri_max_vec.push_back(const_SQP_ri_mpc_);
            }
            
            const_SQP_hi_mpc_  = MPC_Stabilizer_u_mpc_.transpose()*IS_FIPM_SQP_x_phi_max_vec[i]*MPC_Stabilizer_u_mpc_ 
                               + IS_FIPM_SQP_x_pi_max_vec[i].transpose()*MPC_Stabilizer_u_mpc_ 
                               + IS_FIPM_SQP_x_ri_max_vec[i];

            const_A_mpc_.row(constraint_index) = (2*IS_FIPM_SQP_x_phi_max_vec[i]*MPC_Stabilizer_u_mpc_ + IS_FIPM_SQP_x_pi_max_vec[i]).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - 1e+3*MatrixXd::Identity(1,1);
            constraint_index += 1;

            //X min            
            if(s == 0)
            {
                const_SQP_phi_mpc_ = IS_FIPM_SQP_x_phi1_N_stab_mpc_.block(calc_index, 0, input_num, input_num)

                                   + IS_FIPM_SQP_x_phi2_N_stab_mpc_.block(calc_index, 0, input_num, N_stab_mpc)*Sf1_stab_mpc_*SUfx_stab_mpc_*SUf_stab_mpc_;

                const_SQP_phi_mpc_calc_ = 0.5*(const_SQP_phi_mpc_ + const_SQP_phi_mpc_.transpose());

                IS_FIPM_SQP_x_phi_min_vec.push_back(const_SQP_phi_mpc_calc_);

                const_SQP_pi_mpc_  = IS_FIPM_SQP_x_pi1_N_stab_mpc_.block (calc_index, 0, input_num, 3*N_state)*MPC_Stabilizer_state_mpc_

                               + IS_FIPM_SQP_x_pi2_N_stab_mpc_.block(calc_index, 0, input_num, 1)

                               - zmp_min_x_mpc_(i)*IS_FIPM_SQP_x_pi3_N_stab_mpc_.block(calc_index, 0, input_num, 1)
                               
                               + (MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_x_pi4_N_stab_mpc_.block(calc_index2, 0, 3*N_state, N_stab_mpc)*Sf1_stab_mpc_*SUfx_stab_mpc_*SUf_stab_mpc_).transpose()
                                   
                               - (((GRAVITY*b_*b_)*MatrixXd::Identity(1,1)).transpose()*Sf1_stab_mpc_.row(i)*SUfx_stab_mpc_*SUf_stab_mpc_).transpose();

                IS_FIPM_SQP_x_pi_min_vec.push_back(const_SQP_pi_mpc_);

                const_SQP_ri_mpc_ = MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_x_ri1_N_stab_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Stabilizer_state_mpc_

                              + IS_FIPM_SQP_x_ri2_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_

                              - zmp_min_x_mpc_(i)*IS_FIPM_SQP_x_ri3_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_
                          
                              - GRAVITY*b_*b_*zmp_min_x_mpc_(i)*MatrixXd::Identity(1,1);
            
                IS_FIPM_SQP_x_ri_min_vec.push_back(const_SQP_ri_mpc_);
            }

            const_SQP_hi_mpc_ = MPC_Stabilizer_u_mpc_.transpose()*IS_FIPM_SQP_x_phi_min_vec[i]*MPC_Stabilizer_u_mpc_ 
                              + IS_FIPM_SQP_x_pi_min_vec[i].transpose()*MPC_Stabilizer_u_mpc_ 
                              + IS_FIPM_SQP_x_ri_min_vec[i];

            const_A_mpc_.row(constraint_index) = (2*IS_FIPM_SQP_x_phi_min_vec[i]*MPC_Stabilizer_u_mpc_ + IS_FIPM_SQP_x_pi_min_vec[i]).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) =   1e+3*MatrixXd::Identity(1,1);
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            constraint_index += 1;

            //Y max
            if(s == 0)
            {
                const_SQP_phi_mpc_ = IS_FIPM_SQP_y_phi1_N_stab_mpc_.block(calc_index, 0, input_num, input_num)
                
                                   + IS_FIPM_SQP_y_phi2_N_stab_mpc_.block(calc_index, 0, input_num, N_stab_mpc)*Sf1_stab_mpc_*SUfy_stab_mpc_*SUf_stab_mpc_;

                const_SQP_phi_mpc_calc_ = 0.5*(const_SQP_phi_mpc_ + const_SQP_phi_mpc_.transpose());

                IS_FIPM_SQP_y_phi_max_vec.push_back(const_SQP_phi_mpc_calc_);

                const_SQP_pi_mpc_  = IS_FIPM_SQP_y_pi1_N_stab_mpc_.block(calc_index,  0, input_num, 3*N_state)*MPC_Stabilizer_state_mpc_

                                   + IS_FIPM_SQP_y_pi2_N_stab_mpc_.block(calc_index, 0, input_num, 1)

                                   - zmp_max_y_mpc_step_calc(i)*IS_FIPM_SQP_y_pi3_N_stab_mpc_.block(calc_index, 0, input_num, 1)
                                   
                                   + (MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_y_pi4_N_stab_mpc_.block(calc_index2, 0, 3*N_state, N_stab_mpc)*Sf1_stab_mpc_*SUfy_stab_mpc_*SUf_stab_mpc_).transpose()
                                   
                                   - (((GRAVITY*b_*b_)*MatrixXd::Identity(1,1)).transpose()*Sf1_stab_mpc_.row(i)*SUfy_stab_mpc_*SUf_stab_mpc_).transpose();

                IS_FIPM_SQP_y_pi_max_vec.push_back(const_SQP_pi_mpc_);

                const_SQP_ri_mpc_ = MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_y_ri1_N_stab_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Stabilizer_state_mpc_

                              + IS_FIPM_SQP_y_ri2_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_

                              - zmp_max_y_mpc_step_calc(i)*IS_FIPM_SQP_y_ri3_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_
                          
                              - GRAVITY*b_*b_*zmp_max_y_mpc_step_calc(i)*MatrixXd::Identity(1,1);
                              
                IS_FIPM_SQP_y_ri_max_vec.push_back(const_SQP_ri_mpc_);
            }

            const_SQP_hi_mpc_ = MPC_Stabilizer_u_mpc_.transpose()*IS_FIPM_SQP_y_phi_max_vec[i]*MPC_Stabilizer_u_mpc_ 
                              + IS_FIPM_SQP_y_pi_max_vec[i].transpose()*MPC_Stabilizer_u_mpc_ 
                              + IS_FIPM_SQP_y_ri_max_vec[i];

            const_A_mpc_.row(constraint_index) = (2*IS_FIPM_SQP_y_phi_max_vec[i]*MPC_Stabilizer_u_mpc_ + IS_FIPM_SQP_y_pi_max_vec[i]).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - 1e+3*MatrixXd::Identity(1,1);
            constraint_index += 1;

            //Y min
            if(s == 0)
            {
                const_SQP_phi_mpc_ = IS_FIPM_SQP_y_phi1_N_stab_mpc_.block(calc_index, 0, input_num, input_num)
                
                                   + IS_FIPM_SQP_y_phi2_N_stab_mpc_.block(calc_index, 0, input_num, N_stab_mpc)*Sf1_stab_mpc_*SUfy_stab_mpc_*SUf_stab_mpc_;
                
                const_SQP_phi_mpc_calc_ = 0.5*(const_SQP_phi_mpc_ + const_SQP_phi_mpc_.transpose());

                IS_FIPM_SQP_y_phi_min_vec.push_back(const_SQP_phi_mpc_calc_);

                const_SQP_pi_mpc_  = IS_FIPM_SQP_y_pi1_N_stab_mpc_.block (calc_index, 0, input_num, 3*N_state)*MPC_Stabilizer_state_mpc_

                                   + IS_FIPM_SQP_y_pi2_N_stab_mpc_.block(calc_index, 0, input_num, 1)

                                   - zmp_min_y_mpc_step_calc(i)*IS_FIPM_SQP_y_pi3_N_stab_mpc_.block(calc_index, 0, input_num, 1)
                                   
                                   + (MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_y_pi4_N_stab_mpc_.block(calc_index2, 0, 3*N_state, N_stab_mpc)*Sf1_stab_mpc_*SUfy_stab_mpc_*SUf_stab_mpc_).transpose()
                                   
                                   - (((GRAVITY*b_*b_)*MatrixXd::Identity(1,1)).transpose()*Sf1_stab_mpc_.row(i)*SUfy_stab_mpc_*SUf_stab_mpc_).transpose();
                               
                IS_FIPM_SQP_y_pi_min_vec.push_back(const_SQP_pi_mpc_);

                const_SQP_ri_mpc_ = MPC_Stabilizer_state_mpc_.transpose()*IS_FIPM_SQP_y_ri1_N_stab_mpc_.block(calc_index2, 0, 3*N_state, 3*N_state)*MPC_Stabilizer_state_mpc_

                              + IS_FIPM_SQP_y_ri2_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_

                              - zmp_min_y_mpc_step_calc(i)*IS_FIPM_SQP_y_ri3_N_stab_mpc_.block(i, 0, 1, 3*N_state)*MPC_Stabilizer_state_mpc_
                          
                              - GRAVITY*b_*b_*zmp_min_y_mpc_step_calc(i)*MatrixXd::Identity(1,1);
                
                IS_FIPM_SQP_y_ri_min_vec.push_back(const_SQP_ri_mpc_);
            }

            const_SQP_hi_mpc_ = MPC_Stabilizer_u_mpc_.transpose()*IS_FIPM_SQP_y_phi_min_vec[i]*MPC_Stabilizer_u_mpc_ 
                              + IS_FIPM_SQP_y_pi_min_vec[i].transpose()*MPC_Stabilizer_u_mpc_ 
                              + IS_FIPM_SQP_y_ri_min_vec[i];

            const_A_mpc_.row(constraint_index) = (2*IS_FIPM_SQP_y_phi_min_vec[i]*MPC_Stabilizer_u_mpc_ + IS_FIPM_SQP_y_pi_min_vec[i]).transpose();
            const_ub_mpc_.block(constraint_index, 0, 1, 1) =   1e+3*MatrixXd::Identity(1,1);
            const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
            constraint_index += 1;
        
            calc_index  += input_num;
            calc_index2 += 3*N_state;
        }

        //IS equality
        Eigen::MatrixXd b_IS_stab_mpc; b_IS_stab_mpc.resize(N_stab_mpc,1); b_IS_stab_mpc.col(0) = b_IS_plan_mpc_.col(0).segment(0,N_stab_mpc);
        //Pre planned Tail
        Eigen::MatrixXd Const_b_eq; Const_b_eq.setZero(3,1);
        Const_b_eq(0,0) = -(w_/(1 - lambda_is_calc))*(MPC_Stabilizer_state_mpc_(0) + MPC_Stabilizer_state_mpc_(1)/w_ - MPC_Stabilizer_state_mpc_(2));
        Const_b_eq(1,0) = -(w_/(1 - lambda_is_calc))*(MPC_Stabilizer_state_mpc_(3) + MPC_Stabilizer_state_mpc_(4)/w_ - MPC_Stabilizer_state_mpc_(5))
                          +(pow(lambda_is_calc, N_stab_mpc)/(1 + pow(lambda_is_calc, N_step))*(b_IS_step_mpc_.transpose()*Pv_dot_ref_mpc_.col(1))(0,0));
        Const_b_eq(2,0) = -(w_/(1 - lambda_is_calc))*(MPC_Stabilizer_state_mpc_(6) + MPC_Stabilizer_state_mpc_(7)/w_ - MPC_Stabilizer_state_mpc_(8));

        const_SQP_phi_mpc_.setZero(input_num, input_num);
        //const_SQP_pi_mpc_ = (b_IS_stab_mpc.transpose()*SUpx_stab_mpc_*SUp_stab_mpc_).transpose();
        const_SQP_pi_mpc_ = SUpxp_stab_mpc_.transpose()*b_IS_stab_mpc;
        const_SQP_ri_mpc_ = Const_b_eq.block(0, 0, 1, 1);

        //const_SQP_hi_mpc_ = MPC_Stabilizer_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Stabilizer_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Stabilizer_u_mpc_ + const_SQP_ri_mpc_;
        const_SQP_hi_mpc_ = const_SQP_pi_mpc_.transpose()*MPC_Stabilizer_u_mpc_ + const_SQP_ri_mpc_;

        //const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Stabilizer_u_mpc_ + const_SQP_pi_mpc_).transpose();
        const_A_mpc_.row(constraint_index) = const_SQP_pi_mpc_.transpose();
        const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        constraint_index += 1;

        const_SQP_phi_mpc_.setZero(input_num, input_num);
        //const_SQP_pi_mpc_ = (b_IS_stab_mpc.transpose()*SUpy_stab_mpc_*SUp_stab_mpc_).transpose();
        const_SQP_pi_mpc_ = SUpyp_stab_mpc_.transpose()*b_IS_stab_mpc;
        const_SQP_ri_mpc_ = Const_b_eq.block(1, 0, 1, 1);
        //const_SQP_hi_mpc_ = MPC_Stabilizer_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Stabilizer_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Stabilizer_u_mpc_ + const_SQP_ri_mpc_;
        const_SQP_hi_mpc_ = const_SQP_pi_mpc_.transpose()*MPC_Stabilizer_u_mpc_ + const_SQP_ri_mpc_;

        //const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Stabilizer_u_mpc_ + const_SQP_pi_mpc_).transpose();
        const_A_mpc_.row(constraint_index) = const_SQP_pi_mpc_.transpose();
        const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        constraint_index += 1;

        const_SQP_phi_mpc_.setZero(input_num, input_num);
        //const_SQP_pi_mpc_ = (b_IS_stab_mpc.transpose()*SUpz_stab_mpc_*SUp_stab_mpc_).transpose();
        const_SQP_pi_mpc_ = SUpzp_stab_mpc_.transpose()*b_IS_stab_mpc;
        const_SQP_ri_mpc_ = Const_b_eq.block(2, 0, 1, 1);
        //const_SQP_hi_mpc_ = MPC_Stabilizer_u_mpc_.transpose()*const_SQP_phi_mpc_*MPC_Stabilizer_u_mpc_ + const_SQP_pi_mpc_.transpose()*MPC_Stabilizer_u_mpc_ + const_SQP_ri_mpc_;
        const_SQP_hi_mpc_ = const_SQP_pi_mpc_.transpose()*MPC_Stabilizer_u_mpc_ + const_SQP_ri_mpc_;

        //const_A_mpc_.row(constraint_index) = (2*const_SQP_phi_mpc_*MPC_Stabilizer_u_mpc_ + const_SQP_pi_mpc_).transpose();
        const_A_mpc_.row(constraint_index) = const_SQP_pi_mpc_.transpose();
        const_ub_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        const_lb_mpc_.block(constraint_index, 0, 1, 1) = - const_SQP_hi_mpc_;
        constraint_index += 1;

        //delf min max
        double delf_x_max = 0.15, delf_x_min = -0.15;
        double delf_y_max = 0.15, delf_y_min =  0.00;

        double delf_x_max_calc, delf_x_min_calc;
        double delf_y_max_calc, delf_y_min_calc;

        //delf min max
        const_A_mpc_.block (constraint_index, 0, 1, input_num) =   SUfx_stab_mpc_*SUf_stab_mpc_;
        const_ub_mpc_.block(constraint_index, 0, 1, 1)         = - SUfx_stab_mpc_*SUf_stab_mpc_*MPC_Stabilizer_u_mpc_ + delf_x_max*MatrixXd::Ones(1,1);
        const_lb_mpc_.block(constraint_index, 0, 1, 1)         = - SUfx_stab_mpc_*SUf_stab_mpc_*MPC_Stabilizer_u_mpc_ + delf_x_min*MatrixXd::Ones(1,1);
        constraint_index += 1;

        delf_y_max_calc = step_y_norm + (foot_step_mpc_(current_step_num_mpc_, 6)*delf_y_min + (1 - foot_step_mpc_(current_step_num_mpc_, 6))*delf_y_max);
        delf_y_min_calc = step_y_norm - (foot_step_mpc_(current_step_num_mpc_, 6)*delf_y_max + (1 - foot_step_mpc_(current_step_num_mpc_, 6))*delf_y_min);

        const_A_mpc_.block (constraint_index, 0, 1, input_num) =   SUfy_stab_mpc_*SUf_stab_mpc_;
        const_ub_mpc_.block(constraint_index, 0, 1, 1)         = - SUfy_stab_mpc_*SUf_stab_mpc_*MPC_Stabilizer_u_mpc_ + delf_y_max_calc*MatrixXd::Ones(1,1);
        const_lb_mpc_.block(constraint_index, 0, 1, 1)         = - SUfy_stab_mpc_*SUf_stab_mpc_*MPC_Stabilizer_u_mpc_ + delf_y_min_calc*MatrixXd::Ones(1,1);
        constraint_index += 1;
        
        QP_MPC_Stabilizer_.UpdateSubjectToAx(const_A_mpc_, const_lb_mpc_, const_ub_mpc_);

        if(QP_MPC_Stabilizer_.SolveQPoases(100, MPC_Stabilizer_SQP_du_mpc_))
        {
            MPC_Stabilizer_u_mpc_ = MPC_Stabilizer_u_mpc_ + MPC_Stabilizer_SQP_du_mpc_;

            if(print_sec_bool_)
            {
                cout << "IS DCM Stabilizer MPC Solved" << endl;
                if(s == sqp_iter - 1)
                    cout << endl;
            }

            if(s == sqp_iter - 1)
            {
                prev_state.row(0).transpose() = Pcps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + Pcpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(0*N_stab_mpc, N_stab_mpc);
                prev_state.row(1).transpose() = Pcvs_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + Pcvu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(0*N_stab_mpc, N_stab_mpc);
                prev_state.row(2).transpose() = Pvps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + Pvpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(0*N_stab_mpc, N_stab_mpc);

                prev_state.row(3).transpose() = Pcps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + Pcpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(1*N_stab_mpc, N_stab_mpc);
                prev_state.row(4).transpose() = Pcvs_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + Pcvu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(1*N_stab_mpc, N_stab_mpc);
                prev_state.row(5).transpose() = Pvps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + Pvpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(1*N_stab_mpc, N_stab_mpc);

                prev_state.row(6).transpose() = Pcps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + Pcpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(2*N_stab_mpc, N_stab_mpc);
                prev_state.row(7).transpose() = Pcvs_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + Pcvu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(2*N_stab_mpc, N_stab_mpc);
                prev_state.row(8).transpose() = Pvps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + Pvpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(2*N_stab_mpc, N_stab_mpc);

                Eigen::VectorXd lambda_bb_mat, lambda_bb_mat_calc;
                lambda_bb_mat_calc = (prev_state.row(8) - GRAVITY*b_*b_*MatrixXd::Ones(1, N_stab_mpc)).array()/prev_state.row(6).array();
                lambda_bb_mat      = MatrixXd::Ones(N_stab_mpc, 1) - lambda_bb_mat_calc;

                prev_zmp.row(0)    = (prev_state.row(2).array() - (MatrixXd::Ones(1, N_stab_mpc) - lambda_bb_mat.transpose()).array()*prev_state.row(0).array()).array()/lambda_bb_mat.transpose().array();
                prev_zmp.row(1)    = (prev_state.row(5).array() - (MatrixXd::Ones(1, N_stab_mpc) - lambda_bb_mat.transpose()).array()*prev_state.row(3).array()).array()/lambda_bb_mat.transpose().array();

                MPC_Stabilizer_u_mpc_sep_(0) = MPC_Stabilizer_u_mpc_(0*N_stab_mpc);
                MPC_Stabilizer_u_mpc_sep_(1) = MPC_Stabilizer_u_mpc_(1*N_stab_mpc);
                MPC_Stabilizer_u_mpc_sep_(2) = MPC_Stabilizer_u_mpc_(2*N_stab_mpc);

                MPC_Stabilizer_state_mpc_.segment(0,3) = A_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + B_mpc_*MPC_Stabilizer_u_mpc_(0*N_stab_mpc);
                MPC_Stabilizer_state_mpc_.segment(3,3) = A_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + B_mpc_*MPC_Stabilizer_u_mpc_(1*N_stab_mpc);
                MPC_Stabilizer_state_mpc_.segment(6,3) = A_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + B_mpc_*MPC_Stabilizer_u_mpc_(2*N_stab_mpc);

                MPC_Stabilizer_delf_mpc_ = MPC_Stabilizer_u_mpc_.segment(3*N_stab_mpc, 2);
                MPC_Stabilizer_delf_mpc_x_ = MPC_Stabilizer_delf_mpc_.segment(0, 1);
                MPC_Stabilizer_delf_mpc_y_ = MPC_Stabilizer_delf_mpc_.segment(1, 1);

                MPC_Stabilizer_alpha_mpc_  = MPC_Stabilizer_u_mpc_.segment(3*N_stab_mpc + 2, 1*step_time_adj_candidate_num_);

                MPC_Stabilizer_aux_mpc_ = MPC_Stabilizer_u_mpc_.segment(3*N_stab_mpc + 2 + 1*step_time_adj_candidate_num_, 2*step_time_adj_candidate_num_);
                MPC_Stabilizer_aux_mpc_x_ = MPC_Stabilizer_aux_mpc_.segment(0*step_time_adj_candidate_num_, 1*step_time_adj_candidate_num_);
                MPC_Stabilizer_aux_mpc_y_ = MPC_Stabilizer_aux_mpc_.segment(1*step_time_adj_candidate_num_, 1*step_time_adj_candidate_num_);
            }
        }
        else
        {
            cout << "IS FIPM DCM Stabilizer Stepping Not Solved" << endl;
            cout << matlab_tick << endl;

            MPC_Stabilizer_u_mpc_ = MPC_Stabilizer_u_mpc_ + MPC_Stabilizer_SQP_du_mpc_;

            if(s == sqp_iter - 1)
            {
                prev_state.row(0).transpose() = Pcps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + Pcpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(0*N_stab_mpc, N_stab_mpc);
                prev_state.row(1).transpose() = Pcvs_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + Pcvu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(0*N_stab_mpc, N_stab_mpc);
                prev_state.row(2).transpose() = Pvps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + Pvpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(0*N_stab_mpc, N_stab_mpc);

                prev_state.row(3).transpose() = Pcps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + Pcpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(1*N_stab_mpc, N_stab_mpc);
                prev_state.row(4).transpose() = Pcvs_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + Pcvu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(1*N_stab_mpc, N_stab_mpc);
                prev_state.row(5).transpose() = Pvps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + Pvpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(1*N_stab_mpc, N_stab_mpc);

                prev_state.row(6).transpose() = Pcps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + Pcpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(2*N_stab_mpc, N_stab_mpc);
                prev_state.row(7).transpose() = Pcvs_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + Pcvu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(2*N_stab_mpc, N_stab_mpc);
                prev_state.row(8).transpose() = Pvps_stab_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + Pvpu_stab_mpc_*MPC_Stabilizer_u_mpc_.segment(2*N_stab_mpc, N_stab_mpc);

                Eigen::VectorXd lambda_bb_mat, lambda_bb_mat_calc;
                lambda_bb_mat_calc = (prev_state.row(8) - GRAVITY*b_*b_*MatrixXd::Ones(1, N_stab_mpc)).array()/prev_state.row(6).array();
                lambda_bb_mat      = MatrixXd::Ones(N_stab_mpc, 1) - lambda_bb_mat_calc;

                prev_zmp.row(0)    = (prev_state.row(2).array() - (MatrixXd::Ones(1, N_stab_mpc) - lambda_bb_mat.transpose()).array()*prev_state.row(0).array()).array()/lambda_bb_mat.transpose().array();
                prev_zmp.row(1)    = (prev_state.row(5).array() - (MatrixXd::Ones(1, N_stab_mpc) - lambda_bb_mat.transpose()).array()*prev_state.row(3).array()).array()/lambda_bb_mat.transpose().array();

                MPC_Stabilizer_u_mpc_sep_(0) = MPC_Stabilizer_u_mpc_(0*N_stab_mpc);
                MPC_Stabilizer_u_mpc_sep_(1) = MPC_Stabilizer_u_mpc_(1*N_stab_mpc);
                MPC_Stabilizer_u_mpc_sep_(2) = MPC_Stabilizer_u_mpc_(2*N_stab_mpc);

                MPC_Stabilizer_state_mpc_.segment(0,3) = A_mpc_*MPC_Stabilizer_state_mpc_.segment(0,3) + B_mpc_*MPC_Stabilizer_u_mpc_(0*N_stab_mpc);
                MPC_Stabilizer_state_mpc_.segment(3,3) = A_mpc_*MPC_Stabilizer_state_mpc_.segment(3,3) + B_mpc_*MPC_Stabilizer_u_mpc_(1*N_stab_mpc);
                MPC_Stabilizer_state_mpc_.segment(6,3) = A_mpc_*MPC_Stabilizer_state_mpc_.segment(6,3) + B_mpc_*MPC_Stabilizer_u_mpc_(2*N_stab_mpc);

                MPC_Stabilizer_delf_mpc_ = MPC_Stabilizer_u_mpc_.segment(3*N_stab_mpc, 2);
                MPC_Stabilizer_delf_mpc_x_ = MPC_Stabilizer_delf_mpc_.segment(0, 1);
                MPC_Stabilizer_delf_mpc_y_ = MPC_Stabilizer_delf_mpc_.segment(1, 1);

                MPC_Stabilizer_alpha_mpc_  = MPC_Stabilizer_u_mpc_.segment(3*N_stab_mpc + 2, 1*step_time_adj_candidate_num_);

                MPC_Stabilizer_aux_mpc_ = MPC_Stabilizer_u_mpc_.segment(3*N_stab_mpc + 2 + 1*step_time_adj_candidate_num_, 2*step_time_adj_candidate_num_);
                MPC_Stabilizer_aux_mpc_x_ = MPC_Stabilizer_aux_mpc_.segment(0*step_time_adj_candidate_num_, 1*step_time_adj_candidate_num_);
                MPC_Stabilizer_aux_mpc_y_ = MPC_Stabilizer_aux_mpc_.segment(1*step_time_adj_candidate_num_, 1*step_time_adj_candidate_num_);
            }
        }
    }

    if(e_stabilizer_data_txt_[0].is_open())
    {
        e_stabilizer_data_txt_[0] << walking_tick_mpc_            << "," << N_stab_mpc                   << "," << b_                                    << ","
                                  << MPC_Stabilizer_state_mpc_(0) << "," << MPC_Stabilizer_state_mpc_(3) << "," << MPC_Stabilizer_state_mpc_(6)          << ","
                                  << MPC_Stabilizer_state_mpc_(1) << "," << MPC_Stabilizer_state_mpc_(4) << "," << MPC_Stabilizer_state_mpc_(7)          << ","
                                  << MPC_Stabilizer_state_mpc_(2) << "," << MPC_Stabilizer_state_mpc_(5) << "," << MPC_Stabilizer_state_mpc_(8)          << ","
                                  << dcm_refx(0)                  << "," << dcm_refy(0)                  << "," << dcm_refz(0)                           << ","
                                  << com_measured_mpc_(0)         << "," << com_measured_mpc_(1)         << "," << com_measured_mpc_(2)                  << ","
                                  << com_dot_measured_mpc_(0)     << "," << com_dot_measured_mpc_(1)     << "," << com_dot_measured_mpc_(2)              << ","
                                  << step_time_adj_candidate_num_ << "," << step_enable_bool_mpc_        << "," << step_enable_bool_one_tick_before_mpc_ << ","
                                  << MPC_Stabilizer_delf_mpc_x_   << "," << MPC_Stabilizer_delf_mpc_y_   << "," << 0                                     << ","
                                  << step_x_norm                  << "," << step_y_norm                  << "," << 0
                                  << endl;
    }
    else
    {
        cout << "Error: Unable to write to e_stabilizer_data_txt_[0]. File is not open." << endl;
    }

    Eigen::VectorXd data_save_calc;
    data_save_calc.setZero(3*N_stab_mpc);
    data_save_calc << dcm_refx, dcm_refy, dcm_refz;
    e_stabilizer_data_txt_[1] << data_save_calc.transpose() << endl;
    data_save_calc << prev_state.row(0).transpose(), prev_state.row(3).transpose(), Planner_State_Prev_mpc_.row(6).transpose();
    e_stabilizer_data_txt_[2] << data_save_calc.transpose() << endl;
    data_save_calc << prev_state.row(2).transpose(), prev_state.row(5).transpose(), Planner_State_Prev_mpc_.row(8).transpose();
    e_stabilizer_data_txt_[3] << data_save_calc.transpose() << endl;

    Eigen::Map<Eigen::MatrixXd> Sf1_stab_mpc_row(Sf1_stab_mpc_.data(), step_time_adj_candidate_num_*N_stab_mpc, 1);
    e_stabilizer_data_txt_[4] << Sf1_stab_mpc_row.transpose() << endl;

    data_save_calc << zmp_min_x_mpc_step_calc, zmp_max_x_mpc_step_calc, 0*zmp_max_x_mpc_step_calc;
    e_stabilizer_data_txt_[5] << data_save_calc.transpose() << endl;
    data_save_calc << zmp_min_y_mpc_step_calc, zmp_max_y_mpc_step_calc, 0*zmp_max_y_mpc_step_calc;
    e_stabilizer_data_txt_[6] << data_save_calc.transpose() << endl;

    data_save_calc.setZero(1*step_time_adj_candidate_num_);
    data_save_calc << MPC_Stabilizer_alpha_mpc_;
    e_stabilizer_data_txt_[7] << data_save_calc.transpose() << endl;

    data_save_calc.setZero(2*step_time_adj_candidate_num_);
    data_save_calc << MPC_Stabilizer_aux_mpc_x_, MPC_Stabilizer_aux_mpc_y_;
    e_stabilizer_data_txt_[8] << data_save_calc.transpose() << endl;
}