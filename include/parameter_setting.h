#pragma once
#include "wholebody_functions.h"

class parameterSettingConstructor
{
    public:
        parameterSettingConstructor();

        double hz_ = 2000.0;
        double target_x_ = 10.0;
        double target_y_ = 0.0;
        double target_z_ = 0.0;
        double com_height_ = 0.71;
        double target_theta_ = 0.0;
        double step_length_x_ = 0.25;
        double step_length_y_ = 0.0;
        int    is_right_foot_swing_ = 1;
    
        double t_dsp1_        = 0.10 * hz_;
        double t_dsp2_        = 0.10 * hz_;
        double t_total_       = 0.8 * hz_;

        double t_dsp1_const_  = 0.10 * hz_;
        double t_dsp2_const_  = 0.10 * hz_;
        double t_total_const_ = 0.8 * hz_;

        double t_ssp_ = t_total_ - t_dsp1_ - t_dsp2_;
        double foot_height_ = 0.055;

        double t_temp_ = 3.0 * hz_;

        double lfoot_zmp_offset_ = 0.050;
        double rfoot_zmp_offset_ = 0.050;
        double zmp_offset_first_ = 0.015;

        //double pelv_rot_deg_ = 5.0;
        double pelv_rot_deg_ = step_length_x_*100.0/3.0;

};