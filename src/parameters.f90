module parameters
    integer, parameter :: N_total = 32
    integer, parameter :: N_start = 2
    integer, parameter :: N_on_blade = N_total - N_start + 1
    integer, parameter :: N = N_on_blade - 2    
    integer, parameter :: T = 5e5
    integer, parameter :: t_count = 5e2
    integer, parameter :: simulation_type = 1
    integer, parameter :: blade_num = 2
    real, parameter :: radius = 5
    real, parameter :: omega = 72/60*2*3.1415926
    real, parameter :: pitch = 3.0/180*3.1415926
    real, parameter :: total_time = 10
    real, parameter :: r_start = 0.508
    real :: delta_t = total_time/T
    real :: delta_r = radius/(N_on_blade-1)
end module parameters