module structural
    use parameters
    implicit none
    contains

    subroutine bladeConfig(EI,r_location,rho)
        implicit none
        real,dimension(N_on_blade,2,2),intent(out) :: EI
        real,dimension(N_on_blade),intent(out) :: r_location
        real,dimension(N_on_blade) :: rho,twist,EI_f,EI_e
        integer :: i

        OPEN(1, FILE='input/density.csv', FORM='formatted')
        OPEN(2, FILE='input/twist_angle.csv', FORM='formatted')
        OPEN(3, FILE='input/flapwise_stiffness.csv', FORM='formatted')
        OPEN(4, FILE='input/edgewise_stiffness.csv', FORM='formatted')
        OPEN(5, FILE='input/alm_node.csv', FORM='formatted')
    
        Do i = 1, N_total
            if (i >= N_start) then
                read(1, '(E11.4)') rho(i-N_start+1)
                read(2, '(E11.4)') twist(i-N_start+1)
                read(3, '(E11.4)') EI_f(i-N_start+1)
                read(4, '(E11.4)') EI_e(i-N_start+1)
                read(5, '(E11.4)') r_location(i-N_start+1)
            else
                read(1, '(E11.4)')
                read(2, '(E11.4)')
                read(3, '(E11.4)')
                read(4, '(E11.4)')
                read(5, '(E11.4)')
            end if
        End do

        do i = 2, N+1
            EI(i,1,1) = EI_e(i) - (EI_e(i)-EI_f(i))*(cos(twist(i)+pitch)**2)
            EI(i,2,2) = EI_f(i) + (EI_e(i)-EI_f(i))*(cos(twist(i)+pitch)**2)
            EI(i,1,2) = sin(2*(twist(i)+pitch))*((EI_e(i)-EI_f(i))/2)
            EI(i,2,1) = sin(2*(twist(i)+pitch))*((EI_e(i)-EI_f(i))/2)
        end do

        close(1)
        close(2)
        close(3)
        close(4)
    
    end subroutine bladeConfig

    subroutine bladeInit(V_new,M_new,S_new,q_new,V_old,M_old,S_old,q_old)
        implicit none
        real, dimension(N_on_blade,2),intent(out) :: V_new, S_new, M_new, q_new
        real, dimension(N_on_blade,2),intent(out) :: V_old, S_old, M_old, q_old

        integer :: i, j
        ! Get all arrays with initial value 0:
        V_new = 0.0
        V_old = 0.0
        S_new = 0.0
        S_old = 0.0
        M_new = 0.0
        M_old = 0.0
        q_new = 0.0
        q_old = 0.0

        ! For dynamic simulation, load the static simulation result as the initial condition
        if (simulation_type == 1) then
            open(1, file='input/V_init.csv', FORM='formatted')
            open(2, file='input/S_init.csv', FORM='formatted')
            open(3, file='input/M_init.csv', FORM='formatted')
            open(4, file='input/q_init.csv', FORM='formatted')

            do i = 1, N_on_blade
                do j = 1, 2
                    read(1, '(E11.4)') V_old(i,j)
                    read(2, '(E11.4)') S_old(i,j)
                    read(3, '(E11.4)') M_old(i,j)
                    read(4, '(E11.4)') q_old(i,j)
                end do
            end do

            close(1)
            close(2)
            close(3)
            close(4)
        end if
    end subroutine

    subroutine bladeOutput(q_result,M_result,V_new,M_new,S_new,q_new)
        implicit none
        real, dimension(T/t_count+1,N_on_blade,2),intent(in) :: q_result,M_result
        real, dimension(N_on_blade,2),intent(in) :: V_new, S_new, M_new, q_new
    
        integer :: i, j
        open(1, FILE='output/q_0_not_yaw.csv', FORM='formatted')
        open(2, FILE='output/q_1_not_yaw.csv', FORM='formatted')
        open(3, FILE='output/M_0_not_yaw.csv', FORM='formatted')
        open(4, FILE='output/M_1_not_yaw.csv', FORM='formatted')
        Do i = 1, T/t_count+1
            do j = 2,N+1
                write(1, '(E11.4)') q_result(i,j,1)
                write(2, '(E11.4)') q_result(i,j,2)
                write(3, '(E11.4)') M_result(i,j,1)
                write(4, '(E11.4)') M_result(i,j,2)
            end do
        End do
    
        if (simulation_type == 0) then
            open(1, file='input/V_init.csv', FORM='formatted')
            open(2, file='input/S_init.csv', FORM='formatted')
            open(3, file='input/M_init.csv', FORM='formatted')
            open(4, file='input/q_init.csv', FORM='formatted')

            do i = 1,N_on_blade
                write(1, '(E11.4)') V_new(i,:)
                write(2, '(E11.4)') S_new(i,:)
                write(3, '(E11.4)') M_new(i,:)
                write(4, '(E11.4)') q_new(i,:)
            end do
        end if
    end subroutine bladeOutput

    subroutine bladeSolve(V_new,M_new,S_new,q_new,V_old,M_old,S_old,q_old,&
                            rho,EI,r_location,g,F_a,delta_t,delta_r)
        implicit none
        real, dimension(N_on_blade,2),intent(out) :: V_new, S_new, M_new, q_new
        real, dimension(N_on_blade,2),intent(in)  :: V_old, S_old, M_old, q_old
        real,dimension(N_on_blade,2,2),intent(in) :: EI
        real,dimension(N_on_blade),intent(in) :: rho,r_location
        real,intent(in) :: F_a(2), g(2)
        real,intent(in) :: delta_t,delta_r
        real :: damping
        real :: centrifugal_force

        integer :: i, j
        
        ! Set damping factor
        if (simulation_type == 0) then
            damping = 10.0
        else
            damping = 0.1
        end if

        ! Calculation of the velocity new values from the old step time:
        do i = 2, N+1
            ! Interpolate the aerodynamic force
            V_new(i,:) = (1-damping*delta_t)*V_old(i,:) + &
                        delta_t*((-1.0/rho(i))*((M_old(i+1,:)-(2.0*M_old(i,:))+M_old(i-1,:))/delta_r**2)+&
                        ((1.0/rho(i))*((S_old(i+1,:)-S_old(i,:))/delta_r))+(F_a/rho(i))+g)
        end do
    
        ! Boundary conditions at the ROOT position:
        q_new(1:2,:) = 0.0
        V_new(1:2,:) = 0.0

        do i = 2, N+1
            ! Calculation of the bending moment from the new velocity values:
            M_new(i,1) = M_old(i,1) + delta_t * (EI(i,1,1)*((V_new(i+1,1)-(2.0*V_new(i,1))+V_new(i-1,1))/delta_r**2) + &
                                                EI(i,1,2)*((V_new(i+1,2)-(2.0*V_new(i,2))+V_new(i-1,2))/delta_r**2))
            M_new(i,2) = M_old(i,2) + delta_t * (EI(i,2,2)*((V_new(i+1,2)-(2.0*V_new(i,2))+V_new(i-1,2))/delta_r**2) + &
                                                EI(i,2,1)*((V_new(i+1,1)-(2.0*V_new(i,1))+V_new(i-1,1))/delta_r**2))                                     
            ! centrifugal force integration from i to N
            centrifugal_force = 0
            do j = i, N+1
                centrifugal_force = centrifugal_force + rho(j)*omega**2*(r_location(j))*delta_r
            end do

            S_new(i,:) = S_old(i,:) + delta_t*(centrifugal_force*((V_new(i,:)-V_new(i-1,:))/delta_r))
            ! Calculation of the deformation value from the velocity calculated:
            q_new(i,:) = q_old(i,:) + delta_t*((V_old(i,:)))
        end do

        ! Boundary conditions at the TIPS position:
        M_new(N+1:N_on_blade,:) = 0.0
    end subroutine bladeSolve

end module structural