program elasticAcuatorLine
    use structural
    use parameters
    real :: time = 0.0
    real,dimension(2) :: F_a, g
    integer :: j
    integer :: counter = 1
    real,dimension(N+2,2,2) :: EI
    real,dimension(N+2) :: r_location,rho
    real, dimension(N+2,2) :: V_new, V_old, S_new, S_old, M_new, M_old, q_new, q_old
    ! result
    real, dimension(T/t_count+1,N+2,2) :: q_result,M_result
    
    ! Main program

    call bladeConfig(EI,r_location,rho)
    call bladeInit(V_new,M_new,S_new,q_new,V_old,M_old,S_old,q_old)

    do j = 1, T            
        g(1) = 0.0!
        g(2) = -9.8*cos(0-omega*time)

        F_a(1) = 0.0
        F_a(2) = 0.0!-100.0

        call bladeSolve(V_new,M_new,S_new,q_new,V_old,M_old,S_old,q_old,&
                            rho,EI,r_location,g,F_a,delta_t,delta_r)
        ! Upload data to new time-step:
        V_old = V_new 
        M_old = M_new 
        q_old = q_new
        ! Time control:
        time = delta_t*(j-1)

        if (mod(j,t_count) == 0) then
            q_result(counter,:,:) = q_new
            M_result(counter,:,:) = M_new 
            print *, time,M_result(counter,2,2)
            counter = counter + 1
        end if
    end do

    ! Output
    call bladeOutput(q_result,M_result,V_new,M_new,S_new,q_new)

end program elasticAcuatorLine