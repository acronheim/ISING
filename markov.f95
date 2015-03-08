module markov 
  use constants
  implicit none
  private
  public :: gen_config

contains
  
  subroutine gen_config(S,dE,BJ,h)
    integer, intent(inout) :: S(:,:)
    real(dp), intent(out) :: dE
    real(dp), intent(in) :: BJ, h
    integer :: i, j, x(2), S0(L+2,L+2), S_t(L,L)
    real(dp) :: BF, dE_t, r
    
    S_t = S ! initialize trial config
    call random_spin(x)
    
    ! create trial config by flipping 1 spin:
    i = x(1); j = x(2)
    S_t(i,j) = -S_t(i,j)
    
    ! add zero padding
    S0 = 0
    S0(2:L+1,2:L+1) = S_t
    i = i+1; j = j+1 ! adjust indices accordingly 

    ! you can just save the possible energy changes, boltzmann factors in &
    ! an array, may save some computation time. could do this in init_energy
    dE_t = -2._dp*BJ*S0(i,j)*(S0(i-1,j) + S0(i+1,j) + S0(i,j-1) + S0(i,j+1)) -&
      2._dp*h*S0(i,j) ! calculate change in energy

    ! accept trial config 
    if (dE_t < 0._dp) then 
      S = S_t ! if energy decreases always accept
      dE = dE_t
    else ! else accept config with probability of BF
      BF = exp(-dE_t)
      call random_number(r)
      if (r<BF) then
        S = S_t
        dE = dE_t
      else
        S = S
        dE = 0._dp
      endif
    endif
  end subroutine

  subroutine random_spin(x)
    ! returns index of randomly picked spin
    integer, intent(out) :: x(:)
    real(dp) :: u(2)

    call random_number(u)
    u = L*u + 0.5_dp
    x = nint(u) ! index of spin to flip
  end subroutine
end module