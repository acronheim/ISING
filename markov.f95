module markov 
  use constants
  implicit none
  private
  public :: gen_config, sync_branches

contains
  
  subroutine gen_config(S,dE,BJ,h)
    integer, intent(inout) :: S(:,:)
    real(dp), intent(out) :: dE
    real(dp), intent(in) :: BJ, h
    integer, allocatable :: S_t(:,:)
    integer :: i, j, x(2), S_nbrs
    real(dp) :: BF, dE_t, r
    
    allocate(S_t(L+2,L+2))
    
    S_t = S ! initialize trial config
    dE = 0._dp ! init dE

    ! create trial config by flipping 1 spin:
    call random_spin(x) 
    i = x(1)+1; j = x(2)+1 ! adjust for zero padding
    S_t(i,j) = -S_t(i,j)
    
    ! calculate change in BE
    S_nbrs = S_t(i-1,j) + S_t(i+1,j) + S_t(i,j-1) + S_t(i,j+1)
    dE_t = -2._dp*BJ*S_t(i,j)*S_nbrs - 2._dp*h*S_t(i,j) 
   
    if (dE_t < 0._dp) then  
      S = S_t ! if energy decreases always accept
      dE = dE_t
    else ! else accept config with probability of BF
      BF = exp(-dE_t)
      call random_number(r)
      if (r < BF) then
        S = S_t
        dE = dE_t
      endif
    endif

    deallocate(S_t)
  end subroutine

  subroutine random_spin(x)
    ! returns index of randomly picked spin
    integer, intent(out) :: x(:)
    real(dp) :: u(2)

    call random_number(u)
    u = L*u + 0.5_dp
    x = nint(u) ! index of spin to flip
  end subroutine

  subroutine sync_branches(S,BE_branch)
    ! select branch with lowest energy
    integer, intent(inout) :: S(:,:,:)
    real(dp), intent(inout) :: BE_branch(:)
    integer :: i, branch

    branch = minloc(BE_branch,1)
    BE_branch = minval(BE_branch)

    do i = 1,n_br
      S(i,:,:) = S(branch,:,:)       
    enddo
  end subroutine
end module
