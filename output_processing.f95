module output_processing 
  use constants
  implicit none
  private 
  public :: calc_M_chi, calc_corr_function, calc_spec_heat 

contains

  pure subroutine calc_M_chi(L,N_SW,N_SW_2,m,Q,Mag,err_Mag,chi_s,chi_s_err,chi,err_chi,&
      method)
    ! calculates magnetization and susceptibility
    integer, intent(in)      :: L, method
    integer(lng), intent(in) :: m(:), N_SW(:), N_SW_2(:)
    real(dp), intent(out)    :: Q, Mag, err_Mag, chi, chi_s, err_chi, chi_s_err

    real(dp) :: N_SW_mean, N, abs_m_r_block(n_avg), chi_s_block(n_blocks)
    real(dp), allocatable :: m_r(:)
    integer :: i

    allocate(m_r(n_meas))

    ! initialize variables
    N = real(L,dp)**2
    m_r = real(m,dp)
    
    ! calculate magnetization
    Mag = sum(abs(m_r))/(n_meas*N)
    err_Mag = std_err(m_r/N)

    ! calculate unsubtracted susceptibility
    if (method==1) then
      chi = sum(real(N_SW_2,dp)/N**2)/n_meas
      err_chi = std_err(real(N_SW_2,dp)/N**2)
    elseif (method==2) then
      N_SW_mean = sum(real(N_SW,dp))/n_meas
      chi = N_SW_mean/N 
      err_chi = std_err(real(N_SW,dp)/N)
    endif
    
    ! susceptibility
    chi_s = (sum(abs(m_r)**2)/n_meas - sum(abs(m_r)/n_meas)**2)/N
    do i=1,n_blocks ! error for chi_s
      call extract_block(abs(m_r), abs_m_r_block, i)
      chi_s_block(i) = (sum(abs_m_r_block**2)/n_avg - sum(abs_m_r_block/n_avg)**2)/N
    end do
    chi_s_err = std_err(chi_s_block, .true.)

    ! calculate binder cumulant
    Q = 1 - sum((m_r/N)**4/n_meas)/(3._dp*(sum((m_r/N)**2)/n_meas)**2)
    deallocate(m_r)
  end subroutine

  pure subroutine calc_spec_heat(BE,L,Cv,err_Cv)
    ! calculates specific heat, per particle
    real(dp), intent(in)     :: BE(:)
    integer(lng), intent(in) :: L
    real(dp), intent(out)    :: Cv, err_Cv
    
    real(dp)     :: N, mu_BE, mu_BE_2, mu_BE_4

    N = real(L,dp)**2
    mu_BE = sum(BE)/n_meas
    mu_BE_2 = sum((BE-mu_BE)**2)/n_meas
    mu_BE_4 = sum(block_avg(BE-mu_BE)**4)/n_blocks**2

    Cv = mu_BE_2/N
    err_Cv = sqrt(mu_BE_4)/N
  end subroutine

  pure subroutine calc_corr_function(g,r,c_ss_fit,c_ss,eta,err_eta)
    ! calculate correlation function 
    real(dp), intent(in)  :: g(:,:), r(:)
    real(dp), intent(out) :: c_ss(:), c_ss_fit(:), eta, err_eta
    
    real(dp) :: offset

    c_ss = sum(g,1)/n_meas 
    call lin_fit(eta,err_eta,offset,-log(c_ss),log(r))
    c_ss_fit = exp(-offset)*r**(-eta)
  end subroutine 
    
  pure function std_err(A, blockmode)
    ! calculates std error of A from blocked data
    real(dp), intent(in) :: A(:)
    real(dp) :: std_err, sigma_blocks_2, Avg(n_blocks)
    logical, intent(in), optional  :: blockmode

    if (present(blockmode) .and. (blockmode .eqv. .true.)) then
      Avg = A
    else
      Avg = block_avg(A)
    end if
  
    sigma_blocks_2 = sum((Avg - sum(Avg)/n_blocks)**2)/(n_blocks-1)
    std_err = sqrt(sigma_blocks_2/n_blocks)
  end function

  pure subroutine extract_block(input, block, block_nr)
    real(dp), intent(in) :: input(:)
    real(dp), intent(out) :: block(n_avg)
    integer, intent(in):: block_nr   

    block = input(n_avg*(block_nr-1)+1:n_avg*(block_nr))
  end subroutine

  pure function block_avg(A)
    ! returns array containing block average of A
    real(dp), intent(in) :: A(:)
    real(dp) :: block_avg(n_blocks)
    integer :: j

    do j = 0,(n_blocks-1)
      block_avg(j+1) = sum(A(n_avg*j+1:n_avg*(j+1)))/n_avg
    enddo
  end function
  
  pure subroutine lin_fit(slope,err_slope,offset,y,x)
    real(dp), intent(out) :: slope, err_slope, offset
    real(dp), intent(in)  :: y(:), x(:)
    real(dp) :: mu_x, mu_y, ss_yy, ss_xx, ss_yx, s
    ! linear regression
    ! see also: http://mathworld.wolfram.com/LeastSquaresFitting.html

    mu_x = sum(x)/size(x)
    mu_y = sum(y)/size(y)
    ss_yy = sum((y - mu_y)**2)
    ss_xx = sum((x - mu_x)**2)
    ss_yx = sum((x - mu_x)*(y - mu_y))
    
    slope = ss_yx/ss_xx
    offset = mu_y - slope*mu_x
    s = sqrt((ss_yy - slope*ss_yx)/(size(x)-2))

    err_slope = s/(sqrt(ss_xx)*6._dp) 
  end subroutine
end module
