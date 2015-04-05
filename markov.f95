module markov 
  use constants
  use initialize
  use output_processing 
  use io
  use Wolff
  use Swendsen_Wang
  use plotroutines
  implicit none
  private
  public :: singlerun, autorun

contains
  subroutine singlerun(method,calc_css,auto)
    integer, intent(in) :: method
    logical, intent(in) :: calc_css, auto

    real(dp), allocatable :: c_ss(:), r(:), c_ss_fit(:)
    real(dp)              :: BJ, nu, err_nu, chi_s, chi, err_chi, Mag, &
                             err_Mag, Cv, err_Cv
    integer, allocatable  :: S(:,:)
    integer               :: runtime, L, r_max, n_corr

    call user_in(auto,L,r_max,n_corr,BJ)
    allocate(S(L,L),c_ss(r_max),c_ss_fit(r_max),r(r_max))
    call init_lattice(S,L)

    call markov_chain(S,method,auto,r_max,n_corr,BJ,r,Mag,err_Mag,runtime,&
      calc_css,c_ss,c_ss_fit,nu,err_nu,chi_s,chi,err_chi,Cv,err_Cv)
    
    call results_out(BJ,L,r,runtime,calc_css,c_ss,c_ss_fit,nu,err_nu,&
      chi_s,chi,err_chi,Mag,err_Mag,Cv,err_Cv)
    deallocate(S,r,c_ss,c_ss_fit) 
  end subroutine

  subroutine autorun(method,auto)
    integer, intent(in) :: method
    logical, intent(in) :: auto
    
    integer, allocatable  :: S(:,:)
    real(dp), allocatable :: BJ(:), chi_s(:), chi(:), err_chi(:), Mag(:), &
      err_Mag(:), Cv(:), err_Cv(:), Q(:)
    integer  :: i, L, T_s=51, runtime, r_max, n_corr
    real(dp) :: nu, r(1), c_ss(1), c_ss_fit(1), err_nu
    logical  :: calc_css = .false.

    ! initialize
    call user_in(auto,L,r_max,n_corr)
    allocate(S(L,L),BJ(T_s),chi_s(T_s),chi(T_s),err_chi(T_s),Mag(T_s),&
      err_Mag(T_s),Cv(T_s),err_Cv(T_s),Q(T_s))
    call init_lattice(S,L)
    call init_BJ(T_s,BJ)
    
    ! iterate over temperatures, ask temp range from user?
    do i = 1,T_s
      call markov_chain(S,method,auto,r_max,n_corr,BJ(i),r,Mag(i),err_Mag(i),&
        runtime,calc_css,c_ss,c_ss_fit,nu,err_nu,chi_s(i),chi(i),err_chi(i),&
        Cv(i),err_Cv(i))

      ! move this into markov chain routine?
      Q(i) = chi_s(i)/(real(L,dp)**2*Mag(i)**2) + 1._dp
    enddo
    
    call auto_results(L,BJ,Q,Mag,chi_s,Cv)
    deallocate(S,BJ,chi_s,chi,err_chi,Mag,err_Mag,Cv,err_Cv,Q)
  end subroutine

  subroutine markov_chain(S,method,auto,r_max,n_corr,BJ,r,Mag,err_Mag,&
      runtime,calc_css,c_ss,c_ss_fit,nu,err_nu,chi_s,chi,err_chi,Cv,err_Cv)
    integer, intent(inout)  :: S(:,:)
    real(dp), intent(inout) :: BJ
    integer, intent(in)     :: method, r_max, n_corr
    logical, intent(in)     :: calc_css, auto
    integer, intent(out)    :: runtime
    real(dp), intent(out)   :: c_ss(:), r(:), c_ss_fit(:), Mag, err_Mag, nu, &
      err_nu, chi_s, chi, err_chi, Cv, err_Cv

    integer(lng), allocatable :: N_SW(:), N_SW_2(:), m(:)
    real(dp), allocatable     :: g(:,:), BE(:)
    integer(lng) :: start_time, m_tmp, N_SW_tmp, N_SW_2_tmp, end_time
    integer      :: i, j, L
    real(dp)     :: p
    
    ! initialize needed variables
    L = size(S,1)
    j = 0
    N_SW_tmp = 0
    N_SW_2_tmp = 0
    m_tmp = 0
    r = real((/(i,i=1,r_max)/),dp)
    p = 1 - exp(-2._dp*BJ)
    
    ! allocate memory
    allocate(g(n_meas,r_max),BE(n_meas),N_SW(n_meas),N_SW_2(n_meas),m(n_meas))
    
    call animate_lattice()
    call system_clock(start_time)
    
    ! run markov chain
    do i=1,steps
      call gen_config(S,L,m_tmp,N_SW_tmp,N_SW_2_tmp,p,method)

      if (i > meas_start) then
        j = j+1
        m(j) = m_tmp ! record magnetization
        N_SW(j) = N_SW_tmp ! record clustersize
        N_SW_2(j) = N_SW_2_tmp ! squared cluster size

        if (calc_css) call s_corr(g(j,:),S,L,r_max,n_corr)
        call calc_energy(BE(j),S,L,BJ)
      endif

      if ((mod(i,plot_interval) == 0) .and. (.not. auto)) then
        call write_lattice(S,L) ! write lattice to pipe
      endif
    enddo    
    
    call system_clock(end_time)
    call close_lattice_plot()
    
    ! calculate ensemble averages
    call calc_chi(L,N_SW,N_SW_2,m,Mag,err_Mag,chi_s,chi,err_chi,method)
    call calc_spec_heat(BE,L,Cv,err_Cv)
    if (calc_css) call calc_corr_function(g,r,c_ss_fit,c_ss,nu,err_nu)
    
    ! calculate runtime
    runtime = (end_time - start_time)/1000

    deallocate(g,BE,N_SW,N_SW_2,m)
  end subroutine

  subroutine gen_config(S,L,m,N_SW,N_SW_2,p,method)
    integer, intent(inout)    :: S(:,:)
    real(dp), intent(in)      :: p
    integer, intent(in)       :: L, method
    integer(lng), intent(out) :: m, N_SW, N_SW_2 
    
    if (method == 1) then
      call SW_clust(S,L,m,N_SW,N_SW_2,p)
    elseif (method == 2) then
      call Wolff_clust(S,L,m,N_SW,p)
      N_SW_2 = N_SW**2
    endif
  end subroutine
  
  pure subroutine calc_energy(BE,S,L,BJ)
    integer, intent(in)   :: S(:,:), L
    real(dp), intent(in)  :: BJ
    real(dp), intent(out) :: BE
    
    integer :: i, j, k, nn(4,2)

    BE = 0._dp ! initialze energy 

    do i = 1,L
      do j = 1,L
        nn = nn_idx([i,j],L) ! get nearest neighbors of spin i,j
        do k = 1,4
          BE = BE - BJ*S(i,j)*S(nn(k,1),nn(k,2))
        enddo
      enddo
    enddo

    BE = 0.5_dp*BE ! account for double counting of pairs
    !BE = BE - h*sum(S) ! add external field
  end subroutine

  pure subroutine s_corr(g,S,L,r_max,n_corr)
    ! calculate spin-spin correlation function
    integer, intent(in)   :: S(:,:), L, r_max, n_corr
    real(dp), intent(out) :: g(:)
    
    real(dp) :: g_tmp(n_corr,r_max)
    integer  :: i, r_0
    ! intialize starting position             
    r_0 = (L-n_corr)/2
    
    ! iterate over n_corr spins on diagonal, get correlation up to r_max
    do i=1,n_corr
      g_tmp(i,:) = S(i+r_0, i+r_0)*S(i+r_0, i+r_0+1:i+r_0+r_max) 
    enddo

    g = sum(g_tmp,1)/n_corr 
  end subroutine
end module
