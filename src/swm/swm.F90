program swm_agrid
  use module_constants, only : Tdim   ! two dim
  use module_control,   only : imsx,imex,ips,ipe
  use module_control,   only : nip,niE,niV,NPTS,NPTS_MX,ndimFV
  use module_control,   only : nrkl,filenamelen,itsend,dt,ck_stswm
  use module_control,   only : TS0, TS1, TX
  use module_variables, only : FV,FV0,FVrk,AFV,h_tmp
  use module_variables, only : ASV,ASV_ana
  use module_variables, only : lon,lat,lat_v,lon_v,gh_icos,gh_icos_yx,hb_a,fcori
  use module_variables, only : np_ll,g_ll_unif,f_ll,lon_deg_plot,lat_deg_plot,field_plot
  use module_variables, only : L1n,L2n,Linfn,TE0,PE0, Epot0, diag_var, stagecoef
  use readnamelist,     only : stagger,alpha,iswcase,itsbeg,glvl
  use readnamelist,     only : nx_ll,ny_ll,nx_NICAM,mx_NICAM,ck_consv_freq,ncop_freq
  use kinds,            only : rt
  !
!  use icosio    , only: icosio_stop
!  use core_setup, only: iam_nim_task
  implicit none
  !
  integer :: ret, fnlen, LDA, LDB, nstep_ck, nstep_ncop, istage
  character (len=filenamelen) :: filename
  character (len=20)  :: mystring    ! fixed extern func
  character (len=200) :: line, line2
  integer :: ipt, npos, niX          ! niX= nih: position for h field
  integer :: iter, iter_ck, iter_ncop, ip, j, iunit, iv
  integer :: ystart, yend
  real(rt):: Lnorm(3,10)   ! array for 1:3=L1,L2,Linf   1:3=grad,div,vort
  real(rt):: P_sph(3), basis(3,3), invbasis(3,3), X(3), Y(3), temp(3), s
  real(rt):: global_sum(6)
  logical :: cond


  fnlen=filenamelen
  write(6,121)  'bf call start'
  call start         ! --> start_nml_mpi.F90;  call control (readnamelist)
  write(6,121)  'imsx,imex,ndimFV     ', imsx,imex,ndimFV
  write(6,121)  'nip,niE,NPTS,NPTS_MX ', nip,niE,NPTS,NPTS_MX
  call input         ! --> input_grid_metrics.F90


!!  call test_sphere

  call simple_init
  !------------------------------------------------------------
  !   a. check final state of stswm
  !   b. find error norm of A/C grid in time evolution against stswm
  !------------------------------------------------------------
  if ( ck_stswm .EQ. 1 ) then
     call ck_against_stswm
     STOP 'swm stop as instructed by ck_stswm'     
  endif


  if (stagger(1:2).EQ.'Du') then
     call sw_test_init (NPTS_MX,nip,lon_v,lat_v,niV,lon,lat,alpha,hb_a,FVrk,fcori,ndimFV,iswcase,stagger)
  else
     call sw_test_init (NPTS_MX,NPTS,lon_v,lat_v,nip,lon,lat,alpha,hb_a,FVrk,fcori,ndimFV,iswcase,stagger)
  endif
  !
  write(6,101) 'af sw_test_init'


  !!call sw_test_init (NPTS,lon_v,lat_v,alpha,hb_a,FVrk,fcori,ndimFV,iswcase)
  !
  !!write(6, *) 'bf  test_spherical_harmonics (10, 10, 20)'
  !!call test_spherical_harmonics (10, 10, 10)



  !------------------------------------------------------
  ! NICAM test output  Six vector  ASV, ASV_ana
  !------------------------------------------------------

  if (iswcase.EQ.0) then
     if (stagger .EQ. 'A') then
        call test_dyn_A_3func (imsx,imex,ips,ipe,FVrk,ndimFV)
!        STOP 'end compare V(1:3), div, vort'
        ! grad_x
        ystart=1; yend=3
        call error_norm_vector (imsx,imex,ips,ipe,ndimFV,ystart,yend,ASV(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
        Lnorm(1,1)=L1n; Lnorm(2,1)=L2n; Lnorm(3,1)=Linfn
        !
        !
        ! div(V)
        call error_norm_scalar (imsx,imex,ips,ipe,ASV(1,4),ASV_ana(1,4),L1n,L2n,Linfn)
        Lnorm(1,2)=L1n; Lnorm(2,2)=L2n; Lnorm(3,2)=Linfn
        !
        ! vort(V)
        call error_norm_scalar (imsx,imex,ips,ipe,ASV(1,5),ASV_ana(1,5),L1n,L2n,Linfn)
        Lnorm(1,3)=L1n; Lnorm(2,3)=L2n; Lnorm(3,3)=Linfn
        !
        ! div(hV)
        call error_norm_scalar (imsx,imex,ips,ipe,ASV(1,6),ASV_ana(1,6),L1n,L2n,Linfn)
        Lnorm(1,4)=L1n; Lnorm(2,4)=L2n; Lnorm(3,4)=Linfn
        !
        ! ( f + xi ) k \cross V
        ystart=7; yend=9
        call error_norm_vector (imsx,imex,ips,ipe,9,ystart,yend,ASV(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
        Lnorm(1,5)=L1n; Lnorm(2,5)=L2n; Lnorm(3,5)=Linfn
        !
        ! grad(K)
        ystart=10; yend=12
        call error_norm_vector (imsx,imex,ips,ipe,12,ystart,yend,ASV(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
        Lnorm(1,6)=L1n; Lnorm(2,6)=L2n; Lnorm(3,6)=Linfn
        !
        filename='d.NICAM_operator_error'
        open(7, file=filename, status='unknown')
        write(7, '(a)') '# glvl   L1,L2,Linf for grad div vorticity'
        write(7, '(2x,i4,2x,9f12.8)')  glvl, Lnorm(1:3,1), Lnorm(1:3,2), Lnorm(1:3,3)
        write(7, *)
        write(7, '(a)') '# glvl   L1,L2,Linf for div(hV) (f+xi)k \cross V'
        write(7, '(2x,i4,2x,9f12.8)')  glvl, Lnorm(1:3,4), Lnorm(1:3,5)
        close(7)
        !
        !
        STOP 'finish compare with NICAM analytical results ::  A-grid'

     elseif (stagger .EQ. 'C') then
        write(6, *) 'bf C_3func'
        call test_dyn_C_3func (imsx,imex,ips,ipe,FVrk,ndimFV)
        !
        write(6, *) 'af C_3func'
!        STOP 'end compare V(1:3), div, vort'
        ! grad_x
        ystart=1; yend=1
        call error_norm_scalar (imsx,imex,ips,ipe,ASV(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
        Lnorm(1,1)=L1n; Lnorm(2,1)=L2n; Lnorm(3,1)=Linfn
        !
        write(6, *) 'af grad_x'
!!        write(6, 101) 'C-grid: grad_h, L1,2,inf', Lnorm(1:3,1)
        !
        ! div(V)
        call error_norm_scalar (imsx,imex,ips,ipe,ASV(1,4),ASV_ana(1,4),L1n,L2n,Linfn)
        Lnorm(1,2)=L1n; Lnorm(2,2)=L2n; Lnorm(3,2)=Linfn
        !
        ! vort(V)
        ystart=5; yend=5
        call error_norm_vector (imsx,imex,ips,ipe,12,ystart,yend,ASV,ASV_ana,L1n,L2n,Linfn)
        Lnorm(1,3)=L1n; Lnorm(2,3)=L2n; Lnorm(3,3)=Linfn
        !
        ! div(hV)
        call error_norm_scalar (imsx,imex,ips,ipe,ASV(1,6),ASV_ana(1,6),L1n,L2n,Linfn)
        Lnorm(1,4)=L1n; Lnorm(2,4)=L2n; Lnorm(3,4)=Linfn
        !
        ! ( f + xi ) k \cross V
        ystart=7; yend=7
        call error_norm_vector (imsx,imex,ips,ipe,12,ystart,yend,ASV,ASV_ana,L1n,L2n,Linfn)
        Lnorm(1,5)=L1n; Lnorm(2,5)=L2n; Lnorm(3,5)=Linfn

        write(6, 101)  'eta V  calc, ana, diff'
        do ip=1, niE, 1
           write(6, 103)  lon_v(ip), lat_v(ip), &
                ASV(ip,7), ASV_ana(ip,7), ASV(ip,7)-ASV_ana(ip,7)
        enddo
!!        stop 'ck eta V'
        !
        !  U_t^e
        !
        ystart=9; yend=9
        call error_norm_vector (imsx,imex,ips,ipe,12,ystart,yend,ASV,ASV_ana,L1n,L2n,Linfn)
        Lnorm(1,6)=L1n; Lnorm(2,6)=L2n; Lnorm(3,6)=Linfn


        filename='d.NICAM_operator_error'
        open(7, file=filename, status='unknown')
        write(7, '(a)') '# glvl   L1,L2,Linf for grad div vorticity'
        write(7, '(2x,i4,2x,9f12.8)')  glvl, Lnorm(1:3,1), Lnorm(1:3,2), Lnorm(1:3,3)
        write(7, *)
        write(7, '(a)') '# glvl   L1,L2,Linf for div(hV) (f+xi)k \cross V   U_t^e'
        write(7, '(2x,i4,2x,9f12.8)')  glvl, Lnorm(1:3,4), Lnorm(1:3,5), Lnorm(1:3,6)
        close(7)
        !
        !
        STOP 'finish compare with NICAM analytical results  :: C-grid'
     endif
  endif



  if (stagger .EQ. 'A') then
     npos=npts
     filename='d.h_it.0'
     gh_icos(1:npos,1)=lon_v(1:npos); gh_icos(1:npos,2)=lat_v(1:npos); h_tmp(1:npos)=FVrk(1:npos,1)
     call output_icos (nip, nip, Tdim, npos, gh_icos, h_tmp, filename, fnlen)
     ! use the same position at center
     filename='d.vx_it.0'
     call output_icos (nip, nip, Tdim, npos, gh_icos, FVrk(1,2), filename, fnlen)
  elseif (stagger .EQ. 'C' .OR. stagger(1:2) .EQ. 'ZM') then
     ! pos at center
     npos=nip
     gh_icos(1:npos,1)=lon(1:npos); gh_icos(1:npos,2)=lat(1:npos); h_tmp(1:npos)=FVrk(1:npos,1)
     filename='d.h_it.0'
     ! note  ime_E dim for gh_icos, h_tmp in module_allocate.F90 
     !
     call output_icos (NPTS, nip, Tdim, npos, gh_icos, h_tmp, filename, fnlen)
     ! pos at edges, velocity
     npos=npts
     gh_icos(1:npos,1)=lon_v(1:npos); gh_icos(1:npos,2)=lat_v(1:npos)    ! shift grid to edges
     filename='d.vx_it.0'
     call output_icos (NPTS, NPTS, Tdim, npos, gh_icos, FVrk(1,2), filename, fnlen)
  elseif (stagger(1:2) .EQ. 'Du') then
!     !
!     ! pos at center : for vel
!     npos=nip
!     gh_icos(1:npos,1)=lon_v(1:npos); gh_icos(1:npos,2)=lat_v(1:npos); h_tmp(1:npos)=FVrk(1:npos,2)
!     filename='d.vx_it.0'
!     call output_icos (niV, niV, Tdim, npos, gh_icos, h_tmp, filename, fnlen)
!
     ! pos at upp Vertex
     npos=niV
     gh_icos(1:npos,1)=lon(1:npos); gh_icos(1:npos,2)=lat(1:npos); h_tmp(1:npos)=FVrk(1:npos,1)
     filename='d.h_it.0'
     call output_icos (npos, npos, Tdim, npos, gh_icos, h_tmp, filename, fnlen)
  endif

  write(6,101) 'af output_icos'





  !------------------------------------------------------------
  !             use slint to create .nc file
  !------------------------------------------------------------
  call gen_uniform_ll (nx_ll,ny_ll,g_ll_unif, np_ll, Tdim)
  if (stagger .EQ. 'A'  .OR. stagger .EQ. 'C') then
     niX= nip
  elseif (stagger(1:2) .EQ. 'Du') then
     niX= niV
  elseif (stagger(1:2) .EQ. 'ZM') then
     niX= nip
  else
     niX= -1
  endif
  gh_icos_yx(1:niX,1)=lat(1:niX); gh_icos_yx(1:niX,2)=lon(1:niX); filename='dll.h_it.0'
  ! over mountain
  if (iswcase.EQ.5) then
     do ip=1,niX; h_tmp(ip)=FVrk(ip,1)+hb_a(ip); enddo
  endif


  write(6,101) 'bf icos_2_ll'
  write(6,121) 'niX, np_ll, Tdim', niX, np_ll, Tdim
  write(6,121) 'shape(gh_icos_yx),shape(h_tmp),shape(g_ll_unif),shape(f_ll)', &
       shape(gh_icos_yx),shape(h_tmp),shape(g_ll_unif),shape(f_ll)

  call icos_2_ll (niX, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
  write(6,101) 'af icos_2_ll'


  call change_latlon_array_shape &
       (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
  filename='h0_field.nc'
  call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
  write(6,101) 'af nc'



  if (iswcase.EQ.0) then
!!   if (iswcase.EQ.0 .OR. iswcase.EQ.2) then
     if (stagger .EQ. 'A') then
        ! grad(h)
        do ip=1,nip
           X(1)=ASV(ip,1)-ASV_ana(ip,1)
           X(2)=ASV(ip,2)-ASV_ana(ip,2)
           X(3)=ASV(ip,3)-ASV_ana(ip,3)
           call X_2norm(X, s, 3)
           h_tmp(ip)=s
        enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_gradh_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)

        ! vorticity
        do ip=1,nip; h_tmp(ip)=abs(ASV(ip,5)-ASV_ana(ip,5)); enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_vort_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)

        ! div(hV)
        do ip=1,nip; h_tmp(ip)=abs(ASV(ip,6)-ASV_ana(ip,6)); enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_divhV_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)

        ! (f+xi) k \cross V
        do ip=1,nip
           X(1)=ASV(ip,7)-ASV_ana(ip,7)
           X(2)=ASV(ip,8)-ASV_ana(ip,8)
           X(3)=ASV(ip,9)-ASV_ana(ip,9)
           call X_2norm(X, s, 3)
           h_tmp(ip)=s
        enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_etakcrossV_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        !
        !
        ! grad(K) or   total \p V / \p t
        do ip=1,nip
           X(1)=ASV(ip,10)-ASV_ana(ip,10)
           X(2)=ASV(ip,11)-ASV_ana(ip,11)
           X(3)=ASV(ip,12)-ASV_ana(ip,12)
           !
           ! test for kinetic only, fake
           !X(1)=ASV(ip,10)
           !X(2)=ASV(ip,11)
           !X(3)=ASV(ip,12)
           call X_2norm(X, s, 3)
           h_tmp(ip)=s
        enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_grad_K.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        !
        STOP 'the end of nc write for isw_case == 0 or 2 ::  A-grid'
        !
     elseif (stagger .EQ. 'C') then
        ! grad(h)
        do ip=1,niE
           s=ASV(ip,1)-ASV_ana(ip,1)    ! \p h / \p n   directional derivative
           h_tmp(ip)=s
        enddo
        ! repeat
        call icos_2_ll (niE, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_gradh_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)

        ! vorticity
        do ip=1,niE; h_tmp(ip)=abs(ASV(ip,5)-ASV_ana(ip,5)); enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_vort_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)

        ! div(hV)
        do ip=1,nip; h_tmp(ip)=abs(ASV(ip,6)-ASV_ana(ip,6)); enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_divhV_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)

        ! (f+xi) k \cross V
        do ip=1,niE
           X(1)=ASV(ip,7)-ASV_ana(ip,7)
!!           X(2)=ASV(ip,8)-ASV_ana(ip,8)
!!           X(3)=ASV(ip,9)-ASV_ana(ip,9)
!!           call X_2norm(X, s, 3)
           s= abs(X(1))
           h_tmp(ip)= s
        enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_etakcrossV_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)

        ! U_t^e
        do ip=1,niE; h_tmp(ip)=abs(ASV(ip,9)-ASV_ana(ip,9)); enddo
        ! repeat
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)  
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='TC0_ut_diff_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        !
        STOP 'the end of nc write for isw_case == 0 or 2 ::  C-grid'
     endif    ! stagger 
  endif       !


  write(6, 101) 'af iswcase == 0'


  !! if ( iswcase.GT.8 ) STOP 'iswcase .GT. 8, not implemented'
  open (7, file='d.err_vs_time', status='unknown')
  write(7,'(a)') '#  time(day) :  L1,L2,Linf (h_field)      L1,L2,Linf (Velocity)'
  filename='d.TE_PE'
  iunit=8
  open(iunit, file=filename, status='unknown')
  write(iunit, 101) "# time (day)   tot    1:h^*     2:TE      3:PE      4:vort_V      5:div_V     6: dTE      7: dPE     8:Ratio(KE_eff)"



  !------------------------------------------------------------
  !
  !  M a i n          S  W        D y n a m i c s 
  !
  !------------------------------------------------------------  
  !
  !  output every 12 or 6 hours     !  post time = dt * iter
  !
  !
  nstep_ck   = int ( 86400.d0 * ck_consv_freq / dt )  !  default:  ck_consv_freq  0.5  d
  nstep_ncop = int ( 86400.d0 * ncop_freq     / dt )  !  default:  ncop_freq      15.0 d

  write(6, 112)  'ck_consv_freq, ncop_freq', ck_consv_freq, ncop_freq
  write(6, 121)  'nstep_ck, nstep_ncop', nstep_ck, nstep_ncop



  iter=0
  if (iswcase .GE. 5) then     ! initial conserved properties
     if (stagger .EQ. 'A') then
        call TE_PE_calc_A (imsx,imex,ips,ipe,FVrk,global_sum,ndimFV,iter,iunit)
     elseif (stagger .EQ. 'C') then
        call TE_PE_calc_C (imsx,imex,ips,ipe,FVrk,global_sum,ndimFV,iter,iunit)
     endif
     TE0=global_sum(2); PE0=global_sum(3); Epot0=global_sum(6)    ! assign init values
     write(6, 201) 'TE0, PE0, Epot0', TE0, PE0, Epot0
  endif
  
  do iter=itsbeg,itsend
     !
     if (mod(iter, nstep_ck).eq.0) then
        iter_ck = iter
     endif
     if (mod(iter, nstep_ncop).eq.0) then
        iter_ncop = iter
     endif
     !
     !!write(6, 121) 'head, iter=',iter

     call Zerotendencies_sw (imsx,imex,ips,ipe,AFV,ndimFV)
     call Save_State_Var_sw (imsx,imex,ips,ipe,FVrk,FV0,ndimFV) ! updated i-timex RK from previous step
     !
     istage=1

     if (iswcase.EQ.10) then
        tx= (iter-1)*dt +  stagecoef(istage)
        call TC10_hs_t(tx, niX, lon, lat, hb_a(1)) ! T-dep terrain: TC10, covers A/C/Du
     endif

!!     call cpu_time(TS0)
     if (stagger .EQ. 'A') then
        call dyn_A (imsx,imex,ips,ipe,FVrk,AFV,ndimFV,iter)
     elseif (stagger .EQ. 'C') then
        ! call dyn_C (imsx,imex,ips,ipe,FVrk,AFV,ndimFV,iter)
        call dyn_C_combine_Cori_kinetic (imsx,imex,ips,ipe,FVrk,AFV,ndimFV,iter)
        ! call dyn_C_2comp_base (imsx,imex,ips,ipe,FVrk,AFV,ndimFV,iter)
     elseif (stagger(1:2) .EQ. 'Du') then    ! 'DuSPR' or 'DuCVT'
        call dyn_DuCVTSPR (imsx,imex,ips,ipe,FVrk,AFV,ndimFV,iter)
!!        call dyn_DuSPR (imsx,imex,ips,ipe,FVrk,AFV,ndimFV,iter)
     elseif (stagger(1:2) .EQ. 'ZM') then    ! ZM
        call dyn_ZM (imsx,imex,ips,ipe,FVrk,AFV,ndimFV,iter)
     endif

!!     call cpu_time(TS1)
!!     write(6, 101) 'call dyn A/C/Du : ', TS1-TS0


     if (iswcase.EQ.4) then  ! post add
        call add_forcing_TC4 (imsx,imex,ips,ipe,ASV_ana(1,1),AFV,ndimFV,iter,istage)  ! AFV: (inout)
     endif
     call flush(6)
     call RK4th_true_Var_increment (istage,FVrk,AFV,ndimFV,imsx,imex,ips,ipe) 


     if ( mod (iter, 10) .EQ. 1 ) then
        !--print
!        write(6, 121) 'iter, istage =', iter, istage
!        write(6, *)
!        do ipt=1,NPTS, NPTS/10
!           write(line2,'(2x,a,i5,(1x,f10.2),3(1x,f10.3))') &
!                'ipt, FVrk(ipt,1:4), AFV(ipt,1:4)= ', ipt, FVrk(ipt,1:ndimFV)
!           write(line, '(a)') "   |"
!           line2=trim(line2)//trim(line)
!           write(line,'(4(1x,f14.9))') AFV(ipt,1:ndimFV)
!           line2=trim(line2)//trim(line)
!           write(6, '(a)') trim(line2)
!        enddo
     end if



     do istage=2,nrkl

        if (iswcase.EQ.10) then
           tx= (iter-1)*dt +  stagecoef(istage)
           call TC10_hs_t(tx, niX, lon, lat, hb_a(1)) ! T-dep terrain: TC10, covers A/C/Du
        endif

        call RK4th_test_Var_increment (istage,FV0,AFV,FV,ndimFV,imsx,imex,ips,ipe)
        if (stagger .EQ. 'A') then
           call dyn_A (imsx,imex,ips,ipe,FV,AFV,ndimFV,iter)
        elseif (stagger .EQ. 'C') then
           ! call dyn_C (imsx,imex,ips,ipe,FV,AFV,ndimFV,iter)
           call dyn_C_combine_Cori_kinetic (imsx,imex,ips,ipe,FV,AFV,ndimFV,iter)
           ! call dyn_C_2comp_base (imsx,imex,ips,ipe,FV,AFV,ndimFV,iter)
        elseif (stagger(1:2) .EQ. 'Du') then
           call dyn_DuCVTSPR (imsx,imex,ips,ipe,FV,AFV,ndimFV,iter)
!!           call dyn_DuSPR (imsx,imex,ips,ipe,FV,AFV,ndimFV,iter)
        elseif (stagger(1:2) .EQ. 'ZM') then
           call dyn_ZM (imsx,imex,ips,ipe,FV,AFV,ndimFV,iter)
        endif
        if (iswcase.EQ.4) then
           call add_forcing_TC4 (imsx,imex,ips,ipe,ASV_ana(1,1),AFV,ndimFV,iter,istage)  ! AFV: (inout)
        endif
        !add
        call RK4th_true_Var_increment (istage,FVrk,AFV,ndimFV,imsx,imex,ips,ipe)
     enddo
     call flush(6)
     !
     !
!!     if (stagger .EQ. 'A') then     
!!        call remove_radial_component_A (imsx,imex,ips,ipe,FVrk,ndimFV)  ! remove \Rcvec direction
!!     endif



     !
     !   post processing  FVrk(1:nip,4), FVrk(1:niE,2/3)
     !
     !-- post.1   error norm  L1,L2,Linf
     !
!!     cond = (( ((iswcase.GE.2) .OR. (iswcase.LE.4) .OR. (iswcase.EQ.8)) .AND. (iter .EQ. iter_ck)) &
!!          .OR. (iswcase .GE. 5 .AND. iter .EQ. itsend))
!! YGYU, I avoid TC5 STSWM calc. 
     cond = (( ((iswcase.GE.2) .AND. (iswcase.LE.5) .OR. (iswcase.EQ.8)) .AND. (iter .EQ. iter_ck)))
     if (cond) then
        if (stagger .EQ. 'A') then
           !!write(6,101) 'post.1 error norm  L1,L2,Linf'
           !
           ! h
           call error_norm_scalar (imsx,imex,ips,ipe,FVrk(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
           Lnorm(1,1)=L1n; Lnorm(2,1)=L2n; Lnorm(3,1)=Linfn
           !
           ! velocity:  V_x, V_y, V_z
           ystart=2; yend=4
           call error_norm_vector (imsx,imex,ips,ipe,ndimFV,ystart,yend,FVrk,ASV_ana,L1n,L2n,Linfn)
           Lnorm(1,2)=L1n; Lnorm(2,2)=L2n; Lnorm(3,2)=Linfn
           !
           write(7,211)  dt*real(iter)/86400.d0,Lnorm(1:3,1),Lnorm(1:3,2)  ! time_h,L1n,L2n,Linfn, V_norm
           !
        elseif (stagger .EQ. 'C') then
           !
           ! h :   caution
           call error_norm_scalar (imsx,imex,1,nip,FVrk(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
           Lnorm(1,1)=L1n; Lnorm(2,1)=L2n; Lnorm(3,1)=Linfn
           !
           ! velocity:  U_n
           ystart=2; yend=2
           call error_norm_vector (imsx,imex,ips,ipe,ndimFV,ystart,yend,FVrk,ASV_ana,L1n,L2n,Linfn)
           Lnorm(1,2)=L1n; Lnorm(2,2)=L2n; Lnorm(3,2)=Linfn
           !
           write(7,211)  dt*real(iter)/86400.d0,Lnorm(1:3,1),Lnorm(1:3,2)  ! time_h,L1n,L2n,Linfn, V_norm
           !
        elseif (stagger(1:2) .EQ. 'Du') then
           !
           ! h
           call error_norm_scalar (imsx,imex,1,niV,FVrk(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
           Lnorm(1,1)=L1n; Lnorm(2,1)=L2n; Lnorm(3,1)=Linfn
           !
           ! velocity:  V_x, V_y, V_z
           ystart=2; yend=4
           call error_norm_vector (imsx,imex,ips,ipe,ndimFV,ystart,yend,FVrk,ASV_ana,L1n,L2n,Linfn)
           Lnorm(1,2)=L1n; Lnorm(2,2)=L2n; Lnorm(3,2)=Linfn
           !
           write(7,211)  dt*real(iter)/86400.d0,Lnorm(1:3,1),Lnorm(1:3,2)  ! time_h,L1n,L2n,Linfn, V_norm
           !
        elseif (stagger(1:2) .EQ. 'ZM') then
           !
           ! h :   caution
           call error_norm_scalar (imsx,imex,1,nip,FVrk(1,1),ASV_ana(1,1),L1n,L2n,Linfn)
           Lnorm(1,1)=L1n; Lnorm(2,1)=L2n; Lnorm(3,1)=Linfn
           !
           ! velocity:  V_x, V_y, V_z
           ystart=2; yend=4
           call error_norm_vector (imsx,imex,ips,ipe,ndimFV,ystart,yend,FVrk,ASV_ana,L1n,L2n,Linfn)
           Lnorm(1,2)=L1n; Lnorm(2,2)=L2n; Lnorm(3,2)=Linfn
           !
           write(7,211)  dt*real(iter)/86400.d0,Lnorm(1:3,1),Lnorm(1:3,2)  ! time_h,L1n,L2n,Linfn, V_norm
           !
        endif
        !!write(6,101) 'end post.1 error norm  L1,L2,Linf'
     endif
     !
     !-- YGYU ad-hoc add-in for PDF of error at G7 compare A/C grid
!!     if ( iter .EQ. itsend) then
!!        ! h  caution
!!        call error_norm_scalar_collect_bin (imsx,imex,1,nip,FVrk(1,1),ASV_ana(1,1),temp(1),temp(2),temp(3))
!!     endif
!!     write(6,*) 'af post.1'


     !
     !-- post.2  measure conseved quantity
     !
!     if ( (iswcase .GE. 5 .OR. iswcase .EQ. 2).AND. iter .EQ. iter_ck) then
     if ( (iswcase .GE. 5 .AND. iswcase .LE. 9).AND. iter .EQ. iter_ck ) then
        if (stagger .EQ. 'A') then
           call TE_PE_calc_A (imsx,imex,ips,ipe,FVrk,global_sum,ndimFV,iter,iunit)
        elseif (stagger .EQ. 'C') then
           call TE_PE_calc_C (imsx,imex,ips,ipe,FVrk,global_sum,ndimFV,iter,iunit)
        endif
     endif
!!     write(6,*) 'af post.2'



     !
     !-- post.3  output netcdf
     !
     if ( stagger .EQ. 'A' .AND. iter .EQ. iter_ncop ) then
        !------------------------------------------------------------
        !             use slint to create .nc file
        !------------------------------------------------------------
        !
        ! --  output  U, V  fields  from ( h/Vx/Vy/Vz )
        !
        j= int( iter/nstep_ncop )
        do ip=1,nip 
           h_tmp(ip)=FVrk(ip,1)+hb_a(ip)
           !!h_tmp(ip)=FVrk(ip,1)            ! pure h^* 
           !
           P_sph(1)= lon_v(ip)
           P_sph(2)= lat_v(ip) 
           P_sph(3)= 1.d0
           call basis_between_sph_car (P_sph, basis, invbasis, 3)
           X(1)= FVrk(ip,2)
           X(2)= FVrk(ip,3)
           X(3)= FVrk(ip,4)
           call AX_mult(invbasis, X, Y, 3, 3, 3, 3)
           gh_icos(ip,1)= Y(1)      ! u/v/w    re-use the array
           gh_icos(ip,2)= Y(2)      ! u/v/w 
        enddo
        !
        !
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='h'//trim(mystring(j))//'_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        !
        !
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, gh_icos(1,1), g_ll_unif, f_ll)
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='u'//trim(mystring(j))//'_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        !
        !
        if ( iswcase .EQ. 10) then   ! PV for polar vortex
           do ip=1,nip 
              h_tmp(ip)= diag_var(ip,1)
           enddo
           call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
           call change_latlon_array_shape &
                (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
           filename='PV'//trim(mystring(j))//'_field.nc'
           call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        endif


     elseif ( stagger .EQ. 'C' .AND. iter .EQ. iter_ncop ) then
        j= int( iter/nstep_ncop )
        do ip=1,nip 
           h_tmp(ip)=FVrk(ip,1)+hb_a(ip)
        enddo
        !!
        !! need  Nvec, Tvec
        !!do ip=1,NPTS
        !!   !!h_tmp(ip)=FVrk(ip,1)            ! pure h^* 
        !!   !
        !!   P_sph(1)= lon_v(ip)
        !!   P_sph(2)= lat_v(ip) 
        !!   P_sph(3)= 1.d0
        !!   call basis_between_sph_car (P_sph, basis, invbasis, 3)
        !!   X(1)= FVrk(ip,2)
        !!   X(2)= FVrk(ip,3)
        !!   X(3)= FVrk(ip,4)
        !!   call AX_mult(invbasis, X, Y, 3, 3, 3, 3)
        !!   gh_icos(ip,1)= Y(1)      ! u/v/w    re-use the array
        !!   gh_icos(ip,2)= Y(2)      ! u/v/w 
        !!enddo
        !
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        filename='h'//trim(mystring(j))//'_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        !
     elseif ( stagger(1:2) .EQ. 'Du' .AND. iter .EQ. iter_ncop ) then
        j= int( iter/nstep_ncop )
        do iv=1,niV 
           h_tmp(iv)=FVrk(iv,1)+hb_a(iv)
        enddo

        write(6,101) 'a1'
        !
        call icos_2_ll (niV, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
        write(6,101) 'a2'

        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        write(6,101) 'a3'

        filename='h'//trim(mystring(j))//'_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        write(6,101) 'a4'
        !

        if ( iswcase .EQ. 10) then   ! PV for polar vortex
           do ip=1,niX 
              h_tmp(ip)= diag_var(ip,1)
           enddo
           call icos_2_ll (niX, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
           call change_latlon_array_shape &
                (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
           filename='PV'//trim(mystring(j))//'_field.nc'
           call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        endif


     elseif ( stagger(1:2) .EQ. 'ZM' .AND. iter .EQ. iter_ncop ) then
        j= int( iter/nstep_ncop )
        do ip=1,nip 
           h_tmp(ip)=FVrk(ip,1)+hb_a(ip)
        enddo
        call icos_2_ll (nip, np_ll, Tdim, gh_icos_yx, h_tmp, g_ll_unif, f_ll)
        write(6,101) 'a2'
        call change_latlon_array_shape &
             (np_ll, Tdim, nx_ll, ny_ll, g_ll_unif, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
        write(6,101) 'a3'
        !
        filename='h'//trim(mystring(j))//'_field.nc'
        call lon_lat_press_wt_nc4 (Tdim,nx_ll,ny_ll,lon_deg_plot,lat_deg_plot,field_plot,filename)
        write(6,101) 'a4'
        !
     endif
!!     write(6,*) 'af post.3'



     !
     !-- post.4  output  d.h_it.X
     !
     if ( stagger .EQ. 'A' .AND. iter .EQ. iter_ncop ) then
        npos=nip
        gh_icos(1:npos,1)=lon(1:npos); gh_icos(1:npos,2)=lat(1:npos)
        !
        j= int( iter/nstep_ck )
        filename='d.h_it.'//trim(mystring(j))           
        if (iswcase.LE.2 .OR. iswcase.EQ.5 .OR. iswcase.EQ.8) then        
           call output_icos (npts, npts, Tdim, npos, gh_icos, FVrk(1,1), filename, fnlen)
        endif
        !
        ! note  ime_E dim for gh_icos, h_tmp in module_allocate.F90
        ! call output_icos (NPTS, nip, Tdim, npos, gh_icos, FVrk(1,1), filename, fnlen)
        !
     elseif ( stagger(1:2) .EQ. 'Du' .AND. iter .EQ. iter_ncop ) then
        npos=niV
        gh_icos(1:npos,1)=lon(1:npos); gh_icos(1:npos,2)=lat(1:npos)
        !
        j= int( iter/nstep_ck )
        filename='d.h_it.'//trim(mystring(j))           
        if (iswcase.LE.2 .OR. iswcase.EQ.5 .OR. iswcase.EQ.8) then        
           call output_icos (npos, npos, Tdim, npos, gh_icos, FVrk(1,1), filename, fnlen)
        endif
        !
        ! note  ime_E dim for gh_icos, h_tmp in module_allocate.F90
        ! call output_icos (NPTS, nip, Tdim, npos, gh_icos, FVrk(1,1), filename, fnlen)
        !
     endif


!!     write(6,*) 'af post.4'
!!     write(6,*) 'iter', iter

     call flush(6)
     !
     !
  enddo     !  loop iter
  close(7)
  close(iunit)


  write(6,211)  dt*real(iter-1)/86400.d0,Lnorm(1:3,1),Lnorm(1:3,2)  ! time_h,L1n,L2n,Linfn, V_norm
  call datetime
  write (*,'(a)') 'SWM completed.'



  !
  !-----------------------------------------------------
  ! YGYU: testRK4 contains the rest of the code, pass.
  !call testRK4 
  !-----------------------------------------------------
  !call test_sphere
  !


  STOP  'SWM completed.'
  call flush(6)
  include '../swm/myformat.inc'
end program swm_agrid
