!
!- g1: RK4
!- g2: Williamson TC2--TC6 
!- g3: netcdf
!- g4: error norm
!- g5: TE_PE conserve check
!- g6: output lat-lon, h_10.dat
!- g7: diffision, Lap, Lap^2
!
!


!--------------------------------------------------
!    Runge-Kutta
!--------------------------------------------------

!----------------------------------------
! RK 4th order
! dy/dt = f(y, t);  y(t0)=y0
!
! t_n+1 = t_n + h
! y_n+1 = y_n + h/6 * ( k1 + 2k2 + 2k3 + k4)
!
! k1    = f(t_n,       y_n)
! k2    = f(t_n + h/2, y_n + h/2 * k1)
! k3    = f(t_n + h/2, y_n + h/2 * k2)
! k4    = f(t_n + h  , y_n + h   * k3)
!----------------------------------------
!

subroutine RK4th_test_Var_increment (istage,ysv,fj,yi,ndimFV,ims,ime,ips,ipe)
  use kinds, only: rt
  use module_variables,only : stagecoef
  use readnamelist,    only : stagger
  use module_control,  only : nip
  implicit none
  integer, intent(in)  ::  istage,ndimFV,ims,ime,ips,ipe
  real(rt),intent(in)  ::  ysv(ims:ime, ndimFV)  !  save 
  real(rt),intent(in)  ::  fj (ims:ime, ndimFV)  !  j= istage-1   stage
  real(rt),intent(out) ::  yi (ims:ime, ndimFV)  !     istage
  integer  :: ip, j
  !
  if(stagger .EQ. 'A') then
     do ip=ips,ipe
        do j=1, ndimFV
           yi(ip,j)= ysv(ip,j) + stagecoef(istage) * fj(ip,j)
        enddo
     enddo
  elseif(stagger .EQ. 'C') then
     ! h
     j=1
     do ip=1, nip
        yi(ip,j)= ysv(ip,j) + stagecoef(istage) * fj(ip,j)
     enddo
     ! Un
     j=2
     do ip=ips,ipe
        yi(ip,j)= ysv(ip,j) + stagecoef(istage) * fj(ip,j)
     enddo
     !
  elseif(stagger(1:2) .EQ. 'ZM') then
     ! h
     j=1
     do ip=1, nip
        yi(ip,j)= ysv(ip,j) + stagecoef(istage) * fj(ip,j)
     enddo
     ! Vx/Vy/Vz
     do j=2, ndimFV
        do ip=ips,ipe
           yi(ip,j)= ysv(ip,j) + stagecoef(istage) * fj(ip,j)
        enddo
     enddo
     !
  elseif(stagger(1:2) .EQ. 'Du') then
     ! h
     j=1
     do ip=ips,ipe
        yi(ip,j)= ysv(ip,j) + stagecoef(istage) * fj(ip,j)
     enddo
     ! Vx/Vy/Vz
     do j=2, ndimFV
        do ip=1,nip
           yi(ip,j)= ysv(ip,j) + stagecoef(istage) * fj(ip,j)
        enddo
     enddo
  else
     STOP 'stagger not implemented in RK4'
  endif
  return
end subroutine RK4th_test_Var_increment
! RK4th_stage_increment



subroutine RK4th_true_Var_increment (istage,y1,fi,ndimFV,ims,ime,ips,ipe)
  use kinds, only: rt
  use module_variables,only : wtrk
  use readnamelist,    only : stagger
  use module_control,  only : nip
  implicit none
  integer, intent(in)    ::  istage,ndimFV,ims,ime,ips,ipe
  real(rt),intent(inout) ::  y1(ims:ime, ndimFV)        ! y_val
  real(rt),intent(in)    ::  fi(ims:ime, ndimFV)        ! f = dy/dt  at  istage
  integer  :: ip, j
  !
  !write(6,121)  'ips,ipe,ndimFV', ips,ipe,ndimFV

  if(stagger .EQ. 'A') then
     do ip=ips,ipe
        do j=1, ndimFV
           y1(ip,j)= y1(ip,j) + wtrk(istage)*fi(ip,j)
           !        write(6,101)  'y1(ip,j), wtrk(istage), fi(ip,j)', y1(ip,j), wtrk(istage), fi(ip,j)
        enddo
     enddo
  elseif(stagger .EQ. 'C') then
     ! h
     j=1
     do ip=1,nip
        y1(ip,j)= y1(ip,j) + wtrk(istage)*fi(ip,j)
        !        write(6,101)  'y1(ip,j), wtrk(istage), fi(ip,j)', y1(ip,j), wtrk(istage), fi(ip,j)
     enddo
     ! Un
     j=2
     do ip=ips,ipe
        y1(ip,j)= y1(ip,j) + wtrk(istage)*fi(ip,j)
     enddo
  elseif(stagger(1:2) .EQ. 'ZM') then
     ! h
     j=1
     do ip=1,nip
        y1(ip,j)= y1(ip,j) + wtrk(istage)*fi(ip,j)
        !        write(6,101)  'y1(ip,j), wtrk(istage), fi(ip,j)', y1(ip,j), wtrk(istage), fi(ip,j)
     enddo
     ! Vx/Vy/Vz
     do j=2, ndimFV
        do ip=ips,ipe
           y1(ip,j)= y1(ip,j) + wtrk(istage)*fi(ip,j)
        enddo
     enddo
  elseif(stagger(1:2) .EQ. 'Du') then
     ! h
     j=1
     do ip=ips,ipe
        y1(ip,j)= y1(ip,j) + wtrk(istage)*fi(ip,j)
        !        write(6,101)  'y1(ip,j), wtrk(istage), fi(ip,j)', y1(ip,j), wtrk(istage), fi(ip,j)
     enddo
     ! Vx/Vy/Vz
     do j=2, ndimFV
        do ip=ips,ipe
           y1(ip,j)= y1(ip,j) + wtrk(istage)*fi(ip,j)
        enddo
     enddo
  endif

  return
  include '../swm/myformat.inc'
end subroutine RK4th_true_Var_increment
!end subroutine RK4th_time_increment_by_stage




!-- C grid problem
!   h :   ims : ime     ;  I have over-used h(ime+1: imex)=0
!   Un:   imsx: imex
!
subroutine Zerotendencies_sw(ims,ime,ips,ipe,FV,ndimFV)
  use kinds, only: rt
  implicit none
  integer, intent (IN ) :: ims,ime,ips,ipe,ndimFV
  real(rt),intent (OUT) :: FV(ims:ime, ndimFV)    ! forcing of height (h)
  integer               :: ip, j
  do ip=ips,ipe
     do j=1, ndimFV
        FV(ip,j) = 0.d0
     enddo
  enddo
  return
end subroutine Zerotendencies_sw


subroutine Save_State_Var_sw(ims,ime,ips,ipe,FV,FVs,ndimFV)
  use kinds, only: rt
  implicit none
  integer, intent(IN ):: ims,ime,ips,ipe,ndimFV
  real(rt),intent(IN )::  FV   (ims:ime, ndimFV)     ! 
  real(rt),intent(OUT)::  FVs  (ims:ime, ndimFV)     ! save a copy
  integer :: ip, j
  do ip=ips,ipe
     do j=1, ndimFV
        FVs(ip,j)=  FV(ip,j)
     enddo
  enddo
  return
end subroutine Save_State_Var_sw



!----------------------
!
!  dy/dt = -y ,   y(0)=A
!      y = A e^(-t)
!
!  dh /dt = -h
!  dhu/dt = -2 hu
!  dhv/dt = -2 hv
!  h0=u0=v0=1       h, V = exp(-t)
!----------------------

!subroutine testRK4
!  use kinds, only             :  rt
!  use module_control
!  use module_variables
!  use module_constants
!  implicit none
!  integer :: ret, fnlen, ndim2, LDA, LDB, nstep, istage
!  character (len=filenamelen) :: filename
!  character (len=20)  :: mystring    ! fixed extern func
!  !
!  fnlen=filenamelen
!  ndim2=2
!
!  hrk (1:nip)=1.d0
!  u_a (1:nip)=1.d0
!  v_a (1:nip)=1.d0
!  uhrk(1:nip)=1.d0
!  vhrk(1:nip)=1.d0
!
!  call gen_uniform_ll   (nx_ll,ny_ll,g_ll_unif, np_ll, ndim2)
!  g_icos(ims:ime,1)=lon(ims:ime); g_icos(ims:ime,2)=lat(ims:ime); filename='d.h_it.0'
!  !call output_icos_2_ll (nip, np_ll, ndim2, g_icos, hrk, g_ll_unif, f_ll, filename, fnlen)           
!  call output_icos (nip, np_ll, ndim2, g_icos, hrk, filename, fnlen)
!  !                                                                                                  
!  do its=itsbeg,itsend
!     write(6, '(/,2x,a,i6)') 'its=',its
!     call Zerotendencies_sw (ims,ime,ips,ipe,fh_a,fuh_a,fvh_a)
!     call Save_State_Var_sw (ims,ime,ips,ipe,hrk,uhrk,vhrk,h0_a,uh0_a,vh0_a)  !  updated i-time RK   
!     call huv2uv (hrk,uhrk,vhrk,u_a,v_a,ims,ime,ips,ipe)
!
!     filename='d.u_it.'//mystring(its)
!     call output_icos (nip, np_ll, ndim2, g_icos, u_a, filename, fnlen)
!
!     !call iter_icos_forcing_one_third_rule_velocity (ims,ime,ips,ipe,hrk,u_a,v_a,fh_a,fuh_a,fvh_a)
!     call test_exp_forcing (ims,ime,ips,ipe,hrk,u_a,v_a,fh_a,fuh_a,fvh_a)
!
!     !add                                                                                            
!     call RK4th_time_increment_by_stage (1,hrk,uhrk,vhrk,fh_a,fuh_a,fvh_a,ims,ime,ips,ipe)
!     do istage=2,nrkl
!        call RK4th_stage_increment &
!             (istage,h0_a,uh0_a,vh0_a,fh_a,fuh_a,fvh_a,h_a,uh_a,vh_a,ims,ime,ips,ipe)
!        call huv2uv (h_a,uh_a,vh_a,u_a,v_a,ims,ime,ips,ipe)
!!        call iter_icos_forcing_one_third_rule_velocity (ims,ime,ips,ipe,h_a,u_a,v_a,fh_a,fuh_a,fvh_a)
!        call test_exp_forcing (ims,ime,ips,ipe,hrk,u_a,v_a,fh_a,fuh_a,fvh_a)
!        !add                                                                                         
!        call RK4th_time_increment_by_stage (istage,hrk,uhrk,vhrk,fh_a,fuh_a,fvh_a,ims,ime,ips,ipe)
!     enddo
!     !                                                                                               
!     !                                                                                               
!     nstep=1
!     if (mod(its, nstep).eq.0) then
!        filename='d.h_it.'//mystring(its)
!        call output_icos (nip, np_ll, ndim2, g_icos, hrk, filename, fnlen)
!        !call output_icos_2_ll (nip, np_ll, ndim2, g_icos, hrk, g_ll_unif, f_ll, filename, fnlen)     
!     endif
!     !                                                                                               
!     !   Euler one-step                                                                              
!     !   call onestep(ims,ime,ips,ipe,fh_a,fuh_a,fvh_a,h0_a,uh0_a,vh0_a,h_a,uh_a,vh_a,u_a,v_a,dt)    
!  enddo
!
!  return
!end subroutine testRK4
!



!--------------------------------------------------
!
!  initial conditions  +  Williamson test cases
!
!--------------------------------------------------
subroutine simple_init
  use module_constants, only: rt, omega, hbc
  use module_control,   only: glvl, nip, niE, niV, NPTS, dt, ips,ipe,ims,ime
  use module_variables, only: fcori, fcori_V, wtrk, stagecoef, hb_a
  use module_variables, only: lat, lon, latE, lonE, lat_v, lon_v, &
                              wtsph,slen,area,Nvec,Tvec,Rvec,Rcvec,RVvec,dcenter,&
                              Nvec2,Tvec2,Rvec2,areaV_A,dcenter_A,dx_ov_le,lambdas,&
                              Nvec3,Tvec3,Rvec3,slen3,wtTBsph,&
                              RAiv4T_A,RAiv4H_A,wt2sph,wt3sph,&
                              C4E,E4C,wtUt,Eindex,EindexVort,areaV,areaE,areaI,LevdUn,&
                              NoutdUn,arctheta,Edc,Vdc,C4V,V4C,&
                              RAiv4T,RAiv4H,iVtopkey,V_sph_grid,vertex_hex_sph,TE0,PE0
  use module_variables, only: nprox, prox, proxs
  use readnamelist,     only: lgen_sph_var, stagger
  implicit none
  ! loc
  integer :: ip, is, iE, ipt, iV
  integer :: ip2, is2, jp2, js2
  character(len=20) :: filename
  character(len=150) :: line, line2
  real (rt), allocatable :: icos_grid(:,:)
  real (rt), allocatable :: midE_grid(:,:)  


  !
  !-- gr.1  generate sph_var_GX.dat
  !                  overwrite  Amtx.dat data
  !
  if (lgen_sph_var) then
     allocate(icos_grid(2,nip),midE_grid(2,niE))
     write(filename, '(a,i2.2,a)') 'sph_var_G',glvl,'.dat' 
     open (11, file=trim(filename), status='unknown', form='unformatted')
     read (11) &
          icos_grid(1:2,1:nip), wtsph(1:6,1:6,1:nip), slen(1:6,1:nip), area(1:nip),&
          Nvec(1:3,1:6,1:nip), Tvec(1:3,1:6,1:nip), Rvec(1:3,1:6,1:nip), Rcvec(1:3,1:nip), &
          Nvec2(1:3,1:6,1:nip), Tvec2(1:3,1:6,1:nip), Rvec2(1:3,1:6,1:nip),&
          areaV_A(1:6,1:nip),dcenter_A(1:6,1:nip),dx_ov_le(1:niE),lambdas(1:6,1:nip),&
          RAiv4T_A(1:6,1:nip),RAiv4H_A(1:6,1:nip),&
          nprox(1:nip), prox(1:6,1:nip), proxs(1:6,1:nip), &
          midE_grid(1:2,1:niE), dcenter(1:niE), C4E(1:5,1:niE), E4C(1:6,1:nip), &
          wtUt(1:10,1:niE), Eindex(1:10,1:niE), areaV(1:2,1:niE), areaE(1:niE), areaI(1:nip), &
          LevdUn(1:6,1:niE), NoutdUn(1:6,1:nip), EindexVort(1:6,1:niE), arctheta(1:niE), &
          Edc(1:6,1:nip), Vdc(1:6,1:nip), C4V(1:5,1:niV), V4C(1:6,1:nip), &
          RAiv4T(1:6,1:niE), RAiv4H(1:6,1:nip), iVtopkey(1:6,1:nip), &
          V_sph_grid(1:2,1:2,1:niE), vertex_hex_sph(1:2,1:6,1:nip), RVvec(1:3,1:niV),&
          wt2sph(1:3,1:6,1:nip), wt3sph(1:3,1:6,1:nip), &
          Nvec3(1:3,1:3,1:niV),Tvec3(1:3,1:3,1:niV),Rvec3(1:3,1:3,1:niV),&
          slen3(1:3,1:niV),wtTBsph(1:3,1:niV)


     if( stagger.EQ.'A') then
        do ip=1, nip                ! for h, V
           lon(ip)=icos_grid(1,ip)  ! lambda
           lat(ip)=icos_grid(2,ip)  ! theta
           lon_v(ip)=lon(ip)
           lat_v(ip)=lat(ip)
        enddo
     elseif (stagger.EQ.'C') then
        do ip=1, nip                ! for h
           lon(ip)=icos_grid(1,ip)  ! lambda
           lat(ip)=icos_grid(2,ip)  ! theta
        enddo
        do iE=1, NPTS              ! for velocity   C: NPTS=niE
           lon_v(iE)=midE_grid(1,iE)
           lat_v(iE)=midE_grid(2,iE)
        enddo
     elseif (stagger(1:2).EQ.'ZM') then
        do ip=1, nip                ! for h
           lon(ip)=icos_grid(1,ip)  ! lambda
           lat(ip)=icos_grid(2,ip)  ! theta
        enddo
        do iv=1, NPTS              ! for Velocity   ZM: NPTS=niV
           ip=C4V(1,iv)
           is=C4V(2,iv)
           lon_v(iv)=vertex_hex_sph(1,is,ip)  ! lambda
           lat_v(iv)=vertex_hex_sph(2,is,ip)  ! theta   
        enddo
     elseif (stagger(1:2).EQ.'Du') then  ! Du-grid
        do iv=1, niV                ! for h
           ip=C4V(1,iv)
           is=C4V(2,iv)
           lon(iv)=vertex_hex_sph(1,is,ip)  ! lambda
           lat(iv)=vertex_hex_sph(2,is,ip)  ! theta
        enddo
        do ip=1, NPTS                  ! for Velocity  Du: NPTS=nip
           lon_v(ip)=icos_grid(1,ip)  ! lambda
           lat_v(ip)=icos_grid(2,ip)  ! theta
        enddo
     endif



     write(6, 101) 'test Nvec2'
     ipt= nip * 0.65
     write(6, 121) 'ipt=          ', ip
     write(6, 101) 'Nvec(3,5,ipt) ', Nvec(1:3,5,ipt)
     write(6, 101) 'Tvec(3,5,ipt) ', Tvec(1:3,5,ipt)
     write(6, 101) 'Rvec(3,5,ipt) ', Rvec(1:3,5,ipt)
     write(6, 101) 'lambdas( ipt) ', lambdas (1:5,ipt)
     write(6, 101) 'wt  (6,5,ipt) ', wtsph (1:4,5,ipt)
     write(6, 101) 'wt2 (3,5,ipt) ', wt2sph(1:3,5,ipt)
     write(6, 101) 'wt3 (3,5,ipt) ', wt3sph(1:3,5,ipt)

     write(6, 101) 'Nvec2(3,5,ipt)', Nvec2(1:3,5,ipt)
     write(6, 101) 'Tvec2(3,5,ipt)', Tvec2(1:3,5,ipt)
     write(6, 101) 'Rvec2(3,5,ipt)', Rvec2(1:3,5,ipt)
     write(6, 101) 'wt2sph(1:3,5,ipt)', wt2sph(1:3,5,ipt)
     write(6, 101) 'dx_ov_le(.75*niE)', dx_ov_le(0.75*niE)
     write(6, 112) 'RAiv4H_A(1:6,ipt)', RAiv4H_A(1:6,ipt)

     write(6, 101) 'Nvec3(3,3,ipt)', Nvec3(1:3,3,ipt)
     write(6, 101) 'Tvec3(3,3,ipt)', Tvec3(1:3,3,ipt)
     write(6, 101) 'Rvec3(3,3,ipt)', Rvec3(1:3,3,ipt)
     write(6, 101) 'slen3(1:3,ipt)', slen3(1:3,ipt)
     write(6, 101) 'wtTBsph(1:3,ipt)', wtTBsph(1:3,ipt)


     write(6, 121) 'ipt           ', 1, 2
     do ipt=1,2
        write(6, 201) 'dcenter_A(1:6,ipt)',dcenter_A(1:6,ipt)
     enddo
     do ipt=1,2
        write(6, 201) 'areaV_A(1:6,  ipt)',areaV_A(1:6,ipt)
     enddo



!     !
!     !-- check
!     !
!     do ip=1, nip, 1
!        write(6, 121) 'ip                   =', ip
!        write(6, 101) 'lon(ip),  lat(ip)    =', lon(ip), lat(ip)
!        write(6, 111) 'area(ip), slen(1,ip) =', area(ip), slen(1:6,ip)
!        write(6, 101) 'Nvec,Tvec,Rvec(1,1,ip)', Nvec(1,1,ip),Tvec(2,5,ip),Rvec(3,4,ip)
!        write(6, 101) 'Rcvec(2,ip)           ', Rcvec(2,ip)
!        write(6, 121) 'proxs(2,ip)          =', proxs(2,ip)
!        write(6, 121) 'E4C(1:6,ip)          =', E4C(1:6,ip)  
!        write(6, 121) 'nprox(ip)            =', nprox(ip)
!        write(line2, '(2x,a)') 'C4E(1/5,iE=<1,is>)    ='
!        do is=1, nprox(ip)
!           iE=E4C(is,ip)
!           ip2=C4E(1,iE); is2=C4E(5,iE)         ! Left
!           jp2=C4E(2,iE); js2=proxs(is2,ip2)    ! Right
!           if (Edc(is,ip).EQ.0) then            ! 0: Left, 1: Right
!              write(line, '(2(2x,i8))') ip2,is2
!           else
!              write(line, '(2(2x,i8))') jp2,js2
!           endif
!           line2=trim(line2)//trim(line)
!        enddo
!        write(6,'(a)')  trim(line2)
!        write(6,'(/)')
!     enddo
!     write(6,'(//)')
!     !
!     !
!     do iE=1, niE, 1
!        write(6, 121) 'iE                   =', iE
!        write(6, 101) 'lonE(iE), latE(iE)   =', lonE(iE), latE(iE)
!        write(6, 111) 'dcenter(iE)          =', dcenter(iE)
!        write(6, 121) 'C4E(1:5,iE)          =', C4E(1:5,iE)
!        write(6, 112) 'wtUt(1:10,iE)        =', wtUt(1:10,iE)
!        write(6, 121) 'Eindex(1:10,iE)      =', Eindex(1:10,iE)
!        write(6, 111) 'areaE(iE)            =', areaE(iE)
!        write(6, 111) 'LevdUn(1:6,iE)       =', LevdUn(1:6,iE)
!        write(6,'(/)')
!     enddo
!     STOP 'ck with MPAS grid'
!


!     do iV=1, niV, 1
!        write(6, 121) 'iv,ip,is          =', iv, C4V(1,iv), C4V(2,iv)
!        write(6, 121) 'Vdc(is,ip)        =', Vdc(is,ip)
!        write(6, 101) 'lon(iv), lat(iv)  =', lon(iv), lat(iv)
!        write(6,'(/)')
!     enddo


!
!--  from icos.F90
!     write (11) icos_grid(1:2,1:nip)*d2r, wt(1:3,1:6,1:nip), slen(1:6,1:nip), &
!          Nvec(1:3,1:6,1:nip), Tvec(1:3,1:6,1:nip), area(1:nip), &
!          nprox(1:nip), prox(1:6,1:nip), proxs(1:6,1:nip)     
     write(6, 102) 'reading filename='//trim(filename)//' completed'
     write(6, '(//)')
     deallocate(icos_grid,midE_grid)
  endif      ! if (lgen_sph_var)
  !
  !
  if (stagger(1:2).EQ.'Du') then
     do ipt=1, niV                              ! h only on icos centers
        hb_a(ipt)= hbc                          ! background constant
     enddo
  else
     do ipt=1, nip                              ! h only on icos centers
        hb_a(ipt)= hbc                          ! background constant
     enddo
  endif
  !
  !-- Runge Kutta 4th order
  !
  wtrk(1)=dt/6.d0; wtrk(2)=dt/3.d0 
  wtrk(3)=dt/3.d0; wtrk(4)=dt/6.d0
  stagecoef(1)=0.d0; stagecoef(2)=dt/2.d0 
  stagecoef(3)=dt/2.d0; stagecoef(4)=dt
  !
  !
  TE0=-1.d10;  PE0=-1.d20    ! init
  !
  return
  include '../swm/myformat.inc'
end subroutine simple_init






      




!--------------------------------------------------
!  note:    Two set of grid
!           NPTS (lon_v,lat_v) on A/C general grid for velocity field
!           nip  (lon_h, lat_h) on A-grid for h-field
!
!  input:   lon(nip), lat(nip)  from icos after call read_permute
!           \lamda    \theta
!          
!  output:  (h,u,v,hu,hv) on icos grid
!
!  note:     ht :: h true
!  
!            I have also generalized this case with Dual grid layout
!                   watchout  h at vertex niV
!
!--------------------------------------------------


! TC8: a thin layer of TC2
! TC9: a constant thin layer, only evolve h-field
!      no change in velocity   \p V / p t == 0 , fixed
! TC10: Polar vortex erosion
!       Jukes et al. nature 1987
!       Bates and Li (19950  piece wise linear

subroutine sw_test_init (NPTS_MX,NPTS,lon_v,lat_v,nih,lon_h,lat_h,alpha,hb_a,FV,fcori,ndimFV,iswcase,stagger)
  use kinds, only: rt
  use module_constants, only : pi, ae, omega, g, polar_tol, len_mx
  use module_variables, only : omegasw, C4E, Nvec, Tvec, Rcvec, RVvec, Rvec, NoutdUn, nprox, proxs
  use module_variables, only : Eindex, dcenter, slen, wtUt, V_sph_grid, fcori_V
  use module_variables, only : ASV_ana             ! for alpha \grad beta   NICOM test
  use module_variables, only : nlon_st, nlat_st, glon_st, glat_st, h_st, U_st, V_st 

  use module_variables, only : grid_gauss, gh_icos_yx, grid_general  ! spectra  
  use readnamelist,     only : nx_NICAM, mx_NICAM, file_nc_init
  use slint, only       :  slint_init, bilinear_interp
  implicit none
  !
  integer, intent (in)  :: NPTS_MX, NPTS, nih, ndimFV, iswcase
  real(rt),intent (in ) :: lon_v(NPTS), lat_v(NPTS)      ! lon/lat for velocity
  real(rt),intent (out) :: hb_a(NPTS_MX)                 ! hb_a(nih)=surface/background, 4 Dual-grid special
  real(rt),intent (in ) :: lon_h(nih), lat_h(nih)        ! lon/lat for h_field
  real(rt),intent (in ) :: alpha
  real(rt),intent (out) :: FV(NPTS_MX,ndimFV)            ! Four-vector in cartesian
  real(rt),intent (out) :: fcori(NPTS)                   ! coriorios force
  character(len=len_mx), intent(in) :: stagger                ! hardcoded
  real(rt) ::  u(NPTS), v(NPTS)                          ! spherical coord
  !!real(rt) :: h(nih), uh(nih), vh(nih)                 ! spherical coord
  !
  !loc
  real(rt) :: h0, R0    ! specified const
  real(rt) :: u0, u0x, h0x  !  u0x, h0x for TC0
  real(rt) :: lat, lon, latc, lonc
  real(rt) :: r, ht, hs, ut, vt, s, s2, s3, beta, c  ! tmp
  real(rt) :: P_sph(3), P_car(3), basis(3,3), invbasis(3,3), basis3(3,3)
  real(rt) :: Euler(3,3), P1_mat(3,3), P2_mat(3,3)   !  Euler solid angle rotation
  real(rt) :: X(3), Y(3), Y3(3)
  real(rt) :: vec1(3), vec2(3), vec3(3), V_car(4)    ! V_car(4) --> FV(nih,4)
  real(rt) :: u_lat, u_lon, vort
  real(rt) :: gradh_lat, gradh_lon
  real(rt) :: omega_RH,K_RH,R_RH,A_RH,B_RH,C_RH,w,k
  real(rt) :: t, thetap, lambdap, d2r
  real(rt) :: xe, thetab, thetae   ! TC3
  real(rt) :: dx, up, vp, xp   ! TC3
  real(rt) :: yp1, ypn
  real(rt) :: A, B, tx
  real(rt), allocatable :: latp_sample(:), hp(:)     ! TC3 for h
  real(rt), allocatable :: latp(:), lonp(:), wksp(:) ! TC3 for h and velocity
  real(rt), allocatable :: y2(:), work(:), work2(:)  ! spline fit
  real(rt), allocatable :: wksp1(:), wksp2(:), wksp3(:)  ! spline fit
  real(rt), allocatable :: FV_aux(:,:), FV_aux_t(:,:,:,:)
  real(rt), allocatable :: icos_grid_aux(:,:)
  !
  !
  integer  :: ip, ipt, is, nx, ny, i, j, mx, it
  integer  :: m, n, mp, np, ndim
  integer  :: Nsample, Nsample3, ntmp
  integer  :: LDA, LDB, LDT
  integer  :: ims, ime, ips, ipe, itstep, istage
  integer  :: iE, iE2, ipc, isc, is2, ism, jp, js, nb, nvert, niE_loc
  integer, allocatable :: iwksp(:)
  !
  character (len=150) :: filename
  character (len=50)  :: keyword
  logical  :: fexist
  integer  :: ntime, itime
  !
  ! dimension: meter, m/s
  !
  R0=      1.d0
  latc=    0.d0
  lonc=   -pi/2.d0
  d2r=   1.d0/180.d0*pi


  if (iswcase.EQ.2 .OR. iswcase.EQ.8 ) then
     omegasw= 2.d0*pi/(12.d0*86400.d0) ! 2pi/T
     u0=  ae*omegasw                   ! ~ 38 m/s
  elseif (iswcase.EQ.4) then
     u0=  20.d0             ! 20 or 40 m/s 
     omegasw= u0/ae         ! redefine
  elseif (iswcase.EQ.5) then
     u0=  20.d0             ! 20m/s
     omegasw= u0/ae
  else
     omegasw= 2.d0*pi/(12.d0*86400.d0) ! 2pi/T
     u0=  ae*omegasw
  endif
  write(6,201) 'omega, omegasw', omega, omegasw
  !
  !
  ndim=    3                        ! fix
  nx=      nx_NICAM
  mx=      mx_NICAM
  omega_RH=7.848d-6                 ! s^-1
  K_RH    =7.848d-6                 ! s^-1
  R_RH    =4.d0                     !
  thetab=-pi/6.d0                   ! TC3
  thetae= pi/2.d0
  xe= 0.3d0
  !
  !-- s1:  h field defined at center (lat_h,lon_h) for both A and C grid
  !
  do ip=1, nih
     lat= lat_h(ip)
     lon= lon_h(ip)
     if (iswcase.EQ.1) then
        h0=1000.d0
        r  = acos( sin(lat)*sin(latc) + cos(lat)*cos(latc)*cos(lon-lonc) )
        if (r .LT. R0) then
           ht= 0.5d0 * h0 * (1.d0 + cos(pi*r/R0))
        else
           ht= 0.d0          ! original version
           !! ht= 0.5d0 * h0     ! MPAS version
           !ht= 0.1d0 * h0
        endif
     elseif (iswcase .EQ. 2  .OR. iswcase.EQ.3 ) then
        !
        ! h0 is 3000 m in case-2
        !
        h0= 3000.d0
        ! cos\Phi = omega \dot r = cos(delta)
        c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat)   
        ht = h0 - (0.5d0*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g
     elseif (iswcase .EQ. 8 ) then
        !
        ! TC8: a thin layer of TC2
        !
        ! cos\Phi = omega \dot r = cos(delta)
        c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat)   
        hb_a(ip) =  - (0.5d0*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g
        ht =  100.d0    ! 100 meter
!        ht =  1.d0    ! 1 meter
!        ht =  0.01d0   ! 0.01 meter
     elseif (iswcase .EQ. 9 ) then
        !
        ! TC9: a constant thin layer, only evolve h-field
        !      no change in velocity   \p V / p t == 0 , fixed
        !
        ! cos\Phi = omega \dot r = cos(delta)
        hb_a(ip) = 0.d0
!        ht =  100.d0    ! 100 meter
!        ht =  1.d0    ! 1 meter
        ht =  0.01d0   ! 0.01 meter
     elseif (iswcase .EQ. 5) then
        !
        ! h0 is 5960 m in case-5
        !
        h0= 5960.d0
        ! cos\Phi = omega \dot r = cos(delta)
        c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat)
        ht = h0 - (0.5d0*omegasw*omegasw+omegasw*omega)*ae*ae*c*c/g
        !
        !  define surface height
        !
        R0 = pi/9.d0
        lonc= 1.5d0*pi
        latc= pi/6.d0
        r  = min(R0, sqrt((lon-lonc)**2.d0+(lat-latc)**2.d0))
        hs = 2000.d0
        hb_a(ip)= hs * (1.d0 - r/R0)   ! surface height is defined only in TC5
        !!
        !!  strange, donot substract
        ht= ht - hb_a(ip)
        !
     elseif (iswcase .EQ. 6) then
        ! Rossby-Haurwitz
        h0= 8.d3
        w = omega_RH
        k = K_RH
        R = R_RH
!
!        A_RH = 0.5d0 * w * (2.d0 * omega + w) * cos(lat)**2.d0 + & 0.25d0 * K**2.d0 * &
!          cos(lat)**(2.d0*R) * ((R+1.d0)*cos(lat)**2.d0 + 2.d0*R**2.d0 - R - 2.d0 - &
!          2.d0*R**2.d0 * cos(lat)**(-2.d0))
!
        A_RH = 0.5d0 * w * (2.d0 * omega + w) * cos(lat)*cos(lat)  + & 
             0.25d0 * K*K * cos(lat)**(2.d0*R) * ((R+1.d0)*cos(lat)*cos(lat) + 2.d0*R*R - R - 2.d0) - &
             0.5d0  * K*K * cos(lat)**(2.d0*R-2.d0) * R*R
        !
        B_RH = (2.d0*(omega + w)*K / ((R+1.d0)*(R+2.d0))) * cos(lat)**R * ((R*R + &
             2.d0*R + 2.d0) - ((R+1.d0)*cos(lat))**2.0)
        !
        C_RH = 0.25d0 * K*K * cos(lat)**(2.d0*R) * ((R+1.d0)*cos(lat)*cos(lat) - R - 2.d0)
        !
        ht = h0 + ae*ae/g*(A_RH + B_RH*cos(R*lon) + C_RH*cos(2.d0*R*lon))
!        ht = h0 + h0*cos(lat)**2.d0
     else
        ht=0.d0
     endif
     !
     !-- finally specify h field
     !
     FV(ip,1)=ht
     !
  enddo


  !
  !-- set terrain surface field == 0 for other cases
  !
  if (iswcase .NE. 5 .AND. iswcase .NE. 8) then
     do ip=1, nih
        hb_a(ip)= 0.d0
     enddo
  endif
  if (iswcase .EQ. 10) then
     tx=0.d0         ! time 
     call TC10_hs_t(tx, nih, lon_h, lat_h, hb_a(1))  ! T-dep terrain
     !
     !-- test each ingredient in TC10
     !
!     nx=301
!     dx= pi/dble(nx-1)
!     do j=1, nx
!        xp= -pi/2.d0 + dx*dble(j-1)
!        call TC10_func_B(xp, s)
!        call TC10_func_U(xp, s2)
!        write(7, 103) xp, s2, s
!     enddo
!!
!     nx=1001
!     dx= 124.d0/dble(nx-1)
!     do j=1, nx
!        xp=  dx*dble(j-1)
!        call TC10_func_A(xp, s)
!        write(7, 103) xp, s
!     enddo
  endif


  ! 
  !note: in C-grid, we will not touch FV(nih+1:niE,1) --- the h-field
  ! 


  !
  !-- s2:  velocity field
  !
  do ipt=1, NPTS
     lat= lat_v(ipt)     ! lat velocity
     lon= lon_v(ipt)     !
     if (iswcase.EQ.1 .OR. iswcase.EQ.2 .OR. iswcase.EQ.5 &
          .OR. iswcase.EQ.8 .OR. iswcase.EQ.9) then
        !
        ! geostropic flow
        ! V = Omega \cross R
        !
        ut =  u0*(sin(alpha)*sin(lat)*cos(lon) + cos(alpha)*cos(lat))
        vt = -u0*sin(alpha)*sin(lon)
     elseif (iswcase .EQ. 3) then
        ut =  u0*(cos(lat)*cos(alpha) + sin(lat)*sin(alpha)*cos(lon))
        vt = -u0*sin(alpha)*sin(lon)
     elseif (iswcase .EQ. 6) then
        ut =  ae*w*cos(lat) + ae*K*cos(lat)**(R-1.d0)*(R*sin(lat)**2.d0 - cos(lat)**2.d0)*cos(R*lon)
        vt = -ae*K*R*cos(lat)**(R-1.d0)*sin(lat)*sin(R*lon)
     elseif (iswcase .EQ. 10) then  ! Polar vortex
        call TC10_func_U(lat, ut)
        vt = 0.d0
     else
        ut = 0.d0
        vt = 0.d0
     endif
     u(ipt)=ut; v(ipt)=vt
  enddo



  !
  !-- s3: 
  !--      Coriolis (general for both A and C), unified
  !        note velocity and acceleration 
  !        for A-grid at center, C-grid on edge
  !
  do ipt=1, NPTS
     lat= lat_v(ipt)
     lon= lon_v(ipt)
     ! unit omega \dot r = cos(delta)
     c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat) 
     if (iswcase .EQ. 1) then
        fcori(ipt)= 0.d0           ! pure advection
     elseif (iswcase.GE.2 .OR. iswcase.LE.8) then    ! case-3 included
        fcori(ipt)= 2.d0*omega*c   ! 2 * Omega * (\Ome \dot R)
                                   ! global for Earth
     else
        fcori(ipt)= 2.d0*omega*c
     endif
  enddo
  !
  !  C-grid Ringler's symmetric method, fcori_V (1:2, niE)  up and dn
  !
  if (stagger .EQ. 'C') then 
     do ipt=1, NPTS                                             ! iE
        !
        !
        lon= V_sph_grid(1,1,ipt)         ! lam           head point: up
        lat= V_sph_grid(2,1,ipt)         ! theta
        ! unit omega \dot r = cos(delta)
        c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat) 
        fcori_V(1, ipt)= 2.d0*omega*c   !  2 * Omega * (\Ome \dot R)
        !
        !
        lon= V_sph_grid(1,2,ipt)         ! lam           tail point: down
        lat= V_sph_grid(2,2,ipt)         ! theta
        ! unit omega \dot r = cos(delta)
        c  = -sin(alpha)*cos(lat)*cos(lon) + cos(alpha)*sin(lat) 
        fcori_V(2, ipt)= 2.d0*omega*c   !  2 * Omega * (\Ome \dot R)
        !
     enddo
  endif



  !
  !-- s4.1  case3 zonal geostrophic flow with compact support
  !
  !- test u':= up  u_prime function
  !  ntmp=101
  !  dx= pi/dble(ntmp-1)
  !  write(8, '(a)') '# xp, up(xp)'
  !  do i=1, ntmp
  !     xp= -pi/2.d0 + dx*dble(i-1)
  !     call TC3_up(u0, xe, thetab, thetae, xp, up)
  !     write(8, '(2x,2f25.3)')  xp, up 
  !  enddo
  !
  if (iswcase .EQ. 3) then
     !
     ! Euler angle rotation
     !
     !------------------------------------------------------------------------------------
     !     Vector = (e_lam, e_the, e_r)*(U,V,W)^T  = (e_lam', e_the', e_r')*(U', V', W')^T
     !              (ex, ey, ez)  P1 * X  =  (e1', e2', e3') * P2 * X'
     !              (e1, e2, e3)  P1 * X  =  (e1, e2, e3) * Euler * P2 * X'
     !                                 X  =  P1^T * Euler * P2 * X'
     ! 
     !                                 X  =  (u,v,w)^T,    let Y=X'
     !------------------------------------------------------------------------------------
     !
     Euler(1:3,1:3)=0.d0
     Euler(1,1)= cos(alpha)
     Euler(2,2)= 1.d0
     Euler(3,3)= cos(alpha)
     Euler(1,3)=-sin(alpha)
     Euler(3,1)= sin(alpha)
     
     Nsample=1001
     LDA=Nsample; LDB=nih
     allocate (latp_sample(Nsample),hp(Nsample))
     allocate (latp(nih), lonp(nih), iwksp(nih), wksp(nih))  ! for h field
     allocate (y2(LDA), work(LDA), work2(LDB))


     ! prepare table for h field
     !
     call TC3_hp(h0, u0, xe, thetab, thetae, Nsample, latp_sample, hp)

     !  (lat,lon) --> (latp, lonp)
     do ip=1, nih
        lat= lat_h(ip)
        lon= lon_h(ip)        !
        ! zonal flow with compact support
        ! lambda=lon
        ! theta=lat
        ! t:= sin(theta)
        !
        t=sin(lat)*cos(alpha) - cos(lat)*cos(lon)*sin(alpha)
        thetap=asin(t)
        if ( abs(abs(t)-1.d0).LT.1.d-10) then  ! N/S pole
           lambdap=0.d0   ! give any value for lon
        else
           s= asin( cos(lat)*sin(lon) / cos(thetap) )   !    s:= lambda'
           xp= cos(lat)*cos(lon)*cos(alpha)+sin(lat)*sin(alpha)   ! xp:= theta_prime
           if(xp.ge.0.d0) then
              lambdap=s
           else
              lambdap=pi-s
           endif
        endif
        latp(ip)=thetap     ! lat_prime
        lonp(ip)=lambdap    ! lon_prime
        !
        ! ck
        !write(7, 103) lat, lon, latp(ip), lonp(ip)
     enddo



     call indexx (nih,latp,iwksp)   ! sort in ascending order
     do j=1,nih
        wksp(j)=latp(iwksp(j))
     enddo
     !
     ! cube spline interpolation
     !
     yp1=0.d0; ypn=0.d0
     call sple_ab(latp_sample, hp, LDA, yp1, ypn, y2, work)
     call splt_ab(LDA, latp_sample, hp, y2, LDB, wksp, work2)
     !
     !open(7, file='hp.dat', status='unknown')
     !open(8, file='fit_hp.dat', status='unknown')
     !do j=1,LDA
     !   write(7, 103)  latp_sample(j), hp(j), y2(j)
     !enddo
     !do j=1,LDB
     !   write(8, 103)  wksp(j), work2(j)
     !enddo
     !close(7)
     !close(8)
     !!
     !stop 'ck latp lonp'
     !
     do j=1, nih
        ip=iwksp(j)
        FV(ip,1)=work2(j)
     enddo
!     do ip=1, nih
!        write(6,141) 'ip, lat, lon, h', ip,lat_h(ip),lon_h(ip),FV(ip,1)
!     enddo
!     STOP 'ck h field'



     !! pass the test using Euler rotation  lat,lon  <-->  latp,lonp  correct!
     !do ip=1, nih
     !   P_sph(1)=lon_h(ip); P_sph(2)=lat_h(ip); P_sph(3)=1.d0
     !   call vect_from_sph_2_car(P_sph, X, ndim)
     !   P_sph(1)=lonp(ip); P_sph(2)=latp(ip); P_sph(3)=1.d0
     !   call vect_from_sph_2_car(P_sph, Y, ndim)
     !   call AX_mult (Euler, Y, vec1, 3, 3, 3, 3)
     !   X(1:3)=X(1:3)-vec1(1:3)
     !   call X_2norm(X, s, ndim)
     !   write (6, 101) ' ||X- PY||_2 =', s
     !enddo
     deallocate(latp_sample, hp, latp, lonp, wksp, iwksp, y2, work, work2)


     !
     !
     !-- velocity field
     !
     !
     allocate(latp(NPTS), lonp(NPTS), iwksp(NPTS), wksp(NPTS))  ! for h field
     do ip=1, NPTS
        lat= lat_v(ip)
        lon= lon_v(ip)        !
        ! zonal flow with compact support
        ! lambda=lon
        ! theta=lat
        ! t:= sin(theta)
        !
        t=sin(lat)*cos(alpha) - cos(lat)*cos(lon)*sin(alpha)
        thetap=asin(t)
        if ( abs(abs(t)-1.d0).LT.1.d-10) then  ! N/S pole
           lambdap=0.d0   ! give any value for lon
        else
           s= asin( cos(lat)*sin(lon) / cos(thetap) )   !    s:= lambda'
           xp= cos(lat)*cos(lon)*cos(alpha)+sin(lat)*sin(alpha)   ! xp:= theta_prime
           if(xp.ge.0.d0) then
              lambdap=s
           else
              lambdap=pi-s
           endif
        endif
        latp(ip)=thetap     ! lat_prime
        lonp(ip)=lambdap    ! lon_prime
     enddo
     !
     !
     do ip=1, NPTS
        P_sph(1)= lon_v(ip)
        P_sph(2)= lat_v(ip)
        P_sph(3)= 1.d0
        call basis_between_sph_car (P_sph, basis, invbasis, ndim)
        P1_mat(1:3,1:3)=invbasis(1:3,1:3)
        !
        P_sph(1)= lonp(ip)
        P_sph(2)= latp(ip)
        P_sph(3)= 1.d0
        call basis_between_sph_car (P_sph, basis, invbasis, ndim)
        P2_mat(1:3,1:3)=basis(1:3,1:3)
        call AB_mult (P1_mat, Euler, basis,  3, 3, 3)
        call AB_mult (basis,  P2_mat,basis3, 3, 3, 3)
        !
        xp= latp(ip)
        call TC3_up (u0, xe, thetab, thetae, xp, up)
        Y(1)= up     ! U
        Y(2)= 0.d0   ! V
        Y(3)= 0.d0   ! W
        call AX_mult (basis3, Y, X, 3, 3, 3, 3)
        u(ip)=X(1)
        v(ip)=X(2)
     enddo
     !! ck uv
     !call indexx (NPTS,lat_v,iwksp)   ! sort in ascending order
     !open(9, file='uv_TC3.dat', status='unknown')
     !! output "lamb, theta (increasing order), u, v
     !do ip=1, NPTS
     !   j=iwksp(ip)
     !   write(9, 103) lon_v(j), lat_v(j), u(j), v(j)
     !enddo
     !close(9)
     deallocate(latp, lonp, iwksp, wksp)
     !
     !
     !


  elseif (iswcase .EQ. 10) then        !  TC10
     !
     ! Polar vortex erosion
     !
     Nsample = 3001
     Nsample3= 3*Nsample-2
     LDA=Nsample3; LDB=nih
     allocate (latp_sample(Nsample3),hp(Nsample3))   ! 
     allocate (iwksp(nih), wksp(nih))  ! for h field
     allocate (y2(LDA), work(LDA), work2(LDB))
     !
     ! prepare table for h field
     ! 3 segments i
     !

     hp(1)=6000.d0
     thetab=-pi/2.d0; thetae=0.d0   ; h0=hp(1)
     call TC10_h(h0, thetab, thetae, Nsample, latp_sample(1), hp(1))
     thetab=0.d0;     thetae=pi/3.d0; h0=hp(Nsample)
     call TC10_h(h0, thetab, thetae, Nsample, latp_sample(Nsample), hp(Nsample))
     thetab=pi/3.d0;  thetae=pi/2.d0; h0=hp(2*Nsample-1)
     call TC10_h(h0, thetab, thetae, Nsample, latp_sample(2*Nsample-1), hp(2*Nsample-1))
!!     test
!!     do j=1, 3*Nsample-2
!!        write(7,103) latp_sample(j), hp(j)
!!     enddo


     call indexx (nih,lat_h,iwksp)   ! sort in ascending order
     do j=1,nih
        wksp(j)=lat_h(iwksp(j))
     enddo
     !
     ! cube spline interpolation
     !
     yp1=0.d0; ypn=0.d0
     call sple_ab(latp_sample, hp, LDA, yp1, ypn, y2, work)
     call splt_ab(LDA, latp_sample, hp, y2, LDB, wksp, work2)
     

!!     open(7, file='hp.dat', status='unknown')
!!     open(8, file='fit_hp.dat', status='unknown')
!!     do j=1,LDA
!!        write(7, 103)  latp_sample(j), hp(j), y2(j)
!!     enddo
!!     do j=1,LDB
!!        write(8, 103)  wksp(j), work2(j)
!!     enddo
!!     close(7)
!!     close(8)
!!     !
!!     stop 'ck latp lonp'

     
     do j=1, nih
        ip=iwksp(j)
        FV(ip,1)=work2(j)
     enddo

!     do ip=1, nih
!        write(7,141) 'ip, lat, lon, h', ip,lat_h(ip),lon_h(ip),FV(ip,1)
!     enddo
!     STOP 'ck h field'
     deallocate(latp_sample, hp,  wksp, iwksp, y2, work, work2)
     !
     !
     ! end  ::  iswcase == 10
     !
     !

  elseif (iswcase .EQ. 7) THEN
     !
     ! netcdf ERA40_1978Dec_500mb_huv
     !
     it=21                                                !  21-Dec-1978
     filename=trim(file_nc_init)
     write(6,*) 'filename=', filename
     call get_file_dimension(filename, LDA, LDB, LDT)     ! nx,ny,nt
     allocate (lonp(LDA), latp(LDB))
     allocate (FV_aux_t(LDA,LDB,LDT,4))
     write(6, 121) 'nx=nlon, ny=nlat, nt=ntime', LDA, LDB, LDT
     !
     keyword="longitude"
     call get_1d_netcdf(filename, keyword, lonp(1), LDA)
     write(6,*) 'f1'

     keyword="latitude"
     call get_1d_netcdf(filename, keyword, latp(1), LDB)

     keyword="z"
     call get_3d_netcdf(filename, keyword, FV_aux_t(1,1,1,1), LDA, LDB, LDT)
     FV_aux_t(1:LDA,1:LDB,it,1) = FV_aux_t(1:LDA,1:LDB,it,1)/g  ! h:  geopotential / g
     write(6,*) 'f3'

     keyword="u"
     call get_3d_netcdf(filename, keyword, FV_aux_t(1,1,1,2), LDA, LDB, LDT)

     keyword="v"
     call get_3d_netcdf(filename, keyword, FV_aux_t(1,1,1,3), LDA, LDB, LDT)

     !

     !
     !-- interpolate h, u, v from (lonp, latp) to (lon_h,lat_h) or (lon_v,lat_v)
     !
     allocate ( grid_gauss(LDA*LDB, 2 ) )    ! borrow spectra var
     allocate ( wksp1     (LDA*LDB    ) )
     allocate ( wksp2     (LDA*LDB    ) )
     allocate ( wksp3     (LDA*LDB    ) )
     allocate ( grid_general (NPTS, 2 ) )

     ipt=0
     do j=1, LDB        ! lat
        do i=1, LDA     ! lon
           ipt=ipt+1
           grid_gauss(ipt,1) = latp(j) * d2r      ! true lat
           grid_gauss(ipt,2) = lonp(i) * d2r      ! true lon
           wksp1(ipt)=FV_aux_t(i,j,it,1)
           wksp2(ipt)=FV_aux_t(i,j,it,2)      ! U
           wksp3(ipt)=FV_aux_t(i,j,it,3)      ! V
        enddo
     enddo
     !

     !- interpolate h
     gh_icos_yx(1:nih,1)=lat_h(1:nih); gh_icos_yx(1:nih,2)=lon_h(1:nih)
     ipt= LDA * LDB
     call slint_init (grid_gauss,ipt,gh_icos_yx,nih)
     call bilinear_interp (wksp1, FV(1,1))    ! h
     !!
     do ip=1,nih, nih/20
        write(6,101) 'FV(ip,1)',  FV(ip,1)
     enddo
     !
     !
     !- interpolate U,V
     do ipt=1,NPTS
        grid_general(ipt,1)= lat_v(ipt)
        grid_general(ipt,2)= lon_v(ipt)
     enddo
     ipt= LDA * LDB
     call slint_init (grid_gauss,ipt,grid_general,NPTS)
     call bilinear_interp (wksp2, u(1))       ! U  u(1:NPTS)
     call bilinear_interp (wksp3, v(1))       ! V  v(1:NPTS)
     !
     deallocate ( grid_gauss, wksp1, wksp2, wksp3, grid_general )
     !!
     !!do ip=1,NPTS,20
     !!   write(6,101) 'u(ip), v(ip)', u(ip), v(ip)
     !!enddo
     !
     !-- direct output nc4
     !
     filename='h0_field_TC7.nc'
     call lon_lat_press_wt_nc4 (2,LDA,LDB,lonp,latp,FV_aux_t(1,1,it,1),filename)
  ENDIF



  !
  !
  !---  transform velocity (u, v) to cartesian 
  !---    (Four vector:  h/Vx/Vy/Vz)
  !
  !
  IF (iswcase .NE. 4) THEN
!
     do ipt=1, NPTS
        lat= lat_v(ipt)
        lon= lon_v(ipt)
        P_sph(1)=lon; P_sph(2)=lat; P_sph(3)=1.d0
        call basis_between_sph_car (P_sph, basis, invbasis, ndim)
        X(1)=u(ipt); X(2)=v(ipt); X(3)=0.d0
        m=ndim; n=ndim; mp=m; np=n
        call AX_mult(basis, X, Y, m, n, mp, np)       ! Y(1:3) cartesian velocity
        if (stagger .EQ. 'A' .OR. stagger(1:2).EQ.'Du' .OR. stagger(1:2).EQ.'ZM') then
           FV(ipt,2)=Y(1)
           FV(ipt,3)=Y(2)
           FV(ipt,4)=Y(3)
           vec3(1:3)=0.d0    !  aux
        elseif (stagger .EQ. 'C') then
           !
           ! Proj velocity to N T direction for C-grid
           !      Un=(V,N)N,  Ut=(V,T)T
           !  
           ip=C4E(1,ipt); is=C4E(5,ipt)   ! Left:  1:4=i,j,k,Lp, 5=is
           vec1(1:3)= Nvec(1:3,is,ip)     ! for iE;  NoutdUn(is,ip)= + 1 
           call XY_dot (Y,vec1,n,ndim,s)
           vec1(1:3)=  Tvec(1:3,is,ip)   ! 
           call XY_dot (Y,vec1,n,ndim,s2)
           FV(ipt,2)=s             !  Un
           FV(ipt,3)=s2            !  Ut
!           !
!           !-- ck
!           write(6,121) 'iE,is,ip       ', ipt,is,ip
!           write(6,101) 'Nvec(1:3,is,ip)', Nvec(1:3,is,ip)
!           write(6,101) 'Tvec(1:3,is,ip)', Tvec(1:3,is,ip)
!           write(6,101) 'V_xyz(1:3)     ', Y(1:3)
!           write(6,101) 'VdN            ', s
!           write(6,101) 'Vdt            ', s2
           !
           ! backward test  Un* \vec U_n + Ut * \vec U_T = V = omega \cross R
           do j=1, 3
              vec3(j)= s * Nvec(j,is,ip) + s2 * Tvec(j,is,ip)
           enddo
        endif
        !
        ! check with  V = omega \cross R
        !
        vec1(1)=-sin(alpha); vec1(2)=0.d0; vec1(3)=cos(alpha)  ! \vec Omegasw
        call vect_from_sph_2_car (P_sph, P_car, ndim)           ! \vec R
        call XY_cross (vec1, P_car, vec2, ndim)                 !  omega \cross R
        vec2(1:3)=u0*vec2(1:3)
        call X_2norm(vec2, s2, ndim)
        call X_2norm(vec3, s3, ndim)
        call X_2norm(Y, s, ndim)
!        !
!        !-- ck
!        !
!        write(6,121) 'ipt=      ', ipt
!        write(6,101) 'lon, lat  ', lon, lat
!        write(6,101) '\vec Omega', vec1(1:3)
!        write(6,101) '\vec R    ', P_car(1:3)
!        write(6,101) '\vec Omega \cross R', vec2(1:3)
!        write(6,101) 'vec1(1:3)-Y(1:3), norm diff', vec2(1:3)-Y(1:3), s-s2
!        write(6,101) 'vec3(1:3)-Y(1:3), norm diff', vec3(1:3)-Y(1:3), s-s3
!        write(6,101) 'Un, Ut    ', FV(ipt,2)  , FV(ipt,3)
!        write(6,'(//)')
        !
     enddo      ! ipt

!
!
!     !
!     ! check C-grid velocity interpolation scheme
!     !
!     if (stagger .EQ. 'C') then
!        do ipt=1, NPTS
!           iE = ipt
!           ip=C4E(1,iE); jp=C4E(2,iE)
!           is=C4E(5,iE); js=proxs(is,ip)
!           nb= nprox(ip)
!           ism=mod(is+nb-2, nb)+1
!           !
!           nvert = nprox(ip) + nprox(jp) - 2        ! double Benzene
!           s=0.d0
!           do is2=1,nvert                           ! pentagon + hexagon =9; otherwise 10
!              iE2= Eindex(is2,iE)                   ! 
!              ipc=C4E(1,iE2); isc=C4E(5,iE2)        ! ip,is on circle
!              s= s + wtUt(is2,iE)*slen(isc,ipc)*FV(iE2,2)
!           enddo
!           Ut= s/dcenter(iE)
!           write(6,101) 'Ut  vs  Ut_inter', FV(iE,3), Ut
!           write(6,101) 'diff            ', FV(iE,3)-Ut
!        enddo
!     endif
!     !
!     !STOP 'check MPAS Nvec '
!
  ENDIF
  write(6,101) 'u0=', u0




  !
  !-- Single out:  case4 Forced non-linear system with a translating low
  !
  if (iswcase .EQ. 4) then
     ims=1; ime=NPTS
     ips=1; ipe=NPTS      ! these are local variables
     itstep=1
     istage=1
     allocate (FV_aux(NPTS,ndimFV))
     call test_phi0
     !
     ! finish set up  FV/AFV_aux(ip, 1:4) in A, FV/AFV_aux(iE,1:2) in C
     !
     call add_forcing_TC4 (ims,ime,ips,ipe,FV,FV_aux,ndimFV,itstep,istage)
     deallocate (FV_aux)
  endif
  call flush(6)




  !
  !-- assign analytical solutions  (FV_ana for TC4 is done in swm.F90 with time flow)
  !                               
  !--------------------------------------------------------
  !  using NICAM test case
  !--------------------------------------------------------
  !
  !  h= beta=   cos(m lambda) * cos^4(n theta) 
  !     alpha=  sin(  lambda)
  !  V=         alpha \grad beta
  !  gradh=           \grad beta
  !
!
!
! A-grid   ASV(nih, 1:12) == 
!          grad_h,       (div_V, vort, div_hV), (f + xi ) k \cross V, \dot(\vec{V})
! C-grid   ASV(niE, 1:12) ==  
!          (grad_h,0,0), (div_V, vort, div_hV), (eta*Ut,Un,Ut), (\part Un / \part t,0,0)
!                                                               ------------------
!                                                                grad(K) test first
!
!

  if (iswcase .eq. 0) then
     if (mod(nx,2).NE.1) then
        STOP 'wrong input, nx is even, \grad beta(h) diverges as 1/delta'
     endif
     FV(:,:)=-1000.d0   ! initial reset
     !
     write(6,'(//,2x,a)') 'Init for NICAM analytical solution'
     u0x = u0/5.d0        !  note: u0=38.6m/s  end up giving V_max = 200 m/s, too large
     h0x = (u0x*u0x)/(2.d0*g)

     !
     !  First, find all values at center for A-grid and C-grid
     !
     do ipt=1, nih
        lat= lat_h(ipt)
        lon= lon_h(ipt)
        beta= cos(mx*lon) * cos(nx*lat)**4.d0
        u_lon= -u0x * mx*cos(nx*lat)**4.d0/cos(lat)*sin(mx*lon)*sin(lon)
        u_lat= -u0x * 4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(mx*lon)*sin(lon)
        !
        P_sph(1)=lon; P_sph(2)=lat; P_sph(3)=1.d0
        call basis_between_sph_car (P_sph, basis, invbasis, ndim)
        X(1)=u_lon; X(2)=u_lat; X(3)=0.d0
        m=ndim; n=ndim; mp=m; np=n
        call AX_mult(basis, X, Y, m, n, mp, np)
        !
        FV(ipt,1)= h0x * beta      ! h
        if (stagger .EQ. 'A') then
           FV(ipt,2)= Y(1)            ! Vx, y, z
           FV(ipt,3)= Y(2)
           FV(ipt,4)= Y(3)
        endif
        V_car(2) = Y(1)
        V_car(3) = Y(2)
        V_car(4) = Y(3)

!        write(6,121) 'ipt=', ipt
!        write(6, 101) 'lon, lat=', lon, lat
!        write(6, 101) 'u,v, FV(1:4)'
!        write(6, 103) X(1:2)
!        write(6, 105) FV(ipt,1:4)
        !
        ! re-use  u_lon, u_lat  as  grad h \dot (e1, e2)
        !
        if ( abs(pi/2.d0 - abs(lat)) .lt. polar_tol ) then
           gradh_lon=  0.d0
        else
           gradh_lon= - mx*cos(nx*lat)**4.d0/cos(lat)*sin(mx*lon)
        endif
        gradh_lat= -4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(mx*lon)
        X(1)=gradh_lon; X(2)=gradh_lat; X(3)=0.d0
        call AX_mult(basis, X, Y, m, n, mp, np)
        !
        ! analytic grad(h)   --> ASV_ana(1:3)
        ! analytic div, vort --> ASV_ana(4:5)
        !
        ASV_ana(ipt,1)=Y(1)* h0x/ae  ! \grad h _x, y, z
        ASV_ana(ipt,2)=Y(2)* h0x/ae
        ASV_ana(ipt,3)=Y(3)* h0x/ae
        !
        ! div(V), vort,  div(h V)  analytic
        !
        if ( abs(pi/2.d0 - abs(lat)) .lt. polar_tol ) then
           ASV_ana(ipt,4)= 0.d0   ! div
           ASV_ana(ipt,5)= 0.d0   ! vort
           ASV_ana(ipt,6)= 0.d0   ! div ( beta V )
        else
           ASV_ana(ipt,4)= &
                -mx*cos(nx*lat)**4.d0/cos(lat)**2.d0*(cos(lon)*sin(mx*lon)+mx*sin(lon)*cos(mx*lon)) - &
                4.d0*nx*sin(lon)*cos(mx*lon)* &
                ( -sin(lat)*cos(nx*lat)**3.d0*sin(nx*lat)/cos(lat) - &
                   3.d0*nx*cos(nx*lat)**2.d0*sin(nx*lat)**2.d0 + nx*cos(nx*lat)**4.d0 )
           ASV_ana(ipt,5)= &
                -4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(lon)*cos(mx*lon)/cos(lat)
           !
           ! \grad \dot ( h V )  =  \div  ( h0*beta  * v0* alpha* grad(beta) )
           !                     =  h0*u0 *  alpha * [ \grad (beta)^2 + beta  \div (alpha grad(beta)) ]
           ! part-1: alpha * |\grad beta|^2
           s = (gradh_lon*gradh_lon + gradh_lat*gradh_lat) * sin(lon)
           ! part-2:   beta * div (a * grad(b))
           s2=  beta * ASV_ana(ipt,4)       !  b * div (a grad(b))
           ASV_ana(ipt,6) = s + s2
           !
           ASV_ana(ipt,4)= ASV_ana(ipt,4) * u0x / ae
           ASV_ana(ipt,5)= ASV_ana(ipt,5) * u0x / ae           
           ASV_ana(ipt,6)= ASV_ana(ipt,6) * h0x * u0x / ae
        endif
        !
        ! find ( f + xi ) r_unit  \cross  V ;  exact for A-grid, wrong for C-grid due to fcori
        !
        vort = ASV_ana(ipt,5)
        do i=1,3
           j=mod(i,3)+1; k=mod(j,3)+1
!!           vec1(i)= (vort+fcori(ipt))*(Rcvec(j,ipt)*FV(ipt,k+1)-Rcvec(k,ipt)*FV(ipt,j+1))
           vec1(i)= (vort+fcori(ipt))*(Rcvec(j,ipt)*V_car(k+1)-Rcvec(k,ipt)*V_car(j+1))
        enddo
        ASV_ana(ipt,7) = vec1(1)
        ASV_ana(ipt,8) = vec1(2)
        ASV_ana(ipt,9) = vec1(3)
        !
        !
        ! re-use u_lon, u_lat as grad ( 1/2 V*V ) == grad(K)
        !
        if ( abs(pi/2.d0 - abs(lat)) .lt. polar_tol ) then
           A= 0.d0
           B= 0.d0
           u_lon= 0.d0
           u_lat= 0.d0
        else
           A= mx*cos(nx*lat)**4.d0/cos(lat)*sin(lon)*sin(mx*lon)
           B= 4.d0*nx*sin(nx*lat)*cos(nx*lat)**3.d0*sin(lon)*cos(mx*lon)
           u_lon= mx*A*&
                cos(nx*lat)**4.d0/cos(lat)**2.d0* &
                (cos(lon)*sin(mx*lon)+mx*sin(lon)*cos(mx*lon)) + &
                4.d0*nx*B* &
                cos(nx*lat)**3.d0*sin(nx*lat)/cos(lat)* &
                (cos(lon)*cos(mx*lon)-mx*sin(lon)*sin(mx*lon))
           u_lat= mx*A* &
                (-4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(lat) + cos(nx*lat)**4.d0*sin(lat))/ &
                cos(lat)**2.d0 * sin(lon)*sin(mx*lon) + &
                4.d0*nx*B* &
                (-3.d0*nx*cos(nx*lat)**2.d0*sin(nx*lat)**2.d0 + nx*cos(nx*lat)**4.d0)* &
                sin(lon)*cos(mx*lon)
        endif
        X(1)=u_lon; X(2)=u_lat; X(3)=0.d0
        call AX_mult(basis, X, Y, m, n, mp, np)
        !
        ! analytic grad(h)   --> ASV_ana(1:3)
        ! analytic div, vort --> ASV_ana(4:5)
        !
        ASV_ana(ipt,10)=Y(1)* u0x**2.d0/ae  ! \grad(K)  _x, y, z
        ASV_ana(ipt,11)=Y(2)* u0x**2.d0/ae
        ASV_ana(ipt,12)=Y(3)* u0x**2.d0/ae
        !
     enddo


     if (stagger .EQ. 'C') then 
        !
        ! Need to correct the values for C-grid, that are always define at
        ! center of cell:  \grad(h), (f + xi) k \cross V
        ! 
        niE_loc = NPTS
        FV(1:niE_loc,2)=0.d0        !  clear up previous setup
        write(6, 121)    'niE_loc', niE_loc
        do iE=1, niE_loc
           lat= lat_v(iE)
           lon= lon_v(iE)
           beta= cos(mx*lon) * cos(nx*lat)**4.d0
           u_lon= -u0x * mx*cos(nx*lat)**4.d0/cos(lat)*sin(mx*lon)*sin(lon)
           u_lat= -u0x * 4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(mx*lon)*sin(lon)
           !
           P_sph(1)=lon; P_sph(2)=lat; P_sph(3)=1.d0
           call basis_between_sph_car (P_sph, basis, invbasis, ndim)
           X(1)=u_lon; X(2)=u_lat; X(3)=0.d0
           m=ndim; n=ndim; mp=m; np=n
           call AX_mult(basis, X, Y, m, n, mp, np)
           !
           !!  FV(iE,1) ==  h at center:  calculated previously
           !   FV(1:nih,1) -- in fact
           !
           ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)
           vec1(1:3) = Nvec(1:3,is,ip)
           call XY_dot (Y, vec1, 3, 3, s)
           FV(iE,2)= s * NoutdUn(is,ip)      ! U_n^e  \dot V
           !
           vec1(1:3) = Tvec(1:3,is,ip)
           call XY_dot (Y, vec1, 3, 3, s)
           FV(iE,3)= s * NoutdUn(is,ip)      ! U_t^e  \dot V
           ASV_ana(iE,8) = FV(iE,2)          ! U_n
           ASV_ana(iE,9) = FV(iE,3)          ! U_t
           !!write(6, 101) 'C-grid U lon, lat, n, t ', u_lon, u_lat, FV(iE,2), FV(iE,3)
           !
           ! gradient(h) \dot U_n^e  :  on edge
           ! 
           if ( abs(pi/2.d0 - abs(lat)) .lt. polar_tol ) then
              gradh_lon=  0.d0
           else
              gradh_lon=  - mx*cos(nx*lat)**4.d0/cos(lat)*sin(mx*lon)
           endif
           gradh_lat= -4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(mx*lon)
           X(1)=gradh_lon; X(2)=gradh_lat; X(3)=0.d0
           call AX_mult(basis, X, Y, m, n, mp, np)
           !
           ! Y \dot Nvec(is,ip) * NoutdUn
           !
           vec1(1:3) = Nvec(1:3,is,ip)
           call XY_dot (Y, vec1, 3, 3, s)      ! N dot e1
           ASV_ana(iE,1)= ( s * h0x / ae ) * NoutdUn(is,ip)   ! U_n^e * | N_out > < N_out |  \dot \grad h
           ASV_ana(iE,2)= 0.d0
           ASV_ana(iE,3)= 0.d0

           ! 
           ! previous calculated only for C-grid at cell center 
           ! 
           ! div(V)   at center :  previously calculated    ASV_ana(iE,4)
           ! vort     at mid-edge:                          ASV_ana(iE,5)
           ! div(h*V) at center :  previously calculated    ASV_ana(iE,6)
           ! 
           if ( abs(pi/2.d0 - abs(lat)) .lt. polar_tol ) then
              vort = 0.d0   ! vort
           else
              vort  = &
                   -4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(lon)*cos(mx*lon)/cos(lat)
              vort  = vort * u0x / ae
           endif
           ASV_ana(iE,5) = vort
           !
           ! find ( f + xi ) r_unit  \cross  V  for  C-grid due to fcori
           ! 
           !
           vec1(1)= -(vort + 0.5d0*(fcori_V(1,iE)+fcori_V(2,iE))) * u_lat
           vec1(2)=  (vort + 0.5d0*(fcori_V(1,iE)+fcori_V(2,iE))) * u_lon
           vec1(3)=  0.d0
           !
           ! from (lambda, theta, r)  to (x,y,z)
           ! 
           !   X =  P * Y
           !
           call AX_mult(basis, vec1, Y3, m, n, mp, np)
           !
           !    proj to U_n^e direction,    \vec dot N  * NoutdUn
           !
           vec1(1:3)= Nvec(1:3,is,ip)       ! for iE;  NoutdUn(is,ip)= + 1
           call XY_dot (Y3, vec1, 3, 3, s)  !
           ASV_ana(iE,7) = s*NoutdUn(is,ip) ! proj to U_n
           ASV_ana(iE,8) = FV(iE,2)         ! U_n   already appeared in code above
           ASV_ana(iE,9) = FV(iE,3)         ! U_t

           
           !
           !
           ! re-use u_lon, u_lat as grad ( 1/2 V*V ) == grad(K)
           !
           if ( abs(pi/2.d0 - abs(lat)) .lt. polar_tol ) then
              A= 0.d0
              B= 0.d0
              u_lon= 0.d0
              u_lat= 0.d0
           else
              A= mx*cos(nx*lat)**4.d0/cos(lat)*sin(lon)*sin(mx*lon)
              B= 4.d0*nx*sin(nx*lat)*cos(nx*lat)**3.d0*sin(lon)*cos(mx*lon)
              u_lon= mx*A*&
                   cos(nx*lat)**4.d0/cos(lat)**2.d0* &
                   (cos(lon)*sin(mx*lon)+mx*sin(lon)*cos(mx*lon)) + &
                   4.d0*nx*B* &
                   cos(nx*lat)**3.d0*sin(nx*lat)/cos(lat)* &
                   (cos(lon)*cos(mx*lon)-mx*sin(lon)*sin(mx*lon))
              u_lat= mx*A* &
                   (-4.d0*nx*cos(nx*lat)**3.d0*sin(nx*lat)*cos(lat) + cos(nx*lat)**4.d0*sin(lat))/ &
                   cos(lat)**2.d0 * sin(lon)*sin(mx*lon) + &
                   4.d0*nx*B* &
                   (-3.d0*nx*cos(nx*lat)**2.d0*sin(nx*lat)**2.d0 + nx*cos(nx*lat)**4.d0)* &
                   sin(lon)*cos(mx*lon)
           endif
           X(1)=u_lon; X(2)=u_lat; X(3)=0.d0
           call AX_mult(basis, X, Y, m, n, mp, np)
           !
           !    proj to U_n^e direction,    \vec dot N  * NoutdUn
           !
           vec1(1:3)= Nvec(1:3,is,ip)      ! for iE;  NoutdUn(is,ip)= + 1
           call XY_dot (Y, vec1, 3, 3, s)  !
           ASV_ana(iE,10) = s*NoutdUn(is,ip) ! proj to U_n
           !
        enddo     !  loop iE
     endif
!        !
!        write(6, 101) '(hx,y,z), (div,vort, div(hV)), (\yeta k \cross V)'
!        write(6, 103) ASV_ana(iE,1:9)
!
!


  elseif (iswcase.EQ.2 .OR. iswcase.EQ.3 &
       .OR. iswcase.EQ.8 .OR. iswcase.EQ.9) then
     ! steady state let's save initial state to compute Error norms
     ! 1  2  3  4       (10,11,12)
     ! h  Vx Vy Vz        grad(K)
     !
     ! h 

     if (stagger.EQ.'A' .OR. stagger(1:2).EQ.'Du') then
        do ip=1, nih
           ASV_ana(ip, 1)=FV(ip, 1)
           !! write(6, 101) 'FV(ip, 1)', FV(ip, 1)
        enddo
        ! velocity cartesian
        do ipt=1, NPTS
           do j=2,4
              ASV_ana(ipt, j)=FV(ipt, j)
           enddo
           !
           !  grad(K) = - \omega_sw
           !  remove  radial component
           !
           X(1)=0.d0
           X(2)=-sin(alpha)
           X(3)=cos(alpha)
           !
           !   remove radial component, \vec R component  :  \vec F = \vec F  - (F dot R) \vec R
           s=0.d0              ! s= (gradh, Rcvec) not Rvec on edge
           do j=1, 3
              s= s+ X(j)*Rcvec(j,ipt)
           enddo
           do j=1, 3
              X(j)= X(j)-s*Rcvec(j,ipt)    ! grad = grad - (grad,R) R
           enddo
           X(1:3)= -X(1:3)*omegasw
           ASV_ana(ipt,10)=X(1)
           ASV_ana(ipt,11)=X(2)
           ASV_ana(ipt,12)=X(3)
           !
        enddo
     elseif (stagger(1:2).EQ.'ZM') then
        do ip=1, nih
           ASV_ana(ip, 1)=FV(ip, 1)
           !! write(6, 101) 'FV(ip, 1)', FV(ip, 1)
        enddo
        ! velocity cartesian
        do ipt=1, NPTS
           do j=2,4
              ASV_ana(ipt, j)=FV(ipt, j)
           enddo
           !
           !  grad(K) = - \omega_sw
           !  remove  radial component
           !
           X(1)=0.d0
           X(2)=-sin(alpha)
           X(3)=cos(alpha)
           !
           !   remove radial component, \vec R component  :  \vec F = \vec F  - (F dot R) \vec R
           s=0.d0              ! s= (gradh, Rcvec) not Rvec on edge
           do j=1, 3
              s= s+ X(j)*RVvec(j,ipt)
           enddo
           do j=1, 3
              X(j)= X(j)-s*RVvec(j,ipt)    ! grad = grad - (grad,R) R
           enddo
           X(1:3)= -X(1:3)*omegasw
           ASV_ana(ipt,10)=X(1)
           ASV_ana(ipt,11)=X(2)
           ASV_ana(ipt,12)=X(3)
           !
        enddo
     elseif (stagger.EQ.'C') then
        do ip=1, nih
           ASV_ana(ip, 1)=FV(ip, 1)
        enddo
        ! velocity Un 
        do ipt=1, NPTS
           do j=2,3
              ASV_ana(ipt, j)=FV(ipt, j)  !  U_n and U_t  vel
           enddo
        enddo
        do iE=1, NPTS
           !
           !  grad(K) = - \omega_sw
           !  remove  radial component
           !
           ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)
           !
           X(1)=0.d0
           X(2)=-sin(alpha)
           X(3)=cos(alpha)
           !
           !   remove radial component, \vec R component  :  \vec F = \vec F  - (F dot R) \vec R
           s=0.d0              ! s= (gradh, Rcvec) not Rvec on edge
           do j=1, 3
              s= s+ X(j)*Rvec(j,is,ip)
           enddo
           do j=1, 3
              X(j)= X(j)-s*Rvec(j,is,ip)    ! grad = grad - (grad,R) R
           enddo
           vec1(1:3) = Nvec(1:3,is,ip)
           call XY_dot (X, vec1, 3, 3, s)
           ASV_ana(ip,10)= s * NoutdUn(is,ip)      ! U_n^e  \dot V
        enddo
        !
     endif

     !
     ! print out
     !
     write(6, 121) 'ip=', (ip, ip=1, nih, nih/6)
     do ip=1, nih, nih/6
        write(6, 112) 'ASV_ana(ip, 1:4), h/Vx/Vy/Vz, hb_a(ip)', &
             ASV_ana(ip, 1:4), hb_a(ip)
     enddo
     !!
     !! STOP 'test final ASV_ana'
     !!

  elseif (iswcase.EQ.5 .AND. iswcase.LT.7) then    ! test     
!!  elseif (iswcase.EQ. 7) then    ! test     
     !
     !-- read nlon,nlat,glon,glat,h,u,v from stswm solutions
     !
     if (stagger .EQ. 'A') then
        allocate (icos_grid_aux(2,nih))
        INQUIRE (file='stswm_icos_huv.dat',  EXIST=fexist)
        if (fexist) then
           open(unit=11,file='stswm_icos_huv.dat',status='unknown',  &
                form="unformatted", action='read')
           read(11) icos_grid_aux(1:2,1:nih), ASV_ana(1:nih,1), &
                u(1:nih), v(1:nih)
           close(11)
           write(6,*) 'Finish up reading stswm_icos_huv.dat'
           do ipt=1, nih     ! not yet NPTS 
              lon= icos_grid_aux(1,ipt)
              lat= icos_grid_aux(2,ipt)
              P_sph(1)=lon; P_sph(2)=lat; P_sph(3)=1.d0
              call basis_between_sph_car (P_sph, basis, invbasis, ndim)
              X(1)=u(ipt); X(2)=v(ipt); X(3)=0.d0
              m=ndim; n=ndim; mp=m; np=n
              call AX_mult(basis, X, Y, m, n, mp, np)       ! Y(1:3) cartesian velocity
              ASV_ana(ipt,2)=Y(1) 
              ASV_ana(ipt,3)=Y(2) 
              ASV_ana(ipt,4)=Y(3) 
           enddo
        else
           write(6,*) 'Missing stswm_icos_huv.dat, ASV_ana(:,:) == 0.d0, fake data '
           ASV_ana(:,:)=1.d5
        endif
        deallocate(icos_grid_aux)
     else
!!        STOP 'C-grid not implemented'
        write(6,101) 'sw_test_init  TC5 C-grid not implemented'
        ASV_ana(:,:)=1.d5
        ! YGYU a bug at here
     endif
     !
     !
     if (stagger .EQ. 'A') then
        write(6, 101) 'A-grid: from STSWM'
        do ipt=1,NPTS,10**4
           write(6,101) 'ASV_ana(ipt,1:4)', ASV_ana(ipt,1:4)
        enddo
     endif
!!     elseif (stagger .EQ. 'C') then
!!        write(6, 101) 'C-grid: ?? bilinear interpolation of Vx, Vy, Vz from spectra model '
!!        do ipt=1,NPTS
!!           !
!!           ! Proj velocity to N T direction for C-grid
!!           !      Un=(V,N)N,  Ut=(V,T)T
!!           !  
!!           Y(1)= ASV_ana(ipt,2)  
!!           Y(2)= ASV_ana(ipt,3)
!!           Y(3)= ASV_ana(ipt,4)
!!           ip=C4E(1,ipt); is=C4E(5,ipt)   ! Left:  1:4=i,j,k,Lp, 5=is
!!           vec1(1:3)= Nvec(1:3,is,ip)     ! for iE;  NoutdUn(is,ip)= + 1 
!!           call XY_dot (Y,vec1,n,ndim,s)
!!           vec1(1:3)=  Tvec(1:3,is,ip)   ! 
!!           call XY_dot (Y,vec1,n,ndim,s2)
!!           ASV_ana(ipt,2)=s             !  Un
!!           ASV_ana(ipt,3)=s2            !  Ut
!!           !
!!           !-- ck
!!           write(6,121) 'iE,is,ip       ', ipt,is,ip
!!           write(6,101) 'Nvec(1:3,is,ip)', Nvec(1:3,is,ip)
!!           write(6,101) 'Tvec(1:3,is,ip)', Tvec(1:3,is,ip)
!!           write(6,101) 'V_xyz(1:3)     ', Y(1:3)
!!           write(6,101) 'VdN            ', s
!!           write(6,101) 'Vdt            ', s2
!!        enddo
!!        !
!!     endif     ! stagger

  else
     do ipt=1, NPTS
        ASV_ana(ipt,1:6)= 1.d4      ! initial assignment
     enddo
  endif        ! iswcase

  !
  !
  !
  return
  include '../swm/myformat.inc'
end subroutine sw_test_init






subroutine TC3_up(u0, xe, thetab, thetae, thetap, up)
  use kinds, only: rt
  use module_constants, only: tol_exp
  implicit none
  real(rt), intent(in) :: u0, xe, thetab, thetae, thetap
  real(rt), intent(out):: up
  !
  real(rt) :: Tb, Te, T   !  T = theta
  !
  Tb=thetab;  Te=thetae;  T= thetap
  !  tol_exp  due to exponential function
  if (T.LE.(Tb+tol_exp) .OR. T.GE.(Te-tol_exp)) then
     up=0.d0
  else
     up= u0 * exp( (4.d0 - (Te-Tb)**2.d0 / ((T-Tb)*(Te-T))) / xe )
  endif
  return
end subroutine TC3_up


subroutine TC3_hp(h0, u0, xe, thetab, thetae, Npt_tot, thetap, hp)
  use kinds, only: rt
  use module_constants, only  :  g, ae, omega, pi, ndim, tol_exp  ! exp(-1/1.d-4)
  implicit none
  real(rt), intent(in) :: h0, u0, xe, thetab, thetae
  integer,  intent(in) :: Npt_tot
  real(rt), intent(out):: thetap(Npt_tot), hp(Npt_tot)   !  means h_prime
  !
  real(rt) :: Tb, Te, T_L, T      !  T: = theta
  real(rt) :: xl, xh, dx, ddx, xj
  real(rt) :: up, fcori
  real(rt) :: s_int, sum
  real(rt), allocatable :: xpt(:), ypt(:)
  integer  :: Nseg_L, Nseg, Nseg_tot
  integer  :: Npt_L, Npt, Ndiv, LDA
  integer  :: i, j, id
  !
  !  
  ! T_L      Tb                                      Te
  !
  ! ______   _________________________________________>
  !           |         |       |                    |
  !           1         2                           Npt
  !  Left     Xmin      xl      xh                   Xmax
  !             X X X X 
  !               Ndiv
  !   Nseg_L        Nseg
  !
  !         ypt:=  the integrand
  !
  !
  ! tot sequence
  ! 
  ! 1,   Nseg_L+1,                                   Nseg_tot+1
  !
  !
  T_L= -pi/2.d0
  Tb=thetab; Te=thetae
  Nseg_tot=Npt_tot -1
  Nseg_L = (Tb-T_L)/(Te-T_L) * Nseg_tot
  Nseg   = Nseg_tot - Nseg_L
  
  !!write(7, 101) '# thetap, up, integrand_p = fu + u^2*tan/a'
  !-- s.1  the left side trivial
  dx= (Tb-T_L)/dble(Nseg_L)
  do i=1, Nseg_L+1
     T= T_L + dble(i-1)*dx
     thetap(i)= T
     hp(i)= h0
     !!write(7, 103) T, 0.d0, 0.d0, 0.d0, 0.d0
     !!write(9, 103) T, hp(i)
  enddo
  
  Ndiv=251
  LDA= Ndiv
  allocate (xpt(Ndiv), ypt(Ndiv))
  !-- s.2  the right side is an integeral by simpson method
  dx= (Te-Tb)/dble(Nseg)
  ddx= dx/dble(Ndiv-1)
  s_int=h0
  do i=1, Nseg              ! count segment for integration
     id= (Nseg_L+1) + i     ! Right end  ID
     T=  Tb + dble(i-1)*dx  ! the left end
     xl= T 
     xh= T+dx
     thetap(id)= xh         ! Right side
     if (xl.LE.(Tb+tol_exp)) then
        ! integrand (up=0) being zero
        hp(id)= h0
     else
        ! find integrand first
        do j=1, Ndiv
           xj=xl+dble(j-1)*ddx
           xpt(j)=xj
           fcori=2.d0*omega*sin(xj)
           if (xj.GE.(Te-tol_exp)) then  ! to avoid  tg(x), x=pi/2 blow up
              ypt(j)=0.d0
           else
              call TC3_up(u0, xe, thetab, thetae, xj, up)
              ypt(j)=(fcori*ae + up*tan(xj))*up
           endif
           !!write(7, 103) xj, up, fcori*ae*up,  up*up*tan(xj),  ypt(j)
           !!write(7, 103) xj, ypt(j)
        enddo
        ! Simpson 3 pt integration
        call simpson_3pt_int (xpt, ypt, Ndiv, LDA, sum)
        !!write(9, 101) 's_int L, sum, s_int R', s_int, sum, s_int + sum
        !!write(9, 104) id, s_int, -sum/g, s_int - sum/g
        s_int= s_int - sum/g
        hp(id)= s_int
        !!write(9, 103) thetap(id), hp(id)
     endif
     !
  enddo
  deallocate(xpt, ypt)
  !
  return
  include '../swm/myformat.inc'
end subroutine TC3_hp



subroutine TC10_h(h0, thetab, thetae, Npt, theta, h)
  use kinds, only: rt
  use module_constants, only  :  g, ae, omega, pi, ndim, tol_exp  ! exp(-1/1.d-4)
  implicit none
  real(rt), intent(in) :: h0, thetab, thetae
  integer,  intent(in) :: Npt
  real(rt), intent(out):: theta(Npt), h(Npt)   !  means h_prime
  !
  real(rt) :: Tb, Te, T_L, T      !  T: = theta
  real(rt) :: xl, xh, dx, ddx, xj
  real(rt) :: u, fcori
  real(rt) :: s_int, sum
  real(rt), allocatable :: xpt(:), ypt(:)
  integer  :: Ndiv, LDA
  integer  :: i, j, id
  !
  !  
  ! Tb                                      Te
  !
  ! _________________________________________>
  !  |         |       |                    |
  !  1         2                           Npt
  !  Xmin      xl      xh                   Xmax
  !    X X X X 
  !      Ndiv
  !
  Ndiv=251
  LDA= Ndiv
  allocate (xpt(Ndiv), ypt(Ndiv))

  Tb= thetab; Te= thetae
  dx= (Te-Tb)/dble(Npt-1)
  ddx= dx/dble(Ndiv-1)
  do i=1, Npt               ! count segment for integration
     theta(i) = Tb + dble(i-1)*dx
  enddo

  s_int=h0
  h(1)=h0
  do i=1, Npt-1             ! count segment for integration
     xl= theta(i)
     xh= xl + dx
!!     if (xl.LE.(-pi/2.d0 + exp(-1.d0/tol_exp))) then
!!        ! integrand (up=0) being zero
!!        h(i+1)= h0
!!     else
        ! find integrand first
        do j=1, Ndiv
           xj=xl+dble(j-1)*ddx
           xpt(j)=xj
           fcori=2.d0*omega*sin(xj)
           if (xj .LE. (-pi/2.d0 + exp(-1.d0/tol_exp)) .OR.  &
               xj .GE. ( pi/2.d0 - exp(-1.d0/tol_exp))) then  ! to avoid  tg(x), x=pi/2 blow up
              ypt(j)=0.d0
           else
              call TC10_func_U(xj, u)
              ypt(j)=(fcori*ae + u*tan(xj))*u
           endif
           !!write(7, 103) xj, u, fcori*ae*u,  u*u*tan(xj),  ypt(j)
           !!write(7, 103) xj, ypt(j)
        enddo
        ! Simpson 3 pt integration
        call simpson_3pt_int (xpt, ypt, Ndiv, LDA, sum)
        !! write(9, 101) 's_int L, sum, s_int R', s_int, sum, s_int + sum
        !! write(9, 104) id, s_int, -sum/g, s_int - sum/g
        s_int= s_int - sum/g
        h(i+1)= s_int
        !! write(8, 103) theta(i+1), h(i+1)
        !
!!     endif
  enddo
  deallocate(xpt, ypt)
  !
  return
  include '../swm/myformat.inc'
end subroutine TC10_h



subroutine TC10_hs_t(t, nip, lon, lat, hb)
  use kinds, only: rt
  integer, intent(in) :: nip
  real(rt),intent(in) :: t, lon(nip), lat(nip)
  real(rt),intent(out):: hb(nip)       ! output TC10 time depdent terrain
  !
  integer  :: ip
  real(rt) :: t_d, FA, FB, lam, theta
  !
  t_d = t/86400.d0    ! t_d in day; t is in second
  do ip = 1, nip
     lam= lon(ip); theta= lat(ip)
     call TC10_func_A(t_d, FA)
     call TC10_func_B(theta, FB)
     hb(ip) = 720.d0 * FA * FB * sin(lam)
  enddo
  return
end subroutine TC10_hs_t
  



subroutine TC10_func_A(t_d, FA)
  use kinds, only: rt
  use module_constants, only : pi
  implicit none
  real(rt),intent(in) :: t_d
  real(rt),intent(out):: FA
  real(rt) :: t
  !
  ! t_d < 24  only 
  !
  t= t_d - int(t_d/20.d0)*20.d0           ! period of 20 days
  !
  if (t .GE. 0.d0 .AND. t .LT. 4.d0) then
     FA= 0.5d0 * (1.d0 - cos( pi * t / 4.d0 ))
  elseif (t .LE. 16.d0) then
     FA=1.d0
  else
     FA= 0.5d0 * (1.d0 + cos( pi * t / 4.d0 ))     
  endif
  return
end subroutine TC10_func_A
!  elseif (t .GE. 16.d0 .AND. t_d .LE. 20.d0) then     
!  elseif (t .GT. 20.d0 .AND. t .LE. 24.d0) then
!     FA= 0.5d0 * (1.d0 - cos( pi * (t-20.d0) / 4.d0 ))



subroutine TC10_func_B(theta, FB)
  use kinds, only: rt
  use module_constants, only : pi
  implicit none
  real(rt),intent(in) :: theta
  real(rt),intent(out):: FB
  real(rt) :: x
  !
  x= theta
  if (x .LE. 0.d0) then
     FB= 0.d0
  else
     FB= exp( 1.d0 - 1.d0/tan(x)**2.d0 ) / tan(x)**2.d0
  endif
  return
end subroutine TC10_func_B



subroutine TC10_func_U(theta, u)
  use kinds, only: rt
  use module_constants, only : pi
  implicit none
  real(rt),intent(in) :: theta
  real(rt),intent(out):: u   
  real(rt) :: x
  !
  x= theta
  if (x .LT. 0.d0) then
     u= -10.d0 - 20.d0/pi * x
  elseif (x .LT. pi/3.d0) then
     u= -10.d0 + 180.d0/pi * x
  elseif (x .GE. pi/3.d0) then
     u= - 300.d0/pi * ( x - pi/2.d0 )
  endif
  return
end subroutine TC10_func_U




subroutine add_forcing_TC4 (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep,istage)
  use kinds,            only  :  rt
  use readnamelist,     only  :  stagger
  use module_control,   only  :  nip, niE, NPTS, dt
  use module_variables, only  :  lon, lat, lon_v, lat_v, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  omegasw, stagecoef
  use module_variables, only  :  C4E, Nvec, Tvec
  use module_constants, only  :  ae, pi, ndim, omega  ! Earth daily rotation   2 pi / 1 day
  !
  implicit none
  integer, intent(in )::   ims, ime, ips, ipe, ndimFV, itstep, istage
  real(rt),intent(out)::   FV(ims:ime, ndimFV)       ! output TC4 FV_ana(1:4, ndimFV) 
  real(rt),intent(inout)::  AFV(ims:ime, ndimFV)     ! add forcing for Acceleration
  !
  real(rt), external :: F_int_sinx_t_sin_2x_to_14
  real(rt), external :: F_int_tanx_t_sin_2x_to_28 
  !
  ! var
  real(rt) :: F_u, F_v, F_phi, F_h
  real(rt) :: A, B, C1, D, E, F, G
  real(rt) :: u1, u2, v, Phi1, Phi2
  real(rt) :: u1_the
  real(rt) :: u2_lam, u2_the
  real(rt) :: u1_t_tan, u1_d_cos
  real(rt) :: v_lam,  v_the
  real(rt) :: Psi
  real(rt) :: Psi_c, Psi_cc
  real(rt) :: Psi_lam, Psi_the, Psi_the_lam, Psi_the_the, Psi_lam_lam
  real(rt) :: C, C_lam, C_the, C_lam_lam, C_the_lam
  !
  ! const
  real(rt) :: Psi0, h0, g0, Phi0, fcori0
  real(rt) :: lam0, the0
  real(rt) :: u0
  real(rt) :: sigma
  !
  ! local
  real(rt) :: t
  real(rt) :: lam, the, dlam   ! lambda, theta
  real(rt) :: fcori_p
  real(rt) :: P_sph(3), basis(3,3), invbasis(3,3), X(3), Y(3), YV(3)
  real(rt) :: xa, xb, vec1(3), s
  integer  :: ip, m, n, mp, np
  integer  :: ip2, is2
  
  !
  ! only for Shallow Water Test Case 4
  ! consider both A-grid and C-grid 
  !
  g0= 9.80616d0
  u0= omegasw*ae
  h0= 1.d5/g0
  fcori0= 2.d0*omega/sqrt(2.d0)
  Psi0=  -0.03d0*g0*h0/fcori0
  lam0= 0.d0
  the0= pi/4.d0
  sigma= 12.74244d0**2.d0
  t= dt*dble(itstep-1)+stagecoef(istage)
  m=ndim; n=ndim; mp=m; np=n


!!  write(6, 101) 'u0, fcori0, omega, omegasw', u0, fcori0, omega, omegasw



  if (stagger .EQ. 'A') then
     do ip=ips, ipe
        lam= lon_v(ip)     ! velocity grid point   note in A-grid  lat_v == lat  overlap
        the= lat_v(ip)
        !
        dlam= lam - lam0 - omegasw * t
        C= sin(the) * sin(the0) + cos(the) * cos(the0) * cos(dlam)
        Psi=    Psi0 * exp( -sigma * (1.d0 - C) / (1.d0 + C) )

        Psi_c=  2.d0 * sigma / ( 1.d0 + C )**2.d0 * Psi
        Psi_cc=(4.d0 * sigma * sigma / ( 1.d0 + C )**4.d0 - 4.d0 * sigma / ( 1.d0 + C )**3.d0) * Psi
        C_the= cos(the) * sin(the0) - sin(the) * cos(the0) * cos(dlam)
        C_lam=                      - cos(the) * cos(the0) * sin(dlam)
        C_lam_lam=                  - cos(the) * cos(the0) * cos(dlam)
        C_the_lam=                    sin(the) * cos(the0) * sin(dlam)
        Psi_the= Psi_c * C_the
        Psi_lam= Psi_c * C_lam
        Psi_the_the= Psi_cc * C_the * C_the - Psi_c * C
        Psi_the_lam= Psi_cc * C_the * C_lam + Psi_c * C_the_lam
        Psi_lam_lam= Psi_cc * C_lam * C_lam + Psi_c * C_lam_lam

        u1= u0 * (sin( 2.d0 * the ))**14.d0
        u2= - Psi_the / ae
        u1_t_tan = u0 * (sin( 2.d0 * the ))**13.d0 * 2.d0 * sin(the) * sin(the)   ! u1 * tan(the)
        u1_d_cos = u0 * (sin( 2.d0 * the ))**13.d0 * 2.d0 * sin(the)              ! u1 / cos(the)
        v = - Psi_c * cos(the0) * sin(dlam) / ae
        u1_the= 28.d0 * u0 * (sin(2.d0*the))**13.d0 * cos(2.d0*the)
        u2_lam= - Psi_the_lam / ae
        u2_the= - Psi_the_the / ae
        v_lam= - cos(the0) * ( Psi_cc * C_lam * sin(dlam) + Psi_c * cos(dlam) ) / ae
        v_the= - cos(the0) * sin(dlam) * Psi_cc * C_the /ae
        fcori_p= 2.d0 * omega * cos(the)
       

        !
        ! Fphi, h at the same position as velocity
        !
        A=  - fcori(ip) * omegasw * Psi_lam 
        B=  - (u1 + u2) * fcori(ip) * Psi_c * cos(the0) * sin(dlam) / ae
        ! C1= - v * ( fcori * u1 + u1 * u1 * tan(the) / ae )
        C1= - v * ( fcori(ip) +  u1_t_tan / ae) * u1
        D =   v * ( fcori_p * Psi  +  fcori(ip) * Psi_the ) / ae
        F_phi = A + B + C1 + D
        F_h   = F_phi/g0
        ! write(6,101) 'F_phi: A,B,C1,D', A,B,C1,D


        ! F_v
        A= -omegasw * v_lam
        B= u1_d_cos * v_lam / ae
        C1= 2.d0 * u2 * u1_t_tan / ae
        D= v * v_the / ae
        E= fcori(ip) * u2
        F= ( fcori_p * Psi  +  fcori(ip) * Psi_the ) / ae
        G= - u2 * ( - Psi_cc * cos(the0) * cos(the0) * sin(dlam) * sin(dlam) + &
                      Psi_c  *(cos(the)  * cos(the0) * cos(dlam) + sin(the) * sin(the0)) ) / ae / ae
        F_v = A + B + C1 + D + E + F + G


        ! F_u
        A= u1_d_cos * u2_lam / ae
        B= v * ( u1_the + u2_the ) / ae
        C1= - omegasw * u2_lam
        D= - u1_t_tan * v / ae
        E= - Psi_cc * Psi_c * C_the * C_the * cos(the0) * sin(dlam) / ae**3.d0
        F= - fcori(ip) * Psi_c * cos(the0) * sin(dlam) / ae
        G= - fcori(ip) * v
        F_u = A + B + C1 + D + E + F + G


        P_sph(1)=lam; P_sph(2)=the; P_sph(3)=1.d0
        call basis_between_sph_car (P_sph, basis, invbasis, ndim)
        !
        X(1)=F_u; X(2)=F_v; X(3)=0.d0
        call AX_mult(basis, X, Y, m, n, mp, np)
        AFV(ip,1) = AFV(ip,1) + F_h     ! addition upon dynamics
        AFV(ip,2) = AFV(ip,2) + Y(1)    ! p V_x / p t
        AFV(ip,3) = AFV(ip,3) + Y(2)    ! p V_y /
        AFV(ip,4) = AFV(ip,4) + Y(3)    ! p V_z /
        !
        !write(6, 101) 'F: del p (h,Vx,Vy,Vz)/pt', F_h, Y(1:3)


        !
        ! h
        !
        xa= -pi/2.d0; xb= the
        A= F_int_sinx_t_sin_2x_to_14 (xb) - F_int_sinx_t_sin_2x_to_14 (xa)
        B= F_int_tanx_t_sin_2x_to_28 (xb) - F_int_tanx_t_sin_2x_to_28 (xa)
        C1= h0 - ( ae * 2.d0 * omega * u0 * A + u0 * u0 * B ) / g0
        D= fcori(ip) * Psi / g0
        FV(ip,1) = C1 + D


        !
        ! u,v
        !
        X(1)= u1 + u2
        X(2)= v
        X(3)= 0.d0
        call AX_mult(basis, X, Y, m, n, mp, np)
        FV(ip,2) = Y(1)    !  V_x
        FV(ip,3) = Y(2)    !  V_y
        FV(ip,4) = Y(3)    !  V_z


        !     write(6, 121) 'ip=', ip
        !     write(6, 101) 'u1, u2, v', u1, u2, v
        !
        !     write(6, 101) 'lam, the, t, dlam, C, Psi',  lam, the, t, dlam, C, Psi
        !     write(6, 201) 'fcori(ip),C, Psi0, Psi', fcori(ip), C, Psi0, Psi
        !     if (ip .GT. 100)    stop 'ck Psi'
     enddo      ! loop ip
     !
     !
     !
  endif   ! if (stagger .EQ. 'A')



  if (stagger.EQ.'C') then
     !
     !-- step one for h field
     !
     do ip=1, nip
        lam= lon(ip)     ! velocity grid point
        the= lat(ip)
        !
        dlam= lam - lam0 - omegasw * t
        C= sin(the) * sin(the0) + cos(the) * cos(the0) * cos(dlam)
        Psi=    Psi0 * exp( -sigma * (1.d0 - C) / (1.d0 + C) )
        !
        Psi_c=  2.d0 * sigma / ( 1.d0 + C )**2.d0 * Psi
        Psi_cc=(4.d0 * sigma * sigma / ( 1.d0 + C )**4.d0 - 4.d0 * sigma / ( 1.d0 + C )**3.d0) * Psi
        C_the= cos(the) * sin(the0) - sin(the) * cos(the0) * cos(dlam)
        C_lam=                      - cos(the) * cos(the0) * sin(dlam)
        C_lam_lam=                  - cos(the) * cos(the0) * cos(dlam)
        C_the_lam=                    sin(the) * cos(the0) * sin(dlam)
        Psi_the= Psi_c * C_the
        Psi_lam= Psi_c * C_lam
        Psi_the_the= Psi_cc * C_the * C_the - Psi_c * C
        Psi_the_lam= Psi_cc * C_the * C_lam + Psi_c * C_the_lam
        Psi_lam_lam= Psi_cc * C_lam * C_lam + Psi_c * C_lam_lam
        !
        u1= u0 * (sin( 2.d0 * the ))**14.d0
        u2= - Psi_the / ae
        u1_t_tan = u0 * (sin( 2.d0 * the ))**13.d0 * 2.d0 * sin(the) * sin(the)   ! u1 * tan(the)
        u1_d_cos = u0 * (sin( 2.d0 * the ))**13.d0 * 2.d0 * sin(the)              ! u1 / cos(the)
        v = - Psi_c * cos(the0) * sin(dlam) / ae
        u1_the= 28.d0 * u0 * (sin(2.d0*the))**13.d0 * cos(2.d0*the)
        u2_lam= - Psi_the_lam / ae
        u2_the= - Psi_the_the / ae
        v_lam= - cos(the0) * ( Psi_cc * C_lam * sin(dlam) + Psi_c * cos(dlam) ) / ae
        v_the= - cos(the0) * sin(dlam) * Psi_cc * C_the /ae
        fcori_p= 2.d0 * omega * cos(the)
        !
        !
        ! Fphi, h at the same position as velocity
        !
        A=  - fcori(ip) * omegasw * Psi_lam 
        B=  - (u1 + u2) * fcori(ip) * Psi_c * cos(the0) * sin(dlam) / ae
        ! C1= - v * ( fcori * u1 + u1 * u1 * tan(the) / ae )
        C1= - v * ( fcori(ip) +  u1_t_tan / ae) * u1
        D =   v * ( fcori_p * Psi  +  fcori(ip) * Psi_the ) / ae
        F_phi = A + B + C1 + D
        F_h   = F_phi/g0
        AFV(ip,1) = AFV(ip,1) + F_h     ! addition upon dynamics
        !
        ! write(6,101) 'F_phi: A,B,C1,D', A,B,C1,D
        !write(6, 101) 'F: del p (h,Vx,Vy,Vz)/pt', F_h, Y(1:3)


        !
        ! h
        !
        xa= -pi/2.d0; xb= the
        A= F_int_sinx_t_sin_2x_to_14 (xb) - F_int_sinx_t_sin_2x_to_14 (xa)
        B= F_int_tanx_t_sin_2x_to_28 (xb) - F_int_tanx_t_sin_2x_to_28 (xa)
        C1= h0 - ( ae * 2.d0 * omega * u0 * A + u0 * u0 * B ) / g0
        D= fcori(ip) * Psi / g0
        FV(ip,1) = C1 + D

        !     write(6, 121) 'ip=', ip
        !     write(6, 101) 'u1, u2, v', u1, u2, v
        !
        !     write(6, 101) 'lam, the, t, dlam, C, Psi',  lam, the, t, dlam, C, Psi
        !     write(6, 201) 'fcori(ip),C, Psi0, Psi', fcori(ip), C, Psi0, Psi
        !     if (ip .GT. 100)    stop 'ck Psi'
        !
     enddo      ! loop ip=1, nip


     !
     !-- step two for velocity field
     !
     do ip=ips, ipe
        lam= lon_v(ip)     ! velocity grid point
        the= lat_v(ip)
        !
        dlam= lam - lam0 - omegasw * t
        C= sin(the) * sin(the0) + cos(the) * cos(the0) * cos(dlam)
        Psi=    Psi0 * exp( -sigma * (1.d0 - C) / (1.d0 + C) )

        Psi_c=  2.d0 * sigma / ( 1.d0 + C )**2.d0 * Psi
        Psi_cc=(4.d0 * sigma * sigma / ( 1.d0 + C )**4.d0 - 4.d0 * sigma / ( 1.d0 + C )**3.d0) * Psi
        C_the= cos(the) * sin(the0) - sin(the) * cos(the0) * cos(dlam)
        C_lam=                      - cos(the) * cos(the0) * sin(dlam)
        C_lam_lam=                  - cos(the) * cos(the0) * cos(dlam)
        C_the_lam=                    sin(the) * cos(the0) * sin(dlam)
        Psi_the= Psi_c * C_the
        Psi_lam= Psi_c * C_lam
        Psi_the_the= Psi_cc * C_the * C_the - Psi_c * C
        Psi_the_lam= Psi_cc * C_the * C_lam + Psi_c * C_the_lam
        Psi_lam_lam= Psi_cc * C_lam * C_lam + Psi_c * C_lam_lam

        u1= u0 * (sin( 2.d0 * the ))**14.d0
        u2= - Psi_the / ae
        u1_t_tan = u0 * (sin( 2.d0 * the ))**13.d0 * 2.d0 * sin(the) * sin(the)   ! u1 * tan(the)
        u1_d_cos = u0 * (sin( 2.d0 * the ))**13.d0 * 2.d0 * sin(the)              ! u1 / cos(the)
        v = - Psi_c * cos(the0) * sin(dlam) / ae
        u1_the= 28.d0 * u0 * (sin(2.d0*the))**13.d0 * cos(2.d0*the)
        u2_lam= - Psi_the_lam / ae
        u2_the= - Psi_the_the / ae
        v_lam= - cos(the0) * ( Psi_cc * C_lam * sin(dlam) + Psi_c * cos(dlam) ) / ae
        v_the= - cos(the0) * sin(dlam) * Psi_cc * C_the /ae
        fcori_p= 2.d0 * omega * cos(the)
       

        !
        ! Fphi, h at the same position as velocity
        !
        A=  - fcori(ip) * omegasw * Psi_lam 
        B=  - (u1 + u2) * fcori(ip) * Psi_c * cos(the0) * sin(dlam) / ae
        ! C1= - v * ( fcori * u1 + u1 * u1 * tan(the) / ae )
        C1= - v * ( fcori(ip) +  u1_t_tan / ae) * u1
        D =   v * ( fcori_p * Psi  +  fcori(ip) * Psi_the ) / ae
        F_phi = A + B + C1 + D
        F_h   = F_phi/g0
        ! write(6,101) 'F_phi: A,B,C1,D', A,B,C1,D


        ! F_v
        A= -omegasw * v_lam
        B= u1_d_cos * v_lam / ae
        C1= 2.d0 * u2 * u1_t_tan / ae
        D= v * v_the / ae
        E= fcori(ip) * u2
        F= ( fcori_p * Psi  +  fcori(ip) * Psi_the ) / ae
        G= - u2 * ( - Psi_cc * cos(the0) * cos(the0) * sin(dlam) * sin(dlam) + &
                      Psi_c  *(cos(the)  * cos(the0) * cos(dlam) + sin(the) * sin(the0)) ) / ae / ae
        F_v = A + B + C1 + D + E + F + G


        ! F_u
        A= u1_d_cos * u2_lam / ae
        B= v * ( u1_the + u2_the ) / ae
        C1= - omegasw * u2_lam
        D= - u1_t_tan * v / ae
        E= - Psi_cc * Psi_c * C_the * C_the * cos(the0) * sin(dlam) / ae**3.d0
        F= - fcori(ip) * Psi_c * cos(the0) * sin(dlam) / ae
        G= - fcori(ip) * v
        F_u = A + B + C1 + D + E + F + G


        P_sph(1)=lam; P_sph(2)=the; P_sph(3)=1.d0         ! on velocity site
        call basis_between_sph_car (P_sph, basis, invbasis, ndim)
        !
        X(1)=F_u; X(2)=F_v; X(3)=0.d0
        call AX_mult(basis, X, Y, m, n, mp, np)
        !
        ! u,v
        X(1)= u1 + u2
        X(2)= v
        X(3)= 0.d0
        call AX_mult(basis, X, YV, m, n, mp, np)


        !
        ! Proj velocity to N T direction for C-grid
        !      Un=(V,N)N,  Ut=(V,T)T
        !  
        ip2=C4E(1,ip); is2=C4E(5,ip)   ! Left:  1:4=i,j,k,Lp, 5=is
        vec1(1:3)= Nvec(1:3,is2,ip2)     ! for iE;  NoutdUn(is,ip)= + 1 
        call XY_dot (Y,vec1,n,ndim,s)
        AFV(ip,2) = AFV(ip,2) + s        !  F_n proj
        !vec1(1:3)=  Tvec(1:3,is,ip)   ! 
        !call XY_dot (Y,vec1,n,ndim,s2)
        !FV(ipt,3)=s2           !  F_t proj
        call XY_dot (YV,vec1,n,ndim,s)
        FV(ip,2) = s            !  U_n proj
        !
        !
     enddo      ! loop ip = ips,ipe   for both A and C
  endif         ! if (stagger .EQ. 'C')


  return
  include '../swm/myformat.inc'
end subroutine add_forcing_TC4



subroutine test_phi0
  use kinds, only: rt
  use module_constants, only  :  g, ae, pi, omega
  implicit none
  real(rt), external  :: F_int_sinx_t_sin_2x_to_14
  real(rt), external  :: F_int_tanx_t_sin_2x_to_28
  real(rt) :: x, xa, xb, dx, A, B, C1
  real(rt) :: u0, h0
  integer  :: npt, i
  !
  u0=20.d0
  h0=1.d5/g
  !
  xa= -pi/2.d0
  xb=  pi/2.d0
  npt= 101
  dx= (xb-xa)/dble(npt-1)
  !!write(7, 101) '#  theta     int ( sinx * sin(2x)^14 )     int ( tanx * (sin(2x)^28 )    h1(u1, v1=0)'
  do i=1, npt
     x= xa + dx * dble(i-1)
     A= F_int_sinx_t_sin_2x_to_14 (x) - F_int_sinx_t_sin_2x_to_14 (xa)
     B= F_int_tanx_t_sin_2x_to_28 (x) - F_int_tanx_t_sin_2x_to_28 (xa)
     C1= h0 - ( 2.d0 * ae * omega * u0 * A + u0 * u0 * B ) / g       !  phi1/g= h1 field
     !!write(7, 103)  x, A, B, C1
  end do
  return
  !
  include '../swm/myformat.inc'
end subroutine test_phi0



function F_int_sinx_t_sin_2x_to_14 (x) 
  use kinds, only: rt
  implicit none
  real(rt)  :: F_int_sinx_t_sin_2x_to_14
  real(rt), intent (in) :: x
!
!        integral sin(x) sin^14(2 x) dx = (2 cos^15(x) (30558493845 cos(2 x) - 18940553910 cos(4 x) + 8537085415 cos(6 x) - 2746165380 cos(8 x) + 601102677 cos(10 x) - 80604650 cos(12 x) + 5014575 cos(14 x) - 17951149788))/145422675 + constant  
!
!
  F_int_sinx_t_sin_2x_to_14= &
        (2.d0 * cos(x)**15.d0 * (30558493845.d0 * cos(2.d0 *x) - 18940553910.d0 * cos(4.d0 * x) + 8537085415.d0 * cos(6.d0 * x) - 2746165380.d0 * cos(8.d0 * x) + 601102677.d0 * cos(10.d0 * x) - 80604650.d0 * cos(12.d0 * x) + 5014575.d0 * cos(14.d0 * x) - 17951149788.d0))/145422675.d0 
  return
end function F_int_sinx_t_sin_2x_to_14


function F_int_tanx_t_sin_2x_to_28 (x) 
  use kinds, only: rt
  implicit none
  real(rt)  :: F_int_tanx_t_sin_2x_to_28
  real(rt), intent (in) :: x
!
!   integral   tan(x) * sin(2x)**28 =
!   (-5632695584460000 cos(2 x) + 187756519482000 cos(4 x) + 1627223168844000 cos(6 x) - 152552172079125 cos(8 x) - 732250425979800 cos(10 x) + 107683886173500 cos(12 x) + 338435070831000 cos(14 x) - 65806819328250 cos(16 x) - 146237376285000 cos(18 x) + 34635168067500 cos(20 x) + 56675729565000 cos(22 x) - 15585825630375 cos(24 x) - 19182554622000 cos(26 x) + 5937457383000 cos(28 x) + 5541626890800 cos(30 x) - 1889190985500 cos(32 x) - 1333546578000 cos(34 x) + 492832431000 cos(36 x) + 259385490000 cos(38 x) - 102673423125 cos(40 x) - 39113685000 cos(42 x) + 16427747700 cos(44 x) + 4285499400 cos(46 x) - 1895509350 cos(48 x) - 303281496 cos(50 x) + 140408100 cos(52 x) + 10400600 cos(54 x) - 5014575 cos(56 x))/37690512349593600
!
!   
  F_int_tanx_t_sin_2x_to_28 = &
       (-5632695584460000.d0 * cos(2.d0 * x) + 187756519482000.d0 * cos(4.d0 * x) + 1627223168844000.d0 * cos(6.d0 * x) - 152552172079125.d0 * cos(8.d0 * x) - 732250425979800.d0 * cos(10.d0 * x) + 107683886173500.d0 * cos(12.d0 * x) + 338435070831000.d0 * cos(14.d0 * x) - 65806819328250.d0 * cos(16.d0 * x) - 146237376285000.d0 * cos(18.d0 * x) + 34635168067500.d0 * cos(20.d0 * x) + 56675729565000.d0 * cos(22.d0 * x) - 15585825630375.d0 * cos(24.d0 * x) - 19182554622000.d0 * cos(26.d0 * x) + 5937457383000.d0 * cos(28.d0 * x) + 5541626890800.d0 * cos(30.d0 * x) - 1889190985500.d0 * cos(32.d0 * x) - 1333546578000.d0 * cos(34.d0 * x) + 492832431000.d0 * cos(36.d0 * x) + 259385490000.d0 * cos(38.d0 * x) - 102673423125.d0 * cos(40.d0 * x) - 39113685000.d0 * cos(42.d0 * x) + 16427747700.d0 * cos(44.d0 * x) + 4285499400.d0 * cos(46.d0 * x) - 1895509350.d0 * cos(48.d0 * x) - 303281496.d0 * cos(50.d0 * x) + 140408100.d0 * cos(52.d0 * x) + 10400600.d0 * cos(54.d0 * x) - 5014575.d0 * cos(56.d0 * x))/37690512349593600.d0
  !
  !
  return
end function F_int_tanx_t_sin_2x_to_28



subroutine  test_ERA40
  use netcdf
  implicit none
  character(len=100)  ::  filename, key
  integer :: nx, ny, nt
  real, allocatable  :: array(:,:,:)
  integer :: it, ix, iy

  filename = '/Volumes/Software/test/ERA40_1978Dec_500mb_huv.nc'
  key='z'

  call get_file_dimension(filename, nx, ny, nt)
  allocate (array(nx,ny,nt))

  write(6,*) 'nx,ny,nt', nx,ny,nt

  call get_2d_netcdf(filename, key, array, nx, ny, nt)

!  it=20
!  do ix=1, nx, 100
!     do iy=1, ny, 100
!        write(6, '(2x,a,2i6,f14.2)') 'ix,iy,array(ix,iy,it)', ix,iy,array(ix,iy,it)
!     enddo
!  enddo
  deallocate (array)
end subroutine test_ERA40







!--------   NETCDF-4  subroutines -------
!
!!subroutine nc4_r2d_wt (ndim,LDA,LDB,arr,FILE_NAME)
!!  use netcdf
!!  implicit none
!!  integer, intent(in) :: ndim, LDA, LDB
!!  real*8,  intent(in) :: arr(LDA, LDB)
!!  character (len = *), intent(in) :: FILE_NAME
!!  integer :: NX, NY
!!  integer :: ncid, varid, dimids(ndim)
!!  integer :: x_dimid, y_dimid
!!  integer :: chunks(2)
!!  integer :: shuffle, deflate, deflate_level
!!  integer :: x, y
!!  !
!!  NX=LDA; NY=LDB
!!  call cknc( nf90_create(FILE_NAME, nf90_hdf5, ncid) )
!!  call cknc( nf90_def_dim(ncid, "x", NX, x_dimid) )
!!  call cknc( nf90_def_dim(ncid, "y", NY, y_dimid) )
!!  dimids =  (/ x_dimid, y_dimid /)
!!  call cknc( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )
!!  chunks(1) = NX
!!  chunks(2) = NY
!!  call cknc( nf90_def_var_chunking(ncid, varid, 0, chunks) )
!!  shuffle = 1
!!  deflate = 1
!!  deflate_level = 1
!!  call cknc( nf90_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level) )
!!  call cknc( nf90_enddef(ncid) )
!!  !
!!  ! Write the pretend data to the file
!!  call cknc( nf90_put_var(ncid, varid, arr) )
!!  call cknc( nf90_close(ncid) )
!!  print *, '*** SUCCESS writing example file ', FILE_NAME, '!'
!!end subroutine nc4_r2d_wt


subroutine lon_lat_press_wt_nc4 (ndim,LDA,LDB,lons,lats,arr,FILE_NAME)
  use netcdf
  implicit none
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
!  character (len = *), parameter :: LON_UNITS = "radian_east"
!  character (len = *), parameter :: LAT_UNITS = "radian_north"
!  character (len = *), parameter :: PRES_NAME="pressure"
!  character (len = *), parameter :: PRES_UNITS = "hPa"
  character (len = *), parameter :: PRES_NAME="height"
  character (len = *), parameter :: PRES_UNITS = "m"
  character (len = *), parameter :: UNITS = "units"
  !
  integer, intent(in) :: ndim, LDA, LDB
  real*8,  intent(in) :: lons(LDA), lats(LDB), arr(LDA, LDB)
  character (len = *), intent(in) :: FILE_NAME
  integer :: NX, NY         ! alias
  integer :: NLONS, NLATS   ! alias 
  integer :: ncid, varid, dimids(ndim)
  integer :: lat_dimid, lon_dimid
  integer :: pres_varid
  integer :: lat_varid, lon_varid
  integer :: chunks(2)
  integer :: shuffle, deflate, deflate_level
  ! Loop indices
  integer :: lat, lon
  !
  NLONS=LDA; NLATS=LDB
  !! format: lons(1:NLONS); lats(1:NLATS); arr(1:NLONS, 1:NLATS)
  !
  call cknc( nf90_create(FILE_NAME, nf90_clobber, ncid) )
  call cknc( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
  call cknc( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
  call cknc( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
  call cknc( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
  call cknc( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
  call cknc( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
  dimids = (/ lon_dimid, lat_dimid /)
  call cknc( nf90_def_var(ncid, PRES_NAME, NF90_REAL, dimids, pres_varid) )
  call cknc( nf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )
  ! End define mode.
  call cknc( nf90_enddef(ncid) )
  call cknc( nf90_put_var(ncid, lat_varid, lats) )
  call cknc( nf90_put_var(ncid, lon_varid, lons) )
  call cknc( nf90_put_var(ncid, pres_varid, arr) )
  !call cknc( nf90_put_var(ncid, pres_varid, pres_out) )
  ! Close the file.
  call cknc( nf90_close(ncid) )
  print *,"*** SUCCESS writing example file sfc_pres_temp.nc!"
end subroutine lon_lat_press_wt_nc4





! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90
! $Id: pres_temp_4D_rd.f90,v 1.6 2006/12/09 18:44:58 russ Exp $
!!subroutine pres_temp_4D_rd
! for netcdf ERA40_1978Dec_500mb_huv
subroutine lon_lat_press_ERA40_rd_nc4 (ndim,LDA,LDB,ntime,itime,lons,lats,FV,FILE_NAME)
  use netcdf
  implicit none
  !
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
!  character (len = *), parameter :: LON_UNITS = "radian_east"
!  character (len = *), parameter :: LAT_UNITS = "radian_north"
!  character (len = *), parameter :: PRES_NAME="pressure"
!  character (len = *), parameter :: PRES_UNITS = "hPa"
  character (len = *), parameter :: PRES_NAME  ="Geopotential"
  character (len = *), parameter :: PRES_UNITS = "m**2 s**-2"
  character (len = *), parameter :: U_NAME = "U component of wind"
  character (len = *), parameter :: U_UNITS= "m s**-1"
  character (len = *), parameter :: V_NAME = "V component of wind"
  character (len = *), parameter :: V_UNITS= "m s**-1"
  character (len = *), parameter :: UNITS = "units"
  !

  integer, intent(in) :: ndim, LDA, LDB, ntime, itime
  real*8,  intent(out):: lons(LDA), lats(LDB), FV(LDA, LDB, 3)     ! four vector: h,u,v
  character (len = *), intent(in) :: FILE_NAME
  integer :: NX, NY                  ! alias
  integer :: NLONS, NLATS   ! alias 
  integer :: ncid, varid, dimids(ndim)
  integer :: lat_dimid, lon_dimid
  integer :: pres_varid
  integer :: lat_varid, lon_varid
  integer :: U_varid, V_varid
  integer :: chunks(2)
  integer :: shuffle, deflate, deflate_level
  ! Loop indices
  integer :: i, j   ! i for loop lon, lat
  !
  real*8, allocatable  :: pres_in(:,:,:), u_in(:,:,:), v_in(:,:,:)
  real*8  :: z_scale_factor, z_add_offset
  real*8  :: u_scale_factor, u_add_offset
  real*8  :: v_scale_factor, v_add_offset

  NLONS=LDA; NLATS=LDB
  allocate ( pres_in(ntime,LDA,LDB), u_in(ntime,LDA,LDB), v_in(ntime,LDA,LDB) )
  z_scale_factor = 0.198527364076114
  z_add_offset   = 51989.971048818
  u_scale_factor = 0.00160526870951765
  u_add_offset   = 14.6751693505085
  v_scale_factor = 0.00161786332839943
  v_add_offset   = -2.84689352639076


  !- YGYU test
  i = nf90_inq_varid(ncid, LAT_NAME,  lat_varid) 
  write(6, *) 'i=', i, 'finish nf90_inq_varid'
  STOP 'nail -1'



  call cknc( nf90_open(FILE_NAME, nf90_nowrite, ncid) )
  write(6,*) 'p1'
  !
  call cknc( nf90_inq_varid(ncid, LAT_NAME,  lat_varid) )
  write(6,*) 'p2.1'
  call cknc( nf90_inq_varid(ncid, LON_NAME,  lon_varid) )
  write(6,*) 'p2.2'
!  call cknc( nf90_inq_varid(ncid, PRES_NAME, pres_varid) )
!  write(6,*) 'p2.3'
  call cknc( nf90_inq_varid(ncid, U_NAME,    u_varid) )
  write(6,*) 'p2.4'
  call cknc( nf90_inq_varid(ncid, V_NAME,    v_varid) )
  write(6,*) 'p2.X'
  !
  call cknc( nf90_get_var  (ncid, lat_varid, lats) )
  call cknc( nf90_get_var  (ncid, lon_varid, lons) )
  call cknc( nf90_get_var  (ncid, pres_varid,pres_in) )
  call cknc( nf90_get_var  (ncid, u_varid,   u_in) )
  call cknc( nf90_get_var  (ncid, v_varid,   v_in) )
  write(6,*) 'p3'
  ! scale back 
  do i=1, LDA       ! lam
     do j=1, LDB    ! theta
        FV(i,j,1) = z_add_offset + pres_in (itime, j, i) * z_scale_factor
        FV(i,j,2) = u_add_offset + u_in    (itime, j, i) * u_scale_factor
        FV(i,j,3) = v_add_offset + v_in    (itime, j, i) * v_scale_factor
     enddo
  enddo
  !
  ! print
  write(6, *) 'lam'
  do i=1, LDA/20    ! lam
     write(6, *) lons(i)
  enddo
  do j=1, LDB/20    ! theta
     write(6, *) lats(j)
  enddo
  do j=1, LDB/20    ! theta
     do i=1, LDA/20    ! lam
        write(6, '(2x,a,2x,2i7,f12.5)') 'j,i,gh', j, i, pres_in(itime, j, i)         
     enddo
  enddo
  !
  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file.
  call cknc( nf90_close(ncid) )
  !
  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS reading example file ", FILE_NAME, "!"
  return
  !
end subroutine lon_lat_press_ERA40_rd_nc4




subroutine get_file_dimension(fileName, ix, jx, tdim)
  use netcdf
  implicit none
!!#include <netcdf.inc>
  !
  character(len=*), intent(in) :: fileName
  integer, intent(out) :: ix, jx, tdim
  integer  :: ncid , iret, dimid
  !
  iret = nf90_open(trim(fileName), NF90_NOWRITE, ncid)
  !
  iret = nf90_inq_dimid(ncid, "longitude", dimid)
  iret = nf90_inquire_dimension(ncid, dimid, len=ix)
  !
  iret = nf90_inq_dimid(ncid, "latitude", dimid)
  iret = nf90_inquire_dimension(ncid, dimid, len=jx)
  !
  iret = nf90_inq_dimid(ncid, "time", dimid)
  iret = nf90_inquire_dimension(ncid, dimid, len=tdim)
  !
  iret = nf90_close(ncid)
  write(6,*) "ix, jx, tdim = ", ix, jx, tdim
end subroutine get_file_dimension


subroutine get_1d_netcdf(filename, name, array, idim)
  use netcdf
  implicit none
  character(len=*), intent(in) :: name, filename
  integer, intent(in) :: idim
  real, dimension(idim), intent(out) :: array
  integer :: ncid
  !
  integer :: iret, varid
  real    :: scale_factor,   add_offset
  !
  call cknc (  nf90_open (trim(fileName), NF90_NOWRITE, ncid)  )
  call cknc (  nf90_inq_varid (ncid,  name,  varid) )
  call cknc (  nf90_get_var   (ncid, varid,  array) )
  !
  iret = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
  if(iret .eq. 0) array = array * scale_factor
  !
  iret = nf90_get_att(ncid, varid, 'add_offset', add_offset)
  if(iret .eq. 0) array = array + add_offset
  !
  iret = nf90_close(ncid)  
  !
end subroutine get_1d_netcdf


subroutine get_2d_netcdf(filename, name, array, idim, jdim)
  use netcdf
  implicit none
  character(len=*), intent(in) :: name, filename
  integer, intent(in) :: idim, jdim
  real, dimension(jdim,idim), intent(out) :: array
  integer :: ncid
  !
  integer :: iret, varid
  real    :: scale_factor,   add_offset
  !
  call cknc (  nf90_open      (trim(fileName), NF90_NOWRITE, ncid)  )
  call cknc (  nf90_inq_varid (ncid,  name,  varid) )
  call cknc (  nf90_get_var   (ncid, varid,  array) )
  !
  iret = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
  if(iret .eq. 0) array = array * scale_factor
  !
  iret = nf90_get_att(ncid, varid, 'add_offset', add_offset)
  if(iret .eq. 0) array = array + add_offset
  !
  iret = nf90_close(ncid)  
  !
end subroutine get_2d_netcdf


subroutine get_3d_netcdf(filename, name, array_3d, idim, jdim, tdim)
  use netcdf
  implicit none
  !
  character(len=*), intent(in) :: name, filename
  integer, intent(in) :: idim, jdim, tdim
  real, dimension(idim,jdim,tdim), intent(out) :: array_3d      
!!  real, dimension(tdim,jdim,idim)              :: array   ! convection: T,Z,Y,X  
  real, dimension(idim,jdim,tdim)              :: array   
  integer :: ncid
  ! loc
  integer :: iret, varid
  real    :: scale_factor, add_offset
  integer :: ix, iy, it
  !
  call cknc (  nf90_open (trim(fileName), nf90_nowrite, ncid)  )
  call cknc (  nf90_inq_varid (ncid,  name,  varid)  )
  call cknc (  nf90_get_var   (ncid, varid,  array)  )
  iret = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
!!  if(iret .eq. 0) array = array / scale_factor   ! * or /  ?
  if(iret .eq. 0) array = array * scale_factor   ! * or /  ?
  iret = nf90_get_att(ncid, varid, 'add_offset', add_offset)
!!  if(iret .eq. 0) array = array + add_offset
  if(iret .eq. 0) then
     if (trim(name).EQ.'z') then
        array = array + 1.d4 * 9.80616d0       ! 10 km * g
     else
        array = array + add_offset
     endif
  endif

!  do ix=1, idim
!     do iy=1, jdim
!        do it=1, tdim
!           array_3d(ix,iy,it) = array (it,iy,ix)
!        enddo
!     enddo
!  enddo
   array_3d(:,:,:) = array(:,:,:)
  !
  write(6,*) 'scale_factor, add_offset', scale_factor, add_offset
!
!  write(6,'(a)') 'array(:,:,20)'
!  do iy=1,jdim,10
!     write(6,'(6E14.6)') (array(ix,iy,20), ix=1,idim,20)
!     write(6,'(//)')
!  enddo
!
  call cknc (  nf90_close(ncid)  )
  return
end subroutine get_3d_netcdf


subroutine cknc(status)
  use netcdf
  implicit none
  integer, intent (in) :: status
  if(status /= nf90_noerr) then
     write(6, '(2x,a)')  trim(nf90_strerror(status))
     stop 2
  end if
end subroutine cknc






!----------------------------------------
!
!   E r r o r      N o r m
!
!----------------------------------------



!
!-- error between  || X_numeric - X_ana ||__1,2,inf
!
subroutine error_norm_scalar (ims, ime, ips, ipe, X_numeric, X_ana, L1n, L2n, Linfn)
  use kinds,         only: rt
  use readnamelist,  only: stagger
  use module_control,only: nip
  use module_variables,only: area, areaV_A, C4V
  implicit none
  integer, intent(in )::  ims, ime, ips, ipe
  real(rt),intent(in )::  X_numeric (ims:ime)
  real(rt),intent(in )::  X_ana     (ims:ime)
  real(rt),intent(out)::  L1n, L2n, Linfn
  !
  integer  :: ipt, ip, is, istart, iend
  real(rt) :: B1, B2, Binf
  real(rt) :: h0, hx, wt
  !
  L1n=0.d0; L2n=0.d0; Linfn=-1.d0
  B1=0.d0;  B2=0.d0;  Binf=-1.d0
  if (stagger .EQ. 'A') then
     istart=ips; iend=ipe
  elseif (stagger .EQ. 'C' .OR. stagger(1:2) .EQ. 'ZM') then
     istart=1; iend=nip
  elseif (stagger(1:2) .EQ. 'Du') then
     istart=ips; iend=ipe
  endif


  do ipt=istart,iend
     !
     if (stagger.EQ.'A' .OR. stagger.EQ.'C' .OR. stagger(1:2).EQ.'ZM') then
        wt=area(ipt)
     elseif (stagger(1:5) .EQ. 'DuSPR' .OR. stagger(1:5) .EQ. 'DuCVT') then
        ip= C4V(1,ipt)
        is= C4V(2,ipt)
        wt= areaV_A(is,ip)
     endif
     !
     !
     h0=X_ana(ipt)
     hx=X_numeric(ipt)
     L1n=  L1n+abs(hx-h0)*wt
     L2n=  L2n+(hx-h0)**2.d0*wt
     Linfn=max(Linfn, abs(hx-h0))
     !
     B1  = B1 + abs(h0)*wt
     B2  = B2 + h0*h0*wt
     Binf= max(Binf, abs(h0))
  enddo
  L1n=L1n/B1
  L2n=sqrt(L2n/B2)
  Linfn=Linfn/Binf
  return
end subroutine error_norm_scalar




subroutine error_norm_scalar_collect_bin (ims, ime, ips, ipe, X_numeric, X_ana, L1n, L2n, Linfn)
  use kinds,         only: rt
  use readnamelist,  only: stagger
  use module_control,only: nip, glvl
  use module_variables,only: area
  implicit none
  integer, intent(in )::  ims, ime, ips, ipe
  real(rt),intent(in )::  X_numeric (ims:ime)
  real(rt),intent(in )::  X_ana     (ims:ime)
  real(rt),intent(out)::  L1n, L2n, Linfn
  !
  integer  :: ip, istart, iend
  real(rt) :: B1, B2, Binf
  real(rt) :: h0, hx, wt
  !
  integer  :: nbin, ix
  real(rt) :: dx, xmin, xmax, X_err (ims:ime)
  real(rt), allocatable :: freq(:)


  !
  L1n=0.d0; L2n=0.d0; Linfn=-1.d0
  B1=0.d0;  B2=0.d0;  Binf=-1.d0
  if (stagger .EQ. 'A') then
     istart=ips; iend=ipe
  elseif (stagger .EQ. 'C') then
     istart=1; iend=nip
  endif
  !
!!  do ip=ips,ipe
  do ip=istart,iend
     !
     wt=area(ip)
     h0=X_ana(ip)
     hx=X_numeric(ip)
     L1n=  L1n+abs(hx-h0)*wt
     L2n=  L2n+(hx-h0)**2.d0*wt
     Linfn=max(Linfn, abs(hx-h0))
     !
     B1  = B1 + abs(h0)*wt
     B2  = B2 + h0*h0*wt
     Binf= max(Binf, abs(h0))
     !
  enddo
  L1n=L1n/B1
  L2n=sqrt(L2n/B2)
  Linfn=Linfn/Binf


  do ip=istart,iend
     h0=X_ana(ip)
     hx=X_numeric(ip)
     X_err(ip)=abs(hx-h0)/Binf
  enddo

  !
  !  collect bins
  !
  nbin=10*4**(glvl-2)
!  if (stagger .EQ. 'C') then
!     nbin= nbin * 20
!  endif
!  dx = Linfn / dble(nbin)
  dx = 1.d-3/ dble(nbin)
  allocate (freq(nbin))
  freq(1:nbin) = 0.d0
  do ix=1, nbin
     xmin= dble(ix-1)*dx
     xmax= dble(ix  )*dx
     do ip=istart,iend
        if ( X_err(ip).GE.xmin .AND.X_err(ip).LT.xmax ) then
           freq(ix) = freq(ix) + 1.d0
        endif
     enddo
     freq(ix) = freq(ix) / dble(iend-istart)
     write(6, '(2x,4E14.7)')  0.5d0 * (xmin+xmax),  freq(ix)
  enddo
  deallocate (freq)
  return
end subroutine error_norm_scalar_collect_bin





subroutine error_norm_vector (ims, ime, ips, ipe, LDB, ystart, yend, Xvec_numeric, Xvec_ana, L1n, L2n, Linfn)
  use kinds,         only: rt
  use readnamelist,  only: stagger
  use module_control,only: nip
  use module_variables,only: area, areaE, areaV_A, C4V
  integer, intent(in )::  ims, ime, ips, ipe, LDB
  integer, intent(in )::  ystart, yend
  real(rt),intent(in )::  Xvec_numeric (ims:ime,LDB)
  real(rt),intent(in )::  Xvec_ana     (ims:ime,LDB)
  real(rt),intent(out)::  L1n, L2n, Linfn
  !
  integer  :: ip, istart, iend, ipx, isx
  real(rt) :: B1, B2, Binf
  real(rt) :: v0(3), v1(3), dv(3)
  real(rt) :: s1, s2, t1, t2, wt
  !
  !
  L1n=0.d0; L2n=0.d0; Linfn=-1.d0
  B1=0.d0;  B2=0.d0;  Binf=-1.d0



  if (stagger.EQ.'A' .OR. stagger.EQ.'C' .OR. stagger(1:2) .EQ. 'ZM') then
     istart=ips; iend=ipe
  elseif (stagger(1:2) .EQ. 'Du') then
     istart=1; iend=nip
  endif
  !
  if (stagger.EQ.'A' .OR. stagger(1:5).EQ.'DuSPR' .OR. stagger(1:5).EQ.'DuCVT') then
     if (yend-ystart .NE. 2) STOP 'yend-ystart .NE. 2'
     if (LDB .LT. 3) STOP 'ndimFV < 3 in call error_norm_vector, STOP!'
     do ip=istart,iend
        wt = area (ip)
        v0(1:3)=Xvec_ana(ip,ystart:yend)
        v1(1:3)=Xvec_numeric(ip,ystart:yend)
        dv(1:3)=v1(1:3)-v0(1:3)
        !
        s1=sqrt(dv(1)*dv(1)+dv(2)*dv(2)+dv(3)*dv(3))     ! 1 norm  numerator
        s2=    (dv(1)*dv(1)+dv(2)*dv(2)+dv(3)*dv(3))     ! 2 norm  numerator
        !
        t1=sqrt(v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3))     ! 1 norm  dinominator
        t2=    (v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3))     ! 2 norm  dinominator
        !
        L1n=  L1n + s1*wt                                ! 1 norm  accumulator
        L2n=  L2n + s2*wt                                ! 2 norm  accumulator
        Linfn=max(Linfn, abs(s1))                        ! inf norm  accumulator
        !
        B1  = B1 + t1*wt                                 ! 1 norm  accumulator
        B2  = B2 + t2*wt                                 ! 2 norm  accum
        Binf= max(Binf, t1)                              ! inf norm  accum
     enddo
     L1n=L1n/B1
     L2n=sqrt(L2n/B2)
     Linfn=Linfn/Binf
  elseif (stagger(1:2).EQ.'ZM') then
     if (yend-ystart .NE. 2) STOP 'yend-ystart .NE. 2'
     if (LDB .LT. 3) STOP 'ndimFV < 3 in call error_norm_vector, STOP!'
     do ip=istart,iend
        !
        ipx= C4V(1,ip)
        isx= C4V(2,ip)
        wt = areaV_A (isx,ipx)
        !
        v0(1:3)=Xvec_ana(ip,ystart:yend)
        v1(1:3)=Xvec_numeric(ip,ystart:yend)
        dv(1:3)=v1(1:3)-v0(1:3)
        !
        s1=sqrt(dv(1)*dv(1)+dv(2)*dv(2)+dv(3)*dv(3))     ! 1 norm  numerator
        s2=    (dv(1)*dv(1)+dv(2)*dv(2)+dv(3)*dv(3))     ! 2 norm  numerator
        !
        t1=sqrt(v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3))     ! 1 norm  dinominator
        t2=    (v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3))     ! 2 norm  dinominator
        !
        L1n=  L1n + s1*wt                                ! 1 norm  accumulator
        L2n=  L2n + s2*wt                                ! 2 norm  accumulator
        Linfn=max(Linfn, abs(s1))                        ! inf norm  accumulator
        !
        B1  = B1 + t1*wt                                 ! 1 norm  accumulator
        B2  = B2 + t2*wt                                 ! 2 norm  accum
        Binf= max(Binf, t1)                              ! inf norm  accum
     enddo
     L1n=L1n/B1
     L2n=sqrt(L2n/B2)
     Linfn=Linfn/Binf
  elseif (stagger .EQ. 'C') then
     !!if (ystart .NE. 2 .OR. yend.NE.ystart) STOP 'yend-ystart .NE. 0'
     if (yend.NE.ystart) STOP 'yend-ystart .NE. 0'
     do ip=istart,iend
        wt = areaE (ip)
        v0(1:1)=Xvec_ana(ip,ystart:yend)
        v1(1:1)=Xvec_numeric(ip,ystart:yend)
        dv(1)=v1(1)-v0(1)
        !
        s1=sqrt(dv(1)*dv(1))     ! 1 norm  numerator
        s2=    (dv(1)*dv(1))     ! 2 norm  numerator
        !
        t1=sqrt(v0(1)*v0(1))     ! 1 norm  dinominator
        t2=    (v0(1)*v0(1))     ! 2 norm  dinominator
        !
        L1n=  L1n + s1*wt                                ! 1 norm  accumulator
        L2n=  L2n + s2*wt                                ! 2 norm  accumulator
        Linfn=max(Linfn, abs(s1))                        ! inf norm  accumulator
        !
        B1  = B1 + t1*wt                                 ! 1 norm  accumulator
        B2  = B2 + t2*wt                                 ! 2 norm  accum
        Binf= max(Binf, t1)                              ! inf norm  accum
     enddo
     L1n=L1n/B1
     L2n=sqrt(L2n/B2)
     Linfn=Linfn/Binf
  endif
  !
  return
end subroutine error_norm_vector





!-- LDA==nx*ny;  LDB=2
subroutine gen_uniform_ll (nx, ny, g_ll, LDA, LDB)
  use kinds, only: rt
  use module_constants, only : pi
  implicit none
  integer, intent (in)  :: nx, ny, LDA, LDB
  real(rt),intent (out) :: g_ll(LDA, LDB)
  !
  integer  :: i, j, k
  real(rt) :: dx, dy, ymin, ymax
  !
  if (LDA .NE. nx*ny .OR. LDB .NE. 2) STOP 'll(nx,ny) dimension wrong'
  dx=2.d0*pi/dble(nx-1)
  ymin=-pi/2.d0*0.9999d0   ! avoid pole problem   vs. 0.995
  ymax=-ymin
  dy=(ymax-ymin)/dble(ny-1)
  k=0
  do j=1,ny
     do i=1,nx
        k=k+1
        g_ll(k,1)= ymin + dy*dble(j-1)       ! lat
        g_ll(k,2)= dx*dble(i-1)              ! lon
!        g_ll(k,2)= -pi/2.d0 + dy*dble(j-1)
        !! write(6,141) 'k, g_ll(k,1),g_ll(k,2)=', k, g_ll(k,1),g_ll(k,2)
     enddo
  enddo
  include '../swm/myformat.inc'
end subroutine gen_uniform_ll



subroutine change_latlon_array_shape (LDA, LDB, nx, ny, g_ll, f_ll, lon_deg_plot, lat_deg_plot, field_plot)
  use kinds, only: rt
  use module_constants, only : pi
  implicit none
  integer, intent (in)  :: nx, ny, LDA, LDB
  real(rt),intent (out) :: g_ll(LDA, LDB), f_ll(LDA)
  real(rt),intent (out) :: lon_deg_plot(nx), lat_deg_plot(ny), field_plot(nx,ny)
  !
  integer  :: i, j, k
  real(rt) :: dx, dy, d2r
  !
  if (LDA .NE. nx*ny .OR. LDB .NE. 2) STOP 'll(nx,ny) dimension wrong'
!  dx=2.d0*pi/dble(nx-1); dy=pi/dble(ny-1)
!  do i=1,nx
!     lon_plot(i)=dx*dble(i-1)
!  enddo
!  do j=1,ny
!     lat_plot(j)=-pi/2.d0 + dy*dble(j-1)
!  enddo
!
  d2r=1.d0/180.d0*pi
  k=0
  do j=1,ny
     do i=1,nx
        k=k+1
        lon_deg_plot(i)=g_ll(k,2)/d2r     ! longitude
        lat_deg_plot(j)=g_ll(k,1)/d2r     ! latitude
     enddo
  enddo
  !
  k=0
  do j=1,ny
     do i=1,nx
        k=k+1
        field_plot(i,j)= f_ll(k)
     enddo
  enddo
!  do k=1, nx*ny
!     write(6,101) 'll(k,1),ll(k,2)=', ll(k,1),ll(k,2)
!  enddo
  include '../swm/myformat.inc'
end subroutine change_latlon_array_shape


!--- YGYU: using Ning's slint, but curves are not smooth, tolerate it now. 
!
subroutine icos_2_ll (nip, np_ll, Tdim, g_icos, f_icos, g_ll, f_ll)
  use kinds, only: rt
  use slint, only       :  slint_init, bilinear_interp
  implicit none
  integer, intent (in)  :: nip, np_ll, Tdim
  real(rt),intent (in)  :: g_icos(nip,Tdim), g_ll(np_ll,Tdim)
  real(rt),intent (in)  :: f_icos(nip)
  real(rt),intent (out) :: f_ll(np_ll)
  !
  integer :: ip
  !
  call slint_init (g_icos,nip,g_ll,np_ll)
  call bilinear_interp (f_icos, f_ll)
end subroutine icos_2_ll


subroutine output_icos_2_ll (nip, np_ll, Tdim, g_icos, f_icos, g_ll, f_ll, filename, fnlen)
  use kinds, only: rt
  use slint, only       :  slint_init, bilinear_interp
  implicit none
  integer, intent (in)  :: nip, np_ll, Tdim, fnlen
  real(rt),intent (in)  :: g_icos(nip,Tdim), g_ll(np_ll,Tdim)
  real(rt),intent (in)  :: f_icos(nip)
  character (len=fnlen), intent(in) :: filename
  real(rt),intent (out) :: f_ll(np_ll)
  !
  integer :: ip
  !
  call slint_init (g_icos,nip,g_ll,np_ll)  
  call bilinear_interp (f_icos, f_ll)
  open (11, file=trim(filename), status='unknown')
  do ip=1, np_ll
     write(11, 103) g_ll(ip,1),g_ll(ip,2),f_ll(ip)
  enddo
  close(11)
  include '../swm/myformat.inc'
end subroutine output_icos_2_ll




subroutine output_icos (LDA, LDB, Tdim, nip, gh_icos, f_icos, filename, fnlen)
  use kinds, only: rt
  implicit none
  integer, intent (in)  :: LDA, LDB, Tdim
  integer, intent (in)  :: nip, fnlen
  real(rt),intent (in)  :: gh_icos(LDA,Tdim)    ! grid  dim = nip
  real(rt),intent (in)  :: f_icos(LDB)        ! field dim = NPTS
  character (len=fnlen), intent(in) :: filename
  !
  integer :: ip
  !
  open (11, file=trim(filename), status='unknown')
  do ip=1, nip
!     write(11, 103) gh_icos(ip,1),gh_icos(ip,2),f_icos(ip)
     write(11, '(2f10.2,f25.8)') gh_icos(ip,1),gh_icos(ip,2),f_icos(ip)
  enddo
  close(11)
  include '../swm/myformat.inc'
end subroutine output_icos






!
!subroutine huv2uv (h,uh,vh,u,v,ims,ime,ips,ipe)
!  use kinds, only: rt
!  use module_constants,only : hbc
!  implicit none
!  integer, intent(in)  ::  ims,ime,ips,ipe
!  real(rt),intent(in)  ::  h(ims:ime),uh(ims:ime),vh(ims:ime) !
!  real(rt),intent(out) ::  u(ims:ime),v(ims:ime)
!  integer  :: ip
!  !
!  do ip=ips,ipe
!     if ( h(ip) .LT. 0.d0 ) then 
!        write(6,'(2x,a20,2x,i8,f15.1)') 'negative: ip, h(ip)=', ip, h(ip)
!!        STOP 'error: h(ip) < hbc'
!     endif
!     u(ip)= uh(ip)/h(ip)
!     v(ip)= vh(ip)/h(ip)
!  enddo
!  return
!end subroutine huv2uv
!


subroutine onestep(ims,ime,ips,ipe,fh,fuh,fvh,h0,uh0,vh0,h,uh,vh,u,v,dt)
  use kinds, only             :  rt
  use module_variables, only  :  nprox, prox, fcori    ! Coriolis param  f(ip)
  integer, intent(in) :: ims,ime,ips,ipe
  real(rt),intent(in) :: fh (ims:ime), fuh(ims:ime), fvh(ims:ime)
  real(rt),intent(in) :: h0 (ims:ime), uh0(ims:ime), vh0(ims:ime)
  real(rt),intent(in) :: dt
  real(rt),intent(out):: h  (ims:ime), uh (ims:ime), vh (ims:ime)
  real(rt),intent(out):: u  (ims:ime), v  (ims:ime)
  !
  integer             :: ip
  !
  do ip=ips,ipe
     do is=1, nprox(ip)
        h(ip) = h0(ip)  + fh(ip)*dt
        uh(ip)= uh0(ip) + fuh(ip)*dt
        vh(ip)= vh0(ip) + fvh(ip)*dt
        if ( h(ip) .LT. 1.d-10 ) then
           u(ip)=0.d0; v(ip)=0.d0
        else
           u(ip)= uh(ip)/h(ip)
           v(ip)= vh(ip)/h(ip)
        endif
     enddo
  enddo
  return
end subroutine onestep


subroutine test_exp_forcing (ims,ime,ips,ipe,h,u,v,fh,fuh,fvh)
  use kinds, only: rt
  implicit none
  integer, intent(in ):: ims, ime, ips, ipe
  real(rt),intent(in )::  h (ims:ime)    ! h (height) (m) (= H - hb )  wo terrain height
  real(rt),intent(in )::  u (ims:ime)    ! u
  real(rt),intent(in )::  v (ims:ime)    ! 
  real(rt),intent(out):: fh (ims:ime)    ! forcing on h
  real(rt),intent(out):: fuh(ims:ime)    ! uh
  real(rt),intent(out):: fvh(ims:ime)    ! vh
  integer  :: ip
  !
  do ip=ips,ipe
     fh (ip)= -h(ip)
     fuh(ip)= -2.d0*h(ip)*u(ip)
     fvh(ip)= -2.d0*h(ip)*v(ip)
  enddo
  return
end subroutine test_exp_forcing



!--------------------------------------------------
!   sec.4  Spherical Cooridnate test
!--------------------------------------------------
subroutine test_sphere
  use kinds, only: rt
  use module_constants, only : pi, ndim
  implicit none
  real(rt) :: P1(3), P2(3), P3(3)
  real(rt) :: s1, s2, s3
  real(rt) :: basis_sph(3,3), basis_car(3,3)
  real(rt) :: P_car(3), P_sph(3), wt(3)
  real(rt) :: Pmid_car(3), Pmid_sph(3)
  real(rt) :: Nvec(3), Tvec(3), Rvec(3)
  integer  :: i
  !
  !

  P1(1)=0.d0; P1(2)=0.d0;    P1(3)=1.d0
  P2(1)=0.d0; P2(2)=pi/4.d0; P2(3)=1.d0
  P3(1)=pi/3.d0; P3(2)=pi/3.d0; P3(3)=1.d0
  call arc_len(P1, P2, s1, ndim)
  call arc_len(P1, P3, s2, ndim)
  call arc_len(P2, P3, s3, ndim)
  write(6, 101) 'len s1,s2,s3=', s1,s2,s3

  P1(1)=0.d0; P1(2)=pi/2.d0;   P1(3)=1.d0
  !P2(1)=pi/4.d0; P2(2)=atan(1./sqrt(2.)); P2(3)=1.d0
  P2(1)=pi/4.d0; P2(2)=0.d0; P2(3)=1.d0
!  call arc_mid_ntvec (P1, P2, Pmid_car, Pmid_sph, Nvec, Tvec, ndim)
  call arc_mid_ntrvec (P1, P2, Pmid_car, Pmid_sph, Nvec, Tvec, Rvec, ndim)
  
  write(6, 101) 'P1,P2,Pmid_sph', P1,P2,Pmid_sph
  write(6, 101) 'Pmid_car(1:3),Nvec(1:3),Tvec(1:3)', Pmid_car(1:3),Nvec(1:3),Tvec(1:3)



!  call basis_from_sph_2_car(P1, basis_sph, ndim)
!  write(6, 102) 'P1 basis_sph(1:3,i)'
!  do i=1, 3
!     write(6,103) basis_sph(1:3, i)
!  enddo
!
!  call basis_from_sph_2_car(P2, basis_sph, ndim)
!  write(6, 102) 'P2 basis_sph(1:3,i)'
!  do i=1, 3
!     write(6,103) basis_sph(1:3, i)
!  enddo
!
!  call basis_from_sph_2_car(P3, basis_sph, ndim)
!  write(6, 102) 'P3 basis_sph(1:3,i)'
!  do i=1, 3
!     write(6,103) basis_sph(1:3, i)
!  enddo
!
  call vect_from_sph_2_car (P3, P_car, ndim)
  write(6, 102) 'P3 vect_car(1:3)'
  write(6,103) P_car(1:3)

!  call basis_from_car_2_sph(P_car, basis_car, ndim)
!  write(6, 102) 'P3 basis_car(1:3,i)'
!  do i=1, 3
!     write(6,103) basis_car(1:3, i)
!  enddo

  P1(1)=0.d0;    P1(2)=0.d0;   P1(3)=1.d0
  P2(1)=2.d0*pi/(10.d0*2.d0**5.d0); P2(2)=0.d0; P2(3)=1.d0
  P3(1)=0.d0;    P3(2)=P2(1);  P3(3)=1.d0
  call area_sph_triangle     (P1, P2, P3, s1, ndim)
  write(6,201) '1/8 spher area good=', s1
!  write(6,201) '1/8 spher area tan =', s2
!  write(6,201) '1/8 spher area diff=', s2-s1
!  STOP 'test_sphere'


  call arc_len(P1, P2, s1, ndim)
  call arc_len(P2, P3, s2, ndim)
  call arc_len(P3, P1, s3, ndim)
  write(6, 101) 'len s1,s2,s3=', s1,s2,s3
  
  call center_wt_sph_triangle (P1, P2, P3, P_sph, P_car, wt, ndim)
  write(6, 101) 'P_sph, P_car, wt'
  write(6, 103) P_sph
  write(6, 103) P_car
  write(6, 103) wt


  !
  return
  include '../swm/myformat.inc'
end subroutine test_sphere





!
!   Similar to dyn_A
!
subroutine TE_PE_calc_A (ims,ime,ips,ipe,FV,global_sum,ndimFV,itstep,iunit)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a
  use module_variables, only  :  TE0, PE0, Epot0
  use module_constants, only  :  g, ae, pi, ndim
  use module_control,   only  :  dt
  implicit none
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep, iunit
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  global_sum(6)              ! h*, TE, PE, vort, div
  !
  integer  :: j, k
  integer  :: ip, jp, kp
  integer  :: is, js, ks, isp, ism
  integer  :: nb
  !
  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: vort
  real(rt) :: s, t, z
  real(rt) :: div_V, div_hV  ! \div V or hV
  !
  real(rt) :: A1, A2, A3, A4, A5, A6
  
  
  !  A1:  1/A * \int depth h^*                                   Mass
  !  A2:  1/A * \int ( 1/2 h^* V V + 1/2 g [(h* + hs)^2 - hs^2]  Total energy
  !  A3:  1/A * \int (  (Xi + f)^2/ (2h^*)                       Potential enstrophy
  !  A4:  vorticity
  !  A5:  div_V
  !  A6:  Effective KE from eq.(25) in Tomita 
  !       1/A * \int { 1/2 h^* V V + 1/2 g [(h* + hs)^2 - hs^2]  - E_pot_0 }
  !
  !
  !- Williamson Eq.(136-140)
  !
  global_sum(1:6)=0.d0
  !
  do ip=ips,ipe
     nb= nprox(ip)
     ! A1:  h^*
     ! A2:  1/2 * h^* * V^2 + 1/2 * g * (h^2 - hs^2)
     A1= FV(ip,1)   
     A2= 0.5d0*A1*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4)) &
          + 0.5d0*g*((FV(ip,1)+hb_a(ip))*(FV(ip,1)+hb_a(ip)) - hb_a(ip)*hb_a(ip))   
     A6= 0.5d0*g*((FV(ip,1)+hb_a(ip))*(FV(ip,1)+hb_a(ip)) - hb_a(ip)*hb_a(ip))
     !
     !-- s1.  h and velocity from hex center to mid edge
     !
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        !-- weighted h,vx,vy,vz in Cartesian cord for upper corner
        do j=1, ndimFV
           Fc(j,is)=FV(ip,j)*wtsph(1,is,ip)+FV(jp,j)*wtsph(2,is,ip)+FV(kp,j)*wtsph(3,is,ip)  ! h^*,Vx,y,z
        enddo
     enddo
     !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        do j=1, ndimFV
           Fs(j,is)= 0.5d0*(Fc(j,ism) + Fc(j,is))      ! h^*,Vx,Vy,Vz
        enddo
     enddo
     
     !
     !-- s3.1 div(V)
     !-- s3.2 Vorticity = k dot \curl V
     !
     div_V=0.d0; vort=0.d0
     do is= 1, nb
        s=0.d0; t=0.d0
        do j=1,3    !  s= (V dot N); t= (V dot T)
           s= s + Fs(j+1,is)*Nvec(j,is,ip)   ! projection to N
           t= t + Fs(j+1,is)*Tvec(j,is,ip)   ! proj to T
        enddo
        div_V= div_V + s*slen(is,ip)
        vort= vort + t*slen(is,ip)     ! (v_i dot T) * ds_i
     enddo
     div_V= div_V/area(ip)
     vort= vort/area(ip)
     A4= vort
     A5= div_V
     A3= 0.5d0 * ( vort + fcori(ip) ) * ( vort + fcori(ip) ) / FV(ip,1)
     !
     !
     global_sum(1) = global_sum(1) + A1 * area(ip)
     global_sum(2) = global_sum(2) + A2 * area(ip)
     global_sum(3) = global_sum(3) + A3 * area(ip)
     global_sum(4) = global_sum(4) + A4 * area(ip)
     global_sum(5) = global_sum(5) + A5 * area(ip)
     global_sum(6) = global_sum(6) + A6 * area(ip)  !  Epot
     !
     !
  enddo   ! loop ip

  do j=1,6
     global_sum(j) = global_sum(j) / ( 4.d0 * pi * ae * ae )   ! normalize
  enddo

  write(iunit, 211) dt*real(itstep)/86400.d0, global_sum(1:5), &
       (global_sum(2)-TE0)/TE0, (global_sum(3)-PE0)/PE0,  &
       (global_sum(2)-TE0)/(TE0-Epot0)  ! effective Kinet Eng. variation ratio amplitude
                                        ! TE0 - Epot0 is Tomita eq.(25) LHS tot_energy 
  include '../swm/myformat.inc'
end subroutine TE_PE_calc_A


!
!   adopted from dyn_C
!
subroutine TE_PE_calc_C (ims,ime,ips,ipe,FV,global_sum,ndimFV,itstep,iunit)
  use kinds, only: rt
  use module_control,   only  :  dt
  use module_variables, only  :  lon, lat, lon_v, lat_v, nprox, prox, proxs, nip, npp, niE
  use module_variables, only  :  slen, areaI, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a   ! area(ip)=areaI(ip) MPAS
  use module_variables, only  :  areaE, areaV, dcenter, wtUt, RAiv4T, RAiv4H, iVtopkey
  use module_variables, only  :  Eindex, EindexVort, C4E, E4C, arctheta, LevdUn, NoutdUn
  use module_variables, only  :  omegasw, fcori_V
  use module_variables, only  :  TE0, PE0, Epot0
  use module_constants, only  :  g, ae, omega, pi, ndim
  use readnamelist,     only  :  iswcase, iconsv_scheme
  implicit none
  real(rt), parameter  :: third=0.3333333333333333333
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep, iunit
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  global_sum(6)
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  real(rt) :: hfx, divhvvx, divhvvy, phpx, phpy  ! see formula (1)
  integer  :: ip, jp, kp, Lp, iE, iE2, iV
  integer  :: is, isp, ism, is2, js, js2, ks, jsm     ! is plus 1 (is+1)
  integer  :: nb
  integer  :: i, j, k
  integer  :: lda, ldb, mdim
  integer  :: ipc, isc          !  aux (is,ip) circle for two hexagon
  integer  :: nvert, iEh, iph, ish   ! vortex, and Edge head

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, EL, ER, Ut, E1, E2, E3, E4 ! Ener,  E_corner, E0 prime = average Es site
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradEn, gradE_ana(3), gradKn, gradgh    ! Energy gradient, analytic solution
  real(rt) :: gradh(3), gradh_ana(3), gradh_ana_sph(3), gradh_sph(3)  ! h gradient, analytic solution
  real(rt) :: grad2norm                ! Energy gradient  ||  X - X_analytic ||_2  2 norm
  real(rt) :: vort, vort_ana, vortU, vortD    ! analytical sol
  real(rt) :: vortU2, vortD2
  real(rt) :: qeU, qeD
  real(rt) :: h_V_up, h_V_dn, h_V
  real(rt) :: factor_PV, Fe_perp

  real(rt) :: havg              ! analytical sol
  real(rt) :: s, t, sum, ss                ! sum
  real(rt) :: basis(3,3), ibasis(3,3)
  real(rt) :: vec(3), P_sph(3), vel_sph(3) ! 
  real(rt) :: VortV(10)
  real(rt) :: VortE(ims:ime)
  real(rt) :: fcori_av, vort_av, etaUt
  real(rt), allocatable :: qe_diag(:), qe_prog(:), h_e(:), Fe(:)
  real(rt), allocatable :: qe_hex_center(:), qe_hex_vertex(:,:), qe_hex_edge(:)
  real(rt), allocatable :: vort_hex_vertex(:,:), vort_hex_center(:)
  real(rt), allocatable :: Un_dot_gradPV(:), Ut_dot_gradPV(:)

  real(rt) :: A1, A2, A3, A4, A5, A6, areasum

  !--------------  Meta code ----------------
  !  continuity eq.
  !  for each ip
  !      div (hV) =  1/A * he Un dl
  !          he = average LR for edge
  !
  !  for each iE
  !    find ip (Left), jp (Right)
  !      direction dirivatives: grad(K + gh) =  grad (E)
  !          K = weighted Un on 12 edges; h    ??
  !          Ki, h from centers i
  !
  !
  !    KE non-conservative way
  !    find Vort on VcU, VcD   
  !       Vort_U = 1/A Un dot L ds 
  !       \vec L =  DA(1:5, ip)   donor-acceptor relation  \pm 1
  !       there should be a projector angle for R
  !    find local Ut:   linear combination of 12 edges 
  !
  !    KE conservative way
  !    term= eta * U_t^e
  !        = 1/de \sum_e'  W_ee' * l_e' * U_n^e' * 1/2 * (eta_e + eta_e')
  !    eta_e =  1/2 \sum_v  eta_v
  !    eta   =  Vort_U/D + f
  !    E_tot = \sum_e Ae [he U_n^e^2 / 2] + Ai * [ ghi * (1/2 hi * bi)]
  !----------------------------------------
  !    note
  !    FV(ip,1)=h
  !    FV(ip,2)=Un
  !    FV(ip,3)=Ut !  not used 
  !

  allocate (qe_hex_vertex(2,ipe))
  allocate (qe_hex_center(nip))
  allocate (qe_hex_edge  (ipe))
  allocate (vort_hex_vertex(2,ipe))
  allocate (vort_hex_center(nip))
  allocate (qe_diag      (ipe))
  allocate (h_e          (ipe))
  allocate (Fe           (ipe))
  allocate (Un_dot_gradPV(ipe))
  allocate (Ut_dot_gradPV(ipe))


  !
  !- Williamson Eq.(136-140)
  !
  global_sum(1:6)=0.d0
  areasum=0.d0  
  do ip=1,nip
     areasum= areasum + areaI(ip)
  enddo

  
  ! A1: h^*
  t=0.d0
  do ip=1,nip
     A1= FV(ip,1)
     t= t + A1 * areaI(ip)
  enddo
  global_sum(1) = t / areasum

!!
!!  -- A grid definition
!!
!!  ! A2: TE= 1/2 h^* V^2 + 1/2 g (h^2 - hs^2)
!!  t=0.d0
!!  do ip=1,nip
!!     nb= nprox(ip)
!!     s=0.d0
!!     do is= 1, nb          ! loop nn bond
!!        iE=E4C(is,ip)
!!        s= s + FV(iE,2) * FV(iE,2) * areaE(iE)
!!     enddo
!!     E1= s/areaI(ip)
!!     A2= 0.5d0*FV(ip,1)*E1 + 0.5d0*g*((FV(ip,1)+hb_a(ip))*(FV(ip,1)+hb_a(ip)) - hb_a(ip)*hb_a(ip))
!!     t= t + A2 * areaI(ip)
!!  enddo
!!  global_sum(2) = t / areasum
!!


  ! MPAS definition
  !
  ! A2: TE= \sum_e  le * de * (h_e * u_e**2 / 2)  + \sum_i  A_i * [ g * h_i * ( 1/2 h_i + hs) ]
  !                 le * de  = 4 * areaE(iE) 
  t=0.d0
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE)
     s = 0.5d0 * ( FV(ip,1) + FV(jp,1) ) * FV(iE,2) * FV(iE,2) * 0.5d0 * ( areaE(iE)*4.d0 )
     t = t + s
  enddo
  !
  ! --  continue accumulate t
  !                         ss is sum pot ener
  ss=0.d0
  do ip=1,nip
     s = g * FV(ip,1) * ( 0.5d0 * FV(ip,1) + hb_a(ip) ) * areaI(ip)
     t = t + s
     ss=ss + s
  enddo
  global_sum(2) = t / areasum
  global_sum(6) =ss / areasum


  ! A5: div_V
  t=0.d0
  do ip=1,nip
     nb= nprox(ip)
     s=0.d0
     do is= 1, nb          ! loop nn bond
        jp=prox(is,ip)
        iE=E4C(is,ip)
        s= s + FV(iE,2) * slen(is,ip) * NoutdUn(is,ip)        ! NOTE: could use (iE) only 
     enddo
     A5= s / areaI(ip)          ! div_V
     t= t + s
  enddo
  global_sum(5) = t / areasum



  !------------------------------------------------------------
  !  step.2 :
  !  prepare  vorticity = k \dot  curl V :  xiE(iE) by xiU, xiD
  !  output:  qe_diag(iE), h_e(iE), diag= diagnostic
  !
  !  algorithm for PV
  !    vertex --> edge
  !           --> center
  !               grad_n  grad_t 
  !               advect edge values
  !
  !               Add upwiding scheme for PV ,    
  !               d PV /dt =0   -->   \p PV \p t +  V \dot gradPV = 0
  !               output:  qe_diag(iE) correction
  !
  !               Note:  V \dot \grad PV
  !                    =  Un grad_n(PV) * \vec Un \dot \grad_n  +
  !                       Ut grad_t
  !------------------------------------------------------------
  !
  do iE=ips,ipe
     !
     !------  M1:  Q_e^ perp  = Fe^ perp * qe = 1/de \sum_e'  Wee' le' he Un^e' * 1/2 ( q_e + q_e' )
     !                    qe  = 1/2 * ( q_v1 + q_v2),    q_v = ( f + xi ) / h_v
     !        A= qe;   dA/dt=0;  \p qe / \p t = - V \dot \grad qe
     !
     ip=C4E(1,iE); jp=C4E(2,iE); kp=C4E(3,iE); Lp=C4E(4,iE); is=C4E(5,iE)
!!     write(6, 121) 'ip,jp,kp,Lp,is', ip,jp,kp,Lp,is
     h_V_up =  &
          RAiv4T(1,iE)*FV(ip,1) + &
          RAiv4T(2,iE)*FV(jp,1) + &
          RAiv4T(3,iE)*FV(kp,1)
     h_V_dn =  &
          RAiv4T(4,iE)*FV(ip,1) + &
          RAiv4T(5,iE)*FV(Lp,1) + &
          RAiv4T(6,iE)*FV(jp,1)
!
!     write(6, 101) 'RAiv4T(1:6,iE)', RAiv4T(1:6,iE)
!     write(6, 101) 'h_V_up', h_V_up
!     write(6, 101) 'h_V_dn', h_V_dn
!     write(6, *)
     !
     nb= nprox(ip)
     ism= mod(is+nb-2, nb) + 1
     isp= mod(is, nb)+1
     ks=  proxs(isp,ip)
     js=  proxs(is,ip)
     nb=  nprox(jp)
     jsm= mod(js+nb-2, nb) + 1
     !
     vort=0.d0
     do is2=1,3
        iE2= EindexVort(is2,iE)
        vort= vort + FV(iE2,2)*LevdUn(is2,iE)*dcenter(iE2)          ! LevdUn with fixed iE
!!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortU= vort / areaV(1,iE) + fcori_V(1,iE)   ! total rotation quantum for Up Triangle
     vortU2=vort / areaV(1,iE)
     !
     vort=0.d0
     do is2=4,6
        iE2= EindexVort(is2,iE)
        vort= vort + FV(iE2,2)*LevdUn(is2,iE)*dcenter(iE2)
!!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortD= vort / areaV(2,iE) + fcori_V(2,iE)   ! total rotation quantum for Down Triangle
     vortD2=vort / areaV(2,iE)
     !
!     write(6,*)
!     write(6, 121) 'iE', iE
!     write(6, 101) 'vortU, vortD, vort_mid', vortU, vortD, 0.5d0*(vortU+vortD)
!     write(6, 101) 'vortU_D middle analy  ', 2.d0*(omega + omegasw) * sin( lat_v(iE))
!     if (iE.EQ.30)  STOP 'test 2'

!     if ( h_V_up .lt. 1.d-16 .OR. h_V_dn .lt. 1.d-16 ) then
!        write(6, 101)  'FV(ip,1),FV(jp,1),FV(kp,1)', FV(ip,1),FV(jp,1),FV(kp,1)
!        STOP 'h field == 0 in C-grid; MPAS approach fails'
!     endif


     qeU = VortU / h_V_up      !!!  WHAT IF  h == 0 ????, STOP 'failed '
     qeD = VortD / h_V_dn      !!!  WHAT IF  h == 0 ????, STOP 'failed '
     !
     ! -- for each fixed iE, (is,ip)-- up header, 
     !    store point to three (is,ip), (jsm,jp), (ks,kp)
     !
     qe_hex_vertex(1,iE)= qeU
     qe_hex_vertex(2,iE)= qeD
     qe_hex_edge(iE)= 0.5d0 * (qeU + qeD)
     vort_hex_vertex(1,iE)= vortU2
     vort_hex_vertex(2,iE)= vortD2
     !
     !  this is diagnostic q_e  =  PV
     !
!     write(6, 201) 'qe_hex_edge(iE)       ', qe_hex_edge(iE)
!     write(6, 201) 'qe_hex_edge(iE) analy ', 2.d0*(omega + omegasw)*sin(lat_v(iE)) / h_e(iE)
     !
  enddo   ! loop iE


  ! find PV at center 
  do ip=1, nip
     nb= nprox(ip)
     s= 0.d0
     t= 0.d0
     do is=1, nb
        iE= E4C(is,ip)
        iV= iVtopkey(is,ip)    !      flip up or dn
        !  write(6, 121) 'is,ip,iE,iV', is,ip,iE,iV
        !
        ! -- I always have counter clock wise positive for vortex
        !    (is,ip) point to upper,  (iV, iE) should correspond to the same position
        s= s + qe_hex_vertex(iV,iE) * RAiv4H(is,ip)
        t= t + vort_hex_vertex(iV,iE) * RAiv4H(is,ip)
     enddo
     qe_hex_center(ip)= s
     vort_hex_center(ip)=t
     !
  enddo


  !
  !- A3:  1/2 * qe^2 * h^*
  !- A4:  vort
  !
  s= 0.d0
  t= 0.d0
  do ip=1, nip
     A3= 0.5d0 * qe_hex_center(ip) * qe_hex_center(ip) * FV(ip,1)
     s= s + A3 * areaI(ip)
     A4= vort_hex_center(ip)
     t= t + A4 * areaI(ip)
  enddo
  global_sum(3) = s / areasum
  global_sum(4) = t / areasum



  if (TE0.GT.1.d-14 .AND. PE0.GT.1.d-14) then
     write(iunit, 211) dt*real(itstep)/86400.d0, global_sum(1:5), &
          (global_sum(2)-TE0)/TE0, (global_sum(3)-PE0)/PE0, &
          (global_sum(2)-TE0)/(TE0-Epot0)        ! Tomita's def. effective KE
  else
!  forrtl: error (73): floating divide by zero
     write(iunit, 101) 'TE0 or PE0 .LT 1.d-4'
     write(iunit, 211) dt*real(itstep)/86400.d0, global_sum(1:5), &
          (global_sum(2)-TE0)/1.d-14, (global_sum(3)-PE0)/1.d-14, &
          (global_sum(2)-TE0)/1.d-14
  endif


  deallocate ( qe_diag, h_e, Fe )
  deallocate ( qe_hex_center, qe_hex_vertex, qe_hex_edge )
  deallocate ( vort_hex_vertex, vort_hex_center )
  deallocate ( Un_dot_gradPV, Ut_dot_gradPV )

  return

  include '../swm/myformat.inc'
end subroutine TE_PE_calc_C





!--------------------------------------------------
!   sec.6  Test \div  alpha \grad beta as in NICOM
!          for A-grid 
!--------------------------------------------------
!
!
!--- input (h,u,v) at cell center, 
!    output accerlation (fh,fu,fv) by finite volume method
!    test analytic solutions 
!
subroutine test_dyn_A_sph_harmonics (ims,ime,ips,ipe,FV,ndimFV)
  use kinds, only: rt
  use module_variables, only  :  lon, lat, nprox, prox, proxs, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a
  use module_variables, only  :  ASV                                       ! acceler six vectors (grad h, div, vort)
  use module_variables, only  :  omegasw
  use module_constants, only  :  g, ae, pi, ndim
  use readnamelist,     only  :  iswcase

  implicit none
  real(rt), parameter  :: third=0.3333333333333333333
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  real(rt) :: hfx, divhvvx, divhvvy, phpx, phpy  ! see formula (1)
  integer  :: ip, jp, kp
  integer  :: is, js, ks, isp, ism     ! is plus 1 (is+1)
  integer  :: nb
  integer  :: i, j, k
  integer  :: lda, ldb, mdim

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, E0, Ec(npp), E0p         ! Ener,  E_corner, E0 prime = average Es site
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradE(3), gradE_ana(3)   ! Energy gradient, analytic solution
  real(rt) :: gradh(3), gradh_ana(3), gradh_ana_sph(3), gradh_sph(3)  ! h gradient, analytic solution
  real(rt) :: grad2norm                ! Energy gradient  ||  X - X_analytic ||_2  2 norm
  real(rt) :: vort, vort_ana    ! analytical sol
  real(rt) :: havg              ! analytical sol
  real(rt) :: s, t, sum                 ! sum
  real(rt) :: basis(3,3), ibasis(3,3)
  real(rt) :: vec(3), P_sph(3), vel_sph(3) ! 
  real(rt) :: h_base, Phi, Phi2, a1, a2, a3
  real(rt), allocatable :: TestF(:)   ! test field  TF= - omegasw * x  
  real(rt) :: div_V   ! \div V
  
  !--------------  Meta code ----------------
  !  for each ip; specify h, V
  !     check  grad h, div V,  k \dot curl V
  !
  !    FV(ip,1)=h
  !    FV(ip,2/3/4)=Vx/Vy/Vz
  !

  allocate (TestF(ipe))     ! test field F = - omega \dot r
  mdim=ndim; lda=ndim; ldb=ndim
  do ip=ips,ipe
     grad2norm=0.d0
     nb= nprox(ip)
     E0 = 0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*(FV(ip,1)+hb_a(ip))
     h_base = FV(ip,1)
     s=0.d0
     !
     !-- s1.  h and velocity from hex center to mid edge
     !
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        !-- weighted h,vx,vy,vz in Cartesian cord for upper corner
        do j=1, ndimFV
           Fc(j,is)=FV(ip,j)*wtsph(1,is,ip)+FV(jp,j)*wtsph(2,is,ip)+FV(kp,j)*wtsph(3,is,ip)
        enddo
        !-- interpolate E = 1/2*v*v+g(h+h0) to prepare for grad(E)
        Ec(is)= &
             (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*(FV(ip,1)+hb_a(ip)))*wtsph(1,is,ip)+&
             (0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4))+g*(FV(jp,1)+hb_a(jp)))*wtsph(2,is,ip)+&
             (0.5d0*(FV(kp,2)*FV(kp,2)+FV(kp,3)*FV(kp,3)+FV(kp,4)*FV(kp,4))+g*(FV(kp,1)+hb_a(kp)))*wtsph(3,is,ip)
        s=s+Ec(is)
     enddo
     E0p=s/dble(nb)  !  E0p = sum ( E_corner + E_corner_m1 ) / 2
     !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        do j=1, ndimFV
           Fs(j,is)= 0.5d0*(Fc(j,ism) + Fc(j,is))      ! h,Vx,Vy,Vz
        enddo
     enddo


     !
     !-- s2. gradient h
     !
     havg=0.d0
     do is=1,nb
        havg= havg + Fs(1,is)
     enddo
     havg= havg/dble(nb)
!
!     write(6, 121) 'ip=', ip
!     write(6, 101) 'h_base, h_side= Fs(1,1:nb)=', h_base
!     write(6, 105)  Fs(1,1:nb)

     do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds 
        s=0.d0
        do is= 1, nb
           s= s + (Fs(1,is)-h_base)*Nvec(j,is,ip)*slen(is,ip) ! (h-hc)* n ds; 1=h, 2_4=Vx_z
           !s= s + (Fs(1,is)-havg)*Nvec(j,is,ip)*slen(is,ip) ! (h-h_avg)* n ds; 1=h, 2_4=Vx_z
        enddo
        gradh(j)= s/area(ip)
     enddo
     !   remove radial component for gradient
     s=0.d0                      ! s= (gradh, R)
     do j=1, 3
        s= s+ gradh(j)*Rcvec(j,ip)
     enddo
     do j=1, 3
        gradh(j)= gradh(j)-s*Rcvec(j,ip)    ! grad = grad - (grad,R) R
     enddo


     !
     !-- s3. div(V), Vorticity = k dot \curl V
     !
     div_V=0.d0; vort=0.d0
     do is= 1, nb
        s=0.d0; t=0.d0
        do j=1,3    !  s= (V dot N); t= (V dot T)
           s= s + Fs(j+1,is)*Nvec(j,is,ip)   ! projection to N
           t= t + Fs(j+1,is)*Tvec(j,is,ip)   ! proj to T
        enddo
        div_V= div_V + s*slen(is,ip)
        vort= vort + t*slen(is,ip)     ! (v_i dot T) * ds_i
     enddo
     div_V= div_V/area(ip)
     vort= vort/area(ip)
     !
     ASV(ip,1)=gradh(1)
     ASV(ip,2)=gradh(2)
     ASV(ip,3)=gradh(3)
     ASV(ip,4)=div_V
     ASV(ip,5)=vort
     !
  enddo  ! ip  outer loop
  !
  !
  deallocate (TestF)
  return
  include '../swm/myformat.inc'
end subroutine test_dyn_A_sph_harmonics


!
!  h = beta
!  func1 :   grad (h)
!  func2 :   div ( hV )
!  func3 :   (f + xi ) k \cross V
!
!  layout:   ASV(nip, 1:12) ==  grad_h, (div_V, vort, div_hV), (f + xi ) k \cross V, \dot(\vec{V})
!
subroutine test_dyn_A_3func (ims,ime,ips,ipe,FV,ndimFV)
  use kinds, only: rt
  use module_variables, only  :  lon, lat, nprox, prox, proxs, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a
  use module_variables, only  :  ASV                                       ! acceler six vectors (grad h, div, vort)
  use module_variables, only  :  omegasw
  use module_constants, only  :  g, ae, pi, ndim
  use readnamelist,     only  :  iswcase

  implicit none
  real(rt), parameter  :: third=0.3333333333333333333
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  real(rt) :: hfx, divhvvx, divhvvy, phpx, phpy  ! see formula (1)
  integer  :: ip, jp, kp
  integer  :: is, js, ks, isp, ism     ! is plus 1 (is+1)
  integer  :: nb
  integer  :: i, j, k
  integer  :: lda, ldb, mdim

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, E0, Ec(npp), E0p         ! Ener,  E_corner, E0 prime = average Es site
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradE(3), gradE_ana(3)   ! Energy gradient, analytic solution
  real(rt) :: gradh(3), gradh_ana(3), gradh_ana_sph(3), gradh_sph(3)  ! h gradient, analytic solution
  real(rt) :: grad2norm                ! Energy gradient  ||  X - X_analytic ||_2  2 norm
  real(rt) :: vort, vort_ana    ! analytical sol
  real(rt) :: havg              ! analytical sol
  real(rt) :: s, t, sum                 ! sum
  real(rt) :: basis(3,3), ibasis(3,3)
  real(rt) :: vec(3), P_sph(3), vel_sph(3) ! 
  real(rt) :: h_base, Phi, Phi2, a1, a2, a3
  real(rt), allocatable :: TestF(:)   ! test field  TF= - omegasw * x  
  real(rt) :: div_V   ! \div V
  real(rt) :: div_hV   ! \div hV
  
  !--------------  Meta code ----------------
  !  for each ip; specify h, V
  !     check  grad h, div (hV),   (f + xi ) k \cross  V
  !
  !    FV(ip,1)=h
  !    FV(ip,2/3/4)=Vx/Vy/Vz
  !

  allocate (TestF(ipe))     ! test field F = - omega \dot r
  mdim=ndim; lda=ndim; ldb=ndim
  do ip=ips,ipe
     grad2norm=0.d0
     nb= nprox(ip)
     E0 = 0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*(FV(ip,1)+hb_a(ip))
     h_base = FV(ip,1)
     s=0.d0
     !
     !-- s1.  h and velocity from hex center to mid edge
     !
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        !-- weighted h,vx,vy,vz in Cartesian cord for upper corner
        do j=1, ndimFV
           Fc(j,is)=FV(ip,j)*wtsph(1,is,ip)+FV(jp,j)*wtsph(2,is,ip)+FV(kp,j)*wtsph(3,is,ip)
        enddo
        !-- interpolate E = 1/2*v*v+g(h+h0) to prepare for grad(E)
        !Ec(is)= &
        !     (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*(FV(ip,1)+hb_a(ip)))*wtsph(1,is,ip)+&
        !     (0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4))+g*(FV(jp,1)+hb_a(jp)))*wtsph(2,is,ip)+&
        !     (0.5d0*(FV(kp,2)*FV(kp,2)+FV(kp,3)*FV(kp,3)+FV(kp,4)*FV(kp,4))+g*(FV(kp,1)+hb_a(kp)))*wtsph(3,is,ip)
        !
        !-- for Kinetic only, test only, fake
        !
        Ec(is)= &
             (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4)))*wtsph(1,is,ip)+&
             (0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4)))*wtsph(2,is,ip)+&
             (0.5d0*(FV(kp,2)*FV(kp,2)+FV(kp,3)*FV(kp,3)+FV(kp,4)*FV(kp,4)))*wtsph(3,is,ip)
        s=s+Ec(is)
     enddo
     !! E0p=s/dble(nb)  !  E0p = sum ( E_corner + E_corner_m1 ) / 2
     !! use definition for gradient operator, which is easy to analyze theoretically
     !E0p=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*(FV(ip,1)+hb_a(ip))
     !
     !-- for Kinetic only, fake
     !
     E0p=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))


     !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        do j=1, ndimFV
           Fs(j,is)= 0.5d0*(Fc(j,ism) + Fc(j,is))      ! h,Vx,Vy,Vz
        enddo
     enddo


     !
     !-- s2.1 gradient h
     !-- s2.2 Kinetic energy Gradient :  K + g (h + hb_a)
     !        \grad E_j = 1/A  *  \sum_i (Es_i -E0) * n(j,i)  dsi
     !-- s2.3 remove radial component
     !
     !
     !-- use h_base
     !!havg=0.d0
     !!do is=1,nb
     !!   havg= havg + Fs(1,is)
     !!enddo
     !!havg= havg/dble(nb)
!
!     write(6, 121) 'ip=', ip
!     write(6, 101) 'h_base, h_side= Fs(1,1:nb)=', h_base
!     write(6, 105)  Fs(1,1:nb)


     do is= 1, nb
        ism=mod(is+nb-2, nb)+1                
        Es(is)= 0.5d0*(Ec(ism)+Ec(is))    ! KE on six sides
     enddo
     !
     do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds 
        s=0.d0; t=0.d0
        do is= 1, nb
           s= s + (Fs(1,is)-h_base)*Nvec(j,is,ip)*slen(is,ip) ! (h-hc)* n ds; 1=h, 2_4=Vx_z
           !s= s + (Fs(1,is)-havg)*Nvec(j,is,ip)*slen(is,ip) ! (h-h_avg)* n ds; 1=h, 2_4=Vx_z
           t= t + (Es(is) - E0p)*Nvec(j,is,ip)*slen(is,ip) ! E0p: prime, replace E0 by six-side average
        enddo
        gradh(j)= s/area(ip)          ! this is grad(h^*)
        gradE(j)= t/area(ip)          ! grad(1/2*V*V + g(h+hs))
     enddo
     !
     !   remove radial component, \vec R component  :  \vec F = \vec F  - (F dot R) \vec R
     s=0.d0; t=0.d0               ! s= (gradh, Rcvec) not Rvec on edge
     do j=1, 3
        s= s+ gradh(j)*Rcvec(j,ip)
        t= t+ gradE(j)*Rcvec(j,ip)
     enddo
     do j=1, 3
        gradh(j)= gradh(j)-s*Rcvec(j,ip)    ! grad = grad - (grad,R) R
        gradE(j)= gradE(j)-t*Rcvec(j,ip)
     enddo


     !
     !-- s3. div( V), Vorticity = k dot \curl V
     !
     div_V=0.d0; vort=0.d0
     do is= 1, nb
        s=0.d0; t=0.d0
        do j=1,3    !  s= (V dot N); t= (V dot T)
           s= s + Fs(j+1,is)*Nvec(j,is,ip)   ! projection to N
           t= t + Fs(j+1,is)*Tvec(j,is,ip)   ! proj to T
        enddo
        div_V= div_V + s*slen(is,ip)
        vort= vort + t*slen(is,ip)     ! (v_i dot T) * ds_i
     enddo
     div_V= div_V/area(ip)
     vort= vort/area(ip)


     !
     !-- s4.  div (hV) =  \int  (hV \dot N) * ds  / A
     !
     div_hV=0.d0
     do is= 1, nb
        s=0.d0
        do j=1,3    !  s= (V dot N)
           s= s + Fs(j+1,is)*Nvec(j,is,ip)      ! projection to N
        enddo
        div_hV=div_hV + s*Fs(1,is)*slen(is,ip)  ! h^*  (V dot N)
     enddo
     div_hV= div_hV/area(ip)     


     !
     !-- s5.  acceleration for velocity
     !----------------------------------------------------------------------
     ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
     !
     !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
     !----------------------------------------------------------------------
     !
     do i=1,3
        j=mod(i,3)+1; k=mod(j,3)+1
        vec(i)=  (vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))
        ! wo corriorious force
        !        AFV(ip,i+1)=-(vort)*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
     enddo
     !
     !
     ASV(ip,1)=gradh(1)
     ASV(ip,2)=gradh(2)
     ASV(ip,3)=gradh(3)
     ASV(ip,4)=div_V
     ASV(ip,5)=vort
     ASV(ip,6)=div_hV
     ASV(ip,7)=vec(1)
     ASV(ip,8)=vec(2)
     ASV(ip,9)=vec(3)
     !
     ! \dot ( \vec{V} )   full-term
     !
     do i=1,3
        j=mod(i,3)+1; k=mod(j,3)+1
        vec(i)=  -(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
     enddo
     ASV(ip,10)=vec(1)
     ASV(ip,11)=vec(2)
     ASV(ip,12)=vec(3)
     !
     ASV(ip,10)=gradE(1)
     ASV(ip,11)=gradE(2)
     ASV(ip,12)=gradE(3)
     !
  enddo  ! ip  outer loop
  !
  !
  deallocate (TestF)
  return
  include '../swm/myformat.inc'
end subroutine test_dyn_A_3func



!
! adopted from dyn_C
!
! input:  FV
! output: ASV
! layout: ASV(niE, 1:12) ==  
!         (grad_h,0,0), (div_V, vort, div_hV), (eta*Ut,Un,Ut), (\part Un / \part t,0,0)
!
subroutine test_dyn_C_3func (ims,ime,ips,ipe,FV,ndimFV)
  use kinds, only: rt
  use module_control,   only  :  dt
  use module_variables, only  :  lon, lat, lon_v, lat_v, nprox, prox, proxs, nip, npp, niE
  use module_variables, only  :  slen, areaI, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a   ! area(ip)=areaI(ip) MPAS
  use module_variables, only  :  areaE, areaV, dcenter, wtUt, RAiv4T, RAiv4H, iVtopkey
  use module_variables, only  :  Eindex, EindexVort, C4E, E4C, arctheta, LevdUn, NoutdUn
  use module_variables, only  :  omegasw, fcori_V
  use module_constants, only  :  g, ae, omega, pi, ndim
  use module_variables, only  :  ASV                                  ! acceler six vectors (grad h, div, vort)
  use readnamelist,     only  :  iswcase, iconsv_scheme
  implicit none
  real(rt), parameter  :: third=0.3333333333333333333
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  real(rt) :: hfx, divhvvx, divhvvy, phpx, phpy  ! see formula (1)
  integer  :: ip, jp, kp, Lp, iE, iE2, iV
  integer  :: is, isp, ism, is2, js, js2, ks, jsm     ! is plus 1 (is+1)
  integer  :: nb
  integer  :: i, j, k
  integer  :: lda, ldb, mdim
  integer  :: ipc, isc          !  aux (is,ip) circle for two hexagon
  integer  :: nvert, iEh, iph, ish   ! vortex, and Edge head

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, EL, ER, Ut, E1, E2, E3, E4 ! Ener,  E_corner, E0 prime = average Es site
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradEn, gradE_ana(3), gradKn, gradgh    ! Energy gradient, analytic solution
  real(rt) :: gradh(3), gradh_ana(3), gradh_ana_sph(3), gradh_sph(3)  ! h gradient, analytic solution
  real(rt) :: grad2norm                ! Energy gradient  ||  X - X_analytic ||_2  2 norm
  real(rt) :: vort, vort_ana, vortU, vortD    ! analytical sol
  real(rt) :: vortU2, vortD2
  real(rt) :: qeU, qeD
  real(rt) :: h_V_up, h_V_dn, h_V
  real(rt) :: factor_PV, Fe_perp

  real(rt) :: havg              ! analytical sol
  real(rt) :: s, t, sum                 ! sum
  real(rt) :: basis(3,3), ibasis(3,3)
  real(rt) :: vec(3), P_sph(3), vel_sph(3) ! 
  real(rt) :: VortV(10)
  real(rt) :: VortE(ims:ime)
  real(rt) :: fcori_av, vort_av, etaUt
  real(rt), allocatable :: qe_diag(:), qe_prog(:), h_e(:), Fe(:)
  real(rt), allocatable :: qe_hex_center(:), qe_hex_vertex(:,:), qe_hex_edge(:)
  real(rt), allocatable :: vort_hex_vertex(:,:), vort_hex_center(:)
  real(rt), allocatable :: Un_dot_gradPV(:), Ut_dot_gradPV(:)
  real(rt) :: A1, A2, A3, A4, A5, A6, areasum

  !--------------  Meta code ----------------
  !  continuity eq.
  !  for each ip
  !      div (hV) =  1/A * he Un dl
  !          he = average LR for edge
  !
  !  for each iE
  !    find ip (Left), jp (Right)
  !      direction dirivatives: grad(K + gh) =  grad (E)
  !          K = weighted Un on 12 edges; h    ??
  !          Ki, h from centers i
  !
  !
  !    KE non-conservative way
  !    find Vort on VcU, VcD   
  !       Vort_U = 1/A Un dot L ds 
  !       \vec L =  DA(1:5, ip)   donor-acceptor relation  \pm 1
  !       there should be a projector angle for R
  !    find local Ut:   linear combination of 12 edges 
  !
  !    KE conservative way
  !    term= eta * U_t^e
  !        = 1/de \sum_e'  W_ee' * l_e' * U_n^e' * 1/2 * (eta_e + eta_e')
  !    eta_e =  1/2 \sum_v  eta_v
  !    eta   =  Vort_U/D + f
  !    E_tot = \sum_e Ae [he U_n^e^2 / 2] + Ai * [ ghi * (1/2 hi * bi)]
  !----------------------------------------
  !    note
  !    FV(ip,1)=h
  !    FV(ip,2)=Un
  !    FV(ip,3)=Ut !  not used 
  !

  allocate (qe_hex_vertex(2,ipe))
  allocate (qe_hex_center(nip))
  allocate (qe_hex_edge  (ipe))
  allocate (vort_hex_vertex(2,ipe))
  allocate (vort_hex_center(nip))
  allocate (qe_diag      (ipe))
  allocate (h_e          (ipe))
  allocate (Fe           (ipe))
  allocate (Un_dot_gradPV(ipe))
  allocate (Ut_dot_gradPV(ipe))


  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)
     nb= nprox(ip)
     ism=mod(is+nb-2, nb)+1
     !---------------------------------------------
     ! U_n \dot grad E  =  - (E_L - E_R)/dis * (n^out \dot  U_n)
     !       caution:      direction   n^out \dot  U_n
     !       K          =  K_i from n-component only  L, R
     !       E_ip (L,R) =  wt *  En(iE)^2  + gh
     !           areaE  =  1/4 * l_e * d_e
     !---------------------------------------------  
     !
     !-- Left
     !
     E3= FV(ip,1)
     !-- Right
     E4= FV(jp,1)
     gradgh =  (E4 - E3) / dcenter(iE) * NoutdUn(is,ip)    !   \grad(h) * U_n^e
     ASV(iE,1)= gradgh
     ASV(iE,2)= 0.d0
     ASV(iE,3)= 0.d0
  enddo


  !---------------------------------------------
  !  step.1 :  tendency for h field at center of Control Volume
  !  - \p h / \p t =  div_( hV ) = 1/A \int  h_e * Un * le  * (U_n^e \dot N_out)
  !            h_e = (h_L + h_R)/2 
  !---------------------------------------------
  ! 
  ! MPAS Fe = h_e * U_n ^ e
  ! use as single variable
  !
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE); kp=C4E(3,iE); Lp=C4E(4,iE); is=C4E(5,iE)
     Fe(iE) = 0.5d0 * ( FV(ip,1) + FV(jp,1) ) * FV(iE,2)
     h_e(iE)= 0.5d0 * ( FV(ip,1) + FV(jp,1) )
  enddo
  
  do ip=1,nip
     nb= nprox(ip)
     s=0.d0
     t=0.d0
     !write(6, 121) 'ip=', ip
     do is= 1, nb          ! loop nn bond
        jp=prox(is,ip)
        iE=E4C(is,ip)
        !s= s+0.5d0*(FV(ip,1)+FV(jp,1))*FV(iE,2)*slen(is,ip)*NoutdUn(is,ip)
        s= s + Fe(iE) * slen(is,ip) * NoutdUn(is,ip)        ! NOTE: could use (iE) only 
        t= t + FV(iE,2) * slen(is,ip) * NoutdUn(is,ip)        ! NOTE: could use (iE) only
        !
     enddo
     ASV(ip,4)= t / areaI(ip)     ! div_V     center
     ASV(ip,5)= 0.d0              ! vort_V    edge:    MPAS avoid this treatment
     ASV(ip,6)= s / areaI(ip)     ! div_(hV)  center
  enddo



  !------------------------------------------------------------
  !  step.2 :
  !  prepare  vorticity = k \dot  curl V :  xiE(iE) by xiU, xiD
  !  output:  qe_diag(iE), h_e(iE), diag= diagnostic
  !
  !  algorithm for PV
  !    vertex --> edge
  !           --> center
  !               grad_n  grad_t 
  !               advect edge values
  !
  !               Add upwiding scheme for PV ,    
  !               d PV /dt =0   -->   \p PV \p t +  V \dot gradPV = 0
  !               output:  qe_diag(iE) correction
  !
  !               Note:  V \dot \grad PV
  !                    =  Un grad_n(PV) * \vec Un \dot \grad_n  +
  !                       Ut grad_t
  !------------------------------------------------------------
  do iE=ips,ipe
     !
     !------  M1:  Q_e^ perp  = Fe^ perp * qe = 1/de \sum_e'  Wee' le' he Un^e' * 1/2 ( q_e + q_e' )
     !                    qe  = 1/2 * ( q_v1 + q_v2),    q_v = ( f + xi ) / h_v
     !        A= qe;   dA/dt=0;  \p qe / \p t = - V \dot \grad qe
     !
     ip=C4E(1,iE); jp=C4E(2,iE); kp=C4E(3,iE); Lp=C4E(4,iE); is=C4E(5,iE)
!!     write(6, 121) 'ip,jp,kp,Lp,is', ip,jp,kp,Lp,is
     h_V_up =  &
          RAiv4T(1,iE)*FV(ip,1) + &
          RAiv4T(2,iE)*FV(jp,1) + &
          RAiv4T(3,iE)*FV(kp,1)
     h_V_dn =  &
          RAiv4T(4,iE)*FV(ip,1) + &
          RAiv4T(5,iE)*FV(Lp,1) + &
          RAiv4T(6,iE)*FV(jp,1)
!
!     write(6, 101) 'RAiv4T(1:6,iE)', RAiv4T(1:6,iE)
!     write(6, 101) 'h_V_up', h_V_up
!     write(6, 101) 'h_V_dn', h_V_dn
!     write(6, *)
     !
     nb= nprox(ip)
     ism= mod(is+nb-2, nb) + 1
     isp= mod(is, nb)+1
     ks=  proxs(isp,ip)
     js=  proxs(is,ip)
     nb=  nprox(jp)
     jsm= mod(js+nb-2, nb) + 1
     !
     vort=0.d0
     do is2=1,3
        iE2= EindexVort(is2,iE)
        vort= vort + FV(iE2,2)*LevdUn(is2,iE)*dcenter(iE2)          ! LevdUn with fixed iE
!!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortU= vort / areaV(1,iE) + fcori_V(1,iE)   ! total rotation quantum for Up Triangle
     ! add the absolute vorticity
     ASV(iE,5) = vort / areaV(1,iE)    ! Up
     !
     vort=0.d0
     do is2=4,6
        iE2= EindexVort(is2,iE)
        vort= vort + FV(iE2,2)*LevdUn(is2,iE)*dcenter(iE2)
!!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortD= vort / areaV(2,iE) + fcori_V(2,iE)   ! total rotation quantum for Down Triangle
     !
     ! add the absolute vorticity
     ASV(iE,5) = ASV(iE,5) + vort / areaV(2,iE)    ! Down
     ASV(iE,5) = ASV(iE,5) / 2.d0                  ! final  (xi_U + xi_D) / 2
     !
!     write(6,*)
!     write(6, 121) 'iE', iE
!     write(6, 101) 'vortU, vortD, vort_mid', vortU, vortD, 0.5d0*(vortU+vortD)
!     write(6, 101) 'vortU_D middle analy  ', 2.d0*(omega + omegasw) * sin( lat_v(iE))
!     if (iE.EQ.30)  STOP 'test 2'
!


!
!     if ( h_V_up .lt. 1.d-16 .OR. h_V_dn .lt. 1.d-16 ) then
!        write(6, 101)  'FV(ip,1),FV(jp,1),FV(kp,1)', FV(ip,1),FV(jp,1),FV(kp,1)
!        STOP 'h field == 0 in C-grid; MPAS approach fails'
!     endif
!
     qeU = VortU / h_V_up      !!!  WHAT IF  h == 0 ????, STOP 'failed '
     qeD = VortD / h_V_dn      !!!  WHAT IF  h == 0 ????, STOP 'failed '
     !
     ! -- for each fixed iE, (is,ip)-- up header, 
     !    store point to three (is,ip), (jsm,jp), (ks,kp)
     !
     qe_hex_vertex(1,iE)= qeU
     qe_hex_vertex(2,iE)= qeD
     qe_hex_edge(iE)= 0.5d0 * (qeU + qeD)
     !
     !  this is diagnostic q_e  =  PV
     !
!     write(6, 201) 'qe_hex_edge(iE)       ', qe_hex_edge(iE)
!     write(6, 201) 'qe_hex_edge(iE) analy ', 2.d0*(omega + omegasw)*sin(lat_v(iE)) / h_e(iE)
     !
  enddo   ! loop iE


  if ( iconsv_scheme .EQ. 1 .OR. iconsv_scheme .EQ. 2 ) then
     do iE=ips,ipe
        qe_diag(iE)= qe_hex_edge(iE)
     enddo
  else
     STOP 'not implemented here'
  endif


  !------------------------------------------------------------
  !  step.3 : 
  !  output:  (1) grad (K + g(h+b))
  !           (2) Q_e^ perp  symmetrized tendency
  !------------------------------------------------------------

  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)
     nb= nprox(ip)
     ism=mod(is+nb-2, nb)+1
     !---------------------------------------------
     ! U_n \dot grad E  =  - (E_L - E_R)/dis * (n^out \dot  U_n)
     !       caution:      direction   n^out \dot  U_n
     !       K          =  K_i from n-component only  L, R
     !       E_ip (L,R) =  wt *  En(iE)^2  + gh
     !           areaE  =  1/4 * l_e * d_e
     !---------------------------------------------  
     !
     !-- Left
     !
     E=0.d0                          !  Kinetic E.
     do is2=1, nb
        iE2= E4C(is2,ip)
        E= E + FV(iE2,2) * FV(iE2,2) * areaE(iE2)
     enddo
     E1= E/areaI(ip)                 ! test grad(K)
     E3= g*(FV(ip,1) + hb_a(ip))
     EL= E/areaI(ip) + g*(FV(ip,1) + hb_a(ip))  !  E = K_i + gh_i    MPAS def.  areaE/area
     !-- Right
     nb= nprox(jp)
     E=0.d0
     do js2=1, nb
        iE2= E4C(js2,jp)
        E= E + FV(iE2,2) * FV(iE2,2) * areaE(iE2)
     enddo
     E2= E/areaI(jp)                ! test grad(K)
     E4= g*(FV(jp,1) + hb_a(jp))
     ER= E/areaI(jp) + g*(FV(jp,1) + hb_a(jp))  !  E = K_i + gh_i  MPAS def. areaE/area
     gradEn = - (EL - ER) / dcenter(iE) * NoutdUn(is,ip)
     gradKn = - (E1 - E2) / dcenter(iE) * NoutdUn(is,ip)
     gradgh = - (E3 - E4) / dcenter(iE) * NoutdUn(is,ip)


!     !
!     ! test gradient for TC2
!     !
!     if ( iE .LE. 15 ) then
!        s= omegasw**2.d0*ae*sin(lat_v(iE))*Nvec(3,is,ip)    ! (\omega dot \r) \omega \dot N
!        write(6, 101) 'gradKn num', gradKn
!        write(6, 101) 'gradKn ana', - s
!        write(6,*)
!        write(6, 101) 'gradgh num', gradgh
!        write(6, 101) 'gradgh ana', - (1.d0+ 2.d0*omega/omegasw) * s
!        write(6,*)
!     else
!        STOP 'ck gradKn: num vs ana'
!     endif



     !---------------------------------------------
     !   symmetric  (xi + f)_e * Ut^e
     !
     !    KE conservative way
     !    Q_e ^ perp = eta * U_t^e
     !        = 1/de \sum_e' W_ee' * l_e' * he' * U_n^e' * 1/2 * [ qe + qe' ]
     !        =                            ^^^^^^^^^^^^^
     !                                         Fe'
     !    q_e (iE) is prepared in above section
     !---------------------------------------------
     !
     nvert = nprox(ip) + nprox(jp) - 2        ! double Benzene
     if (iconsv_scheme .EQ. 1 .OR. iconsv_scheme .EQ. 3) then
        !
        !-- K_i conserve, qe diagnostic 
        !
        s=0.d0
        t=0.d0
        do is2=1,nvert                           ! pentagon + hexagon =9; otherwise 10
           iE2= Eindex(is2,iE)                   ! 
           ipc= C4E(1,iE2); isc= C4E(5,iE2)        ! ip,is on circle
!!           s= s + wtUt(is2,iE)*slen(isc,ipc)*FV(iE2,2)*h_e(iE2)*0.5d0*(qe_diag(iE) + qe_diag(iE2)) 
           s= s + wtUt(is2,iE) * slen(isc,ipc) * Fe(iE2) * 0.5d0 * ( qe_diag(iE) + qe_diag(iE2) )
           t= t + wtUt(is2,iE) * slen(isc,ipc) * FV(iE2,2)
           !
           !
!           if ( iE .LE. 20 ) then
!              write(6, 121) 'iE,is2, iE2,ipc,isc', iE,is2, iE2,ipc,isc
!              write(6, 101) 'wtUt(is2,iE) ', wtUt(is2,iE)
!              write(6, 101) 'slen(isc,ipc)', slen(isc,ipc)
!              write(6, 101) 'Fe(iE2), h_e(iE), Un(iE)', Fe(iE2), h_e(iE2), FV(iE2,2)
!              write(6, 201) 'qe_diag(iE)           ', qe_diag(iE)
!              write(6, 201) 'qe_diag(iE2)          ', qe_diag(iE2)
!              write(6, 201) 'qe_hex_edge(iE) analy ', 2.d0*(omega + omegasw)*sin(lat_v(iE)) / h_e(iE)
!              write(6,*)
!              !
!           else
!              STOP 'stop check Ut'
!           endif
           !
           !
        enddo
        etaUt= s / dcenter(iE)
        Ut=    t / dcenter(iE)
     endif
     !
     ! \P Un / \pt - (f+vort)*U_t  + U_n^e \dot grad(E)= 0
     !             as \U_n \dot (k \cross U_n^e) = -1
     !
     !!!AFV(iE,2)=  etaUt - gradEn
     !!!AFV(iE,3)=  0.d0   ! in case defined
     !
     ! -- previous dated  12/28/2017
     ! AFV(iE,2)= -(fcori(iE)+vort)*Ut-gradEn          !  coriolios force is not included
     !
     ASV(iE,7) = -etaUt      ! U_n^e  direction
     ASV(iE,8) = FV(iE,2)    ! U_n
     ASV(iE,9) = Ut          ! U_t
     ASV(iE,10)= etaUt - gradEn   ! total tendency
     ASV(iE,11)= 0.d0
  enddo       ! iE

!
!  write(iunit, 211) dt*real(itstep)/86400.d0, global_sum(1:5), &
!       (global_sum(2)-TE0)/TE0, (global_sum(3)-PE0)/PE0
!
  deallocate ( qe_diag, h_e, Fe )
  deallocate ( qe_hex_center, qe_hex_vertex, qe_hex_edge )
  deallocate ( vort_hex_vertex, vort_hex_center )
  deallocate ( Un_dot_gradPV, Ut_dot_gradPV )
  return

  include '../swm/myformat.inc'
end subroutine test_dyn_C_3func






!
!--  YGYU:  29-Oct-2018
!           this code shows NR77 sub for spherical harmonics fail at N=16
!           enough test
!
subroutine test_spherical_harmonics (nx, ny, MM)
  use module_constants, only  :  pi
  implicit none
  integer, intent (in) :: nx, ny, MM
  complex :: HSC(0:MM,0:MM)
  COMPLEX, external :: Y_lm_sph_harm     !  External function
  !
  integer :: ix, iy, L
  real    :: theta(nx), phi(ny), grid(nx*ny,2), dx, dy
  real    :: lat, lon, H_r
  integer :: ipt, NPTS
  integer :: NN, m, n
  complex :: Y_lm, Y_lmm, sum


  if ( nx.LT.1 .OR. ny.LT.1 ) STOP 'nx,ny < 1 test_spherical_harmonics fails'
  dx = 2.d0 * pi / (nx -1)
  dy = 1.d0 * pi / (ny -1)
  !
  ipt = 0
  do ix = 1, nx
     do iy = 1, ny
        ipt = ipt + 1
        theta (ix) = dx * ( ix - 1 )
        phi   (iy) = dy * ( ix - 1 )
!!        HSC(ix, iy) = ( cos(theta(ix)) ,  sin(phi(iy)) )
        grid(ipt, 1)=  phi(iy)           ! lon
        grid(ipt, 2)=  theta(ix)         ! lat
     enddo
  enddo
  NPTS = ipt

  write(6,*) 'pt -1'

  !
  !  Guo Shuo Hong: ED:  page 95:  1/r =  \sum_n 1/2**n Pn(cos theta)
  !
  NN = MM
  HSC(0:NN, 0:MM) = CMPLX(0.d0, 0.d0)
  do n=0, NN
     m=0
     HSC(n, m) = CMPLX ( 2.d0**(-n), 0.d0 )
  enddo
  !
  !  HSC(n, -m) = (-1)**m * CONJG ( HSC(n, m) )
  !
  write(6,*) 'pt 2'


  !
  !--  Spec trans of real function
  !
  do ipt=1, NPTS
     lon = grid(ipt, 1)
     lat = grid(ipt, 2)
     sum = CMPLX(0.d0, 0.d0)
     
     write(6, 101) 'lon,lat,sum', lon,lat,sum

     do n=0, MM
        m = 0            ! order m = zero
        Y_lm = Y_lm_sph_harm (n,m,lat,lon)
        write(6, 101) 'n, m, Y_lm', 1.0*n, 1.0*m, Y_lm
        call flush(6)

        sum = sum + HSC(n,m) * Y_lm
        do m= 1, n       ! order m > 0
           Y_lm  = Y_lm_sph_harm (n,m,lat,lon)
           write(6, 101) 'Y_lm', Y_lm; call flush(6)

           Y_lmm = CMPLX( (-1.d0)**dble(m), 0.d0) * CONJG(Y_lm)     ! NR77 book
           write(6, 101) 'Y_lmm', Y_lmm; call flush(6)

           sum   = sum + HSC(n,m) * Y_lm
           write(6, 101) 'sum', sum ; call flush(6)

        enddo
     enddo

     H_r = REAL (sum)

     write(6, 101) 'lon, lat, Fvalue', lon, lat, H_r
     call flush(6) 


  enddo
  !
  STOP 't2'
  !
  !
  include '../swm/myformat.inc'
end subroutine test_spherical_harmonics



subroutine ck_against_stswm
  use module_constants, only : pi
  use module_control,   only : nip,niE,NPTS,ndimFV,filenamelen
  use module_variables, only : lon_deg_plot,lat_deg_plot,field_plot
  use readnamelist,     only : glvl
  use readnamelist,     only : nx_ll,ny_ll
  implicit none
  !
  real    :: d2r, temp(3)
  integer :: iunit, nx, ny, i, j, k
  character (len=filenamelen) :: filename, filename_stsw, filename_icos

  d2r=1.d0/180.d0*pi
  iunit=11
  open(unit=iunit,file='stswm_reg_latlon.dat',status='unknown',&
       form="unformatted", action='read')
  read(iunit) nx, ny
  
  k=0
  do j=1,ny
     do i=1,nx
        k=k+1
        read(iunit) lon_deg_plot(i), lat_deg_plot(j), field_plot(i,j), temp(2:3)
        !
        !!   from stswm
        !!write(11) RLON, RLAT, h_A2, u_A2, v_A2
        !
        lon_deg_plot(i)= lon_deg_plot(i)/d2r
        lat_deg_plot(j)= lat_deg_plot(j)/d2r
        !
        if ( mod( k, nx*ny/4) .EQ. 1) then
           write(6, 103) lon_deg_plot(i), lat_deg_plot(j), field_plot(i,j), temp(2:3)
        endif
        !
     enddo
  enddo

  filename='h_field_stswm.nc'
  call lon_lat_press_wt_nc4 (2,nx,ny,lon_deg_plot,lat_deg_plot,field_plot,filename)
  write(6,*) 'nail 6'
  !
  include '../swm/myformat.inc'  
end subroutine ck_against_stswm






!----------------------------------------
!
!   E r r o r      N o r m
!
!----------------------------------------



!
!               kp 
!               .      ______
!               .  .  /   
!               .    / .P2   
!             P3___V/    . jp
!            /  .   \ . 
!           /   .  . \ *    P1 (dcenter_A, Nvec2)
!    ______/    .     \_*_____
!          \   ip        *
!           \             * Nvec2
!            \
!
!   meta code:
!    for (is,ip)
!      find gradFX on V
!    find div( gradFX ) via FV, using gradFV(mid, edge)
!
subroutine cal_Laplace_scalar_A (npt, FX, diff_term)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, proxs, nip
  use module_variables, only  :  slen, area, Nvec, Nvec2, areaV_A, dcenter_A
  use module_control,   only  :  hypervis
  implicit none
  integer, intent(in ):: npt
  real(rt),intent(in ):: FX(npt)    ! scalar
  real(rt),intent(out):: diff_term(npt)  ! scalar
  !
  integer  :: ip, jp, kp, nb, nb2, j
  integer  :: is, js, ks, isp, ism, jsm     ! is plus 1 (is+1)
  real(rt) :: E1, E2, E3
  real(rt) :: s1, s2, s3
  real(rt) :: N1(3), N2(3), N3(3)
  real(rt) :: s,  t
  real(rt) :: gradFX_V(3,6), Lap        !  Laplace term
  
  IF (npt .NE. nip) STOP 'npt .NE. nip in subroutine cal_grad_A'

  !
  ! s1. find grad(FX) on vertex_U
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        js=proxs(is,ip)
        nb2=nprox(jp)
        jsm=mod(js+nb2-2, nb2) + 1
        !
        E1= 0.5d0*(FX(ip) + FX(jp))  ! i-j
        E2= 0.5d0*(FX(jp) + FX(kp))  ! j-k
        E3= 0.5d0*(FX(kp) + FX(ip))  ! k-i
        !
        s1= dcenter_A(is, ip)
        s2= dcenter_A(jsm,jp)
        s3= dcenter_A(isp,ip)
        ! N_out_vec for upper Triangle
        N1(1:3)= Nvec2(1:3,is, ip)    ! watch out this is the true Nvec2 on A-grid
        N2(1:3)= Nvec2(1:3,jsm,jp)    ! note: length, Nvec2, areaV_A all exactly calc.
        N3(1:3)=-Nvec2(1:3,isp,ip)    ! 
        !
        do j=1,3  ! 3-component
           s= E1*N1(j)*s1 + E2*N2(j)*s2 + E3*N3(j)*s3   ! FX * Nvec * ds
           gradFX_V(j,is)= s/areaV_A(is,ip)
        enddo
        !--
        !!write(6,121) 'is,ip,jp,kp         ', is,ip,jp,kp
        !!write(6,101) 'FX(ip),FX(jp),FX(kp)', FX(ip),FX(jp),FX(kp)
        !!write(6,101) 'gradFX_V(1:3,is)    ', gradFX_V(1:3,is)
     enddo    ! loop is 


     !
     ! s2. find div ( grad(FX) ) from mid-point
     !
     s=0.d0      ! div( grad(FX) )
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        !
        ! V \dot N
        t=0.d0
        do j=1,3
           t= t + 0.5d0 * (gradFX_V(j,is) + gradFX_V(j,ism)) * Nvec(j,is,ip)
        enddo
        s= s + t*slen(is,ip)
     enddo
     Lap = s/area(ip)
     diff_term(ip)= s/area(ip)*hypervis     ! 2nd order \viscosity * Lap 
     !
     !! write(6,241) 'ip, Lap _X1', ip, Lap 
  enddo     !  ip

  return
  include '../swm/myformat.inc'
end subroutine cal_Laplace_scalar_A





!
!               kp 
!               .      ______
!               .  .  /   
!               .    / .P2   
!             P3___V/    . jp
!            /  .   \ . 
!           /   .  . \ *    P1 (dcenter_A, Nvec2)
!    ______/    .     \_*_____
!          \   ip        *
!           \             * Nvec2
!            \
!
!   meta code:
!    for (is,ip)
!      find gradFX on V
!    find div( gradFX ) via FV, using gradFV(mid, edge)
!

subroutine cal_Laplace_square_A (npt, FX, diff_term)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, proxs, nip
  use module_variables, only  :  wtsph, slen, area, Nvec, Nvec2, areaV_A, dcenter_A
  use module_control,   only  :  hypervis
  implicit none
  integer, intent(in ):: npt
  real(rt),intent(in ):: FX(npt)    ! scalar
  real(rt),intent(out):: diff_term(npt)  ! scalar    Laplace^2 
  !
  integer  :: ip, jp, kp, nb, nb2, j
  integer  :: is, js, ks, isp, ism, jsm     ! is plus 1 (is+1)
  real(rt) :: E1, E2, E3, EV
  real(rt) :: s1, s2, s3
  real(rt) :: N1(3), N2(3), N3(3)
  real(rt) :: s,  t
  real(rt) :: gradFX_V(3,6), gradFY_V(3,6), Lap 
  real(rt), allocatable :: FY(:)  ! scalar

  
  IF (npt .NE. nip) STOP 'npt .NE. nip in subroutine cal_grad_A'
  allocate (FY(npt))
  !
  ! s1. find grad(FX) on vertex_U
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        js=proxs(is,ip)
        nb2=nprox(jp)
        jsm=mod(js+nb2-2, nb2) + 1
        !
        do j=1,3  ! 3-component
           gradFX_V(j,is)= &
               ((FX(ip) + FX(jp)) * dcenter_A(is, ip) * Nvec2(j,is, ip) + &
                (FX(jp) + FX(kp)) * dcenter_A(jsm,jp) * Nvec2(j,jsm,jp) - &
                (FX(kp) + FX(ip)) * dcenter_A(isp,ip) * Nvec2(j,isp,ip)   )/ &
                areaV_A(is,ip) * 0.5d0
        enddo                                         ! no need to eliminate Radial 
        !                                             ! bc it is orthogonal to Nvec
     enddo    ! loop is 

     !
     ! s2. find div ( grad(FX) ) from mid-point
     !
     s=0.d0      ! div( grad(FX) )
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        !
        ! V \dot N
        s=   s +  0.5d0 * ( &
             (gradFX_V(1,is) + gradFX_V(1,ism)) * Nvec(1,is,ip) + &
             (gradFX_V(2,is) + gradFX_V(2,ism)) * Nvec(2,is,ip) + &
             (gradFX_V(3,is) + gradFX_V(3,ism)) * Nvec(3,is,ip) ) * slen(is,ip)
     enddo
     Lap   = s/area(ip)
     FY(ip)= s/area(ip)      !  extra copy
!!     write(6,241) 'ip, Lap _X2', ip, Lap 
  enddo     !  ip


  !
  !-- repeat the grad(FX), div( grad(FX) ) again; rename as FY
  !
  ! s3. find grad(FY) on vertex_U
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        js=proxs(is,ip)
        nb2=nprox(jp)
        jsm=mod(js+nb2-2, nb2) + 1
        !
        do j=1,3  ! 3-component
           gradFY_V(j,is)= &
               ((FY(ip) + FY(jp)) * dcenter_A(is, ip) * Nvec2(j,is, ip) + &
                (FY(jp) + FY(kp)) * dcenter_A(jsm,jp) * Nvec2(j,jsm,jp) - &
                (FY(kp) + FY(ip)) * dcenter_A(isp,ip) * Nvec2(j,isp,ip)   )/ &
                areaV_A(is,ip) * 0.5d0
        enddo                                         ! no need to eliminate Radial 
        !                                             ! bc it is orthogonal to Nvec
     enddo    ! loop is 


     !
     ! s4. find div ( grad(FY) ) from mid-point
     !
     s=0.d0      ! div( grad(FY) )
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        !
        ! V \dot N
        s=   s +  0.5d0 * ( &
             (gradFY_V(1,is) + gradFY_V(1,ism)) * Nvec(1,is,ip) + &
             (gradFY_V(2,is) + gradFY_V(2,ism)) * Nvec(2,is,ip) + &
             (gradFY_V(3,is) + gradFY_V(3,ism)) * Nvec(3,is,ip) ) * slen(is,ip)
     enddo
     Lap    = s/area(ip)                     ! 4th order  = - \viscosity * Lap^2 
     diff_term(ip)= -s/area(ip)*hypervis     ! 4th order  = - \viscosity * Lap^2 
     !
!!      write(6,241) 'ip, Lap^2, diff_term ', ip, Lap, diff_term(ip)
     !
  enddo     !  ip

!!  stop 'ck subroutine cal_Laplace_square_A'

  deallocate (FY)
  return
  include '../swm/myformat.inc'
end subroutine cal_Laplace_square_A




!!  old version save a copy
!!
!!subroutine cal_Laplace_square_A_sv (npt, FX, diff_term)
!!  use kinds, only: rt
!!  use module_variables, only  :  nprox, prox, proxs, nip
!!  use module_variables, only  :  wtsph, slen, area, Nvec, Nvec2, areaV_A, dcenter_A
!!  use module_control,   only  :  hypervis
!!  implicit none
!!  integer, intent(in ):: npt
!!  real(rt),intent(in ):: FX(npt)    ! scalar
!!  real(rt),intent(out):: diff_term(npt)  ! scalar    Laplace^2 
!!  !
!!  integer  :: ip, jp, kp, nb, nb2, j
!!  integer  :: is, js, ks, isp, ism, jsm     ! is plus 1 (is+1)
!!  real(rt) :: E1, E2, E3, EV
!!  real(rt) :: s1, s2, s3
!!  real(rt) :: N1(3), N2(3), N3(3)
!!  real(rt) :: s,  t
!!  real(rt) :: gradFX_V(3,6), gradFY_V(3,6), Lap 
!!  real(rt), allocatable :: FY(:)  ! scalar
!!
!!  
!!  IF (npt .NE. nip) STOP 'npt .NE. nip in subroutine cal_grad_A'
!!  allocate (FY(npt))
!!  !
!!  ! s1. find grad(FX) on vertex_U
!!  !
!!  do ip=1, nip
!!     nb= nprox(ip)
!!     !
!!     do is= 1, nb
!!        isp=mod(is, nb)+1
!!        jp=prox(is,ip)
!!        kp=prox(isp,ip)
!!        js=proxs(is,ip)
!!        nb2=nprox(jp)
!!        jsm=mod(js+nb2-2, nb2) + 1
!!        !
!!        E1= 0.5d0*(FX(ip) + FX(jp))  ! i-j
!!        E2= 0.5d0*(FX(jp) + FX(kp))  ! j-k
!!        E3= 0.5d0*(FX(kp) + FX(ip))  ! k-i
!!        EV= FX(ip)*wtsph(1,is,ip)+FX(jp)*wtsph(2,is,ip)+FX(kp)*wtsph(3,is,ip)  ! on vertex_U
!!        !
!!        s1= dcenter_A(is, ip)
!!        s2= dcenter_A(jsm,jp)
!!        s3= dcenter_A(isp,ip)
!!        ! N_out_vec for upper Triangle
!!        N1(1:3)= Nvec2(1:3,is, ip)    ! watch out this is the true Nvec2 on A-grid
!!        N2(1:3)= Nvec2(1:3,jsm,jp)    ! note: length, Nvec2, areaV_A all exactly calc.
!!        N3(1:3)=-Nvec2(1:3,isp,ip)    ! 
!!        !
!!        do j=1,3  ! 3-component
!!           s= E1*N1(j)*s1 + E2*N2(j)*s2 + E3*N3(j)*s3 ! FX * Nvec * ds
!!           t= N1(j)*s1 + N2(j)*s2 + N3(j)*s3          ! Nvec * ds
!!!!           gradFX_V(j,is)= s/areaV_A(is,ip)         ! Tomita
!!           gradFX_V(j,is)= (s - EV*t)/areaV_A(is,ip)  ! remove constant field
!!        enddo                                         ! no need to eliminate Radial 
!!        !                                             ! bc it is orthogonal to Nvec
!!        !
!!!!        write(6,121) 'is,ip,jp,kp         ', is,ip,jp,kp
!!!!        write(6,101) 'FX(ip),FX(jp),FX(kp)', FX(ip),FX(jp),FX(kp)
!!!!        write(6,101) 'gradFX_V(1:3,is)    ', gradFX_V(1:3,is)
!!     enddo    ! loop is 
!!
!!
!!     !
!!     ! s2. find div ( grad(FX) ) from mid-point
!!     !
!!     s=0.d0      ! div( grad(FX) )
!!     do is= 1, nb
!!        ism=mod(is+nb-2, nb)+1
!!        !
!!        ! V \dot N
!!        t=0.d0
!!        do j=1,3   !  ! shit!!!!!!!
!!           t= 0.5d0 * (gradFX_V(j,is) + gradFX_V(j,ism)) * Nvec(j,is,ip)
!!        enddo
!!        s= s + t*slen(is,ip)
!!     enddo
!!     Lap   = s/area(ip)
!!     FY(ip)= s/area(ip)      !  extra copy
!!  enddo     !  ip
!!
!!
!!  !
!!  !-- repeat the grad(FX), div( grad(FX) ) again; rename as FY
!!  !
!!!  do ip=1,nip
!!!     FY(ip)= Lap(ip)
!!!  enddo
!!  !
!!  ! s3. find grad(FY) on vertex_U
!!  !
!!  do ip=1, nip
!!     nb= nprox(ip)
!!     !
!!     do is= 1, nb
!!        isp=mod(is, nb)+1
!!        jp=prox(is,ip)
!!        kp=prox(isp,ip)
!!        js=proxs(is,ip)
!!        nb2=nprox(jp)
!!        jsm=mod(js+nb2-2, nb2) + 1
!!        !
!!        E1= 0.5d0*(FY(ip) + FY(jp))  ! i-j
!!        E2= 0.5d0*(FY(jp) + FY(kp))  ! j-k
!!        E3= 0.5d0*(FY(kp) + FY(ip))  ! k-i
!!        !
!!        s1= dcenter_A(is, ip)
!!        s2= dcenter_A(jsm,jp)
!!        s3= dcenter_A(isp,ip)
!!        ! N_out_vec for upper Triangle
!!        N1(1:3)= Nvec2(1:3,is, ip)    ! watch out this is the true Nvec2 on A-grid
!!        N2(1:3)= Nvec2(1:3,jsm,jp)    ! note: length, Nvec2, areaV_A all exactly calc.
!!        N3(1:3)=-Nvec2(1:3,isp,ip)    ! 
!!        !
!!        do j=1,3  ! 3-component
!!           s= E1*N1(j)*s1 + E2*N2(j)*s2 + E3*N3(j)*s3   ! FY * Nvec * ds
!!           gradFY_V(j,is)= s/areaV_A(is,ip)
!!        enddo
!!     enddo    ! loop is
!!     
!!     !
!!     ! s4. find div ( grad(FY) ) from mid-point
!!     !
!!     s=0.d0      ! div( grad(FY) )
!!     do is= 1, nb
!!        ism=mod(is+nb-2, nb)+1
!!        !
!!        ! V \dot N
!!        t=0.d0
!!        do j=1,3
!!           t= 0.5d0 * (gradFY_V(j,is) + gradFY_V(j,ism)) * Nvec(j,is,ip)
!!        enddo
!!        s= s + t*slen(is,ip)
!!     enddo
!!     Lap    = s/area(ip)                     ! 4th order  = - \viscosity * Lap^2 
!!     diff_term(ip)= -s/area(ip)*hypervis     ! 4th order  = - \viscosity * Lap^2 
!!  enddo     !  ip
!!  !
!!  !
!!
!!  deallocate (FY)
!!  return
!!  include '../swm/myformat.inc'
!!end subroutine cal_Laplace_square_A_sv
!!







!
!  aux_v =  grad ( aux_s )
!
subroutine cal_grad_A (npt, aux_s, aux_v)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, proxs, nip
  use module_variables, only  :  slen, area, Nvec, wtsph
  implicit none
  integer, intent(in ):: npt
  real(rt),intent(in ):: aux_s(npt)    ! scalar
  real(rt),intent(out):: aux_v(3,npt)  ! vector
  !
  integer  :: ip, jp, kp, nb, j
  integer  :: is, js, ks, isp, ism     ! is plus 1 (is+1)
  real(rt) :: Ec(6), Es(6), s, E0
  !
  !
  IF (npt .NE. nip) STOP 'npt .NE. nip in subroutine cal_grad_A'
  !
  do ip=1, nip
     nb= nprox(ip)
     E0= aux_s(ip)
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        ! E at corner
        Ec(is)= aux_s(ip)*wtsph(1,is,ip) + &
             aux_s(jp)*wtsph(2,is,ip)+ aux_s(kp)*wtsph(3,is,ip)
     enddo
     !
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        Es(is)= 0.5d0*(Ec(ism)+Ec(is))    ! KE on six sides
     enddo
     !
     do j=1, 3
        s=0.d0
        do is= 1, nb
           ism=mod(is+nb-2, nb)+1
           s = s + (Es(is) - E0)*Nvec(j,is,ip)*slen(is,ip) 
        enddo
        aux_v(j,ip)= s/area(ip)       ! this is grad(h^*)
     enddo
  enddo
  return
end subroutine cal_grad_A

!
!  aux_s =  div ( aux_v )
!
subroutine cal_diverg_A (npt, aux_v, aux_s)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, proxs, nip
  use module_variables, only  :  slen, area, Nvec, wtsph
  implicit none
  integer, intent(in ):: npt
  real(rt),intent(in ):: aux_v(3,npt)  ! vector
  real(rt),intent(out):: aux_s(npt)    ! scalar
  !
  integer  :: ip, jp, kp, nb, j
  integer  :: is, js, ks, isp, ism     ! is plus 1 (is+1)
  real(rt) :: Fc(3,6), Fs(3,6), s, t               ! Fc : vector at corner
  !
  IF (npt .NE. nip) STOP 'npt .NE. nip in subroutine cal_grad_A'
  !
  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        ! vector at corner
        do j=1, 3
           Fc(j,is)= aux_v(j,ip)*wtsph(1,is,ip) + &
             aux_v(j,jp)*wtsph(2,is,ip)+ aux_v(j,kp)*wtsph(3,is,ip)
        enddo
     enddo
     !
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        do j=1, 3
           Fs(j,is)= 0.5d0*(Fc(j,ism) + Fc(j,is))      ! h^*,Vx,Vy,Vz
        enddo
     enddo
     !
     s=0.d0
     do is= 1, nb
        t=0.d0
        do j=1, 3
           t = t + Fs(j,is)*Nvec(j,is,ip)    !  V \dot N
        enddo
        s= s + t*slen(is,ip)
     enddo
     aux_s(ip)= s/area(ip)          ! this is  div(V)
  enddo
  return
end subroutine cal_diverg_A



!
!  aux_s =  Un(1:niE)
!
!  U_n \dot \Lap V = Un \dot \grad (div V) - Un \dot n_out \part t ( vort )
!
subroutine  cal_Laplace_Un_C (npt, aux_s, diff_term)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, proxs, niE
  use module_variables, only  :  slen, areaI, Nvec
  use module_variables, only  :  areaE, areaV, dcenter, wtUt
  use module_variables, only  :  Eindex, EindexVort, C4E, E4C, LevdUn, NoutdUn
  use module_control,   only  :  hypervis
  implicit none
  integer, intent(in ):: npt
  real(rt),intent(out):: aux_s(npt)     ! scalar, Un
  real(rt),intent(out):: diff_term(npt)
  !
  real(rt) :: laplace_Un     ! scalar,  \Laplace Un
  integer  :: iE, iE2, ip, jp, kp, nb, j, nb2
  integer  :: is, js, ks, isp, ism, is2     ! is plus 1 (is+1)
  real(rt) :: Vc(3,6), s, t               ! Vc : vector at corner
  real(rt) :: Lap, div_L, div_R, vort, vortU, vortD
  !
  !
  IF (npt .NE. niE) STOP 'npt .NE. niE in subroutine  cal_Laplace_Un_C'
  !
  do iE=1, niE
     ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)  !! ; kp=C4E(3,iE); Lp=C4E(4,iE); 
     nb= nprox(ip)
     ism= mod(is+nb-2, nb) + 1
     !
     ! Left
     !
     s=0.d0
     do is2= 1, nb          ! loop nn bond
        !        write(6, 121) 'is=', is
        iE2=E4C(is2,ip)
        s= s+aux_s(iE2)*slen(is2,ip)*NoutdUn(is2,ip)
     enddo
     div_L= s / areaI(ip)
     !
     ! Right
     !
     nb2= nprox(jp)
     s=0.d0
     do js= 1, nb2        ! loop nn bond
        !        write(6, 121) 'is=', is
        iE2=E4C(js,jp)
        s= s+aux_s(iE2)*slen(js,jp)*NoutdUn(js,jp)
     enddo
     div_R= s / areaI(jp)
     !
     !
     vort=0.d0
     do is2=1,3
        iE2= EindexVort(is2,iE)
        vort= vort + aux_s(iE2)*LevdUn(is2,iE)*dcenter(iE2)          ! LevdUn with fixed iE
!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortU= vort / areaV(1,iE)     ! total rotation quantum for Up Triangle
     !
     !
     vort=0.d0
     do is2=4,6
        iE2= EindexVort(is2,iE)
        vort= vort + aux_s(iE2)*LevdUn(is2,iE)*dcenter(iE2)
!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortD= vort / areaV(2,iE)    ! total rotation quantum for Down Triangle
     !
     !
!     write(6, 121) 'iE,ip,is=', iE,ip,is
!
     Lap = (div_R - div_L)/dcenter(iE)*NoutdUn(is,ip) -  &
           (vortU - vortD)/slen(is,ip)*NoutdUn(is,ip)
     laplace_Un = Lap
     diff_term(iE)= &
          ((div_R - div_L)/dcenter(iE)*NoutdUn(is,ip) -  &
           (vortU - vortD)/slen(is,ip)*NoutdUn(is,ip))*hypervis   
    ! 2nd ord =  \nu * Lap
!
!     write(6, 101) 'div_R, div_L', div_R, div_L
!     write(6, 101) 'vortU, vortD', vortU, vortD
!     write(6, 201) 'Lap', Lap
!
  enddo
  return
  include '../swm/myformat.inc'
end subroutine cal_Laplace_Un_C



!
!  aux_s =  h(1:niE)    ! e.g.
!
!  \Lap Psi =  1/A_I \sum_e  ( Psi_L - Psi_R ) / de * Le
!
subroutine  cal_Laplace_scalar_C (npt, aux_s, laplace_s)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, nip, niE
  use module_variables, only  :  slen, areaI, dcenter
  use module_variables, only  :  E4C
  implicit none
  integer, intent(in ):: npt
  real(rt),intent(in ):: aux_s(npt)     ! scalar, Un
  real(rt),intent(out):: laplace_s(npt) ! scalar,  \Laplace Un
  !
  integer  :: ip, jp, is, nb, iE
  real(rt) :: s
  !
  !
  IF (npt .NE. nip) STOP 'npt .NE. nip in subroutine  cal_Laplace_scalar'
  !
  do ip=1, nip
     nb= nprox(ip)
     s=0.d0
     do is= 1, nb          ! loop nn bond
        jp=prox(is,ip)
        iE=E4C(is,ip)
        s= s + (aux_s(jp)-aux_s(ip))/dcenter(iE)*slen(is,ip)
     enddo
     laplace_s(ip)= s / areaI(ip)
  enddo
  return
  include '../swm/myformat.inc'
end subroutine cal_Laplace_scalar_C



!
!  aux_s =  Un(1:niE)
!  Lap   =  \grad \dot \grad
!  Lap^2 =  \grad^4
!  U_n \dot \grad^4 V = Un \dot \grad [ Lap (div V) ] - Un \dot n_out \part t  Lap ( vort )
!
subroutine  cal_Laplace_square_Un_C (npt, aux_s, diff_term)
  use kinds, only: rt
  use module_variables, only  :  nprox, prox, proxs, nip, niE
  use module_variables, only  :  slen, areaI
  use module_variables, only  :  areaV, dcenter
  use module_variables, only  :  Eindex, EindexVort, C4E, E4C, LevdUn, NoutdUn, RAiv4T, Edc
  use module_control,   only  :  hypervis
  implicit none
  integer, intent(in ):: npt
  real(rt),intent(in ):: aux_s(npt)     ! scalar, Un
  real(rt),intent(out):: diff_term(npt) ! scalar,  \Laplace Un
  !
  real(rt) :: laplace_square_Un(npt) ! scalar,  \Laplace Un
  !
  integer  :: ip, jp, kp, Lp, nb, nb2, j, iE, iE2, iE3
  integer  :: is, js, ks, isp, ism, jsm, ksp, is2     ! is plus 1 (is+1)
  real(rt) :: Vc(3,6), s, t               ! Vc : vector at corner
  real(rt) :: Lap, div_L, div_R, vort, vortU, vortD
  real(rt) :: Lap_U(6), V1, V2, V3, V4
  real(rt), allocatable :: div(:), vort_I(:), laplace_s1(:), laplace_s2(:)
  !
  !
  IF (npt .NE. niE)  STOP 'npt .NE. niE in subroutine  cal_Laplace_Un_C'
  allocate (div(nip), vort_I(nip), laplace_s1(nip), laplace_s2(nip))
  !
  do ip=1, nip
     s=0.d0
     nb= nprox(ip)
     do is= 1, nb          ! loop nn bond
        iE2=E4C(is,ip)
        s= s+aux_s(iE2)*slen(is,ip)*NoutdUn(is,ip)
     enddo
     div(ip)= s / areaI(ip)
!!     write(6,201) 'div(ip)', div(ip)
  enddo
  !
  !  calc Lap at cell center for Div only
  !  accumulate Lap
  call cal_Laplace_scalar_C (nip, div, laplace_s1)
  do iE=1, niE
     ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)     
     div_L= laplace_s1(ip)
     div_R= laplace_s1(jp)
     Lap = (div_R - div_L)/dcenter(iE)*NoutdUn(is,ip)    !  the div part
     laplace_square_Un(iE)= Lap
!!     write(6,143) 'iE, laplace_s1(ip),laplace_square_Un(iE)', &
!!          iE, laplace_s1(ip),laplace_square_Un(iE)
  enddo     !  iE



!
!
!                        jsm2
!                 . V3____
!                 .  /.    \
!            V4 ____/V1 .jp \
!              /  iE\ . 
!             /   .  \____
!                ip  /V2
!                   /
!
! meta-code:
!  for each iE --> (is,ip)
!    for is cal vortU(V1), vortD(V2), 
!    for isp cal vortU(V4)
!    for (jsm2,jp) vortU(V3)
!    \Lap zeta =  \loop triangle (V2-V1)/le*de 
!    \hat Un \dot n_out  \part t ( Lap (vort) )
!
!    major: find Lap(zeta) at V1
!

  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb          ! loop nn bond
        iE=E4C(is,ip)
        !
        vort=0.d0
        do is2=1,3
           iE2= EindexVort(is2,iE)
           vort= vort + aux_s(iE2)*LevdUn(is2,iE)*dcenter(iE2)          ! LevdUn with fixed iE
           !!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
        enddo
        vortU= vort / areaV(1,iE)
        !
        vort=0.d0
        do is2=4,6
           iE2= EindexVort(is2,iE)
           vort= vort + aux_s(iE2)*LevdUn(is2,iE)*dcenter(iE2)
        enddo
        vortD= vort / areaV(2,iE) 
        !
        if (Edc(is,ip).EQ.0) then   ! first encounter
           V1=vortU
           V2=VortD
        else
           V1=vortD
           V2=VortU
        endif
        !
        !-- V4
        !
        isp=mod(is, nb)+1
        iE2=E4C(isp,ip)
        vort=0.d0
        do is2=1,3
           iE3= EindexVort(is2,iE2)
           vort= vort + aux_s(iE3)*LevdUn(is2,iE2)*dcenter(iE3)          ! LevdUn with fixed iE
        enddo
        vortU= vort / areaV(1,iE2)
        vort=0.d0
        do is2=4,6
           iE3= EindexVort(is2,iE2)
           vort= vort + aux_s(iE3)*LevdUn(is2,iE2)*dcenter(iE3)          ! LevdUn with fixed iE
        enddo
        vortD= vort / areaV(2,iE2)
        if (Edc(isp,ip).EQ.0) then   ! first encounter        
           V4=VortU
        else
           V4=VortD
        endif
        !
        !-- V3, check if it is (jsm,jp) or (ksp,kp)
        !
        jp= prox(is,ip)
        js= proxs(is,ip)
        nb2=nprox(jp)
        jsm=mod(js+nb2-2, nb2)+1
        iE2=E4C(jsm,jp)
        vort=0.d0
        do is2=1,3
           iE3= EindexVort(is2,iE2)
           vort= vort + aux_s(iE3)*LevdUn(is2,iE2)*dcenter(iE3)          ! LevdUn with fixed iE
        enddo
        vortU= vort / areaV(1,iE2)
        vort=0.d0
        do is2=4,6
           iE3= EindexVort(is2,iE2)
           vort= vort + aux_s(iE3)*LevdUn(is2,iE2)*dcenter(iE3)          ! LevdUn with fixed iE
        enddo
        vortD= vort / areaV(2,iE2)
        if (Edc(jsm,jp).EQ.0) then   ! first encounter        
           V3=VortD
        else
           V3=VortU
        endif
        !
        !
        !-- zeta1=V1=(is,ip), V2= V(ism1,ip), V3=(jsm,jp), V4=(isp,ip)
        !
        ! Lap_U, Lap_D for the triangle, FVM to get Lap
        !
        iE2= E4C(jsm,jp)
        iE3= E4C(isp,ip)
        s= (V2-V1)/slen(is,ip)*dcenter(iE)+ &
             (V3-V1)/slen(jsm,jp)*dcenter(iE2)+ &
             (V4-V1)/slen(isp,ip)*dcenter(iE3)
        Lap_U(is)= s/areaV(1,iE)
     enddo    ! is

     !
     !-- avoid double count for each edge,  left and right twice
     !
     do is= 1, nb          ! loop nn bond
        if (Edc(is,ip).EQ.0) then  ! first encounter
           iE=E4C(is,ip)
           ism=mod(is+nb-2, nb)+1
           vortU = Lap_U(is)
           vortD = Lap_U(ism)
           Lap = (vortU - vortD)/slen(is,ip)*NoutdUn(is,ip)
           laplace_square_Un(iE)= laplace_square_Un(iE) - Lap
           diff_term(iE)= -hypervis * laplace_square_Un(iE)   ! 4th ord  - \nu\Lap^2
           
!!           write(6,143) 'is, Lap(zeta)_Up, Lap(zeta)_D, Un dot Pt P Lap', &
!!                is, VortU, vortD, Lap
        endif
     enddo
  enddo   ! loop  ip


  deallocate (div, vort_I, laplace_s1, laplace_s2)
  return
  include '../swm/myformat.inc'
end subroutine cal_Laplace_square_Un_C




subroutine remove_radial_component_A (ims,ime,ips,ipe,FV,ndimFV)
  use kinds, only: rt
  use module_variables, only  :  Rcvec
  integer, intent(in   )::  ims, ime, ips, ipe, ndimFV
  real(rt),intent(inout)::  FV (ims:ime, ndimFV)       ! Four Vector
  !
  real(rt) :: s
  integer  :: ip, j
  !
  !-- s6. remove radial component for \part V / \part t
  !
  do ip=ips,ipe
     s=0.d0               ! s= (Vector, Rcvec) not Rvec on edge
     do j=1, 3
        s= s+ FV(ip,j+1)*Rcvec(j,ip)
     enddo
     do j=1, 3
        FV(ip,j+1)= FV(ip,j+1)-s*Rcvec(j,ip)    ! Vec = Vec - (grad,R) R
     enddo
!     if ( mod(ip,23).EQ.1 ) then
!        write(6,143) 'ip, Radial value', ip, s    ! output G4: [1.d-17, 2.d-15]
!     endif
  enddo
  include '../swm/myformat.inc'
end subroutine remove_radial_component_A


subroutine test_double_precision
  use kinds, only: rt
  use module_constants, only  :  pi, ae
  IMPLICIT NONE
  real(rt) :: Lap(4), x
  integer  :: i

  x=pi
  do i=1, 8
     Lap(i)= x
     x= x/ae
     write(6,203) i, Lap(i)
  enddo
  write(6,204) 'pi', pi
  write(6,204) '4.d0*atan(1.d0)',  4.d0*atan(1.d0)
  write(6,204) 'pi-4.d0*atan(1.d0)',  pi- 4.d0*atan(1.d0)

  stop 'ck  test_double_precision'

203  format (2x, i8, 4x, 10E100.50)
204  format (2x, a, 4x, 10E100.50)

end subroutine test_double_precision



subroutine test_mach_precision
  use kinds, only: rt
  IMPLICIT NONE
  real(rt) :: DEPS, DA
  REAL EPS, A
  INTEGER N, DN
 
      EPS = 1.0
      DEPS = 1.0D0
      N = 0
      DN = 0
 
! determination of epsilon for reals (single precision) *
 
1     EPS = EPS / 2.0
      A = EPS + 1.0
!      PRINT *, A
!      write(6, 203) A
      IF (A .GT. 1.0) THEN
         N = N + 1
         GO TO 1
      ELSE
         EPS = ABS(2.0 * EPS)
      END IF

! determination of epsilon for double precision  *
 
2     DEPS = DEPS / 2.0D0
      DA = DEPS + 1.0D0
!      PRINT *, DA
      write(6, 203) DA
      IF (DA .GT. 1.0D0) THEN
         DN = DN + 1
         GO TO 2
      ELSE
         DEPS = DABS(2.0D0 * DEPS)
      END IF
 
      IF(N .eq. DN)THEN
           PRINT*,'Your computer probably has a math coprocessor.'
           PRINT*,'Uncomment the two print statements and rerun.'
           STOP
      ENDIF

      PRINT*,'                              ',-N
      WRITE(*,10) eps
10    FORMAT(' Single precision =',e16.8,' or 2')
      PRINT*,'                              ',-DN
      WRITE(*,11) deps
11    FORMAT(' Double precision =',e16.8,' or 2')

      STOP
203  format (2x, 10E100.50)
    end subroutine test_mach_precision



