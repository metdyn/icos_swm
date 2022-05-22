
subroutine dyn_DuCVTSPR (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_constants, only  :  g, ae, pi
  use readnamelist,     only  :  iswcase, idamp, imethod_DuCVT, jmethod_DuCVT
  use module_variables, only  :  lon, lat, nprox, prox, proxs, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, wtsph, wt2sph
  use module_variables, only  :  Nvec2, hb_a, areaV_A, dcenter_A, RAiv4H_A
  use module_variables, only  :  omegasw, lambdas, Vdc, Edc, V4C, C4V, diag_var
  use module_control,   only  :  nip, niV, hypervis    ! 
  implicit none
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  integer  :: ip, jp, kp, iV
  integer  :: is, ism, isp, ispx, js, jsm, jsm2, ks, ksp
  integer  :: nb, nb2, is2, iV1, iV2
  integer  :: i, j, k
  integer  :: ipx, jpx, kpx, is2m, ia

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, E0, Kc(npp), K_s(npp)     ! Ener, E0 at ip,  E_corner
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradE(3), gradE_ana(3)      ! Energy gradient, analytic solution
  real(rt) :: vort                        !
  real(rt) :: vort_V, abvort_V            ! vorticity on Vertex 
  real(rt) :: s, t, z, s1, s2, s3, t2, r      ! sum, r=lambdas(is,ip)=ratio
  real(rt) :: div_V, div_hV               ! \div V or hV
  real(rt) :: hi, hj, hk                  !
  real(rt) :: hs, hs1, hs2, hs3, V1(3), V2(3), V3(3)
  !
  real(rt), allocatable :: diff_term(:,:), div_hV_vert(:,:)
  integer  :: imethod, jmethod
  integer  :: key(3)
  real(rt) :: hv_value(3), vort_value(3), abvort_value(3)  !  div(hV), vorticity on Vertex for Du

  !--------------  Meta code ----------------
  !
  !  imethod=1: intentionally confuse GC with circumcenter,  use 1/2 (hi + hj)
  !  imethod=2: intentionally confuse GC with circumcenter,  use  X = (x1 + lam*x2) / (1+lam) 
  !  imethod=3: linear interpolate  div(hV) from ip to (vertex)  use 1/2(V1+V2)
  !  imethod=4: linear interpolate  div(hV) from ip to (vertex)  use  wt(1:3,is,ip)
  !  imethod=5: compute div(hV), grad(div(hV)) from ip, then correct div(hV) by gradient
  !
  !----------------------------------------
  !    note
  !    FV(iV,1)=h^* ,  hb_a is used over mountain
  !    FV(ip,2)=Vx
  !    FV(ip,3)=Vy
  !    FV(ip,4)=Vz
  !    ndimFV=4
  !

  imethod= imethod_DuCVT    !  default: 3: 1/2(V1+V2);   4: Vs(is)*wt(1:3,is,ip)
  jmethod= jmethod_DuCVT    !  default: 1: 1/2(V1+V2);   2: Tomita's method

  allocate (diff_term(nip,ndimFV))
  allocate (div_hV_vert(6,nip))
  !--------------------------------------------------------------
  !  pre-calc. hypervisocosity = - eta * \nebala^4 V_x,V_y,V_z, h
  !--------------------------------------------------------------
  do j=2, ndimFV   !  no_h!,vx,vy,vz in Cartesian cord
     if (idamp .LT. 0) then
        do ip=1,nip
           diff_term(ip,j)=0.d0
        enddo
     elseif (idamp .EQ. 2) then
        call cal_Laplace_scalar_A (nip, FV(1,j), diff_term(1,j)) ! 2nd ord term =  \vis * Lap
     elseif (idamp .EQ. 4) then
        call cal_Laplace_square_A (nip, FV(1,j), diff_term(1,j)) ! 4th ord term = -\vis * Lap^2
     else
        STOP  'idamp not implemented'
     endif
  enddo


  !
  !-- part 1  h(iV) = FV(iV, 1)  
  !-- Triangle due grid
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- s1.  hs, Vs(1:3) for each triangle
     !
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq 
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
           jp=prox(is,ip);    js=proxs(is,ip)
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           nb2=  nprox(jp);   jsm= mod(js+nb2-2, nb2) + 1
           iV = V4C(is,ip)
           !
          if (imethod .EQ. 3) then
              !
              ! find div(hV) for ip, jp, kp
              !
              key(1)=ip; key(2)=jp; key(3)=kp
              !
              ! find hv_value(1:3)
              !
              z=0.d0
              do ia=1, 3                          ! loop 3-points
                 ipx=key(ia)
                 nb2=nprox(ipx)
                 s=0.d0
                 do is2=1, nb2
                    jpx=prox(is2,ipx)
                    is2m=mod(is2+nb2-2, nb2)+1    ! is2 minus
                    !
                    iV1= V4C(is2,ipx)
                    iV2= V4C(is2m,ipx)
                    hs = 0.5d0*(FV(iV1,1)+FV(iV2,1))
                    t=0.d0
                    do j=1,3
                       V1(j)= 0.5d0*(FV(ipx,j+1)+FV(jpx,j+1))    ! velocity
                       t= t + V1(j)*Nvec(j,is2,ipx)       ! V \dot N
                    enddo
                    s= s + t * hs * slen(is2,ipx)
                 enddo
                 hv_value(ia)= s / area(ipx)
              enddo
              div_hV= hv_value(1)*wtsph(1,is,ip) + &
                   hv_value(2)*wtsph(2,is,ip) + &
                   hv_value(3)*wtsph(3,is,ip)
              !
           elseif (imethod .EQ. 4) then           !  use Tomita' method
              !
              ! find div(hV) for ip, jp, kp
              !
              key(1)=ip; key(2)=jp; key(3)=kp
!              Vx(1:3,1)=FV(ip,2:4)
!              Vx(1:3,2)=FV(ip,2:4)
!              Vx(1:3,3)=FV(ip,2:4)
              !
              ! find hv_value(1:3)
              !
              z=0.d0
              do ia=1, 3                          ! loop 3-points
                 ipx=key(ia)
                 nb2=nprox(ipx)
                 s=0.d0; s2=0.d0
                 do is2=1, nb2
                    jpx=prox(is2,ipx)
                    ispx=mod(is2, nb2)+1
                    kpx=prox(ispx,ipx)
                    !
                    do j=2,ndimFV
                       Fc(j,is2)= FV(ipx,j)*wtsph(1,is2,ipx)+FV(jpx,j)*wtsph(2,is2,ipx)+FV(kpx,j)*wtsph(3,is2,ipx)
                    enddo
                 enddo
                 !
                 do is2=1, nb2
                    is2m=mod(is2+nb2-2, nb2)+1    ! is2 minus
                    do j=2,ndimFV
                       Fs(j,is2)= 0.5d0*(Fc(j,is2m) + Fc(j,is2))   !  Vx,Vy,Vz on six sides, Fs = component on side
                    enddo
                 enddo
                 !
                 do is2=1, nb2
                    jpx=prox(is2,ipx)
                    is2m=mod(is2+nb2-2, nb2)+1    ! is2 minus
                    !
                    iV1= V4C(is2,ipx)
                    iV2= V4C(is2m,ipx)
                    hs = 0.5d0*(FV(iV1,1)+FV(iV2,1))
                    !
                    t=0.d0; t2=0.d0
                    do j=1,3
                       t=  t + Fs(j+1,is2)*Nvec(j,is2,ipx)       ! V \dot N
                       t2= t2+ Fs(j+1,is2)*Tvec(j,is2,ipx)       ! V \dot T
                    enddo
                    s=  s  + t  * hs * slen(is2,ipx)
                    s2= s2 + t2      * slen(is2,ipx)
                 enddo
                 hv_value(ia)=   s / area(ipx)
                 vort_value(ia)= s2/ area(ipx)
                 abvort_value(ia) = s2/ area(ipx) + fcori(ipx)
              enddo      !  ia=1,3
              div_hV= hv_value(1)*wtsph(1,is,ip) + &
                   hv_value(2)*wtsph(2,is,ip) + &
                   hv_value(3)*wtsph(3,is,ip)
              !
              !
              abvort_V= abvort_value(1)*wtsph(1,is,ip) + &
                   abvort_value(2)*wtsph(2,is,ip) + &
                   abvort_value(3)*wtsph(3,is,ip)
              !
              !
              !
           endif                 !  imethod
           !
           !
           !
           AFV(iV,1)= -div_hV    !  acc for h(iV)
           !
           ! TC10: polar vortex
           if (iswcase.EQ.10) then
              !-- add PV
              diag_var(iV,1)= abvort_V/g/FV(iV,1)    ! \eta / h:  Bates Li, Omega/Phi
           endif
           !
        endif                    !  if ( Vdc == 0 )
     enddo                       !  is
  enddo                          !  ip


!
!
!---  Very bad algorithm on CVT grid
!
!  !
!  !-- single out imethod = 1 or 2
!  !
!  do ip=1, nip
!     nb= nprox(ip)
!     ! 
!     !- loop 6 sides
!     do is=1, nb
!        isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
!        jp=prox(is,ip);    js=proxs(is,ip)
!        kp=prox(isp,ip);   ks=proxs(isp,ip)
!        nb2=nprox(jp);     jsm=mod(js+nb2-2, nb2) + 1
!        iV=V4C(is,ip)
!        !
!        do j=1,3       ! vec j-comp
!           V1(j)= 0.5d0*(FV(ip,j+1)+FV(jp,j+1))
!           V2(j)= 0.5d0*(FV(jp,j+1)+FV(kp,j+1))
!           V3(j)= 0.5d0*(FV(kp,j+1)+FV(ip,j+1))
!        enddo
!        
!        if (imethod .EQ. 1) then                   !  CVT  method-1
!           !
!           !-- apporx. with one-half
!           !
!           iV2= V4C(ism,ip)
!           hs1= 0.5d0*(FV(iV,1)+FV(iV2,1))
!           !
!           nb2=nprox(kp);ksp=mod(ks,nb2)+1;iV2=V4C(ksp,kp)
!           hs2= 0.5d0*(FV(iV,1)+FV(iV2,1))
!           !
!           iV2= V4C(isp,ip)
!           hs3= 0.5d0*(FV(iV,1)+FV(iV2,1))
!           !
!           s1=   V1(1)*Nvec2(1,is,ip) +V1(2)*Nvec2(2,is,ip) +V1(3)*Nvec2(3,is,ip)
!           s2=   V2(1)*Nvec2(1,jsm,jp)+V2(2)*Nvec2(2,jsm,jp)+V2(3)*Nvec2(3,jsm,jp)
!           s3= -(V3(1)*Nvec2(1,isp,ip)+V3(2)*Nvec2(2,isp,ip)+V3(3)*Nvec2(3,isp,ip))! Nvec2(ks,kp)=-(isp,ip)
!           !
!           s=   hs1*s1*dcenter_A(is,ip) + &
!                hs2*s2*dcenter_A(jsm,jp) + &
!                hs3*s3*dcenter_A(isp,ip)
!           div_hV_vert(is,ip)  = s/areaV_A(is,ip)
!           div_hV_vert(jsm,jp) = div_hV_vert(is,ip)
!           div_hV_vert(ks, kp) = div_hV_vert(is,ip)
!!!              if (ip.EQ.11) then
!!!                 write(6,101) 'V1, V2, V3', V1, V2, V3
!!!                 write(6,201) 'hs1, hs2, hs3', hs1, hs2, hs3
!!!                 write(6,201) 'div_hV', div_hV
!!!                 write(6,146) 'ip,is,dcenter_A(is,ip),areaV_A(is,ip)', &
!!!                      ip,is,dcenter_A(is,ip),areaV_A(is,ip)
!!!              endif
!!!              if (ip.GT.30) STOP 'ck 1'
!
!
!        elseif (imethod .EQ. 2) then  
!           !
!           !-- exact Lambda method
!           !
!           iV2= V4C(ism,ip)
!           r=lambdas(is,ip)/(1.d0-lambdas(is,ip))
!           hs1= (FV(iV,1) + r*FV(iV2,1))/(1.d0 + r)
!           !
!           nb2=nprox(jp);jsm2=mod(js+nb2-3,nb2)+1;iV2= V4C(jsm2,jp)
!           r=lambdas(jsm,jp)/(1.d0-lambdas(jsm,jp))
!           hs2= (FV(iV,1) + r*FV(iV2,1))/(1.d0 + r)
!           !
!           iV2= V4C(isp,ip)
!           r=lambdas(ks,kp)/(1.d0-lambdas(ks,kp))
!           hs3= (FV(iV,1) + r*FV(iV2,1))/(1.d0 + r)
!           !
!           s1=  V1(1)*Nvec2(1,is,ip) +V1(2)*Nvec2(2,is,ip) +V1(3)*Nvec2(3,is,ip)
!           s2=  V2(1)*Nvec2(1,jsm,jp)+V2(2)*Nvec2(2,jsm,jp)+V2(3)*Nvec2(3,jsm,jp)
!           s3=-(V3(1)*Nvec2(1,isp,ip)+V3(2)*Nvec2(2,isp,ip)+V3(3)*Nvec2(3,isp,ip))! Nvec2(ks,kp)=-(isp,ip)
!           !
!           s=   hs1*s1*dcenter_A(is,ip) + &
!                hs2*s2*dcenter_A(jsm,jp) + &
!                hs3*s3*dcenter_A(isp,ip)
!           !
!           !
!           div_hV_vert(is,ip)= s/areaV_A(is,ip)
!           div_hV_vert(jsm,jp) = div_hV_vert(is,ip)
!           div_hV_vert(ks, kp) = div_hV_vert(is,ip)
!
!           !!              if (ip.EQ.11) then
!           !!                 write(6,101) 'V1, V2, V3', V1, V2, V3
!           !!                 write(6,201) 'hs1, hs2, hs3', hs1, hs2, hs3
!           !!                 write(6,201) 'div_hV', div_hV
!           !!              endif
!           !              if (ip.GT.3) STOP 'ck 1'
!        endif
!     enddo   ! is
!  enddo      ! ip
!
!
!
!  IF (imethod.EQ.1 .OR. imethod.EQ.2) THEN
!  do ip=1, nip
!     nb= nprox(ip)
!     do is= 1, nb                       ! each vertex
!        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq 
!           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
!           kp=prox(isp,ip);   ks=proxs(isp,ip)
!           nb2=nprox(kp);     ksp=mod(ks,nb2)+1
!           iV = V4C(is,ip)
!           s1= div_hV_vert(ism,ip)
!           s2= div_hV_vert(ksp,kp)
!           s3= div_hV_vert(isp,ip)
!           div_hV= s1*wt2sph(1,is,ip)+s2*wt2sph(2,is,ip)+s3*wt2sph(3,is,ip)
!           AFV(iV,1)= -div_hV    !  acc for h(iV)
!
!!           write(6,147) 'iV, ip, is, AFV(iV,1), wt2sph(1:3,is,ip)', &
!!                iV, ip, is, AFV(iV,1), wt2sph(1:3,is,ip)
!
!        endif  ! if _ Vdc
!     enddo     ! is
!  enddo        ! ip
!  ENDIF
!
!!  stop 'ck 1'
!




  !
  !-- part 2  V(ip) = FV(ip, 2:4)  velocity Vx/Vy/Vz  
  !-- hex grid
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- s1. prep  h for mid edge:  exact 1/2
     !
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        jp=prox(is,ip); js=proxs(is,ip)
        j=V4C(is,ip); k=V4C(js,jp)    ! iV1= upper corner, iV2= lower
        Fs(1,is)= 0.5d0*(FV(j,1)+hb_a(j) + FV(k,1)+hb_a(k)) !  average (h^*+hs)  up + down
     enddo


     !
     !-- s2. prep  Fs(2:4,is) and Es(is)
     !
     if (jmethod.EQ.1) then             !  velocity interp using 1/2*(L + R)
        !
        !-- V, Kin, E_tot for mid edge
        !
        do is= 1, nb
           jp=prox(is,ip)        
           !
           do j=2,ndimFV
              Fs(j,is)= 0.5d0*(FV(ip,j)+FV(jp,j))  !  average V  ip-jp pair
           enddo
           !
           K_s(is)= 0.5d0 * &
                (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4)) + &
                0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4)))
           !-- E = K + g*ht
           Es(is) = K_s(is) + g*Fs(1,is)  ! (KE + PE) on six sides,  Es = energy on side
           !
        enddo

     elseif (jmethod.EQ.2) then    !  velocity interp using Tomita's method
        do is= 1, nb
           isp=mod(is, nb)+1
           jp=prox(is,ip)
           kp=prox(isp,ip)
           do j=2, ndimFV
              Fc(j,is)=FV(ip,j)*wtsph(1,is,ip)+FV(jp,j)*wtsph(2,is,ip)+FV(kp,j)*wtsph(3,is,ip)  ! Vx,y,z
           enddo
           !-- interpolate E = 1/2*v*v+g(h+h0) to prepare for grad(E)
           Kc(is)= &
                (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4)))*wtsph(1,is,ip)+&
                (0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4)))*wtsph(2,is,ip)+&
                (0.5d0*(FV(kp,2)*FV(kp,2)+FV(kp,3)*FV(kp,3)+FV(kp,4)*FV(kp,4)))*wtsph(3,is,ip)
        enddo
        !
        !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
        do is= 1, nb
           ism=mod(is+nb-2, nb)+1
           do j=2, ndimFV
              Fs(j,is)= 0.5d0*(Fc(j,ism) + Fc(j,is))   ! h^*,Vx,Vy,Vz on six sides, Fs = component on side
           enddo
           Es(is)     = 0.5d0*(Kc(ism)   + Kc(is)  ) + g*Fs(1,is)  ! (KE + PE) on six sides,  Es = energy on side
        enddo
        !
     endif          !  jmethod


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
        vort = vort  + t*slen(is,ip)     ! (v_i dot T) * ds_i
     enddo
     div_V= div_V/area(ip)
     vort = vort /area(ip)

     
     !-- s4:  gradE = grad[ K + g (h + hb_a) ]
     !
     s=0.d0
     nb2= nprox(ip)
     do is2=1, nb2
        iV2= V4C(is2,ip)
        s= s + (FV(iV2,1)+hb_a(iV2))*RAiv4H_A(is2,ip)
     enddo
     hi=s
     E0=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*hi
     !!E0=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))
     !
     do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds
        t=0.d0; z=0.d0
        do is= 1, nb
           t= t + Es(is)  *Nvec(j,is,ip)*slen(is,ip) ! E0 center, E0p: prime, replace E0 by six-side average
           z= z +          Nvec(j,is,ip)*slen(is,ip) !       n ds
        enddo
        gradE(j)= (t - z*E0)/area(ip)          ! grad( 1/2*V*V + g(h+hs))
     enddo


!     write(6, 101) 'div_V, vort, fcori(ip), gradE(1:3)', &
!          div_V, vort, fcori(ip), gradE(1:3)
!     stop 'ck 1'


     !
     !-- s5.  acceleration for velocity
     !----------------------------------------------------------------------
     ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
     !
     !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
     !----------------------------------------------------------------------
     !
     do i=1,3
        j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
        if (idamp .GT. 0) then
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)+diff_term(ip,i+1)
        else
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
           ! wo Coriolis force
           !        AFV(ip,i+1)=-(vort)*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
        endif
     enddo
     !
     !     AFV(ip,1)= dh/dt, AFV(ip,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
     !


     !
     !-- s6. remove radial component for \part V / \part t
     !
     s=0.d0               ! s= (Vector, Rcvec) not Rvec on edge
     do j=1, 3
        s= s+ AFV(ip,j+1)*Rcvec(j,ip)
     enddo
     do j=1, 3
        AFV(ip,j+1)= AFV(ip,j+1)-s*Rcvec(j,ip)    ! Vec = Vec - (grad,R) R
     enddo
     !
     !

     
     !
     !-- s7. reset to zero if sw_test_case_1 or advection
     !
     if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
        AFV(ip,2:4)=0.d0
     endif
     !
     !
  enddo  ! ip  1:nip



  deallocate (diff_term, div_hV_vert)
  return
  include '../swm/myformat.inc'
end subroutine dyn_DuCVTSPR





subroutine dyn_C_2comp_base (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_variables, only  :  lon, lat, lon_v, lat_v, nprox, prox, proxs, nip, npp, niE
  use module_variables, only  :  slen, areaI, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a   ! area(ip)=areaI(ip) MPAS
  use module_variables, only  :  areaE, areaV, dcenter, wtUt, RAiv4T, RAiv4H, iVtopkey
  use module_variables, only  :  Eindex, EindexVort, C4E, E4C, arctheta, LevdUn, NoutdUn
  use module_variables, only  :  omegasw, fcori_V
  use module_constants, only  :  g, ae, omega, pi, ndim
  use readnamelist,     only  :  iswcase, iconsv_scheme
  implicit none
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  integer  :: ip, jp, kp, Lp, iE, iE2, iV
  integer  :: is, isp, ism, is2, js, js2, ks, jsm     ! is plus 1 (is+1)
  integer  :: nb
  integer  :: i, j, k
  integer  :: ipc, isc          !  aux (is,ip) circle for two hexagon
  integer  :: nvert, iEh, iph, ish   ! vortex, and Edge head

  real(rt) :: E, EL, ER, Ut, E1, E2, E3, E4 ! Ener,  E_corner, E0 prime = average Es site
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradEn, gradE_ana(3), gradKn, gradgh    ! Energy gradient, analytic solution
  real(rt) :: vort, vort_ana, vortU, vortD    ! analytical sol
  real(rt) :: qeU, qeD
  real(rt) :: h_V_up, h_V_dn, h_V, etaUt

  real(rt) :: s, t, sum                 ! sum
  real(rt), allocatable :: qe_diag(:), Fe(:), h(:), Un(:), Ah(:), AUn(:)

  !    note
  !    FV(ip,1)=h
  !    FV(ip,2)=Un

  allocate (qe_diag      (ipe))
  allocate (Fe           (ipe))
  allocate (h            (nip))
  allocate (Ah           (nip))
  allocate (Un           (ipe))
  allocate (AUn          (ipe))


  do ip=1, nip
     h(ip)= FV(ip,1)
  enddo
  do iE=1, niE
     Un(iE)= FV(iE,2)
  enddo


  !---------------------------------------------
  !  step.1 :  tendency for h field at center of Control Volume
  !  - \p h / \p t =  div_( hV ) = 1/A \int  h_e * Un * le  * (U_n^e \dot N_out)
  !            h_e = (h_L + h_R)/2
  !---------------------------------------------
  ! MPAS Fe = h_e * U_n ^ e
  !
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE)
     Fe(iE) = 0.5d0 * ( h(ip) + h(jp) ) * Un(iE)
  enddo
  !
  do ip=1,nip
     nb= nprox(ip)
     s=0.d0
     do is= 1, nb          ! loop nn bond
        jp=prox(is,ip)
        iE=E4C(is,ip)
        s= s + (h(ip)+h(jp)) * Un(iE) * slen(is,ip) * NoutdUn(is,ip)  !  1/2*(h(ip) + h(jp))
!!        write(6, 141) 'is, Un(iE)', is, Un(iE)
     enddo
     Ah(ip)= -0.5d0 * s / areaI(ip)  ! acc for h ; Ringler paper Eq. (64) use linear areaI; 
                                     ! potential error due to areaI, \sum 1/4 * de * le ; 
                                     ! 
!!     write(6, 141) 'ip, h(ip), Ah(ip)', ip, h(ip), Ah(ip)
  enddo



  !---------------------------------------------
  !   step. 2 
  !   symmetric  (xi + f)_e * Ut^e
  !
  !    KE conservative way
  !    Q_e ^ perp = eta * U_t^e
  !        = 1/de \sum_e' W_ee' * l_e' * he' * U_n^e' * 1/2 * [ qe + qe' ]
  !        =                            ^^^^^^^^^^^^^
  !                                         Fe'
  !    q_e (iE) = 1/2 * ( q_v1 + q_v2)    
  !         q_v = ( f + xi ) / h_v
  !          xi = 
  !---------------------------------------------
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE); kp=C4E(3,iE); Lp=C4E(4,iE); is=C4E(5,iE)
     h_V_up =  &
          RAiv4T(1,iE)*h(ip) + &
          RAiv4T(2,iE)*h(jp) + &
          RAiv4T(3,iE)*h(kp)
     h_V_dn =  &
          RAiv4T(4,iE)*h(ip) + &
          RAiv4T(5,iE)*h(Lp) + &
          RAiv4T(6,iE)*h(jp)
     !
     vort=0.d0
     do is2=1,3
        iE2= EindexVort(is2,iE)
        vort= vort + Un(iE2)*LevdUn(is2,iE)*dcenter(iE2)          ! LevdUn with fixed iE
!!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortU= vort / areaV(1,iE) + fcori_V(1,iE)   ! total rotation quantum for Up Triangle
     !
     vort=0.d0
     do is2=4,6
        iE2= EindexVort(is2,iE)
        vort= vort + Un(iE2)*LevdUn(is2,iE)*dcenter(iE2)
     enddo
     vortD= vort / areaV(2,iE) + fcori_V(2,iE)   ! total rotation quantum for Down Triangle
     !
     qeU = VortU / h_V_up      !!!  WHAT IF  h == 0 ????, STOP 'failed '
     qeD = VortD / h_V_dn      !!!  WHAT IF  h == 0 ????, STOP 'failed '
     qe_diag(iE)= 0.5d0 * (qeU + qeD)      !  qe_hex_edge(iE)
     !
!!     write(6, 241) 'iE, qe_hex_edge(iE)', iE, qe_diag(iE)
     !
  enddo   ! loop iE



  !------------------------------------------------------------
  !  step.3 : 
  !  output:  grad (K + g(h+b))  +  Coriolis term
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
        E= E + Un(iE2) * Un(iE2) * areaE(iE2)
     enddo
     EL= E/areaI(ip) + g*(h(ip) + hb_a(ip))  !  E = K_i + gh_i    MPAS def.  areaE/area
     !-- Right
     nb= nprox(jp)
     E=0.d0
     do js2=1, nb
        iE2= E4C(js2,jp)
        E= E + Un(iE2) * Un(iE2) * areaE(iE2)
     enddo
     ER= E/areaI(jp) + g*(h(jp) + hb_a(jp))  !  E = K_i + gh_i  MPAS def. areaE/area
     gradEn = (ER - EL) / dcenter(iE) * NoutdUn(is,ip)

     !
     !  Coriolis term
     !
     nvert = nprox(ip) + nprox(jp) - 2        ! double Benzene
     s=0.d0
     t=0.d0
     do is2=1,nvert                           ! pentagon + hexagon =9; otherwise 10
        iE2= Eindex(is2,iE)                   ! 
        ipc= C4E(1,iE2); isc= C4E(5,iE2)        ! ip,is on circle
!!           s= s + wtUt(is2,iE)*slen(isc,ipc)*FV(iE2,2)*h_e(iE2)*0.5d0*(qe_diag(iE) + qe_diag(iE2)) 
        s= s + wtUt(is2,iE) * slen(isc,ipc) * Fe(iE2) * 0.5d0 * ( qe_diag(iE) + qe_diag(iE2) )
     enddo
     etaUt= s / dcenter(iE)
     AUn(iE)=  etaUt - gradEn

!!     write(6, 241) 'iE, gradEn, etaUt ', iE, gradEn, etaUt
     ! 
     ! as 
     ! \P Un / \pt - (f+vort)*U_t  + U_n^e \dot grad(E)= 0
     !             as \U_n \dot (k \cross U_n^e) = -1
     !
  enddo       ! iE



  do ip=1, nip
     AFV(ip,1)= Ah(ip)
  enddo
  do iE=ips,ipe     
     AFV(iE,2)= AUn(iE)
  enddo
  if ( iswcase .EQ. 1 ) then
     do iE=ips,ipe     
        AFV(iE,2)=0.d0    
     enddo
  endif

  deallocate ( qe_diag,Fe,h,Un,Ah,AUn )

  return
  include '../swm/myformat.inc'
end subroutine dyn_C_2comp_base




!----------------------------------------
!          Diagram  'A' grid
!
!             |
!             |
!  ip         |       kp   <--- (1) average h,v at vertex on sphere
!          /     \            map to GST
!     is /         \ 
!      /   js        \  (2) (hs,us,vs) side for flux
!    /                 \                      /
!  /         jp          \                  /
!                           - - - - - - - -
!                 (1) (hc,uc,vc) corner    \
!                                             \
!-----------------------------------------------


!--------------------------------------------------
!  sec.1  Dynamic equation
!--------------------------------------------------
!
!--- input (h,u,v) at cell center, 
!    output accerlation (fh,fu,fv) by finite volume method
subroutine dyn_A (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_variables, only  :  lon, lat, nprox, prox, proxs, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a, areaV_A, dcenter_A
  use module_variables, only  :  omegasw, diag_var
  use module_constants, only  :  g, ae, pi, ndim
  use readnamelist,     only  :  iswcase, idamp, iwachs
  use module_control,   only  :  hypervis
  use module_control,   only  :  nip          ! newly added
  implicit none
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  real(rt) :: hfx, divhvvx, divhvvy, phpx, phpy  ! see formula (1)
  integer  :: ip, jp, kp, Lp
  integer  :: is, js, ks, isp, ism      ! is plus 1 (is+1)
  integer  :: nb
  integer  :: i, j, k, ia, iap, iam, key(6)
  integer  :: lda, ldb, mdim
  integer  :: nc                ! counter

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, E0, Ec(npp)    ! Ener, E0 at ip,  E_corner
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: Fx(6,npp)         ! corner point;  1:4 h*/Vxyz; 5:htot; 6:Etot 
  real(rt) :: wtx(6)            ! weight corner point
  
  real(rt) :: Gc (ndimFV, npp)  ! filed at center      : c;  note index changed  ! h\vec V
  real(rt) :: Gs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: gradE(3), gradE_ana(3)   ! Energy gradient, analytic solution
  real(rt) :: gradh(3), gradh_ana(3), gradh_ana_sph(3), gradh_sph(3)  ! h gradient, analytic solution
  real(rt) :: vort, vort_ana    ! analytical sol
  real(rt) :: havg, h_base, ht_base   ! h base at ip, htot base at ip
  real(rt) :: s, t, z, sum                 ! sum
  real(rt) :: basis(3,3), ibasis(3,3)
  real(rt) :: vec(3), P_sph(3), vel_sph(3) ! 
  real(rt) :: lonc, latc, h0, Phi, Phi2, a1, a2, a3
  real(rt) :: div_V, div_hV  ! \div V or hV

  real(rt):: grad4_vec(4)                 ! \grad^4 (FV: h,Vx,Vy,Vz)
  real(rt), allocatable :: diff_term(:,:)
  integer :: ilevel, ncount, nloop,  npt_mx

  
  !
  !--------------  Meta code ----------------
  !
  !  for each ip 
  !   loop is  --- corner
  !   use isp
  !    find ip/jp/kp corner ip_sph(3) ...
  !    center_sph(1:3,is) for ijk
  !    cal weight: wt(1:3,is) for i/j/k pts
  !    cal u/v/h   by wt
  !   loop is  --- edge
  !    use ism          
  !    center (3,is) & center (3,ism)
  !    arc length:  slen
  !    nvec
  !    tvec
  !    h/u/v at mid
  !    div hV:  flux
  !    div 1/2 V*V
  !    vort_k
  !    grad(E)
  !----------------------------------------
  !    note
  !    FV(ip,1)=h^*    ,   hb_a is used over mountain 
  !    FV(ip,2)=Vx
  !    FV(ip,3)=Vy
  !    FV(ip,4)=Vz
  !
  !


  ! for idamp
  allocate (diff_term(nip,ndimFV))
  lonc=-pi/2.d0; latc=0.d0; h0=1000.d0    ! sw-test-1
  mdim=ndim; lda=ndim; ldb=ndim


  !--------------------------------------------------------------
  !  pre-calc. hypervisocosity = - eta * \nebala^4 V_x,V_y,V_z, h
  !--------------------------------------------------------------
  do j=2, ndimFV   !  no_h!,vx,vy,vz in Cartesian cord
     if (idamp .LT. 0) then
        do ip=1,nip
           diff_term(ip,j)=0.d0
        enddo
     elseif (idamp .EQ. 2) then
        call cal_Laplace_scalar_A (nip, FV(1,j), diff_term(1,j)) ! 2nd ord term =  \vis * Lap
     elseif (idamp .EQ. 4) then
        call cal_Laplace_square_A (nip, FV(1,j), diff_term(1,j)) ! 4th ord term = -\vis * Lap^2
     else
        STOP  'idamp not implemented'
     endif
  enddo

  Gc(1:ndimFV, 1:npp)=0.d0; Gs(1:ndimFV, 1:npp)=0.d0   !  Gvec = h \vec V  full mom. vector
  do ip=ips,ipe
     nb= nprox(ip)
     E0=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*(FV(ip,1)+hb_a(ip))
     h_base = FV(ip,1)
     !
     !-- s1.  h and velocity from hex center to mid edge
     !
     if (iwachs.EQ.0) then  ! use  Tomita ip-jp-kp -> vertex --> 1/2
        do is= 1, nb
           isp=mod(is, nb)+1
           jp=prox(is,ip)
           kp=prox(isp,ip)
           !-- weighted h,vx,vy,vz in Cartesian cord for upper corner
           do j=1, ndimFV
              Fc(j,is)=FV(ip,j)*wtsph(1,is,ip)+FV(jp,j)*wtsph(2,is,ip)+FV(kp,j)*wtsph(3,is,ip)  ! h^*,Vx,y,z
           enddo
           !-- interpolate E = 1/2*v*v+g(h+h0) to prepare for grad(E)
           Ec(is)= &
                (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*(FV(ip,1)+hb_a(ip)))*wtsph(1,is,ip)+&
                (0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4))+g*(FV(jp,1)+hb_a(jp)))*wtsph(2,is,ip)+&
                (0.5d0*(FV(kp,2)*FV(kp,2)+FV(kp,3)*FV(kp,3)+FV(kp,4)*FV(kp,4))+g*(FV(kp,1)+hb_a(kp)))*wtsph(3,is,ip)
           !!!-- weighted h \vec V
           !!do j=2, ndimFV   !  hVx, hVy, hVz   (2:4)
           !!   Gc(j,is)=FV(ip,1)*FV(ip,j)*wtsph(1,is,ip)+&
           !!        FV(jp,1)*FV(jp,j)*wtsph(2,is,ip)+&
           !!        FV(kp,1)*FV(kp,j)*wtsph(3,is,ip)         ! (h^* Vx  h^* Vy, h^* Vz)
           !!enddo
        enddo
        !
        !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
        do is= 1, nb
           ism=mod(is+nb-2, nb)+1
           do j=1, ndimFV
              Fs(j,is)= 0.5d0*(Fc(j,ism) + Fc(j,is))   ! h^*,Vx,Vy,Vz on six sides, Fs = component on side
              Es(is)  = 0.5d0*(Ec(ism)   + Ec(is)  )   ! (KE + PE) on six sides,  Es = energy on side
              !! Gs(j,is)= 0.5d0*(Gc(j,ism) + Gc(j,is))   ! h^* (Vx,Vy,Vz) vector:  j=1 unused
           enddo
        enddo
        !
        !
     elseif (iwachs.EQ.1) then  ! use 4-pt interp. ip/Lp/jp/kp  counter-clock count
        !
        !
        do is= 1, nb
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
           jp=prox(is, ip)
           kp=prox(isp,ip)
           Lp=prox(ism,ip)
           key(1)=ip
           key(2)=Lp
           key(3)=jp
           key(4)=kp
           !
           ! loop component
           do ia=1,4                   ! 4-corner pts
              do j=1, ndimFV           ! h*, Vxyz
                 Fx(j,ia)= FV(key(ia),j)
              enddo
           enddo
           !
           ! j = 5  htot
           j=5
           do ia=1, 4
              Fx(j,ia)= Fx(1,ia) + hb_a(key(ia))
           enddo
           !
           ! j = 6  Etot = K + g*(h + hb)
           j=6
           do ia=1, 4
              s=0.d0
              do k=2,ndimFV
                 s= s + 0.5d0 * Fx(k,ia) * Fx(k,ia)
              enddo
              Fx(j,ia)= s + g*Fx(5,ia)
           enddo
           !
           ! weighted linear interp
           ! 
!           Fs(1:4,ia)=0.d0
!           Es(ia)=0.d0
!

!           write(6,121) 'ip, is', ip, is
           do j=1, ndimFV
!              write(6,201) 'bf wachs, Fs(j,is)=', Fs(j,is)
              s=0.d0
              do ia=1, 4
                 s= s+Fx(j,ia)*wtsph(ia,is,ip)
              enddo
              Fs(j,is)=s
!              write(6,201) 'af wachs, Fs(j,is)=', Fs(j,is)
           enddo

           j=6
!           write(6,201) 'bf wachs, Es(is)=', Es(is)
           t=0.d0
           do ia=1, 4
              t= t+Fx(j,ia)*wtsph(ia,is,ip)
           enddo
           Es(is)=t
!           write(6,201) 'af wachs, Es(is)=', Es(is) 

        enddo    ! is
     endif       ! iwachs
!     stop 'ck 1'



     !
     !-- s2 gradient h; gradE = grad[ K + g (h + hb_a) ]
     !      note: remove \Rcvec component at final step for AFV(ip,2:4), not at here
     !
     do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds 
        s=0.d0; t=0.d0; z=0.d0
        do is= 1, nb
!           s= s + (Fs(1,is)-h_base)*Nvec(j,is,ip)*slen(is,ip) ! (h-hc)* n ds; 1=h, 2_4=Vx_z
!           t= t + (Es(is) - E0)*Nvec(j,is,ip)*slen(is,ip) ! E0 center, E0p: prime, replace E0 by six-side average
!
           s= s + Fs(1,is)*Nvec(j,is,ip)*slen(is,ip) ! (h) * n ds; 1=h, 2_4=Vx_z
           t= t + Es(is)  *Nvec(j,is,ip)*slen(is,ip) ! E0 center, E0p: prime, replace E0 by six-side average
           z= z +          Nvec(j,is,ip)*slen(is,ip) !       n ds
        enddo
        gradh(j)= (s - z*h_base)/area(ip)          ! this is grad(h^*)
        gradE(j)= (t - z*E0    )/area(ip)          ! grad( 1/2*V*V + g(h+hs))
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


     !
     !-- s4.  div (hV) =  \int  (hV \dot N) * ds  / A
     !
!     !-- method-2:  proven worse than method-1 by TC5_A_G6_90d_nodamp and ef32
!     div_hV=0.d0
!     do is= 1, nb
!        s=0.d0
!        do j=1,3    !  s= ( hV dot N)
!          s= s + Gs(j+1,is)*Nvec(j,is,ip)      ! projection to N
!        enddo
!        div_hV=div_hV + s*slen(is,ip)  ! \int (h^*V dot N) ds
!     enddo
!     div_hV= div_hV/area(ip)
!     AFV(ip,1)= -div_hV


!     !-- method-1  standard
     div_hV=0.d0
     do is= 1, nb
        s=0.d0
        do j=1,3    !  s= (V dot N)
          s= s + Fs(j+1,is)*Nvec(j,is,ip)      ! projection to N
        enddo
        div_hV=div_hV + s*Fs(1,is)*slen(is,ip)  ! h^*  (V dot N)
     enddo
     div_hV= div_hV/area(ip)
     AFV(ip,1)= -div_hV


!-- YGYU check difference
!     write(6,121) 'ip', ip
!     write(6,201) 'diff div_hV(M2) -M1',&
!          -AFV(ip,1), div_hV, -AFV(ip,1)-div_hV
!
!     if (ip.GT.11) then
!        write(6101) 'stop at method-2 - m1'
!        STOP 
!     endif


     !-- method-3
!!     ---  YGYU: 26-March-2020 :
!!         Above code use full FV
!!         Code below in logic is more accurate, 
!!         but crashes at 38d-- 56d in TC2 G5  ( L2 = 0.1E-1 )
!!     !
!!     !-- s4.  div (hV) =  < V(center), \grad (h) >  + h(center)* div_V
!!     !        
!!     !        discard inferior:  1/A  \int (h V)_mid \dot n  ds ; bc twice interpolation
!!     !
!!     s=0.d0    ! s = (V, g(h))
!!     do j=1,3
!!        s= s + gradh(j)*FV(ip,j+1)  ! V is in plane, dot  gradh(1:3), eliminate \R vec
!!     enddo
!!     AFV(ip,1)= -(FV(ip,1)*div_V + s)




     !
     !-- s5.  acceleration for velocity
     !----------------------------------------------------------------------
     ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
     !
     !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
     !----------------------------------------------------------------------
     !
     do i=1,3
        j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
        if (idamp .GT. 0) then
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)+diff_term(ip,i+1)
        else
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
           ! wo Coriolis force
           !        AFV(ip,i+1)=-(vort)*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
        endif
     enddo
     !
     !     AFV(ip,1)= dh/dt, AFV(ip,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
     !



     !
     !-- s6. remove radial component for \part V / \part t
     !
     s=0.d0               ! s= (Vector, Rcvec) not Rvec on edge
     do j=1, 3
        s= s+ AFV(ip,j+1)*Rcvec(j,ip)
     enddo
     do j=1, 3
        AFV(ip,j+1)= AFV(ip,j+1)-s*Rcvec(j,ip)    ! Vec = Vec - (grad,R) R
     enddo
     !
     !

     
     
     !
     !-- s7. reset to zero if sw_test_case_1 or advection
     !
     if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
        AFV(ip,2:4)=0.d0
     elseif (iswcase.EQ.10) then
        !-- add PV
        diag_var(ip,1)= (vort+fcori(ip))/g/FV(ip,1)    ! \eta / h
     endif
     !
     !
  enddo  ! ip  outer loop



  deallocate (diff_term)
  return
  include '../swm/myformat.inc'
end subroutine dyn_A





!--------------------------------------------------
!   sec.5  C-grid dynamics
!--------------------------------------------------
!
!---------------------------------------------------------------------
!          Diagram 'C' grid
!                                    |
!                                  /    \       
!                               /         \
!                            /              \
!                          /                  \
!                        /                       \
!               t_ev   |                          |
!               L_ev   |           (ks,kp)        |
!         .   <-----   | . .  . . .  . .          |
!          .           |              .  .        |
!           .          | V2          .    .       | 
!            .       /  \           .      .    /   \
!             .    /      \        .         /         \
!                /          \    .          /.            \
!              /  .           \           /   .  LevdUn      \      LevdUn = L_e(Vort) \dot U_n^e
!            /      .         . \   V1  /       .              \
!       V3 /         .      .     \   /           .              \
!         |           .   .         |               .             |
!         |             .           | t_ev           .            |
!         |                         | L_ev                        |
!         |    (is,ip) iE . . . . . | ----->        (js,jp)       | --->
!         |                         |                             | NoutdUn = n^out \dot U_n^e
!       V4|                         |                             |
!          \              (4)      / \   (5)  (1/2/3/4/5) LevdUn  |
!            \                   /     \                        /
!               \              /         \                    /
!                 \         /     (Ls,Lp)   \              /
!                   \    /                     \        /
!                     V5                          \   /
!
!
!---------------------------------------------------------------------
! Note:
!       L_ev(is,ip) / n_out ,   n_out / n_ev [is, ip] =  n_out /  NoutdUn(is,ip)
!
subroutine dyn_C (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_control,   only  :  dt
  use module_variables, only  :  lon, lat, lon_v, lat_v, nprox, prox, proxs, nip, npp, niE
  use module_variables, only  :  slen, areaI, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a   ! area(ip)=areaI(ip) MPAS
  use module_variables, only  :  areaE, areaV, dcenter, wtUt, RAiv4T, RAiv4H, iVtopkey
  use module_variables, only  :  Eindex, EindexVort, C4E, E4C, arctheta, LevdUn, NoutdUn
  use module_variables, only  :  omegasw, fcori_V
  use module_constants, only  :  g, ae, omega, pi, ndim
  use readnamelist,     only  :  iswcase, iconsv_scheme
  implicit none
  real(rt), parameter  :: third=0.3333333333333333333
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
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
  real(rt) :: qeU, qeD
  real(rt) :: h_V_up, h_V_dn, h_V
  real(rt) :: factor_PV, Fe_perp

  real(rt) :: havg              ! analytical sol
  real(rt) :: s, t, sum                 ! sum
  real(rt) :: basis(3,3), ibasis(3,3)
  real(rt) :: vec(3), P_sph(3), vel_sph(3) ! 
  real(rt) :: lonc, latc, h0, Phi, Phi2, a1, a2, a3
  real(rt) :: VortV(10)
  real(rt) :: VortE(ims:ime)
  real(rt) :: fcori_av, vort_av, etaUt
  real(rt), allocatable :: qe_diag(:), qe_prog(:), h_e(:), Fe(:)
  real(rt), allocatable :: qe_hex_center(:), qe_hex_vertex(:,:), qe_hex_edge(:) 
  real(rt), allocatable :: Un_dot_gradPV(:), Ut_dot_gradPV(:)


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


!  allocate (TestF(ipe))     ! test field F = - omega \dot r
!  allocate (div_V(ipe))    
!  allocate (qe_prog(ipe))
  allocate (qe_hex_vertex(2,ipe))
  allocate (qe_hex_center(nip))
  allocate (qe_hex_edge  (ipe))
  allocate (qe_diag      (ipe))
  allocate (h_e          (ipe))
  allocate (Fe           (ipe))
  allocate (Un_dot_gradPV(ipe))
  allocate (Ut_dot_gradPV(ipe))


  lonc=-pi/2.d0; latc=0.d0; h0=1000.d0    ! sw-test-1
  mdim=ndim; lda=ndim; ldb=ndim
  factor_PV= 0.5d0
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
     !write(6, 121) 'ip=', ip
     do is= 1, nb          ! loop nn bond
        jp=prox(is,ip)
        iE=E4C(is,ip)
        !s= s+0.5d0*(FV(ip,1)+FV(jp,1))*FV(iE,2)*slen(is,ip)*NoutdUn(is,ip)
        s= s + Fe(iE) * slen(is,ip) * NoutdUn(is,ip)        ! NOTE: could use (iE) only 
        !write(6, 121) 'is,jp,iE', is,jp,iE
        !write(6, 101) 'FV(ip,1),FV(jp,1),FV(iE,2),slen(is,ip)', &
        !     FV(ip,1),FV(jp,1),FV(iE,2),slen(is,ip)
        !
     enddo
     AFV(ip,1)= -s / areaI(ip)  ! acc for h ; Ringler paper Eq. (64) use linear areaI; 
                                ! potential error due to areaI, \sum 1/4 * de * le ; 
  enddo

!!  write(6,101) 'dyn C pt 2'



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
     !
     vort=0.d0
     do is2=4,6
        iE2= EindexVort(is2,iE)
        vort= vort + FV(iE2,2)*LevdUn(is2,iE)*dcenter(iE2)
!!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortD= vort / areaV(2,iE) + fcori_V(2,iE)   ! total rotation quantum for Down Triangle
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
     do is=1, nb
        iE= E4C(is,ip)
        iV= iVtopkey(is,ip)    !      flip up or dn
        !  write(6, 121) 'is,ip,iE,iV', is,ip,iE,iV
        !
        ! -- I always have counter clock wise positive for vortex
        !    (is,ip) point to upper,  (iV, iE) should correspond to the same position
        s= s + qe_hex_vertex(iV,iE) * RAiv4H(is,ip)
     enddo
     qe_hex_center(ip)= s
     !
     !--  For test TC2 only
     !    passing this test means iVtopkey is correct
     !
!     write(6, 201) 'qe_hex_center(ip)       ', qe_hex_center(ip)
!     write(6, 201) 'qe_hex_center(ip) analy ', 2.d0*(omega + omegasw)*sin(lat(ip)) / FV(ip,1)
     !
     ! comment: this interpolation has a bit larger error
     !-0.000000085926  vs.  qe_hex_center(ip) analy    -0.000000087323
     !
  enddo


  !
  ! gradient PV
  ! N dir
  !
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)
     Un_dot_gradPV(iE) = &
          - ( qe_hex_center(ip) - qe_hex_center(jp) ) / dcenter(iE) * FV(iE,2) * NoutdUn(is,ip)
  enddo
  !
  ! T dir  ( ?? still need to verify )
  !
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE)
     is=C4E(5,iE); js=proxs(is,ip)
     nb= nprox(ip)
     ism=mod(is+nb-2, nb)+1
     !
     nvert = nprox(ip) + nprox(jp) - 2        ! double Benzene
     s=0.d0
     do is2=1,nvert                           ! pentagon + hexagon =9; otherwise 10
        iE2= Eindex(is2,iE)                   ! 
        ipc=C4E(1,iE2); isc=C4E(5,iE2)        ! ip,is on circle
        s= s + wtUt(is2,iE)*slen(isc,ipc)*FV(iE2,2)
     enddo
     Ut= s/dcenter(iE)
     Ut_dot_gradPV(iE) = &
          Ut * ( qe_hex_vertex(2,iE) - qe_hex_vertex(1,iE) ) / slen(is,ip) * LevdUn(1,iE)
     !
     !    from Formula on notebook 9
     !               down_vertex          up_vertex
     !
  enddo
  !
  ! PV advection, update the edge values
  !
  if ( iconsv_scheme .EQ. 1 .OR. iconsv_scheme .EQ. 2 ) then
     do iE=ips,ipe
        qe_diag(iE)= qe_hex_edge(iE)
     enddo
  elseif (iconsv_scheme .EQ. 3) then
     do iE=ips,ipe
        qe_diag(iE)= qe_hex_edge(iE) &
             - factor_PV * dt * ( Un_dot_gradPV(iE) + Ut_dot_gradPV(iE) )
     enddo
  endif


!!  write(6,101) 'dyn C pt 3'


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
        !-- K_i conserve,  qe diagnostic 
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


!        if ( iE .LE. 20 ) then
!           write(6, 101) 'FV(iE,3), Ut, diff', FV(iE,3), Ut, FV(iE,3) - Ut
!        else
!           STOP 'stop check Ut'
!        endif


     elseif (iconsv_scheme .EQ. 2) then
        !
        !-- Potential Estrophy conserve:  \int xi^2 dV
        !
        !   eta^e * Ut^e  = ( f + xi )^e * U_t^e  = Fe^perp * qe
        !---------------------------------------------
        !   Ut      = 1 / de * \sum_e'  wtUt (e, e') * l_e' * Un_e'
        !   Fe^perp = 1 / de * \sum_e'  wtUt (e, e') * l_e' * Un_e' * h_e'
        !   Qe^perp = Fe^perp * q_e
        !---------------------------------------------
        s=0.d0
        do is2=1,nvert                           !  pentagon + hexagon =9; otherwise 10
           iE2= Eindex(is2,iE)
           ipc=C4E(1,iE2); isc=C4E(5,iE2)        ! ip,is on circle       
           s= s+ wtUt(is2,iE)*slen(isc,ipc)*Fe(iE2)
        enddo
        Fe_perp=s/dcenter(iE)
        etaUt = Fe_perp * qe_diag(iE)
     endif


     !
     ! \P Un / \pt - (f+vort)*U_t  + U_n^e \dot grad(E)= 0
     !             as \U_n \dot (k \cross U_n^e) = -1
     !
     AFV(iE,2)=  etaUt - gradEn
     AFV(iE,3)=  0.d0   ! in case defined
     !
     ! -- previous dated  12/28/2017
     ! AFV(iE,2)= -(fcori(iE)+vort)*Ut-gradEn          !  coriolios force is not included
     !
     !
  enddo       ! iE


!!  write(6,101) 'dyn C pt 4'



  !------------------------------------------------------------
  !   reset acceleration if sw_test_case_1
  !         passive advection only
  !------------------------------------------------------------
  !
  if ( iswcase .EQ. 1 ) then
     do iE=ips,ipe     
        AFV(iE,2)=0.d0    
     enddo
!  elseif( iswcase .GT. 2 ) then
!     STOP 'iswcase .GT. 2,  STOP'
  endif

  deallocate ( qe_diag, h_e, Fe )
  deallocate ( qe_hex_center, qe_hex_vertex, qe_hex_edge )
  deallocate ( Un_dot_gradPV, Ut_dot_gradPV )

  return
  include '../swm/myformat.inc'
end subroutine dyn_C



subroutine dyn_C_combine_Cori_kinetic (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_control,   only  :  dt
  use module_variables, only  :  lon, lat, lon_v, lat_v, nprox, prox, proxs, nip, npp, niE
  use module_variables, only  :  slen, areaI, Nvec, Tvec, Rvec, Rcvec, wtsph, hb_a   ! area(ip)=areaI(ip) MPAS
  use module_variables, only  :  Nvec2, Tvec2, Rvec2, dx_ov_le
  use module_variables, only  :  areaE, areaV, dcenter, wtUt, RAiv4T, RAiv4H, iVtopkey
  use module_variables, only  :  Eindex, EindexVort, C4E, E4C, arctheta, LevdUn, NoutdUn
  use module_variables, only  :  omegasw, fcori_V
  use module_constants, only  :  g, ae, omega, pi, ndim
  use module_control,   only  :  hypervis
  use readnamelist,     only  :  iswcase, iconsv_scheme, imodify_C, idamp
  implicit none
  real(rt), parameter  :: third=0.3333333333333333333
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  real(rt) :: h1, h2, h3
  real(rt) :: hfx, divhvvx, divhvvy, phpx, phpy  ! see formula (1)
  integer  :: ip, jp, kp, Lp, iE, iE1, iE2, iE3, iV
  integer  :: is, isp, ism, is2, js, js2, ks, jsm, jsp, Ls, Lsm     ! is plus 1 (is+1)
  integer  :: jpk, jpL, jsk, jsL
  integer  :: nb, nb2
  integer  :: i, j, k
  integer  :: lda, ldb, mdim
  integer  :: ipc, isc          !  aux (is,ip) circle for two hexagon
  integer  :: nvert, iEh, iph, ish   ! vortex, and Edge head

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, EL, ER, Ut, E1, E2, E3, E4 ! Ener,  E_corner, E0 prime = average Es site
  real(rt) :: Ki, Kj, Kk, KL
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradEn, gradEn_FD, gradE_ana(3), gradKn, gradgh    ! Energy gradient, analytic solution
  real(rt) :: gradh(3), gradh_ana(3), gradh_ana_sph(3), gradh_sph(3)  ! h gradient, analytic solution
  real(rt) :: grad2norm                ! Energy gradient  ||  X - X_analytic ||_2  2 norm
  real(rt) :: vort, vort_ana, vortU, vortD    ! analytical sol
  real(rt) :: qeU, qeD
  real(rt) :: h_V_up, h_V_dn, h_V, K_V_up, K_V_dn
  real(rt) :: factor_PV, Fe_perp

  real(rt) :: havg              ! analytical sol
  real(rt) :: s, t, sum, s2, s3                 ! sum
  real(rt) :: grad_interp  
  real(rt) :: basis(3,3), ibasis(3,3)
  real(rt) :: vec1(3), vec2(3), vec3(3), P_sph(3), vel_sph(3) ! 
  real(rt) :: lonc, latc, h0, Phi, Phi2, a1, a2, a3
  real(rt) :: VortV(10)
  real(rt) :: VortE(ims:ime)
  real(rt) :: fcori_av, vort_av, etaUt
  real(rt), allocatable :: qe_diag(:), qe_prog(:),h_e(:), h_e_4pt(:), Fe(:)
  real(rt), allocatable :: qe_hex_center(:), qe_hex_vertex(:,:), qe_hex_edge(:) 
  real(rt), allocatable :: Un_dot_gradPV(:), Ut_dot_gradPV(:)
  real(rt), allocatable :: gradK_d_Ut(:)
  real(rt), allocatable :: eta_edge(:)
  real(rt), allocatable :: K_center(:)
  real(rt), allocatable :: diff_term(:,:)
                           !aux_s(:)     ! Laplace damping term


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


!  allocate (TestF(ipe))     ! test field F = - omega \dot r
!  allocate (div_V(ipe))    
!  allocate (qe_prog(ipe))
  allocate (qe_hex_vertex(2,ipe))
  allocate (qe_hex_center(nip))
  allocate (qe_hex_edge  (ipe))
  allocate (qe_diag      (ipe))
  allocate (h_e          (ipe))
  allocate (h_e_4pt      (ipe))
  allocate (Fe           (ipe))
  allocate (Un_dot_gradPV(ipe))
  allocate (Ut_dot_gradPV(ipe))
  allocate (gradK_d_Ut   (ipe))
  allocate (eta_edge     (ipe))
  allocate (K_center     (nip))
  allocate (diff_term (niE,ndimFV))
  

  lonc=-pi/2.d0; latc=0.d0; h0=1000.d0    ! sw-test-1
  mdim=ndim; lda=ndim; ldb=ndim
  factor_PV= 0.5d0




  !--------------------------------------------------------------
  !--  pre-calc.  hypervisocosity = + eta * \nebala^2 U_n  only 
  !--------------------------------------------------------------
  do j=2, 2   !  no_h!,  Un only   C-grid  ndimFV should be 2 
     if (idamp .LT. 0) then
        do iE=1,niE
           diff_term(iE,j)=0.d0
        enddo
     elseif (idamp .EQ. 2) then
        call cal_Laplace_Un_C   (niE, FV(1,j), diff_term(1,j))  ! 2nd ord term =  \vis * Lap
     elseif (idamp .EQ. 4) then
           call cal_Laplace_square_Un_C (niE, FV(1,j), diff_term(1,j)) ! 4th ord term = -\vis * Lap^2
     else
        STOP  'idamp not implemented'
     endif
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
     nb= nprox(ip)
     ism= mod(is+nb-2, nb) + 1
     Fe(iE) = 0.5d0 * ( FV(ip,1) + FV(jp,1) ) * FV(iE,2)
     h_e(iE)= 0.5d0 * ( FV(ip,1) + FV(jp,1) )          !  questionalbe
     h_V_up = FV(ip,1)*wtsph(1,is,ip)+FV(jp,1)*wtsph(2,is,ip)+FV(kp,1)*wtsph(3,is,ip)  ! h
     h_V_dn = FV(ip,1)*wtsph(1,ism,ip)+FV(Lp,1)*wtsph(2,ism,ip)+FV(jp,1)*wtsph(3,ism,ip)  ! h
     h_e_4pt(iE)= 0.5d0 * ( h_V_up +  h_V_dn )

!     write(6, 101) 'wtsph(1,is,ip)',  wtsph(1,is,ip)
!     write(6, 101) 'wtsph(1,ism,ip)',  wtsph(1,ism,ip)
!     write(6, 101) 'h_e(iE) - h_e_4pt(iE)', h_e(iE) - h_e_4pt(iE) 
!     write(6, 101) 'h_e(iE)', h_e(iE)
  enddo



  do ip=1,nip
     nb= nprox(ip)
     s=0.d0
!     write(6, 121) 'ip=', ip
     do is= 1, nb          ! loop nn bond
!        write(6, 121) 'is=', is
        jp=prox(is,ip)
        iE=E4C(is,ip)
        !s= s+0.5d0*(FV(ip,1)+FV(jp,1))*FV(iE,2)*slen(is,ip)*NoutdUn(is,ip)
        if (imodify_C.EQ.0 .OR. imodify_C.EQ.1)  then !  TRSK or (Cori+Kinetic)
           s= s + Fe(iE) * slen(is,ip) * NoutdUn(is,ip)        ! NOTE: could use (iE) only 
        elseif (imodify_C.EQ.2 .OR. imodify_C.EQ.3)  then !  Linear H_corr  or (Cori+Kinetic)+H_corr
           s= s + h_e_4pt(iE) * FV(iE,2) * slen(is,ip) * NoutdUn(is,ip)
        endif
        !write(6, 121) 'is,jp,iE', is,jp,iE
        !write(6, 101) 'FV(ip,1),FV(jp,1),FV(iE,2),slen(is,ip)', &
        !     FV(ip,1),FV(jp,1),FV(iE,2),slen(is,ip)
        !
     enddo
     AFV(ip,1)= -s / areaI(ip)  ! acc for h ; Ringler paper Eq. (64) use linear areaI; 
                                ! potential error due to areaI, \sum 1/4 * de * le ; 
  enddo

!  write(6, 101)  'nail 1'


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
     !
     vort=0.d0
     do is2=4,6
        iE2= EindexVort(is2,iE)
        vort= vort + FV(iE2,2)*LevdUn(is2,iE)*dcenter(iE2)
!!        write(6, 121) 'iE,is2,EindexVort(is2,iE)', iE,is2,iE2
     enddo
     vortD= vort / areaV(2,iE) + fcori_V(2,iE)   ! total rotation quantum for Down Triangle

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
     eta_edge(iE)= 0.5d0 * (vortU + vortD)
     !
     !  this is diagnostic q_e  =  PV
     !
!     write(6, 201) 'qe_hex_edge(iE)       ', qe_hex_edge(iE)
!     write(6, 201) 'qe_hex_edge(iE) analy ', 2.d0*(omega + omegasw)*sin(lat_v(iE)) / h_e(iE)
     !
  enddo   ! loop iE


!  write(6, 101)  'nail 2'



  ! find PV at center 
  do ip=1, nip
     nb= nprox(ip)
     s= 0.d0
     do is=1, nb
        iE= E4C(is,ip)
        iV= iVtopkey(is,ip)    !      flip up or dn
        !  write(6, 121) 'is,ip,iE,iV', is,ip,iE,iV
        !
        ! -- I always have counter clock wise positive for vortex
        !    (is,ip) point to upper,  (iV, iE) should correspond to the same position
        s= s + qe_hex_vertex(iV,iE) * RAiv4H(is,ip)
     enddo
     qe_hex_center(ip)= s
     !
     !--  For test TC2 only
     !    passing this test means iVtopkey is correct
     !
!     write(6, 201) 'qe_hex_center(ip)       ', qe_hex_center(ip)
!     write(6, 201) 'qe_hex_center(ip) analy ', 2.d0*(omega + omegasw)*sin(lat(ip)) / FV(ip,1)
     !
     ! comment: this interpolation has a bit larger error
     !-0.000000085926  vs.  qe_hex_center(ip) analy    -0.000000087323
     !
  enddo


!  write(6, 101)  'nail 3'



  !
  ! gradient PV
  ! N dir
  !
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)
     Un_dot_gradPV(iE) = &
          - ( qe_hex_center(ip) - qe_hex_center(jp) ) / dcenter(iE) * FV(iE,2) * NoutdUn(is,ip)
  enddo
  !
  ! T dir  ( ?? still need to verify )
  !
  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE)
     is=C4E(5,iE); js=proxs(is,ip)
     nb= nprox(ip)
     ism=mod(is+nb-2, nb)+1
     !
     nvert = nprox(ip) + nprox(jp) - 2        ! double Benzene
     s=0.d0
     do is2=1,nvert                           ! pentagon + hexagon =9; otherwise 10
        iE2= Eindex(is2,iE)                   ! 
        ipc=C4E(1,iE2); isc=C4E(5,iE2)        ! ip,is on circle
        s= s + wtUt(is2,iE)*slen(isc,ipc)*FV(iE2,2)
     enddo
     Ut= s/dcenter(iE)
     Ut_dot_gradPV(iE) = &
          Ut * ( qe_hex_vertex(2,iE) - qe_hex_vertex(1,iE) ) / slen(is,ip) * LevdUn(1,iE)
     !
     !    from Formula on notebook 9
     !               down_vertex          up_vertex
     !
  enddo
  !
  ! PV advection, update the edge values
  !
  if ( iconsv_scheme .EQ. 1 .OR. iconsv_scheme .EQ. 2 ) then
     do iE=ips,ipe
        qe_diag(iE)= qe_hex_edge(iE)
     enddo
  elseif (iconsv_scheme .EQ. 3) then
     do iE=ips,ipe
        qe_diag(iE)= qe_hex_edge(iE) &
             - factor_PV * dt * ( Un_dot_gradPV(iE) + Ut_dot_gradPV(iE) )
     enddo
  endif


!  write(6, 101)  'nail 4'
!  stop 'nail 4'



  !------------------------------------------------------------
  !  step.3 : 
  !  -- K_center (ip)
  !  -- gradK_d_Ut (iE)
  !  -- gradEn, gradgh
  !  -- etaUt
  !
  !  output:  (1) U_t \dot grad (K)
  !           (2) Q_e^ perp  symmetrized tendency
  !            
  !------------------------------------------------------------

  ! find K_center(ip)
  do ip=1, nip
     nb= nprox(ip)
     E= 0.d0
     do is=1, nb
        iE= E4C(is,ip)
        E= E + FV(iE,2) * FV(iE,2) * areaE(iE)
     enddo
     K_center(ip)= E/areaI(ip)
  enddo


!  write(6, 101)  'nail 5'
!  stop 'nail 5'



  do iE=ips,ipe
     ip=C4E(1,iE); jp=C4E(2,iE); is=C4E(5,iE)
     kp=C4E(3,iE); Lp=C4E(4,iE)     
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
     !-- Kinetic energy for ip, jp, kp, Lp
     !
     Ki = K_center(ip)
     Kj = K_center(jp)
     Kk = K_center(kp)
     KL = K_center(Lp)
     K_V_up = Ki*wtsph(1,is,ip)  + Kj*wtsph(2,is,ip)  + Kk*wtsph(3,is,ip) 
     K_V_dn = Ki*wtsph(1,ism,ip) + KL*wtsph(2,ism,ip) + Kj*wtsph(3,ism,ip) 
     gradK_d_Ut(iE) = ( K_V_up - K_V_dn ) / slen(is,ip) * NoutdUn(is,ip)
!     gradK_d_Ut(iE) = ( K_V_up - K_V_dn ) / slen(is,ip) 
  enddo 


!  write(6, 101)  'nail 6'
!  stop 'nail 6'


!
!
!     E3= g*(FV(ip,1) + hb_a(ip))
!     EL= E/areaI(ip) + g*(FV(ip,1) + hb_a(ip))  !  E = K_i + gh_i    MPAS def.  areaE/area
!     !-- Right
!     nb= nprox(jp)
!     E=0.d0
!     do js2=1, nb
!        iE2= E4C(js2,jp)
!        E= E + FV(iE2,2) * FV(iE2,2) * areaE(iE2)
!     enddo
!     E2= E/areaI(jp)                ! test grad(K)
!     E4= g*(FV(jp,1) + hb_a(jp))
!     ER= E/areaI(jp) + g*(FV(jp,1) + hb_a(jp))  !  E = K_i + gh_i  MPAS def. areaE/area
!     gradEn = - (EL - ER) / dcenter(iE) * NoutdUn(is,ip)
!     gradKn = - (E1 - E2) / dcenter(iE) * NoutdUn(is,ip)
!     gradgh = - (E3 - E4) / dcenter(iE) * NoutdUn(is,ip)
!
  do iE=ips,ipe
     !---------------------------------------------
     !    (xi + f)_e * Ut^e  - U_n \dot  grad K
     !  = 1/de \sum_e' W_ee' * l_e' * ( \eta * U_n +  U_t \dot grad K )
     !    KE conservative way
     !    Q_e ^ perp = eta * U_t^e
     !        = 1/de \sum_e' W_ee' * l_e' * he' * U_n^e' * 1/2 * [ qe + qe' ]
     !        =                            ^^^^^^^^^^^^^
     !                                         Fe'
     !    q_e (iE) is prepared in above section
     !---------------------------------------------
     !
     !-- grad (g h)
     !
     ip=C4E(1,iE); jp=C4E(2,iE); kp=C4E(3,iE); Lp=C4E(4,iE)
     is=C4E(5,iE); js=proxs(is,ip)
     nb=  nprox(ip)
     ism= mod(is+nb-2, nb)+1
     isp= mod(is, nb)+1
     ks=  proxs(isp,ip)
     nb2= nprox(jp)
     jsm= mod(js+nb2-2, nb2) + 1
     jsp= mod(js, nb2)+1
!!     nb2= nprox(Lp)
!!     Ls=  proxs(ism,ip)
!!     Lsm= mod(Ls+nb2-2, nb2) + 1
!!
     do is2=1,nb2
        jpk=prox(is2,jp)
        jpL=prox(is2,jp)
        if (jpk.EQ.kp) then
           jsk=is2
        endif
        if (jpL.EQ.Lp) then
           jsL=is2
        endif
     enddo


!     write(6, 121)  'iE, is, ip, js, jp', iE, is, ip, js, jp
!     write(6, 121)  'is, ism, isp', is, ism, isp
!     write(6, 121)  'js, jsm, jsp, jsk, jsL', js, jsm, jsp, jsk, jsL


     !
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
     nb2= nprox(jp)
     E=0.d0
     do js2=1, nb2
        iE2= E4C(js2,jp)
        E= E + FV(iE2,2) * FV(iE2,2) * areaE(iE2)
     enddo
     E2= E/areaI(jp)                ! test grad(K)
     E4= g*(FV(jp,1) + hb_a(jp))
     ER= E/areaI(jp) + g*(FV(jp,1) + hb_a(jp))  !  E = K_i + gh_i  MPAS def. areaE/area
     gradEn = - (EL - ER) / dcenter(iE) * NoutdUn(is,ip)
     gradKn = - (E1 - E2) / dcenter(iE) * NoutdUn(is,ip)
     gradgh = - (E3 - E4) / dcenter(iE) * NoutdUn(is,ip)
     !
     !-- I give a general name gradEn_FD to refer to all choices
     !
     if (imodify_C .EQ. 0) then     ! TRSK     
        gradEn_FD = gradEn
     elseif (imodify_C .EQ. 1) then     ! (cori + Kin) 
        gradEn_FD = gradgh
     elseif (imodify_C .EQ. 2) then     ! (linear H_corr)
        gradEn_FD = gradEn
     elseif (imodify_C .EQ. 3) then     !  (cori + Kin) + (linear H_corr)
        gradEn_FD = gradgh
     endif


     
     !
     !  Linear correction on gradient force
     !
     !  (linear H_corr) or (cori+Kin)+H_corr
     !
     if (imodify_C .EQ. 2  .OR. imodify_C .EQ. 3) then  
        !
        !new linear correction
        !
        ! line integral for grad(h) at Vertex_up
        !
        !--  Vertex Up
        !
        iE1 = EindexVort(1,iE)    ! iE
        iE2 = EindexVort(2,iE)
        iE3 = EindexVort(3,iE)
        if (imodify_C .EQ. 2) then                         ! linear H_corr
           h1 = 0.5d0 * ( ( FV(ip,1) + FV(jp,1) + hb_a(ip) + hb_a(jp) )*g +  K_center(ip) + K_center(jp))
           h2 = 0.5d0 * ( ( FV(jp,1) + FV(kp,1) + hb_a(jp) + hb_a(kp) )*g +  K_center(jp) + K_center(kp))
           h3 = 0.5d0 * ( ( FV(kp,1) + FV(ip,1) + hb_a(kp) + hb_a(ip) )*g +  K_center(kp) + K_center(ip))
        elseif  (imodify_C .EQ. 3) then                    ! linear H_corr +  (Cori+Kin)
           h1 = 0.5d0 *   ( FV(ip,1) + FV(jp,1) + hb_a(ip) + hb_a(jp) )*g 
           h2 = 0.5d0 *   ( FV(jp,1) + FV(kp,1) + hb_a(jp) + hb_a(kp) )*g 
           h3 = 0.5d0 *   ( FV(kp,1) + FV(ip,1) + hb_a(kp) + hb_a(ip) )*g 
        endif
        !
        ! line integral for grad(h) at Vertex_dn
        !
        do j=1, 3
           vec1(j) = h1 * (- Tvec2(j, is, ip) ) * dcenter (iE1) + &
                h2 * (- Tvec2(j, jsk, jp) ) * dcenter (iE2) + &
                h3 * (  Tvec2(j, isp, ip) ) * dcenter (iE3) 
        enddo
        vec1(1:3) = vec1(1:3) / areaV(1,iE)


!        write(6, 101) 'Tvec2(1:3, is, ip)', Tvec2(1:3, is, ip)
!        write(6, 101) 'h1,h2,h3', h1,h2,h3
!        write(6, 101) 'dcenter(iE1),dcenter(iE2),dcenter(iE3)', dcenter(iE1),dcenter(iE2),dcenter(iE3)
!        write(6, 121) 'jsm, jsk', jsm, jsk
!        write(6, 121) 'is, isp', is, isp
!        write(6, 101) 'vec1(1:3)', vec1(1:3)


        !
        !-- Vertex dn
        !
        iE1 = EindexVort(4,iE)
        iE2 = EindexVort(5,iE)
        iE3 = EindexVort(6,iE)    ! iE
        !
        if (imodify_C .EQ. 2) then                         ! linear H_corr
           h1 = 0.5d0 * ( ( FV(ip,1) + FV(Lp,1) + hb_a(ip) + hb_a(Lp) )*g +  K_center(ip) + K_center(Lp))
           h2 = 0.5d0 * ( ( FV(jp,1) + FV(Lp,1) + hb_a(jp) + hb_a(Lp) )*g +  K_center(jp) + K_center(Lp))
           h3 = 0.5d0 * ( ( FV(jp,1) + FV(ip,1) + hb_a(jp) + hb_a(ip) )*g +  K_center(jp) + K_center(ip))
        elseif  (imodify_C .EQ. 3) then                    ! linear H_corr +  (Cori+Kin)
           h1 = 0.5d0 *   ( FV(ip,1) + FV(Lp,1) + hb_a(ip) + hb_a(Lp) )*g 
           h2 = 0.5d0 *   ( FV(jp,1) + FV(Lp,1) + hb_a(jp) + hb_a(Lp) )*g 
           h3 = 0.5d0 *   ( FV(jp,1) + FV(ip,1) + hb_a(jp) + hb_a(ip) )*g 
        endif
        !
        ! line integral for grad(h) at Vertex_dn
        do j=1, 3
           vec2(j) = h1 * (- Tvec2(j, ism, ip) ) * dcenter (iE1) + &
                h2 * (  Tvec2(j, jsL, jp) ) * dcenter (iE2) + &
                h3 * (  Tvec2(j, is,  ip) ) * dcenter (iE3) 
        enddo
        vec2(1:3) = vec2(1:3) / areaV(2,iE)


!        write(6, 101) 'h1,h2,h3', h1,h2,h3
!        write(6, 101) 'vec2(1:3)', vec2(1:3)
!        write(6, 111) 'areaV(1:2,iE)', areaV(1:2,iE)



        !
        !  use mid point is the best choice
        !
        vec3(1:3) = 0.5d0 * ( vec1(1:3) + vec2(1:3) )
        !
        ! grad(h) \dot Nout  \dot Un
        !
        s=0.d0
        do j=1,3
           s= s + vec3(j) * Nvec2(j,is,ip)
        enddo
        s= s * NoutdUn(is,ip)      ! this is mid point gradient from linear interp
        !
        !-- Linear correction on grad (E) due to gradient from T direction
        !
        s2=0.d0
        do j=1,3
           s2= s2 + ( vec1(j) - vec2(j)) * Nvec2(j,is,ip)     ! grad (Up - Dn) proj to Nvec2
        enddo
        grad_interp = gradEn_FD +  s2 * dx_ov_le(iE) * NoutdUn(is,ip)
        !
        !
!        write(6,101) 'grad_interp, gradEn_FD, diff, percent', &
!             grad_interp, gradEn_FD, grad_interp - gradEn_FD, (grad_interp - gradEn_FD)/gradEn_FD
        !
        !
     endif
     !
     !! if (iE .GT. 30) STOP 'end ck grad_interp'
     !
     !
     nvert = nprox(ip) + nprox(jp) - 2        ! double Benzene
     if (iconsv_scheme .EQ. 1 .OR. iconsv_scheme .EQ. 3) then
        !
        !-- K_i conserve,  qe diagnostic 
        !
        s=0.d0
        t=0.d0
        do is2=1,nvert                           ! pentagon + hexagon =9; otherwise 10
           iE2= Eindex(is2,iE)                   ! 
           ipc= C4E(1,iE2); isc= C4E(5,iE2)        ! ip,is on circle
           !
           !
           if (imodify_C .EQ. 0  .OR. imodify_C .EQ. 2) then   ! TRSK  or  Linear H_corr
              s= s + wtUt(is2,iE) * slen(isc,ipc) * Fe(iE2) * 0.5d0 * ( qe_diag(iE) + qe_diag(iE2) )
           elseif (imodify_C .EQ. 1  .OR. imodify_C .EQ. 3) then !  YGYU:  (Cori + Kinetic) OR (Cori+Kinetic)+H_linear
              s= s + wtUt(is2,iE) * slen(isc,ipc) * ( eta_edge(iE2)*FV(iE2,2) +  gradK_d_Ut(iE2) ) 
           endif
           !
           !
           t= t + wtUt(is2,iE) * slen(isc,ipc) * FV(iE2,2)
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


!        if ( iE .LE. 20 ) then
!           write(6, 101) 'FV(iE,3), Ut, diff', FV(iE,3), Ut, FV(iE,3) - Ut
!        else
!           STOP 'stop check Ut'
!        endif



     elseif (iconsv_scheme .EQ. 2) then
        !
        !-- Potential Estrophy conserve:  \int xi^2 dV
        !
        !   eta^e * Ut^e  = ( f + xi )^e * U_t^e  = Fe^perp * qe
        !---------------------------------------------
        !   Ut      = 1 / de * \sum_e'  wtUt (e, e') * l_e' * Un_e'
        !   Fe^perp = 1 / de * \sum_e'  wtUt (e, e') * l_e' * Un_e' * h_e'
        !   Qe^perp = Fe^perp * q_e
        !---------------------------------------------
        s=0.d0
        do is2=1,nvert                           !  pentagon + hexagon =9; otherwise 10
           iE2= Eindex(is2,iE)
           ipc=C4E(1,iE2); isc=C4E(5,iE2)        ! ip,is on circle       
           s= s+ wtUt(is2,iE)*slen(isc,ipc)*Fe(iE2)
        enddo
        Fe_perp=s/dcenter(iE)
        etaUt = Fe_perp * qe_diag(iE)
     endif


     !
     ! \P Un / \pt - (f+vort)*U_t  + U_n^e \dot grad(E)= 0
     !             as \U_n \dot (k \cross U_n^e) = -1
     !
     if (imodify_C .EQ. 0) then   ! TRSK  
        AFV(iE,2)=  etaUt - gradEn
     elseif (imodify_C .EQ. 1) then   !   Cori + Kinet
        AFV(iE,2)=  etaUt - gradgh
     elseif (imodify_C .EQ. 2) then   !   Linear H_corr
        AFV(iE,2)=  etaUt - grad_interp
     elseif (imodify_C .EQ. 3) then   !   Cori + Kinet  +  Linear H_corr
        AFV(iE,2)=  etaUt - grad_interp
     endif

!!  YGYU  why ???

     if (idamp .GT. 0) then
        AFV(iE,2)=  etaUt - gradEn_FD + diff_term(iE,2)   !  revert efforts for linear corr in grad(h)
     else
        AFV(iE,2)=  etaUt - gradEn_FD                     !  revert efforts for linear corr in grad(h)
     endif
     AFV(iE,3)=  0.d0   ! in case defined
     !
     ! -- previous dated  12/28/2017
     ! AFV(iE,2)= -(fcori(iE)+vort)*Ut-gradEn          !  coriolios force is not included
     !
     !
  enddo       ! iE



  !------------------------------------------------------------
  !   reset acceleration if sw_test_case_1
  !         passive advection only
  !------------------------------------------------------------
  !
  if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
     do iE=ips,ipe     
        AFV(iE,2)=0.d0
     enddo
  endif
  


  deallocate ( h_e, h_e_4pt, Fe )
  deallocate ( qe_diag, qe_hex_center, qe_hex_vertex, qe_hex_edge )
  deallocate ( Un_dot_gradPV, Ut_dot_gradPV )
  deallocate ( eta_edge, K_center, gradK_d_Ut, diff_term)

  return
  include '../swm/myformat.inc'
end subroutine dyn_C_combine_Cori_kinetic





!
!
!--  This code shows
!      (i) the consequence of not following strict Gravity center from FV theory
!
!
subroutine dyn_DuSPR (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_constants, only  :  g, ae, pi, ndim
  use readnamelist,     only  :  iswcase, idamp
  use module_variables, only  :  lon, lat, nprox, prox, proxs, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, wtsph, wt2sph
  use module_variables, only  :  Nvec2, hb_a, areaV_A, dcenter_A, Vdc, V4C, RAiv4H_A
  use module_variables, only  :  omegasw
  use module_control,   only  :  nip, niV, hypervis    ! 
  implicit none
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  integer  :: ip, jp, kp, iV
  integer  :: is, ism, isp, js, jsm, ks, ksp
  integer  :: nb, nb2, is2, iV2
  integer  :: i, j, k

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, E0, Kc(npp)    ! Ener, E0 at ip,  E_corner
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradE(3), gradE_ana(3)      ! Energy gradient, analytic solution
  real(rt) :: vort                        !
  real(rt) :: s, t, z, s1, s2, s3       ! sum
  real(rt) :: div_V, div_hV               ! \div V or hV
  real(rt) :: hi, hj, hk                  !
  real(rt) :: hs1, hs2, hs3, V1(3), V2(3), V3(3)
  !
  real(rt), allocatable :: diff_term(:,:), div_hV_vert(:,:)

  
  !
  !--------------  Meta code ----------------
  !
  !  for each ip 
  !   loop is  --- Vertex   
  !    use Vdc(is,ip) --> VTC: Vertex Triple Count  
  !      RAiv4H_A(is,ip), 
  !      note 3 vertex:  (is,ip), (jsm,jp), (ks,kp)
  !
  !  for each ip
  !    A-grid Velocity 
  !
  !----------------------------------------
  !    note
  !    FV(iV,1)=h^* ,  hb_a is used over mountain
  !    FV(ip,2)=Vx
  !    FV(ip,3)=Vy
  !    FV(ip,4)=Vz
  !    ndimFV=4
  !


!!#define debug_GC
#undef debug_GC

  allocate (diff_term(nip,ndimFV))
  allocate (div_hV_vert(6,nip))

  !--------------------------------------------------------------
  !  pre-calc. hypervisocosity = - eta * \nebala^4 V_x,V_y,V_z, h
  !--------------------------------------------------------------
  do j=2, ndimFV   !  no_h!,vx,vy,vz in Cartesian cord
     if (idamp .LT. 0) then
        do ip=1,nip
           diff_term(ip,j)=0.d0
        enddo
     elseif (idamp .EQ. 2) then
        call cal_Laplace_scalar_A (nip, FV(1,j), diff_term(1,j)) ! 2nd ord term =  \vis * Lap
     elseif (idamp .EQ. 4) then
        call cal_Laplace_square_A (nip, FV(1,j), diff_term(1,j)) ! 4th ord term = -\vis * Lap^2
     else
        STOP  'idamp not implemented'
     endif
  enddo


  !
  !-- s. 1  h(iV) = FV(iV, 1)  
  !-- Triangle due grid
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- s1.  hs, Vs(1:3) for each triangle
     !
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq 
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1  
           jp=prox(is,ip);    js=proxs(is,ip)
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           nb2=  nprox(jp);   jsm= mod(js+nb2-2, nb2) + 1
           iV = V4C(is,ip)      

           !- ip
           s=0.d0
           nb2= nprox(ip)
           do is2=1, nb2
              iV2= V4C(is2,ip)
              s= s + FV(iV2,1)*RAiv4H_A(is2,ip)
           enddo
           hi=s

           !- jp
           s=0.d0
           nb2= nprox(jp)
           do is2=1, nb2
              iV2= V4C(is2,jp)
              s= s + FV(iV2,1)*RAiv4H_A(is2,jp)
           enddo
           hj=s

           !- kp
           s=0.d0
           nb2= nprox(kp)
           do is2=1, nb2
              iV2= V4C(is2,kp)
              s= s + FV(iV2,1)*RAiv4H_A(is2,kp)
           enddo
           hk=s

           hs1= 0.5d0*(hi+hj)
           hs2= 0.5d0*(hj+hk)
           hs3= 0.5d0*(hk+hi)

           do j=1,3      ! vec j-comp 
              V1(j)= 0.5d0*(FV(ip,j+1)+FV(jp,j+1))
              V2(j)= 0.5d0*(FV(jp,j+1)+FV(kp,j+1))
              V3(j)= 0.5d0*(FV(kp,j+1)+FV(ip,j+1))
           enddo

           s1=  V1(1)*Nvec2(1,is,ip) +V1(2)*Nvec2(2,is,ip) +V1(3)*Nvec2(3,is,ip)
           s2=  V2(1)*Nvec2(1,jsm,jp)+V2(2)*Nvec2(2,jsm,jp)+V2(3)*Nvec2(3,jsm,jp)
           s3=-(V3(1)*Nvec2(1,isp,ip)+V3(2)*Nvec2(2,isp,ip)+V3(3)*Nvec2(3,isp,ip))! Nvec2(ks,kp)=-(isp,ip)           
           s=s+ hs1*s1*dcenter_A(is,ip) + &
                hs2*s2*dcenter_A(jsm,jp) + &
                hs3*s3*dcenter_A(isp,ip)


#ifdef debug_GC
           !-- test
           div_hV_vert(is, ip) = s/areaV_A(is,ip)
           div_hV_vert(jsm,jp) = div_hV_vert(is,ip)
           div_hV_vert(ks, kp) = div_hV_vert(is,ip)

!           write(6,101) 'inside debug_GC 1'
#else
           div_hV= s/areaV_A(is,ip)
           AFV(iV,1)= -div_hV         !  acc for h(iV)
#endif
        endif         !  if Vdc(is,ip)=0
     enddo   ! is
  enddo      ! ip



#ifdef debug_GC
  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq 
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           nb2=nprox(kp);     ksp=mod(ks,nb2)+1
           iV = V4C(is,ip)
           s1= div_hV_vert(ism,ip)
           s2= div_hV_vert(ksp,kp)
           s3= div_hV_vert(isp,ip)
           div_hV= s1*wt2sph(1,is,ip)+s2*wt2sph(2,is,ip)+s3*wt2sph(3,is,ip)
           AFV(iV,1)= -div_hV    !  acc for h(iV)
!           write(6, 101) 'AFV(iV,1)', AFV(iV,1)
!           write(6, 101) 'wt2sph(1:3,is,ip)', wt2sph(1:3,is,ip)
           !
        endif  ! if _ Vdc
     enddo     ! is
  enddo        ! ip
!  write(6,101) 'inside debug_GC 2'
#endif



  !
  !-- s. 2  V(ip) = FV(ip, 2:4)  velocity Vx/Vy/Vz  
  !-- hex grid
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- h for mid edge
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        j=V4C(is,ip); k=V4C(ism,ip)    ! iV1= upper corner, iV2= lower 
        Fs(1,is)= 0.5d0*(FV(j,1)+hb_a(j) + FV(k,1)+hb_a(k)) !  average (h^*+hs) on side
     enddo
     !
     !-- velocity from hex center to mid edge
     !
     do is= 1, nb
        isp=mod(is, nb)+1
        jp=prox(is,ip)
        kp=prox(isp,ip)
        !-- weighted h,vx,vy,vz in Cartesian cord for upper corner
        do j=2, ndimFV
           Fc(j,is)=FV(ip,j)*wtsph(1,is,ip)+FV(jp,j)*wtsph(2,is,ip)+FV(kp,j)*wtsph(3,is,ip)  ! h^*,Vx,y,z
        enddo
        !-- interpolate E = 1/2*v*v to prepare for grad(E)
        Kc(is)= &
             (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4)))*wtsph(1,is,ip)+&
             (0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4)))*wtsph(2,is,ip)+&
             (0.5d0*(FV(kp,2)*FV(kp,2)+FV(kp,3)*FV(kp,3)+FV(kp,4)*FV(kp,4)))*wtsph(3,is,ip)
     enddo
     !
     !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        do j=2, ndimFV
           Fs(j,is)= 0.5d0*(Fc(j,ism) + Fc(j,is))   ! h^*,Vx,Vy,Vz on six sides, Fs = component on side
        enddo
        Es(is) = 0.5d0*(Kc(ism)+Kc(is)) + g*Fs(1,is)  ! (KE + PE) on six sides,  Es = energy on side
     enddo


     s=0.d0
     nb2= nprox(ip)
     do is2=1, nb2
        iV2= V4C(is2,ip)
        s= s + (FV(iV2,1)+hb_a(iV2))*RAiv4H_A(is2,ip)
     enddo
     hi=s
     E0=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*hi
     !!hi=0.d0       ! donot correct grad(ht)
     !!E0=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))
     !
     !-- gradient h; gradE = grad[ K + g (h + hb_a) ]
     !      note: remove \Rcvec component at final step for AFV(ip,2:4), not at here
     !
     do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds 
        t=0.d0; z=0.d0
        do is= 1, nb
           t= t + Es(is)  *Nvec(j,is,ip)*slen(is,ip) ! E0 center, E0p: prime, replace E0 by six-side average
           z= z +          Nvec(j,is,ip)*slen(is,ip) !       n ds
        enddo
        gradE(j)= (t - z*E0    )/area(ip)          ! grad( 1/2*V*V + g(h+hs))
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


     !
     !-- s5.  acceleration for velocity
     !----------------------------------------------------------------------
     ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
     !
     !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
     !----------------------------------------------------------------------
     !
     do i=1,3
        j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
        if (idamp .GT. 0) then
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)+diff_term(ip,i+1)
        else
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
           ! wo Coriolis force
           !        AFV(ip,i+1)=-(vort)*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
        endif
     enddo
     !
     !     AFV(ip,1)= dh/dt, AFV(ip,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
     !


     !
     !-- s6. remove radial component for \part V / \part t
     !
     s=0.d0               ! s= (Vector, Rcvec) not Rvec on edge
     do j=1, 3
        s= s+ AFV(ip,j+1)*Rcvec(j,ip)
     enddo
     do j=1, 3
        AFV(ip,j+1)= AFV(ip,j+1)-s*Rcvec(j,ip)    ! Vec = Vec - (grad,R) R
     enddo
     !
     !
     
     
     !
     !-- s7. reset to zero if sw_test_case_1 or advection
     !
     if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
        AFV(ip,2:4)=0.d0
     endif
     !
     !
  enddo  ! ip  outer loop


  deallocate (diff_term, div_hV_vert)
  return
  include '../swm/myformat.inc'
end subroutine dyn_DuSPR




!
!-- Dual grid SPR  version 2:  average h, V by 1/2 at mid-point
!-- conclusion: this code on SPR-dy-grid does not work
!
subroutine dyn_DuSPRv2 (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_constants, only  :  g, ae, pi, ndim
  use readnamelist,     only  :  iswcase, idamp
  use module_variables, only  :  lon, lat, nprox, prox, proxs, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, wtsph
  use module_variables, only  :  Nvec2, hb_a, areaV_A, dcenter_A, Vdc, E4C, RAiv4H_A
  use module_variables, only  :  omegasw
  use module_control,   only  :  nip, niV, hypervis    ! 
  implicit none
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
  !
  real(rt) :: hc, uc, vc                         ! u on sphere at corner
  integer  :: ip, jp, kp, iV
  integer  :: is, ism, isp, js, jsm, jsm2, ks
  integer  :: nb, nb2, is2, iV2
  integer  :: i, j, k

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: E, E0, Kc(npp), K_s(npp)     ! Ener, E0 at ip,  E_corner
  real(rt) :: Es(npp)                     ! Energy on side
  real(rt) :: gradE(3), gradE_ana(3)      ! Energy gradient, analytic solution
  real(rt) :: vort                        !
  real(rt) :: s, t, z, s1, s2, s3       ! sum
  real(rt) :: div_V, div_hV               ! \div V or hV
  real(rt) :: hi, hj, hk                  !
  real(rt) :: hs1, hs2, hs3, V1(3), V2(3), V3(3)
  !
  real(rt), allocatable :: diff_term(:,:)

  
  !
  !--------------  Meta code ----------------
  !
  !  for each ip 
  !   loop is  --- Vertex   
  !    use Vdc(is,ip) --> VTC: Vertex Triple Count  
  !      RAiv4H_A(is,ip), 
  !      note 3 vertex:  (is,ip), (jsm,jp), (ks,kp)
  !
  !  for each ip
  !    A-grid Velocity 
  !
  !----------------------------------------
  !    note
  !    FV(iV,1)=h^* ,  hb_a is used over mountain
  !    FV(ip,2)=Vx
  !    FV(ip,3)=Vy
  !    FV(ip,4)=Vz
  !    ndimFV=4
  !




  allocate (diff_term(nip,ndimFV))
  !--------------------------------------------------------------
  !  pre-calc. hypervisocosity = - eta * \nebala^4 V_x,V_y,V_z, h
  !--------------------------------------------------------------
  do j=2, ndimFV   !  no_h!,vx,vy,vz in Cartesian cord
     if (idamp .LT. 0) then
        do ip=1,nip
           diff_term(ip,j)=0.d0
        enddo
     elseif (idamp .EQ. 2) then
        call cal_Laplace_scalar_A (nip, FV(1,j), diff_term(1,j)) ! 2nd ord term =  \vis * Lap
     elseif (idamp .EQ. 4) then
        call cal_Laplace_square_A (nip, FV(1,j), diff_term(1,j)) ! 4th ord term = -\vis * Lap^2
     else
        STOP  'idamp not implemented'
     endif
  enddo


  !
  !-- s. 1  h(iV) = FV(iV, 1)  
  !-- Triangle due grid
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- s1.  hs, Vs(1:3) for each triangle
     !
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq 
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1  
           jp=prox(is,ip);    js=proxs(is,ip)
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           nb2=  nprox(jp);   jsm= mod(js+nb2-2, nb2) + 1
           iV = E4C(is,ip)      

!
!-- apporx. with one-half
!
           iV2= E4C(ism,ip)
           hs1= 0.5d0*(FV(iV,1)+FV(iV2,1))
           nb2=nprox(jp);jsm2=mod(js+nb2-3,nb2)+1;iV2= E4C(jsm2,jp)
           hs2= 0.5d0*(FV(iV,1)+FV(iV2,1))
           iV2= E4C(isp,ip)
           hs3= 0.5d0*(FV(iV,1)+FV(iV2,1))

           do j=1,3      ! vec j-comp 
              V1(j)= 0.5d0*(FV(ip,j+1)+FV(jp,j+1))
              V2(j)= 0.5d0*(FV(jp,j+1)+FV(kp,j+1))
              V3(j)= 0.5d0*(FV(kp,j+1)+FV(ip,j+1))
           enddo

           s1=  V1(1)*Nvec2(1,is,ip) +V1(2)*Nvec2(2,is,ip) +V1(3)*Nvec2(3,is,ip)
           s2=  V2(1)*Nvec2(1,jsm,jp)+V2(2)*Nvec2(2,jsm,jp)+V2(3)*Nvec2(3,jsm,jp)
           s3=-(V3(1)*Nvec2(1,isp,ip)+V3(2)*Nvec2(2,isp,ip)+V3(3)*Nvec2(3,isp,ip))! Nvec2(ks,kp)=-(isp,ip)
           
           s=s+ hs1*s1*dcenter_A(is,ip) + &
                hs2*s2*dcenter_A(jsm,jp) + &
                hs3*s3*dcenter_A(isp,ip)

           div_hV= s/areaV_A(is,ip)
           AFV(iV,1)= -div_hV         !  acc for h(iV)
        endif
     enddo   ! is
  enddo      ! ip



  !
  !-- s. 2  V(ip) = FV(ip, 2:4)  velocity Vx/Vy/Vz  
  !-- hex grid
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- h for mid edge
     do is= 1, nb
        ism=mod(is+nb-2, nb)+1
        j=E4C(is,ip); k=E4C(ism,ip)    ! iV1= upper corner, iV2= lower
        Fs(1,is)= 0.5d0*(FV(j,1)+hb_a(j) + FV(k,1)+hb_a(k)) !  average (h^*+hs)  up + down
     enddo
     !
     !-- V for mid edge
     do is= 1, nb
        jp=prox(is,ip)        
        do j=2,ndim
           Fs(j,is)= 0.5d0*(FV(ip,j)+FV(jp,j))  !  average V  ip-jp pair
        enddo
     enddo
     !
     !-- Kinetic
     do is= 1, nb
        jp=prox(is,ip)
        K_s(is)= 0.5d0 * &
             (0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4)) + &
              0.5d0*(FV(jp,2)*FV(jp,2)+FV(jp,3)*FV(jp,3)+FV(jp,4)*FV(jp,4)))
     enddo
     !
     !-- E = K + g*ht
     do is= 1, nb
        Es(is) = K_s(is) + g*Fs(1,is)  ! (KE + PE) on six sides,  Es = energy on side
     enddo



     s=0.d0
     nb2= nprox(ip)
     do is2=1, nb2
        iV2= E4C(is2,ip)
        s= s + (FV(iV2,1)+hb_a(iV2))*RAiv4H_A(is2,ip)
     enddo
     hi=s
     E0=0.5d0*(FV(ip,2)*FV(ip,2)+FV(ip,3)*FV(ip,3)+FV(ip,4)*FV(ip,4))+g*hi
     !
     !-- gradient h; gradE = grad[ K + g (h + hb_a) ]
     !      note: remove \Rcvec component at final step for AFV(ip,2:4), not at here
     !
     do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds
        t=0.d0; z=0.d0
        do is= 1, nb
           t= t + Es(is)  *Nvec(j,is,ip)*slen(is,ip) ! E0 center, E0p: prime, replace E0 by six-side average
           z= z +          Nvec(j,is,ip)*slen(is,ip) !       n ds
        enddo
        gradE(j)= (t - z*E0    )/area(ip)          ! grad( 1/2*V*V + g(h+hs))
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


     !
     !-- s5.  acceleration for velocity
     !----------------------------------------------------------------------
     ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
     !
     !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
     !----------------------------------------------------------------------
     !
     do i=1,3
        j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
        if (idamp .GT. 0) then
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)+diff_term(ip,i+1)
        else
           AFV(ip,i+1)=-(vort+fcori(ip))*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
           ! wo Coriolis force
           !        AFV(ip,i+1)=-(vort)*(Rcvec(j,ip)*FV(ip,k+1)-Rcvec(k,ip)*FV(ip,j+1))-gradE(i)
        endif
     enddo
     !
     !     AFV(ip,1)= dh/dt, AFV(ip,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
     !


     !
     !-- s6. remove radial component for \part V / \part t
     !
     s=0.d0               ! s= (Vector, Rcvec) not Rvec on edge
     do j=1, 3
        s= s+ AFV(ip,j+1)*Rcvec(j,ip)
     enddo
     do j=1, 3
        AFV(ip,j+1)= AFV(ip,j+1)-s*Rcvec(j,ip)    ! Vec = Vec - (grad,R) R
     enddo
     !
     !
     
     
     !
     !-- s7. reset to zero if sw_test_case_1 or advection
     !
     if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
        AFV(ip,2:4)=0.d0
     endif
     !
     !
  enddo  ! ip  outer loop


  deallocate (diff_term)
  return
  include '../swm/myformat.inc'
end subroutine dyn_DuSPRv2





subroutine dyn_ZM (ims,ime,ips,ipe,FV,AFV,ndimFV,itstep)
  use kinds, only: rt
  use module_constants, only  :  g, ae, pi
  use readnamelist,     only  :  iswcase, idamp, imethod_ZM
  use module_variables, only  :  lon, lat, nprox, prox, proxs, npp, fcori  ! Coriolis param  f(ip)
  use module_variables, only  :  slen, area, Nvec, Tvec, Rvec, Rcvec, RVvec, wtsph, wt2sph
  use module_variables, only  :  Nvec2, Tvec2, hb_a, areaV_A, dcenter_A, RAiv4H_A
  use module_variables, only  :  Nvec3, Tvec3, slen3, wtTBsph, areaTB
  use module_variables, only  :  omegasw, lambdas, Vdc, Edc, V4C, C4V
  use module_control,   only  :  nip, niV, hypervis    ! 
  implicit none
  !
  integer, intent(in )::  ims, ime, ips, ipe, ndimFV, itstep
  real(rt),intent(in )::  FV (ims:ime, ndimFV)       ! Four Vector
  real(rt),intent(out)::  AFV(ims:ime, ndimFV)       ! Acceleration
  !
  integer  :: ip, jp, kp, Lp, Mp, Np, iV
  integer  :: is, ism, isp, isp2, ispx, js, jsm, jsm2, ks, ksp, Ls
  integer  :: nb, nb2, is2, iV1, iV2, iV3, iVx
  integer  :: i, j, k
  integer  :: ipx, jpx, kpx, is2m, ia, ib

  real(rt) :: Fc (ndimFV, npp)  ! filed at center      : c;  note index changed
  real(rt) :: Fs (ndimFV, npp)  ! field at mid of side : s        index changed
  real(rt) :: Gc(ndimFV,npp), Gs(ndimFV,npp)  !  Gc=  h * V  (2:4) 
  real(rt) :: hv_c (3, npp), hv_s(3, npp)  ! 
  real(rt) :: E, E0, E0_base, Kc(npp), K_s(npp)     ! Ener, E0 at ip,  E_corner
  real(rt) :: Es(npp), Ec(npp)             ! Energy on side
  real(rt) :: gradE(3), gradE_ana(3)      ! Energy gradient, analytic solution
  real(rt) :: gradE_val(3,npp)
  real(rt) :: gradh(3), gradK(3)
  real(rt) :: vort, vort_val(npp)                       !
  real(rt) :: s, t, z, s1, s2, s3, r      ! sum, r=lambdas(is,ip)=ratio
  real(rt) :: div_V, div_hV, div_hV_val(npp)       ! \div V or hV
  real(rt) :: hi, hj, hk
  real(rt) :: h1, h2, h3, h4, h5, h6, htc(3)
  real(rt) :: hc, Phi, hs, hs1, hs2, hs3, V1(3), V2(3), V3(3)
  real(rt) :: Nvecs(3,3), Tvecs(3,3), length(3)
  !
  real(rt), allocatable :: diff_term(:,:), div_hV_vert(:,:)
  integer  :: imethod, jmethod
  integer  :: key(3)
  real(rt) :: hv_value(3)


  !--------------  Meta code ----------------
  !
  !  imethod=1: Tomita
  !  imethod=2: intentionally confuse GC with circumcenter,  use  X = (x1 + lam*x2) / (1+lam) 
  !
  !----------------------------------------
  !    note
  !    FV(iV,1)=h^* ,  hb_a is used over mountain
  !    FV(ip,2)=Vx
  !    FV(ip,3)=Vy
  !    FV(ip,4)=Vz
  !    ndimFV=4


  imethod= imethod_ZM    !  default: 1= Tomita, TC2 fails within 12 days
                         !  Try 2:  Triangle;  should put GC at Triangular grid


  allocate (diff_term(nip,ndimFV))
  allocate (div_hV_vert(6,nip))
  !--------------------------------------------------------------
  !  pre-calc. hypervisocosity = - eta * \nebala^4 V_x,V_y,V_z, h
  !--------------------------------------------------------------
  do j=2, ndimFV   !  no_h!,vx,vy,vz in Cartesian cord
     if (idamp .LT. 0) then
        do ip=1,nip
           diff_term(ip,j)=0.d0
        enddo
     elseif (idamp .EQ. 2) then
        call cal_Laplace_scalar_A (nip, FV(1,j), diff_term(1,j)) ! 2nd ord term =  \vis * Lap
     elseif (idamp .EQ. 4) then
        call cal_Laplace_square_A (nip, FV(1,j), diff_term(1,j)) ! 4th ord term = -\vis * Lap^2
     else
        STOP  'idamp not implemented'
     endif
  enddo




  IF (imethod .EQ. 1) THEN
  !
  !-- part 1   DIV( hV ) _ ip
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- s1.  div(hV)
     !
     do is= 1, nb                       ! each vertex
        isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
        jp=prox(is,ip);    js=proxs(is,ip)
        kp=prox(isp,ip);   ks=proxs(isp,ip)
        j=1
        hc= FV(ip,j)*wtsph(1,is,ip)+FV(jp,j)*wtsph(2,is,ip)+FV(kp,j)*wtsph(3,is,ip)  ! h^*
        iV = V4C(is,ip)
        do j=1, 3
           hv_c(j,is)=  hc * FV(iV,j+1)
        enddo
     enddo
     !
     do is= 1, nb                       ! each vertex
        ism=mod(is+nb-2, nb)+1
        do j=1, 3
           hv_s(j,is)= 0.5d0 * ( hv_c(j,is) + hv_c(j,ism) )
        enddo
     enddo
     !
     !
     div_hV=0.d0
     do is= 1, nb
        s=0.d0
        do j=1,3    !  s= (hV dot N)
          s= s +  hv_s(j,is) *Nvec(j,is,ip)      ! projection to N
        enddo
        div_hV=div_hV + s*slen(is,ip)  ! (hV dot N)  ds
     enddo
     div_hV= div_hV/area(ip)
     AFV(ip,1)= -div_hV
!     write(6, 143) 'ip, AFV(ip,1)', ip, AFV(ip,1)
!     stop 'ck 1'
  enddo ! ip


  !
  !-- part 2:   dyn for velocity
  !
  !
  !  for each (is,ip)
  !  Loop three hex (ipx,jpx,kpx) 
  !  corner (Fc) -->  edge (Fs) --> FV alg.  
  !
  do ip=1, nip
     nb= nprox(ip)
     !
     !-- s1.  hs, Vs(1:3) for each triangle
     !
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq 
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
           jp=prox(is,ip);    js=proxs(is,ip)
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           nb2=  nprox(jp);   jsm= mod(js+nb2-2, nb2) + 1
           iV = V4C(is,ip)           
           !
           ! find vort,  grad(K+gh) for ip, jp, kp
           !
           key(1)=ip; key(2)=jp; key(3)=kp
           !          Vx(1:3,1)=FV(ip,2:4)
           !          Vx(1:3,2)=FV(ip,2:4)
           !          Vx(1:3,3)=FV(ip,2:4)
           !
           !  find hv_value(1:3)
           !
           z=0.d0
           do ia=1, 3                          ! loop 3-points
              ipx=key(ia)
              nb2=nprox(ipx)
              s=0.d0
              do is2=1, nb2
                 jpx=prox(is2,ipx)
                 ispx=mod(is2, nb2)+1
                 kpx=prox(ispx,ipx)
                 is2m=mod(is2+nb2-2, nb2)+1    ! is2 minus
                 iVx= V4C(is2,ipx)
                 !
                 Phi= g*(FV(ipx,1)+hb_a(ipx))*wtsph(1,is2,ipx)+&
                      g*(FV(jpx,1)+hb_a(jpx))*wtsph(2,is2,ipx)+&
                      g*(FV(kpx,1)+hb_a(kpx))*wtsph(3,is2,ipx)
                 Ec(is2)= 0.5d0*(FV(iVx,2)*FV(iVx,2)+FV(iVx,3)*FV(iVx,3)+FV(iVx,4)*FV(iVx,4)) + Phi
                 j=1
                 hc= FV(ipx,j)*wtsph(1,is2,ipx)+FV(jpx,j)*wtsph(2,is2,ipx)+FV(kpx,j)*wtsph(3,is2,ipx)  ! h^*
                 Fc(1,is2)= hc
                 do j=2,ndimFV
                    Fc(j,is2)=FV(iVx,j)             ! vel
                 enddo
              enddo
              !
              !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
              do is2=1, nb2
                 is2m=mod(is2+nb2-2, nb2)+1    ! is2 minus
                 do j=1, ndimFV
                    Fs(j,is2)= 0.5d0*(Fc(j,is2m) + Fc(j,is2))   ! h^*,Vx,Vy,Vz on six sides, Fs = component on side
                 enddo
                 Es(  is2)= 0.5d0*(Ec(  is2m) + Ec(  is2))   ! (KE + PE) on six sides,  Es = energy on side
              enddo
              !
              E0= 0.d0
              do is2=1,nb2
                 iVx= V4C(is2,ipx)
                 s=0.d0
                 do j=2,ndimFV
                    s= s + 0.5d0*FV(iVx,j)*FV(iVx,j)
                 enddo
                 E0= E0 + s*RAiv4H_A(is2,ipx)
              enddo
!!              E0_base= E0 + g*(FV(ipx,1)+ hb_a(ipx))
              E0_base= g*(FV(ipx,1)+ hb_a(ipx))            ! only gh
              !-- s3.1 div(V)
              !-- s3.2 Vorticity = k dot \curl V
              !
              div_V=0.d0; vort=0.d0
              do is2= 1, nb2
                 s=0.d0; t=0.d0
                 do j=1,3    !  s= (V dot N); t= (V dot T)
                    s= s + Fs(j+1,is2)*Nvec(j,is2,ipx)   ! projection to N
                    t= t + Fs(j+1,is2)*Tvec(j,is2,ipx)   ! proj to T
                 enddo
                 div_V= div_V + s*slen(is2,ipx)
                 vort = vort  + t*slen(is2,ipx)     ! (v_i dot T) * ds_i
              enddo
              div_V= div_V/area(ipx)
              vort_val(ia) = vort /area(ipx)
              !
              !-- s4:  gradE = grad[ K + g (h + hb_a) ]
              !
              do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds
                 t=0.d0; z=0.d0
                 do is2= 1, nb2
                    t= t + Es(is2)  *Nvec(j,is2,ipx)*slen(is2,ipx) ! E0 cent,E0p: prime,replace E0 by 6-side averg
                    z= z +           Nvec(j,is2,ipx)*slen(is2,ipx) !       n ds
                 enddo
                 gradE_val(j,ia)= (t - z*E0_base)/area(ip)          ! grad( 1/2*V*V + g(h+hs))
              enddo
           enddo    !  loop ia
           !
           ! --  weighted average
           !
           vort=0.d0
           do j=1, 3      !  ia
              vort = vort +  vort_val(j) * wtsph(j,is,ip)
           enddo
           !
           ! remove Rcvec component in gradE_val
           !
           do ia=1, 3                          ! loop 3-points
              ipx=key(ia)
              s=0.d0
              do j=1, 3
                 s= s+ gradE_val(j,ia)*Rcvec(j,ipx)
              enddo
              do j=1, 3
                 gradE_val(j,ia)= gradE_val(j,ia) - s*Rcvec(j,ipx)    ! Vec = Vec - (grad,R) R
              enddo
           enddo
           !
           gradE(1:3)=0.d0
           do k=1, 3      !  component
              do j=1, 3   !  ia
                 gradE(k) = gradE(k) + gradE_val(k,j) * wtsph(j,is,ip)
              enddo
           enddo

!           write(6, 201) 'div_V, vort, fcori(ip), gradE(1:3)', &
!                div_V, vort, fcori(iV), gradE(1:3)
!           stop 'ck 1'


           !
           !-- s5.  acceleration for velocity
           !----------------------------------------------------------------------
           ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
           !
           !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
           !----------------------------------------------------------------------
           !
           do i=1,3
              j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
              if (idamp .GT. 0) then
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)+diff_term(iV,i+1)
              else
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
                 ! wo Coriolis force
                 !        AFV(iV,i+1)=-(vort)*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
              endif
           enddo
           !
           !     AFV(iV,1)= dh/dt, AFV(iV,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
           !
                      
           
           !
           !-- s6. remove radial component for \part V / \part t
           !
           s=0.d0               ! s= (Vector, RVvec) not Rvec on edge
           do j=1, 3
              s= s+ AFV(iV,j+1)*RVvec(j,iV)
           enddo
           do j=1, 3
              AFV(iV,j+1)= AFV(iV,j+1)-s*RVvec(j,iV)    ! Vec = Vec - (grad,R) R
           enddo


           !
           !-- s7. reset to zero if sw_test_case_1 or advection
           !
           if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
              AFV(iV,2:4)=0.d0
           endif
           !
           !
        endif    ! Vdc
     enddo       ! is
  enddo          ! ip 


  ELSEIF (imethod.EQ.2) THEN
     !
     !  GC at Triangle, base is on Triangle (ip/jp/kp)
     !  maj:  vel interp from vertex to center
     !  Every quantity on Triangle
     !
  !  Part.1   Velocity
  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq     
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
           jp=prox(is,ip);    js=proxs(is,ip)
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           nb2=nprox(jp);     jsm=mod(js+nb2-2, nb2) + 1
           iV = V4C(is,ip)           
           key(1)=ip; key(2)=jp; key(3)=kp
           !
           ! interp  velocity
           !
           z=0.d0
           do ia=1, 3                          ! loop 3-points
              ipx=key(ia)
              nb2=nprox(ipx)              
              Fc(1,ia)= FV(ipx,1)
              t=0.d0
              do j=2,ndimFV                    ! vel
                 s=0.d0
                 do is2=1, nb2
                    iV2= V4C(is2,ipx)
                    s= s + FV(iV2,j)*RAiv4H_A(is2,ipx)
                 enddo
                 Fc(j,ia)= s           ! velocity field on corner
              enddo
              !
              ! Kinetic
              t=0.d0
              do is2=1, nb2
                 iV2= V4C(is2,ipx)
                 s=0.d0
                 do j=2,ndimFV
                    s= s + 0.5d0*FV(iV2,j)*FV(iV2,j)
                 enddo
                 t= t + s * RAiv4H_A(is2,ipx)
              enddo
              Kc(ia)= t
              Ec(ia)= t + g*(FV(ipx,1)+hb_a(ipx))
           enddo        !  ia three corners
           !
           !  3 -sides
           do ia=1, 3           ! loop 3-points
              ib=mod(ia,3)+1    ! ia plus 1
              do j=1,ndimFV
                 Fs(j,ia)= 0.5d0*(Fc(j,ia)+Fc(j,ib))     ! h^*,Vx/Vy/Vz
              enddo
              Es(ia)= 0.5d0*(Ec(ia)+Ec(ib))
           enddo
           E0=0.5d0*(FV(iV,2)*FV(iV,2)+FV(iV,3)*FV(iV,3)+FV(iV,4)*FV(iV,4))
           s= (FV(ip,1)+hb_a(ip))*wtsph(1,is,ip)+&
              (FV(jp,1)+hb_a(jp))*wtsph(2,is,ip)+&
              (FV(kp,1)+hb_a(kp))*wtsph(3,is,ip)
           E0_base= E0 + s
           !
           ! vort,  gradE
           ! V \dot T  dl 
           !
           Tvecs(1:3,1)=  Tvec2(1:3,is, ip)
           Tvecs(1:3,2)=  Tvec2(1:3,jsm,jp)
           Tvecs(1:3,3)= -Tvec2(1:3,isp,ip)
           !
           Nvecs(1:3,1)=  Nvec2(1:3,is, ip)
           Nvecs(1:3,2)=  Nvec2(1:3,jsm,jp)
           Nvecs(1:3,3)= -Nvec2(1:3,isp,ip)
           !
           length(1)= dcenter_A(is, ip)
           length(2)= dcenter_A(jsm,jp)
           length(3)= dcenter_A(isp,ip)


!           write(6,101) 'Tvecs(1:3,1:3)', Tvecs(1:3,1:3)
!           write(6,101) 'Nvecs(1:3,1:3)', Nvecs(1:3,1:3)
!           write(6,101) 'length(1:3)   ', length(1:3)
!           stop 'ck 1'

           !
           div_V=0.d0; vort=0.d0
           do ia=1, 3     ! 3 sides
              s=0.d0; t=0.d0
              do j=1,3
                 s= s + Fs(j+1,ia)*Nvecs(j,ia) !  V \dot N
                 t= t + Fs(j+1,ia)*Tvecs(j,ia) !  V \dot T
              enddo
              div_V = div_V  + s*length(ia)     ! (v_i dot N) * ds_i
              vort  = vort   + t*length(ia)     ! (v_i dot T) * ds_i
           enddo
           div_V= div_V/areaV_A(is,ip)
           vort=  vort /areaV_A(is,ip)
           !
           ! grad
           !
           do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds 
              t=0.d0; z=0.d0
              do ia= 1, 3
                 t= t + Es(ia) * Nvecs(j,ia)*length(ia)
                 z= z +          Nvecs(j,ia)*length(ia)  ! N ds
              enddo
              gradE(j)= (t - z*E0_base)/area(ip)         ! grad( 1/2*V*V + g(h+hs))
           enddo
           !

!           write(6, 201) 'div_V, vort, fcori(ip), gradE(1:3)', &
!                div_V, vort, fcori(iV), gradE(1:3)
!           stop 'ck 1'


           !-- s5.  acceleration for velocity
           !----------------------------------------------------------------------
           ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
           !
           !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
           !----------------------------------------------------------------------
           !
           do i=1,3
              j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
              if (idamp .GT. 0) then
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)+diff_term(iV,i+1)
              else
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
                 ! wo Coriolis force
                 !        AFV(iV,i+1)=-(vort)*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
              endif
           enddo
           !
           !     AFV(iV,1)= dh/dt, AFV(iV,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
                      
           
           !-- s6. remove radial component for \part V / \part t
           !
           s=0.d0               ! s= (Vector, RVvec) not Rvec on edge
           do j=1, 3
              s= s+ AFV(iV,j+1)*RVvec(j,iV)
           enddo
           do j=1, 3
              AFV(iV,j+1)= AFV(iV,j+1)-s*RVvec(j,iV)    ! Vec = Vec - (grad,R) R
           enddo


           !
           !-- s7. reset to zero if sw_test_case_1 or advection
           !
           if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
              AFV(iV,2:4)=0.d0
           endif
           !
           !
        endif    ! Vdc
     enddo       ! is
  enddo          ! ip 



  !-- part 2.  div(hV)
  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb                       ! each vertex
        isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
        jp=prox(is,ip);    js=proxs(is,ip)
        kp=prox(isp,ip);   ks=proxs(isp,ip)
        nb2=nprox(jp);     jsm=mod(js+nb2-2, nb2) + 1
        iV = V4C(is,ip)           
        key(1)=ip; key(2)=jp; key(3)=kp
        !
        ! interp  velocity
        !
        z=0.d0
        do ia=1, 3                          ! loop 3-points
           ipx=key(ia)
           nb2=nprox(ipx)              
           Fc(1,ia)= FV(ipx,1)
           t=0.d0
           do j=2,ndimFV                    ! vel
              s=0.d0
              do is2=1, nb2
                 iV2= V4C(is2,ipx)
                 s= s + FV(iV2,j)*RAiv4H_A(is2,ipx)
              enddo
              Fc(j,ia)= s           ! velocity field on corner
              Gc(j,ia)= Fc(1,ia)*s  ! hV    on corner   Gc(1,ia) undefined
              t= t + 0.5d0*s*s               ! Kc
           enddo
        enddo        !  ia three corners
        !
        !  3 -sides
        do ia=1, 3           ! loop 3-points
           ib=mod(ia,3)+1    ! ia plus 1
           do j=1,ndimFV
              Gs(j,ia)= 0.5d0*(Gc(j,ia)+Gc(j,ib))     !     (h*V_j)
           enddo
        enddo

        !
        ! vort,  gradE
        ! V \dot T  dl 
        !
        Tvecs(1:3,1)=  Tvec2(1:3,is, ip)
        Tvecs(1:3,2)=  Tvec2(1:3,jsm,jp)
        Tvecs(1:3,3)= -Tvec2(1:3,isp,ip)
        !
        Nvecs(1:3,1)=  Nvec2(1:3,is, ip)
        Nvecs(1:3,2)=  Nvec2(1:3,jsm,jp)
        Nvecs(1:3,3)= -Nvec2(1:3,isp,ip)
        !
        length(1)= dcenter_A(is, ip)
        length(2)= dcenter_A(jsm,jp)
        length(3)= dcenter_A(isp,ip)
        
        
!       write(6,101) 'Tvecs(1:3,1:3)', Tvecs(1:3,1:3)
!       write(6,101) 'Nvecs(1:3,1:3)', Nvecs(1:3,1:3)
!       write(6,101) 'length(1:3)   ', length(1:3)
!       stop 'ck 1'
        !
        !
        div_hV=0.d0
        do ia=1, 3     ! 3 sides
           r=0.d0
           do j=1,3
              r= r + Gs(j+1,ia)*Nvecs(j,ia) !  hV \dot N
           enddo
           div_hV= div_hV + r*length(ia)     ! (vh_i dot N) * ds_i
        enddo
        div_hV_val(is)= div_hV/areaV_A(is,ip)
     enddo   !  is 
     !
     s=0.d0
     do is=1, nb
        s= s +  div_hV_val(is) * RAiv4H_A(is,ip)
     enddo
     AFV(ip,1)= -s
  enddo       ! ip
        

  ELSEIF (imethod.EQ.3) THEN
     !
     !  Mix gradient ( Hex + Tri-A )
     !  GC at Triangle, base is on Triangle (ip/jp/kp)
     !  grad( gh) on Tri-A
     !  grad ( 1/2 V*V),   \xi on hexagon
     !

  !  Part.1   grad h
  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq     
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
           jp=prox(is,ip);    js=proxs(is,ip)
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           nb2=nprox(jp);     jsm=mod(js+nb2-2, nb2) + 1
           iV = V4C(is,ip)           
           key(1)=ip; key(2)=jp; key(3)=kp
           !
           !--  
           z=0.d0
           do ia=1, 3                          ! loop 3-points
              ipx=key(ia)
              Ec(ia)= g*(FV(ipx,1)+hb_a(ipx))
           enddo        !  ia three corners
           !
           !  3 -sides
           do ia=1, 3           ! loop 3-points
              ib=mod(ia,3)+1    ! ia plus 1
              Es(ia)= 0.5d0*(Ec(ia)+Ec(ib))
           enddo
           s= (FV(ip,1)+hb_a(ip))*wtsph(1,is,ip)+&
              (FV(jp,1)+hb_a(jp))*wtsph(2,is,ip)+&
              (FV(kp,1)+hb_a(kp))*wtsph(3,is,ip)
           E0_base= s
           !
           ! grad(h)
           !
           Nvecs(1:3,1)=  Nvec2(1:3,is, ip)
           Nvecs(1:3,2)=  Nvec2(1:3,jsm,jp)
           Nvecs(1:3,3)= -Nvec2(1:3,isp,ip)
           !
           length(1)= dcenter_A(is, ip)
           length(2)= dcenter_A(jsm,jp)
           length(3)= dcenter_A(isp,ip)
           !
           do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds 
              t=0.d0; z=0.d0
              do ia= 1, 3
                 t= t + Es(ia) * Nvecs(j,ia)*length(ia)
                 z= z +          Nvecs(j,ia)*length(ia)  ! N ds
              enddo
              gradh(j)= (t - z*E0_base)/area(ip)         ! grad( g(h+hs))  only 
           enddo
           !



           !
           !-- grad( K ),   Vort :  interpolate from Hex
           !
           key(1)=ip; key(2)=jp; key(3)=kp
           !          Vx(1:3,1)=FV(ip,2:4)
           !          Vx(1:3,2)=FV(ip,2:4)
           !          Vx(1:3,3)=FV(ip,2:4)
           !
           !  find hv_value(1:3)
           !
           z=0.d0
           do ia=1, 3                          ! loop 3-points
              ipx=key(ia)
              nb2=nprox(ipx)
              s=0.d0
              do is2=1, nb2
                 jpx=prox(is2,ipx)
                 ispx=mod(is2, nb2)+1
                 kpx=prox(ispx,ipx)
                 is2m=mod(is2+nb2-2, nb2)+1    ! is2 minus
                 iVx= V4C(is2,ipx)
                 !
                 Ec(is2)= 0.5d0*(FV(iVx,2)*FV(iVx,2)+FV(iVx,3)*FV(iVx,3)+FV(iVx,4)*FV(iVx,4))  ! K only
                 j=1
                 hc= FV(ipx,j)*wtsph(1,is2,ipx)+FV(jpx,j)*wtsph(2,is2,ipx)+FV(kpx,j)*wtsph(3,is2,ipx)  ! h^*
                 Fc(1,is2)= hc
                 do j=2,ndimFV
                    Fc(j,is2)=FV(iVx,j)             ! vel
                 enddo
              enddo
              !
              !--  One-half rule for mid of edge: 1/2 * (X_ims + X_im)
              do is2=1, nb2
                 is2m=mod(is2+nb2-2, nb2)+1    ! is2 minus
                 do j=1, ndimFV
                    Fs(j,is2)= 0.5d0*(Fc(j,is2m) + Fc(j,is2))   ! h^*,Vx,Vy,Vz on six sides, Fs = component on side
                 enddo
                 Es(  is2)= 0.5d0*(Ec(  is2m) + Ec(  is2))   ! (KE + PE) on six sides,  Es = energy on side
              enddo
              !
              E0= 0.d0
              do is2=1,nb2
                 iVx= V4C(is2,ipx)
                 s=0.d0
                 do j=2,ndimFV
                    s= s + 0.5d0*FV(iVx,j)*FV(iVx,j)
                 enddo
                 E0= E0 + s*RAiv4H_A(is2,ipx)
              enddo
              E0_base= E0 
!!              E0_base= E0 + g*(FV(ipx,1)+ hb_a(ipx))            ! only gh
              !-- s3.1 div(V)
              !-- s3.2 Vorticity = k dot \curl V
              !
              div_V=0.d0; vort=0.d0
              do is2= 1, nb2
                 s=0.d0; t=0.d0
                 do j=1,3    !  s= (V dot N); t= (V dot T)
                    s= s + Fs(j+1,is2)*Nvec(j,is2,ipx)   ! projection to N
                    t= t + Fs(j+1,is2)*Tvec(j,is2,ipx)   ! proj to T
                 enddo
                 div_V= div_V + s*slen(is2,ipx)
                 vort = vort  + t*slen(is2,ipx)     ! (v_i dot T) * ds_i
              enddo
              div_V= div_V/area(ipx)
              vort_val(ia) = vort /area(ipx)
              !
              !-- s4:  gradE = grad[ K + g (h + hb_a) ]
              !
              do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds
                 t=0.d0; z=0.d0
                 do is2= 1, nb2
                    t= t + Es(is2)  *Nvec(j,is2,ipx)*slen(is2,ipx) ! E0 cent,E0p: prime,replace E0 by 6-side averg
                    z= z +           Nvec(j,is2,ipx)*slen(is2,ipx) !       n ds
                 enddo
                 gradE_val(j,ia)= (t - z*E0_base)/area(ip)          ! grad( 1/2*V*V + g(h+hs))
              enddo
           enddo    !  loop ia
           !
           ! --  weighted average
           !
           vort=0.d0
           do j=1, 3      !  ia
              vort = vort +  vort_val(j) * wtsph(j,is,ip)
           enddo
           !
           ! remove Rcvec component in gradE_val
           !
           do ia=1, 3                          ! loop 3-points
              ipx=key(ia)
              s=0.d0
              do j=1, 3
                 s= s+ gradE_val(j,ia)*Rcvec(j,ipx)
              enddo
              do j=1, 3
                 gradE_val(j,ia)= gradE_val(j,ia) - s*Rcvec(j,ipx)    ! Vec = Vec - (grad,R) R
              enddo
           enddo
           !
           gradK(1:3)=0.d0
           do k=1, 3      !  component
              do j=1, 3   !  ia
                 gradK(k) = gradK(k) + gradE_val(k,j) * wtsph(j,is,ip)
              enddo
           enddo


           !
           !-- reassemble
           !
           do j=1,3
              gradE(j)= gradh(j) + gradK(j)
           enddo



!           write(6, 201) 'div_V, vort, fcori(ip), gradE(1:3)', &
!                div_V, vort, fcori(iV), gradE(1:3)
!           stop 'ck 1'


           !-- s5.  acceleration for velocity
           !----------------------------------------------------------------------
           ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
           !
           !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
           !----------------------------------------------------------------------
           !
           do i=1,3
              j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
              if (idamp .GT. 0) then
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)+diff_term(iV,i+1)
              else
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
                 ! wo Coriolis force
                 !        AFV(iV,i+1)=-(vort)*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
              endif
           enddo
           !
           !     AFV(iV,1)= dh/dt, AFV(iV,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
                      
           
           !-- s6. remove radial component for \part V / \part t
           !
           s=0.d0               ! s= (Vector, RVvec) not Rvec on edge
           do j=1, 3
              s= s+ AFV(iV,j+1)*RVvec(j,iV)
           enddo
           do j=1, 3
              AFV(iV,j+1)= AFV(iV,j+1)-s*RVvec(j,iV)    ! Vec = Vec - (grad,R) R
           enddo


           !
           !-- s7. reset to zero if sw_test_case_1 or advection
           !
           if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
              AFV(iV,2:4)=0.d0
           endif
           !
           !
        endif    ! Vdc
     enddo       ! is
  enddo          ! ip 



  !-- part 2.  div(hV)
  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb                       ! each vertex
        isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
        jp=prox(is,ip);    js=proxs(is,ip)
        kp=prox(isp,ip);   ks=proxs(isp,ip)
        nb2=nprox(jp);     jsm=mod(js+nb2-2, nb2) + 1
        iV = V4C(is,ip)           
        key(1)=ip; key(2)=jp; key(3)=kp
        !
        ! interp  velocity
        !
        z=0.d0
        do ia=1, 3                          ! loop 3-points
           ipx=key(ia)
           nb2=nprox(ipx)              
           Fc(1,ia)= FV(ipx,1)
           t=0.d0
           do j=2,ndimFV                    ! vel
              s=0.d0
              do is2=1, nb2
                 iV2= V4C(is2,ipx)
                 s= s + FV(iV2,j)*RAiv4H_A(is2,ipx)
              enddo
              Fc(j,ia)= s           ! velocity field on corner
              Gc(j,ia)= Fc(1,ia)*s  ! hV    on corner   Gc(1,ia) undefined
              t= t + 0.5d0*s*s               ! Kc
           enddo
        enddo        !  ia three corners
        !
        !  3 -sides
        do ia=1, 3           ! loop 3-points
           ib=mod(ia,3)+1    ! ia plus 1
           do j=1,ndimFV
              Gs(j,ia)= 0.5d0*(Gc(j,ia)+Gc(j,ib))     !     (h*V_j)
           enddo
        enddo

        !
        ! vort,  gradE
        ! V \dot T  dl 
        !
        Tvecs(1:3,1)=  Tvec2(1:3,is, ip)
        Tvecs(1:3,2)=  Tvec2(1:3,jsm,jp)
        Tvecs(1:3,3)= -Tvec2(1:3,isp,ip)
        !
        Nvecs(1:3,1)=  Nvec2(1:3,is, ip)
        Nvecs(1:3,2)=  Nvec2(1:3,jsm,jp)
        Nvecs(1:3,3)= -Nvec2(1:3,isp,ip)
        !
        length(1)= dcenter_A(is, ip)
        length(2)= dcenter_A(jsm,jp)
        length(3)= dcenter_A(isp,ip)
        
        
!       write(6,101) 'Tvecs(1:3,1:3)', Tvecs(1:3,1:3)
!       write(6,101) 'Nvecs(1:3,1:3)', Nvecs(1:3,1:3)
!       write(6,101) 'length(1:3)   ', length(1:3)
!       stop 'ck 1'
        !
        !
        div_hV=0.d0
        do ia=1, 3     ! 3 sides
           r=0.d0
           do j=1,3
              r= r + Gs(j+1,ia)*Nvecs(j,ia) !  hV \dot N
           enddo
           div_hV= div_hV + r*length(ia)     ! (vh_i dot N) * ds_i
        enddo
        div_hV_val(is)= div_hV/areaV_A(is,ip)
     enddo   !  is 
     !
     s=0.d0
     do is=1, nb
        s= s +  div_hV_val(is) * RAiv4H_A(is,ip)
     enddo
     AFV(ip,1)= -s
  enddo       ! ip
        


  ELSEIF (imethod.EQ.4) THEN
     !
     !  Vel-dyn at Tri-B
     !
     !
  !  Part.1   Velocity

     write(6,101) 'inside dyn_ZM  imethod==4'

  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb                       ! each vertex
        if ( Vdc(is,ip) .EQ. 0 ) then   ! uniq     
           isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1; isp2=mod(is+1, nb)+1
           jp=prox(is,ip);    js=proxs(is,ip)
           kp=prox(isp,ip);   ks=proxs(isp,ip)
           Lp=prox(ism,ip);   Ls=proxs(ism,ip)
           nb2=nprox(jp)
           jsm=mod(js+nb2-2, nb2)+1 
           jsm2=mod(js+nb2-3,nb2)+1
           Mp=prox(jsm2,jp)
           Np=prox(isp2,ip)

           iV = V4C(is,  ip)
           iV1= V4C(ism, ip)
           iV2= V4C(jsm2,jp)
           iV3= V4C(isp, ip)
           key(1)=ip; key(2)=jp; key(3)=kp
           !
           ! obtain velocity, interp hc
           !
           h1= FV(ip,1); h2= FV(jp,1); h3= FV(kp,1) 
           h4= FV(Lp,1); h5= FV(Mp,1); h6= FV(Np,1)
           Fc(1,1)= &               !  hx(1) pt_1
                h1*wtsph(1,ism,ip)+ &
                h4*wtsph(2,ism,ip)+ &
                h2*wtsph(3,ism,ip)
           Fc(1,2)= &               !  hx(2) pt_2
                h2*wtsph(1,jsm2,jp)+&
                h5*wtsph(2,jsm2,jp)+&
                h3*wtsph(3,jsm2,jp)
           Fc(1,3)= &              !  hx(3) pt_3
                h1*wtsph(1,isp,ip)+&
                h3*wtsph(2,isp,ip)+&
                h6*wtsph(3,isp,ip)

           h1= FV(ip,1)+hb_a(ip); h2= FV(jp,1)+hb_a(jp); h3= FV(kp,1)+hb_a(kp)
           h4= FV(Lp,1)+hb_a(Lp); h5= FV(Mp,1)+hb_a(Mp); h6= FV(Np,1)+hb_a(Np)
           htc(1)= &               !  htx(1) pt_1
                h1*wtsph(1,ism,ip)+ &
                h4*wtsph(2,ism,ip)+ &
                h2*wtsph(3,ism,ip)
           htc(2)= &               !  htx(2) pt_2
                h2*wtsph(1,jsm2,jp)+&
                h5*wtsph(2,jsm2,jp)+&
                h3*wtsph(3,jsm2,jp)
           htc(3)= &              !  htx(3) pt_3
                h1*wtsph(1,isp,ip)+&
                h3*wtsph(2,isp,ip)+&
                h6*wtsph(3,isp,ip)
           do j=2,ndimFV
              Fc(j,1)= FV(iV1,j)
              Fc(j,2)= FV(iV2,j)
              Fc(j,3)= FV(iV3,j)
           enddo
           Ec(1:3)=0.d0
           do ia=1,3     ! 3-pt  loop corner
              do j=1, 3  ! Vx/y/z
                 Ec(ia)= Ec(ia) + 0.5d0*Fc(j+1,ia)*Fc(j+1,ia)
              enddo
              Ec(ia)= Ec(ia) + g*htc(ia)
           enddo
           !
           !  3 -sides
           do ia=1, 3           ! loop 3-points
              ib=mod(ia,3)+1    ! ia plus 1
              do j=1,ndimFV
                 Fs(j,ia)= 0.5d0*(Fc(j,ia)+Fc(j,ib))     ! h^*,Vx/Vy/Vz
              enddo
              Es(ia)= 0.5d0*(Ec(ia)+Ec(ib))
           enddo
           E0=0.5d0*(FV(iV,2)*FV(iV,2)+FV(iV,3)*FV(iV,3)+FV(iV,4)*FV(iV,4))
           ! reuse h1/h2/h3 tot energy
           ! if iV is GC of Triangle-A
           s= h1*wtsph(1,is,ip)+h2*wtsph(2,is,ip)+h3*wtsph(3,is,ip)
!-- if iV is GC of Triangle-B
!           s= htc(1)*wtTBsph(1,iV)+&
!              htc(2)*wtTBsph(2,iV)+&
!              htc(3)*wtTBsph(3,iV)
           E0_base= E0 + s
           !
           ! vort,  gradE
           ! V \dot T  dl 
           !
           div_V=0.d0; vort=0.d0
           do ia=1, 3     ! 3 sides
              s=0.d0; t=0.d0
              do j=1,3
                 s= s + Fs(j+1,ia)*Nvec3(j,ia,iV) !  V \dot N
                 t= t + Fs(j+1,ia)*Tvec3(j,ia,iV) !  V \dot T
              enddo
              div_V = div_V  + s*slen3(ia,iV)     ! (v_i dot N) * ds_i
              vort  = vort   + t*slen3(ia,iV)     ! (v_i dot T) * ds_i
           enddo
           div_V= div_V/areaTB(iV)
           vort=  vort /areaTB(iV)
           !
           ! grad
           !
           do j=1,3    !  loop x, y, z   gradh(j) = 1/A * \sum_is h * N(j,is) ds 
              t=0.d0; z=0.d0
              do ia= 1, 3
                 t= t + Es(ia) * Nvec3(j,ia,iV)*slen3(ia,iV)
                 z= z +          Nvec3(j,ia,iV)*slen3(ia,iV)  ! N ds
              enddo
              gradE(j)= (t - z*E0_base)/areaTB(iV)         ! grad( 1/2*V*V + g(h+hs))
           enddo

           write(6, 201) 'div_V, vort, fcori(ip), gradE(1:3)', &
                div_V, vort, fcori(iV), gradE(1:3)
           stop 'ck 1'


           !-- s5.  acceleration for velocity
           !----------------------------------------------------------------------
           ! \p V / \pt + ( xi[eta] + f) * r_unit \cross V = - grad( E = K+gh )
           !
           !                          (r_unit \cross V)_i = \Epsil_ijk Rj Vk
           !----------------------------------------------------------------------
           !
           do i=1,3
              j=mod(i,3)+1; k=mod(j,3)+1     ! i-j-k permutate
              if (idamp .GT. 0) then
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)+diff_term(iV,i+1)
              else
                 AFV(iV,i+1)=-(vort+fcori(iV))*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
                 ! wo Coriolis force
                 !        AFV(iV,i+1)=-(vort)*(RVvec(j,iV)*FV(iV,k+1)-RVvec(k,iV)*FV(iV,j+1))-gradE(i)
              endif
           enddo
           !
           !     AFV(iV,1)= dh/dt, AFV(iV,2:4)= dVx/dt, dVy/dt, dVz/dt; d denote partial/partial
                      
           
           !-- s6. remove radial component for \part V / \part t
           !
           s=0.d0               ! s= (Vector, RVvec) not Rvec on edge
           do j=1, 3
              s= s+ AFV(iV,j+1)*RVvec(j,iV)
           enddo
           do j=1, 3
              AFV(iV,j+1)= AFV(iV,j+1)-s*RVvec(j,iV)    ! Vec = Vec - (grad,R) R
           enddo


           !
           !-- s7. reset to zero if sw_test_case_1 or advection
           !
           if ( iswcase .EQ. 1  .OR. iswcase .EQ. 9 ) then
              AFV(iV,2:4)=0.d0
           endif
           !
           !
        endif    ! Vdc
     enddo       ! is
  enddo          ! ip 



  !-- part 2.  div(hV)
  do ip=1, nip
     nb= nprox(ip)
     do is= 1, nb                       ! each vertex
        isp=mod(is, nb)+1; ism=mod(is+nb-2, nb)+1
        jp=prox(is,ip);    js=proxs(is,ip)
        kp=prox(isp,ip);   ks=proxs(isp,ip)
        nb2=nprox(jp);     jsm=mod(js+nb2-2, nb2) + 1
        iV = V4C(is,ip)           
        key(1)=ip; key(2)=jp; key(3)=kp
        !
        ! interp  velocity
        !
        z=0.d0
        do ia=1, 3                          ! loop 3-points
           ipx=key(ia)
           nb2=nprox(ipx)              
           Fc(1,ia)= FV(ipx,1)
           t=0.d0
           do j=2,ndimFV                    ! vel
              s=0.d0
              do is2=1, nb2
                 iV2= V4C(is2,ipx)
                 s= s + FV(iV2,j)*RAiv4H_A(is2,ipx)
              enddo
              Fc(j,ia)= s           ! velocity field on corner
              Gc(j,ia)= Fc(1,ia)*s  ! hV    on corner   Gc(1,ia) undefined
              t= t + 0.5d0*s*s               ! Kc
           enddo
        enddo        !  ia three corners
        !
        !  3 -sides
        do ia=1, 3           ! loop 3-points
           ib=mod(ia,3)+1    ! ia plus 1
           do j=1,ndimFV
              Gs(j,ia)= 0.5d0*(Gc(j,ia)+Gc(j,ib))     !     (h*V_j)
           enddo
        enddo

        !
        ! vort,  gradE
        ! V \dot T  dl 
        !
        Tvecs(1:3,1)=  Tvec2(1:3,is, ip)
        Tvecs(1:3,2)=  Tvec2(1:3,jsm,jp)
        Tvecs(1:3,3)= -Tvec2(1:3,isp,ip)
        !
        Nvecs(1:3,1)=  Nvec2(1:3,is, ip)
        Nvecs(1:3,2)=  Nvec2(1:3,jsm,jp)
        Nvecs(1:3,3)= -Nvec2(1:3,isp,ip)
        !
        length(1)= dcenter_A(is, ip)
        length(2)= dcenter_A(jsm,jp)
        length(3)= dcenter_A(isp,ip)
        
        
!       write(6,101) 'Tvecs(1:3,1:3)', Tvecs(1:3,1:3)
!       write(6,101) 'Nvecs(1:3,1:3)', Nvecs(1:3,1:3)
!       write(6,101) 'length(1:3)   ', length(1:3)
!       stop 'ck 1'
        !
        !
        div_hV=0.d0
        do ia=1, 3     ! 3 sides
           r=0.d0
           do j=1,3
              r= r + Gs(j+1,ia)*Nvecs(j,ia) !  hV \dot N
           enddo
           div_hV= div_hV + r*length(ia)     ! (vh_i dot N) * ds_i
        enddo
        div_hV_val(is)= div_hV/areaV_A(is,ip)
     enddo   !  is 
     !
     s=0.d0
     do is=1, nb
        s= s +  div_hV_val(is) * RAiv4H_A(is,ip)
     enddo
     AFV(ip,1)= -s
  enddo       ! ip
        


        


  ENDIF  !  imethod 


  deallocate (diff_term, div_hV_vert)
  return
  include '../swm/myformat.inc'
end subroutine dyn_ZM

