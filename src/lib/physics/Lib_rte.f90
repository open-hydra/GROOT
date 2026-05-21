!> @brief Contains RTE integration subroutines
module GROOT_Lib_DTM
  use, intrinsic :: iso_fortran_env,only: I4 => int32, R8 => real64
  implicit none

  private
  public :: co_flux 
  public :: co_source

  contains

!> @brief Subroutine commputes radiative heat flux on all boundary cells with DTM
subroutine co_flux(grid)
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  use GROOT_Global_m, only: Nr, Ngg, pi, sigma, itermax, tol, imax
  use GROOT_Lib_Geometry, only: getP0fb, distance     !  Library containing grid manipulation procedures
  use GROOT_Lib_Raytracing, only: getray, march        !  Library containing ray tracing procedures
  implicit none
  type(GROOT_domain_type),  intent(inout) :: grid  !< Block-level data [1:Nb]

  real(kind=R8) :: R(1:3,1:3)        ! rotation matrix R from the face frame to the global frame
  real(kind=R8) :: P0(1:3)           ! ray starting point
  real(kind=R8) :: s(1:3)            ! ray vector
  real(kind=R8) :: dphi, dthe        
  real(kind=R8) :: phi, the
  real(kind=R8) :: qnew, q0, err, qtot
  real(kind=R8), allocatable :: Intensity(:) ! intensity at ray end (1:Ngg/Nbands)
  real(kind=R8), allocatable :: flux(:,:,:)  ! impinging heat flux (i1,i2,Ngg/Nbands)
  real(kind=R8), allocatable :: ftot(:)      ! impinging heat flux on single face (1:Ngg/Nbands)
  integer  :: iface,  i, j, k, ib       ! block, cell and face index
  integer  :: iter                      ! DTM iteration
  integer  :: ithe, iphi, Nthe, Nphi
  integer  :: bcm(1:5)                  ! array with block/cell indexes of P0
  integer  :: cnt                       ! counter of lost rays
  logical  :: rayflag                   ! ray error flag
  integer  :: pos(1:imax,1:4)           ! array with indexes of cells traversed by ray
  integer  :: nstep                     ! number of ray steps
  integer  :: face                      ! 
  real(kind=R8) :: ds(1:imax)                ! array with distances of ray steps
  integer  :: isp(1:6)                  ! species pointers

  cnt = 0
  ! initialize angles
  Nphi = 2*Nr**(0.5) 
  Nthe = 0.5*Nr**(0.5)
  dphi = 2*pi/Nphi
  dthe = sin((pi/2)/Nthe)

  ! allocate variables
  allocate(Intensity(0:Ngg), flux(1:Nthe,1:Nphi,0:Ngg),ftot(0:Ngg))
 
  iter = 1
  err = 2e3

  ! if the walls are not black it is necessary to iterate the computation until convergence
  ! the variable used for convergence is the total heat flux
  do while ((err>=tol).and.(iter<=itermax)) 
    ! compute total heat flux
    call co_q0(grid, q0)

    print *, 'iter',iter

    ! for each point on each face of each block
    do ib = 1, grid%nb
      do iface = 1, 6
        do j = 1, grid%blk(ib)%face(iface)%n(2)
          do i = 1, grid%blk(ib)%face(iface)%n(1)
            ! if BC = 101, 200, 300 the cell is not a wall but a connection/periodic/symmetry
            if (grid%blk(ib)%face(iface)%bc(i,j).gt.300) then
              ! initialize impinging flux at 0
              flux(1:Nthe,1:Nphi,0:Ngg) = 0.d0
              ftot(0:Ngg) = 0.d0

              ! compute face center P0, its i-j-k indexes (bcm) 
              ! and the rotation matrix R from the face frame to the global frame
              call getP0fb(grid,ib,iface,i,j,P0,R,bcm)
              !$omp parallel private(face,ithe,iphi,Intensity,rayflag,the,phi,s,nstep,pos,ds)
              !$omp do collapse(2)
              do ithe = 1, Nthe
                do iphi = 1, Nphi

                  ! compute ray direction
                  ! azimuth angle goes from 0 to 2*pi
                  phi = (pi/Nphi) + (iphi - 1)*dphi
                  ! elevation angle goes from 0 to pi/2
                  the = pi/(4*Nthe) + (ithe-1)*dthe

                  call getray(phi, the, R, s)

                  ! for each angle trace the ray until boundary is met, than integrate RTE backwards
              
                  call march(grid, s, P0, bcm, pos, ds, nstep, face, rayflag)
                  
                  call intRTE(2*Nr,s, grid, pos, ds, nstep, face, rayflag, Intensity)

                  if (.not.rayflag) then
                    cnt = cnt+1
                    write(*,"(6(I4,1x))") ib, iface, j, i, ithe, iphi
                  endif
                  ! update impinging flux

                  flux(ithe,iphi,0:Ngg) = Intensity(0:Ngg)*cos(pi/2-the)*sin(pi/2-the)
                enddo
              enddo
              !$omp end do      
              !$omp end parallel  

              ftot(0:Ngg) = 0.d0
              do ithe = 1, Nthe
                do iphi =1, Nphi
                  ftot(0:Ngg) = ftot(0:Ngg)+flux(ithe,iphi,0:Ngg)
                enddo
              enddo

              grid%blk(ib)%face(iface)%Q(i,j,0:Ngg) = ftot(0:Ngg)*dphi*dthe
              grid%blk(ib)%face(iface)%qtot(i,j) = sum(grid%blk(ib)%face(iface)%Q(i,j,0:Ngg))   
            endif 
          enddo
        enddo     
      enddo
   enddo
   
   ! update the total wall heat flux
   call co_q0(grid, qnew)
   err = abs(qnew-q0)
   if (iter==1) write(*,"(a13,i5)")' rays lost = ',cnt
   print *, 'error = ', err
   iter = iter + 1

  enddo
 
  ! compute net total power at wall to check for imbalances (useful if source is computed)
  call co_fluxtot(grid, qtot)
  write(*,"(a19,E15.6,a2)") ' total heat power = ',qtot,' W'

  call co_sourcetot(grid, Qtot)
  write(*,"(a21,E15.6,a2)")' total heat source = ',Qtot,' W'

end subroutine co_flux

!> @brief Integrates heat flux on all faces
subroutine co_q0(grid, q0)
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  implicit none
  type(GROOT_domain_type), intent(in) :: grid      !< Block-level data [1:Nb]
  real(kind=R8),           intent(out):: q0        !< Total heat flux

  integer :: ib, iface, i, j

  q0 = 0.d0
  do ib = 1, grid%nb
    do iface = 1, 6
      do j = 1, grid%blk(ib)%face(iface)%n(2)
        do i = 1, grid%blk(ib)%face(iface)%n(1)
          q0 = q0 + grid%blk(ib)%face(iface)%Qtot(i,j)
        enddo
      enddo
    enddo
  enddo

end subroutine co_q0

!> @brief Subroutine computes divergence of the radiative heat flux (source term) for each grid point
subroutine co_source(grid, Nrs)
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  use GROOT_Lib_Geometry, only: getP0cb    !  Library containing grid manipulation procedures
  use GROOT_Lib_Raytracing, only: getray, march        !  Library containing ray tracing procedures
  use GROOT_Global_m, only: Ngg, sigma, imax, optically_thin, pi
  implicit none
  type(GROOT_domain_type),  intent(inout) :: grid  !< Block-level data [1:Nb]
  integer,                  intent(in)    :: Nrs       !< Number of rays

  real(kind=R8) :: R(1:3,1:3)              ! rotation matrix from cell face frame to global frame
  real(kind=R8) :: P0(1:3)                 ! ray starting point
  real(kind=R8) :: s(1:3)                  ! ray vector
  real(kind=R8) :: ds(1:imax)              ! array with distances travelled by ray
  real(kind=R8) :: div(0:Ngg)              ! divergence of heat flux
  real(kind=R8) :: Gvec(0:Ngg)             ! spectral/graygas irradiance
  real(kind=R8) :: phi, the, dphi, dthe
  real(kind=R8) :: Gtot(0:Ngg)                    ! incident radiation (integrated over solid angle)
  real(kind=R8) :: Intensity(0:Ngg)               ! intensity at ray end
  real(kind=R8) :: Qtot                    ! integral of source term
  real(kind=R8) :: deta                    ! narrow band width
  real(kind=R8), allocatable :: G(:,:,:)     ! impinging intensity (1:Nthe,1:Nphi)
  integer  :: bcm(1:5)                ! array with block/cell indexes of P0
  integer  :: pos(1:imax,1:4)         ! array with indexes of cells traversed by ray
  integer  :: Nphi, Nthe, ithe, iphi
  integer  :: i, j, k, ib
  integer  :: ig                      ! gray gas/narrow band index
  integer  :: cnt                     ! counter of lost rays
  integer  :: nstep                   ! number of ray steps
  integer  :: face
  logical  :: rayflag                 ! ray error flag

  cnt = 0
  ! initialize angles   
  Nphi = 2*(Nrs/2.d0)**(0.5) 
  Nthe = (Nrs/2.d0)**(0.5)
  dphi = 2*pi/Nphi
  dthe = sin(pi/Nthe)

  deta = 25.d0*100.d0

  allocate(G(1:Nthe,1:Nphi,0:Ngg))

  ! for each cell
  do ib = 1, grid%nb
    do k = 1, grid%blk(ib)%dim(3)
      do j = 1, grid%blk(ib)%dim(2)
        do i = 1, grid%blk(ib)%dim(1)
          
          if (optically_thin) then 
            div(0:Ngg) = 0.d0
            do ig=0,Ngg

                div(ig) = grid%blk(ib)%ka(i,j,k,ig)*(4*pi*grid%blk(ib)%a(i,j,k,ig)*grid%blk(ib)%Ib(i,j,k)-0.d0)
              
            enddo
          else 
            ! initialize impinging flux at 0
            G(1:Nthe,1:Nphi,0:Ngg) = 0.d0
            ! get cell center P0 and the rotation matrix R
            ! from the cell frame to the global frame
            call getP0cb(grid, ib, i, j, k, P0, R)
            ! save point indexes
            bcm = (/ib, 0, i , j, k/)
            Gtot(0:Ngg) = 0.d0 
            !$omp parallel private(ithe,iphi,the,phi,s,ig,nstep,pos,ds,face,rayflag,Intensity)
            !$omp do collapse(2)
            do ithe = 1, Nthe
              do iphi = 1, Nphi
                ! azimuth angle is in [0, 2*pi]
                phi = (pi/Nphi) + (iphi - 1)*dphi
                ! elevation angle is in [-pi/2, pi/2]
                the = pi/(2*Nthe) + (ithe-1)*dthe - pi/2
                call getray(phi, the, R, s)

                ! for each angle trace the ray until boundary is met, than integrate RTE backwards
                call march(grid, s, P0, bcm, pos, ds, nstep, face, rayflag)

                call intRTE(Nrs, s, grid, pos, ds, nstep, face, rayflag, Intensity)

                if (.not.rayflag) cnt = cnt+1
                G(ithe,iphi,0:Ngg) = Intensity*sin(pi/2-the)*dphi*dthe !integrated over solid angle
              enddo
            enddo
            !$omp end do      
            !$omp end parallel

            do ithe=1,Nthe 
              do iphi=1,Nphi 
                Gtot(0:Ngg)=Gtot(0:Ngg)+G(ithe,iphi,0:Ngg)
              enddo
            enddo
            div(0:Ngg) = 0.d0
            do ig=0,Ngg
                div(ig) = grid%blk(ib)%ka(i,j,k,ig)*(4*pi*grid%blk(ib)%a(i,j,k,ig)*grid%blk(ib)%Ib(i,j,k)-Gtot(ig))
            enddo

          endif
          grid%blk(ib)%G(i,j,k) = sum(Gtot(0:Ngg)) ! per ora non moltiplico lui per Deta, solo la fine
          grid%blk(ib)%source(i,j,k) = sum(div(0:Ngg))!4*pi*grid%blk(ib)%ib(i,j,k)*coeff-Gk
        enddo
      enddo
    enddo
  enddo

  ! compute heat power 
  call co_sourcetot(grid, Qtot)
  write(*,"(a13,i5)")' rays lost = ',cnt
  write(*,"(a21,E15.6,a2)")' total heat source = ',Qtot,' W'
end subroutine co_source


!> @brief Subroutine computes integral of heat flux on enclosure walls
subroutine co_fluxtot(grid, qtot)
  use GROOT_Lib_Geometry, only : co_area !  Function to compute face area
  use GROOT_Global_m, only: sigma
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  implicit none
  type(GROOT_domain_type),  intent(in) :: grid     !< Block-level data [1:Nb]
  real(kind=R8),            intent(out):: qtot     !< Net heat power

  real(kind=R8):: A, Atot, qrad
  integer :: ib, ifa, i, j
  
  qtot = 0.d0
  Atot = 0.d0
  do ib = 1, grid%nb
    do ifa = 1, 6
      do j = 1, grid%blk(ib)%face(ifa)%n(2)
        do i = 1, grid%blk(ib)%face(ifa)%n(1)
          if (grid%blk(ib)%face(ifa)%bc(i,j).gt.3) then
            A = co_area(face = grid%blk(ib)%face(ifa), i1 = i, i2 = j)
            Atot = Atot+A
            qrad = grid%blk(ib)%face(ifa)%eps(i,j)*grid%blk(ib)%face(ifa)%Qtot(i,j) &
                -sigma*grid%blk(ib)%face(ifa)%eps(i,j)*grid%blk(ib)%face(ifa)%T(i,j)**4
            qtot = qtot + A*qrad
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine co_fluxtot

!> @brief Subroutine computes the integral of the divergence of the radiative heat flux
subroutine co_sourcetot(grid, Qtot)
  use GROOT_Lib_Geometry, only : co_volume !  Function computing cell volume
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  implicit none
  type(GROOT_domain_type),  intent(in) :: grid  !< Block-level data [1:Nb]
  real(kind=R8),            intent(out):: Qtot     !< Total heat power

  real(kind=R8) :: V, Vtot
  integer  :: i, j, k, ib

  Qtot = 0.d0
  Vtot = 0.d0
  do ib = 1, grid%nb
    do k = 1, grid%blk(ib)%dim(3)
      do j = 1, grid%blk(ib)%dim(2)
        do i = 1, grid%blk(ib)%dim(1)
          V = co_volume(blck = grid%blk(ib), i = i, j = j, k = k)
          Vtot = Vtot + V
          Qtot = Qtot + V*grid%blk(ib)%source(i, j, k)
        enddo
      enddo
    enddo
  enddo

end subroutine co_sourcetot

!> @brief Integrates RTE for all gray gases / narrow bands
subroutine intRTE(Nrs,ray, grid, pos, ds, nstep, ifa1, rayflag, Intensity)
  use GROOT_Lib_Geometry, only : co_area, co_normal
  use GROOT_Lib_Math,     only : scalar
  use GROOT_Mod_Face,     only: get_ind
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  use GROOT_Global_m, only: Ngg, sigma, imax, pi
  implicit none
  integer,                  intent(in)    :: Nrs              !< Number of rays
  real(kind=R8),            intent(in)    :: ray(1:3)         !< current ray versor
  type(GROOT_domain_type),  intent(inout) :: grid             !< block level data
  integer,                  intent(in)    :: pos(:,:)         !< indexes of cells traversed by ray
  real(kind=R8),            intent(in)    :: ds(:)            !< distances travelled by ray
  integer,                  intent(in)    :: nstep            !< number of ray steps
  integer,                  intent(in)    :: ifa1             !< index of starting face
  logical,                  intent(in)    :: rayflag          !< ray error flag
  real(kind=R8),            intent(out)   :: Intensity(0:Ngg) !< intensity at ray end

  integer :: ib,ibw, i, j, k    ! cell indexes
  integer :: ind1, ind2     ! face indexes
  integer :: jj             ! ray step counter
  integer :: ig             ! gray gas/narrow band index
  real(kind=R8):: Ivec(1:imax)   ! intensities alogn ray steps
  real(kind=R8):: ka, ks         ! abs and scatter coefficient
  real(kind=R8):: ibb            ! black body intensity
  real(kind=R8):: I0             ! intensity at ray start
  real(kind=R8):: S, om          ! appoggio

  Intensity = 0.d0
  ! get face indexes of final point
  ibw = pos(nstep,1); i = pos(nstep,2); j = pos(nstep,3); k = pos(nstep,4)
  call get_ind(i, j, k, ifa1, ind1, ind2)
  do ig = 0, Ngg
    if (rayflag) then
      ! set initial condition for integration
      associate(f => grid%blk(ibw)%face(ifa1))
        I0 = (f%a(ind1,ind2,ig) * f%eps(ind1,ind2) * sigma * f%T(ind1,ind2)**4 &
           + (1 - f%eps(ind1,ind2)) * f%Q(ind1,ind2,ig)) / pi
      end associate
    else
      I0 = 0.d0
    endif
    ib = pos(nstep,1); i = pos(nstep,2); j = pos(nstep,3); k = pos(nstep,4)
    Ivec(1) = I0
    do jj = 1, nstep
      ib = pos(nstep+1-jj,1); i = pos(nstep+1-jj,2); j = pos(nstep+1-jj,3); k = pos(nstep+1-jj,4)
      ka = grid%blk(ib)%ka(i,j,k,ig)
      ibb = grid%blk(ib)%Ib(i,j,k)

      om = 1/(ka+1.d-40)
      S = grid%blk(ib)%a(i,j,k,ig)*ibb

      Ivec(jj+1) = Ivec(jj)*exp(-ka*ds(nstep+1-jj))+S*(1-exp(-ka*ds(nstep+1-jj))) 
    enddo   
    Intensity(ig) = Ivec(nstep+1)
  enddo 

end subroutine intRTE


end module GROOT_Lib_DTM
