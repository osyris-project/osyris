!=============================================================================
! Modified version of AMR2CELL.f90 from the RAMSES source
!=============================================================================
subroutine ramses_data(infile,lmax2,xcenter,ycenter,zcenter,deltax,deltay,deltaz,lscale,&
                     & data_array,data_names,ncells,ncpu,ndim,levelmin,levelmax,nstep,boxsize,t,ud,ul,ut)

  implicit none
  
  ! Array dimensions
  integer, parameter :: nmax=10000000
  integer, parameter :: nvarmax=50

  ! Subroutine arguments
  character(LEN=*)              , intent(in ) :: infile
  integer                       , intent(in ) :: lmax2
  real*8                        , intent(in ) :: xcenter,ycenter,zcenter,deltax,deltay,deltaz,lscale
  
  integer                       , intent(out) :: ncells,ncpu,ndim,levelmin,levelmax,nstep
  real*8                        , intent(out) :: boxsize,t,ud,ul,ut
  real*8,dimension(nmax,nvarmax), intent(out) :: data_array
  character(len=500)            , intent(out) :: data_names
  
  ! Variables
  integer :: i,j,k,twotondim,ivar,nboundary,ngrid_current,nvar_tot
  integer :: nx,ny,nz,ilevel,idim,icell,ngrp,nlevelmax,lmax,ind,ipos,ngrida
  integer :: ngridmax,icpu,ncpu_read,imin,imax,jmin,jmax,kmin,kmax,nvarh
  integer :: nx_full,ny_full,nz_full,ix,iy,iz,s_start,s_end
  integer, dimension(:  ), allocatable :: cpu_list
  integer, dimension(:,:), allocatable :: son,ngridfile,ngridlevel,ngridbound
  
  real*8 :: xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0,zmin=0.0,zmax=1.0
  real*8 :: xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dx2,gamma,mu,boxlen
!   real*8 :: xi,yi,zi,dxi,rhoi,ui,vi,wi,tempi,bxi,byi,bzi,bi,veli,ul,ud,ut
  real*8, dimension(:,:    ), allocatable :: x,xg
  real*8, dimension(:,:,:  ), allocatable :: var
  real*8, dimension(1:8,1:3)              :: xc
  real*8, dimension(    1:3)              :: xbound=(/0.0d0,0.0d0,0.0d0/)
  
  character(LEN=5  ) :: nchar,ncharcpu
  character(LEN=50 ) :: string
  character(LEN=128) :: nomfich,repository
  
  logical                            :: ok,ok_cell
  logical, dimension(:), allocatable :: ref
  
  type level
     integer :: ilevel
     integer :: ngrid
     integer :: imin
     integer :: imax
     integer :: jmin
     integer :: jmax
     integer :: kmin
     integer :: kmax
  end type level

  type(level),dimension(1:100) :: grid

!=============================================================================

  lmax = lmax2
  repository = trim(infile)

  !-----------------------------------------------
  ! Reading RAMSES data
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read (10) ncpu
  read (10) ndim
  read (10) nx,ny,nz
  read (10) nlevelmax
  read (10) ngridmax
  read (10) nboundary
  read (10) ngrid_current
  read (10) boxlen
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
  
  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  if(ndim==2)then
     write(*,*)'Output file contains 2D data'
     write(*,*)'Aborting'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  open(unit=10,file=nomfich,form='formatted',status='old')
  read (10,'(13x,I11)') ncpu
  read (10,'(13x,I11)') ndim
  read (10,'(13x,I11)') levelmin
  read (10,'(13x,I11)') levelmax
  read (10,*)
  read (10,'(13x,I11)') nstep
  read (10,*)
  
  read (10,'(13x,E23.15)') boxlen
  read (10,'(13x,E23.15)') t
  read (10,*)
  read (10,*)
  read (10,*)
  read (10,*)
  read (10,*)
  read (10,*)
  read (10,'(13x,E23.15)') ul
  read (10,'(13x,E23.15)') ud
  read (10,'(13x,E23.15)') ut
  read (10,'(13x,E23.15)') mu
  read (10,'(13x,I11)') ngrp
  close(10)
  
  ! Open hydro file descriptor
  open(unit=11,file=TRIM(repository)//'/hydro_file_descriptor.txt',form='formatted',status='old')
  read(11,'(13x,I11)') nvarh
  do i = 1,nvarh
     read(11,'(14x,a)') string
     data_names = trim(data_names)//' '//string
  enddo
  close(11)
  ! The total number of variables returned to python will
  ! be nvarh + 5: 3 coordinates, dx and ilevel
  nvar_tot = nvarh+5
  data_names = trim(data_names)//' '//"level"
  data_names = trim(data_names)//' '//"x"
  data_names = trim(data_names)//' '//"y"
  data_names = trim(data_names)//' '//"z"
  data_names = trim(data_names)//' '//"dx"
  
  allocate(cpu_list(1:ncpu))
  
  if(deltax > 0.0d0)then
      xmin = xcenter - 0.5d0*deltax*lscale/(boxlen*ul)
      xmax = xcenter + 0.5d0*deltax*lscale/(boxlen*ul)
  endif
  if(deltax > 0.0d0)then
     ymin = ycenter - 0.5d0*deltay*lscale/(boxlen*ul)
     ymax = ycenter + 0.5d0*deltay*lscale/(boxlen*ul)
  endif
  if(deltax > 0.0d0)then
     zmin = zcenter - 0.5d0*deltaz*lscale/(boxlen*ul)
     zmax = zcenter + 0.5d0*deltaz*lscale/(boxlen*ul)
  endif
  
  !-----------------------
  ! Map parameters
  !-----------------------
  if(lmax==0)then
     lmax=nlevelmax
  endif
  write(*,*)'time=',t
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax
  ncpu_read=ncpu
  do j=1,ncpu
     cpu_list(j)=j
  end do

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1
     kmax=int(zzmax*dble(nz_full))+1
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin
     grid(ilevel)%kmax=kmax
  end do

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------
  icell = 0

  ! Loop over processor files
  do k=1,ncpu_read
     icpu=cpu_list(k)
     write(ncharcpu,'(i5.5)') icpu

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=10,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     do i=1,21
        read(10)
     end do
     ! Read grid numbers
     read(10)ngridlevel
     ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
     read(10)
     if(nboundary>0)then
        do i=1,2
           read(10)
        end do
        read(10)ngridbound
        ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
     endif
     read(10)
     read(10)
     read(10)
     read(10)
     read(10)
     read(10)

     ! Open HYDRO file and skip header
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11) nvarh
     read(11)
     read(11)
     read(11)
     read(11) gamma

     ! Loop over levels
     do ilevel=1,lmax

        ! Geometry
        dx=0.5d0**ilevel
        dx2=0.5d0*dx
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5d0)*dx
           xc(ind,2)=(dble(iy)-0.5d0)*dx
           xc(ind,3)=(dble(iz)-0.5d0)*dx
        end do

        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(xg (1:ngrida,1:ndim))
           allocate(son(1:ngrida,1:twotondim))
           allocate(var(1:ngrida,1:twotondim,1:nvar_tot))
           allocate(x  (1:ngrida,1:ndim))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

           ! Read AMR data
           if(ngridfile(j,ilevel)>0)then
              read(10) ! Skip grid index
              read(10) ! Skip next index
              read(10) ! Skip prev index
              ! Read grid center
              do idim=1,ndim
                 if(j.eq.icpu)then
                    read(10)xg(:,idim)
                 else
                    read(10)
                 endif
              end do
              read(10) ! Skip father index
              do ind=1,2*ndim
                 read(10) ! Skip nbor index
              end do
              ! Read son index
              do ind=1,twotondim
                 if(j.eq.icpu)then
                    read(10)son(:,ind)
                 else
                    read(10)
                 end if
              end do
              ! Skip cpu map
              do ind=1,twotondim
                 read(10)
              end do
              ! Skip refinement map
              do ind=1,twotondim
                 read(10)
              end do
           endif

           ! Read HYDRO data
           read(11)
           read(11)
           if(ngridfile(j,ilevel)>0)then
              ! Read hydro variables
              do ind=1,twotondim
                 do ivar=1,nvarh
                    if(j.eq.icpu)then
                       read(11)var(:,ind,ivar)
                    else
                       read(11)
                    end if
                 end do
              end do
           end if
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

              ! Compute cell center
              do i=1,ngrida
                 x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                 x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                 x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
              end do
              ! Check if cell is refined
              do i=1,ngrida
                 ref(i)=son(i,ind)>0.and.ilevel<lmax
              end do
              ! Store data cube
              do i=1,ngrida
                ok_cell= .not.ref(i).and. &
                     & (x(i,1)+dx2)>=xmin.and.&
                     & (x(i,2)+dx2)>=ymin.and.&
                     & (x(i,3)+dx2)>=zmin.and.&
                     & (x(i,1)-dx2)<=xmax.and.&
                     & (x(i,2)-dx2)<=ymax.and.&
                     & (x(i,3)-dx2)<=zmax

                 if(ok_cell)then
                 
                    icell = icell + 1
                 
                    var(i,ind,nvarh+1) = real(ilevel)
                    var(i,ind,nvarh+2) = x(i,1)*boxlen
                    var(i,ind,nvarh+3) = x(i,2)*boxlen
                    var(i,ind,nvarh+4) = x(i,3)*boxlen
                    var(i,ind,nvarh+5) = dx*boxlen
                 
                    do ivar = 1,nvar_tot                    
                       data_array(icell,ivar) = var(i,ind,ivar)
                    enddo
                 
!                     rhoi = var(i,ind,1)*ud
! 
!                     bi   = sqrt(0.25d0*((var(i,ind,5)+var(i,ind,8))**2 &
!                          &         +(var(i,ind,6)+var(i,ind,9))**2 &
!                          &         +(var(i,ind,7)+var(i,ind,10))**2) &
!                          &         *4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)
! 
!                     bxi = 0.5d0*(var(i,ind,5)+var(i,ind, 8))*sqrt(4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)
!                     byi = 0.5d0*(var(i,ind,6)+var(i,ind, 9))*sqrt(4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)
!                     bzi = 0.5d0*(var(i,ind,7)+var(i,ind,10))*sqrt(4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)
! 
!                     tempi = var(i,ind,17)
! 
!                     xi = x(i,1)*ul*boxlen
!                     yi = x(i,2)*ul*boxlen
!                     zi = x(i,3)*ul*boxlen
!                     dxi = dx*ul*boxlen
! 
!                     ui = var(i,ind,2)*ul/ut
!                     vi = var(i,ind,3)*ul/ut
!                     wi = var(i,ind,4)*ul/ut
!                     
!                     veli = sqrt(ui*ui+vi*vi+wi*wi)
! 
!                     icell = icell + 1
!                     data_array(icell, 1) = real(ilevel)
!                     data_array(icell, 2) = xi
!                     data_array(icell, 3) = yi
!                     data_array(icell, 4) = zi
!                     data_array(icell, 5) = dxi
!                     data_array(icell, 6) = rhoi
!                     data_array(icell, 7) = veli
!                     data_array(icell, 8) = tempi
!                     data_array(icell, 9) = bi
!                     data_array(icell,10) = ui
!                     data_array(icell,11) = vi
!                     data_array(icell,12) = wi
!                     data_array(icell,13) = bxi
!                     data_array(icell,14) = byi
!                     data_array(icell,15) = bzi
                                        
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,var,ref,x)
        endif

     end do
     ! End loop over levels

     close(10)
     close(11)

  end do
  ! End loop over cpus
  
  ncells = icell
  write(*,*) 'Read ',ncells,' cells'
  
  boxsize = boxlen*ul
  
!   write(*,*) data_names

end subroutine ramses_data
