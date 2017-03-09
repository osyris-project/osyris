subroutine ramses_data(infile,xcenter,ycenter,zcenter,lmax2,data_array,ncells)
  !--------------------------------------------------------------------------
  ! Ce programme calcule le cube cartesien pour les
  ! variables hydro d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
!   integer, intent(in):: nmax
  integer, parameter:: nmax=10000000,nvarmax=15
  real(KIND=8),dimension(nmax,nvarmax), intent(out)::data_array
  character(LEN=*), intent(in)::infile
  real(KIND=8), intent(in) :: xcenter,ycenter,zcenter
  integer, intent(in)::lmax2
  integer, intent(out):: ncells
  
  integer::ndim,n,i,j,k,twotondim,ncoarse,type=0,domax=0
  integer::ivar,nvar,ncpu,ncpuh,nboundary,ngrid_current
  integer::nx,ny,nz,ilevel,idim,jdim,kdim,icell
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout,lmax
  integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer::ngridmax,nstep_coarse,icpu,ncpu_read
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  real(KIND=8)::gamma,smallr,smallc,gammah
  real(KIND=8)::boxlen,boxlen2,ul,ud,ut,mu
  real(KIND=8)::t,aexp,hexp,t2,aexp2,hexp2
  real(KIND=8)::omega_m,omega_l,omega_k,omega_b
  real(KIND=8)::scale_l,scale_d,scale_t
  real(KIND=8)::omega_m2,omega_l2,omega_k2,omega_b2,maxrho

  integer::nx_sample=0,ny_sample=0,nz_sample=0,eos=0
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin,levelmin
  integer::ix,iy,iz,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,dmax,dummy
  real(KIND=8), dimension(1) :: order_min
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dx2
  real(KIND=8),dimension(:,:),allocatable::x,xg
  real(KIND=8),dimension(:,:,:),allocatable::var
  real(KIND=8),dimension(:),allocatable::varmax
  real(KIND=4),dimension(:,:,:),allocatable::toto
  real(KIND=8),dimension(:)  ,allocatable::rho
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:)  ,allocatable::isp
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  logical::ok,ok_part,ok_cell
  real(KIND=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  character(LEN=1)::proj='z'

  type level
     integer::ilevel
     integer::ngrid
!      real(KIND=4),dimension(:,:,:),pointer::cube
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level

  type(level),dimension(1:100)::grid

  ! Temporary space for reading labels from the info file.
  character(LEN=128)::temp_label
  
  real(kind=8) :: xi,yi,zi,dxi,rhoi,ui,vi,wi,tempi,bxi,byi,bzi,bi,veli

  lmax = lmax2
    
  repository = trim(infile)
  

  maxrho = 0.0d0
  
!   write(*,*) 'Writing file: '//TRIM(outfich)
!   open(unit=20,file=TRIM(outfich),form='formatted')
!   write(20,*) '# ilevel               x (cm)               y (cm)               z (cm)'//&
!        &'               dx (cm)              rho (g/cm3)          u (cm/s)            '//&
!        &' v (cm/s)             w (cm/s)             T (K)                Bx (G)       '//&
!        &'        By (G)               Bz (G)               B (G)'
!   close(20)
!   open(unit=20,file=TRIM(outfich),form='formatted',access='append')

 !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
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
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)ngrid_current
  read(10)boxlen
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
  read(10,*)
  read(10,*)
  read(10,'(13x,I11)')levelmin
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  
  read(10,'(13x,E23.15)')boxlen
  read(10,'(13x,E23.15)')t
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'(13x,E23.15)')ul
  read(10,'(13x,E23.15)')ud
  read(10,'(13x,E23.15)')ut
  read(10,'(13x,E23.15)')mu
  read(10,*) ! ngrp
  read(10,*)
  read(10,'(14x,A80)') ordering
!   write(*,*) 'ORDERING IS',ordering
  ordering = '1'
!   write(*,'(" ordering type=",A20)')TRIM(ordering)
  read(10,*)
  
!   nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
!   open(unit=10,file=nomfich,form='formatted',status='old')
!   read(10,*)
!   read(10,*)
!   read(10,'(A13,I11)')temp_label,levelmin
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
! 
!   read(10,*)
!   read(10,'(A13,E23.15)')temp_label,t
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
!   read(10,*)
! 
!   read(10,'(A14,A80)')temp_label,ordering
!   write(*,'(XA14,A20)')temp_label,TRIM(ordering)
!   read(10,*)
  allocate(cpu_list(1:ncpu))
!   if(TRIM(ordering).eq.'hilbert')then
!      allocate(bound_key(0:ncpu))
!      allocate(cpu_read(1:ncpu))
!      cpu_read=.false.
!      do impi=1,ncpu
!         read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
!      end do
!   endif
  close(10)

  !-----------------------
  ! Map parameters
  !-----------------------
  if(lmax==0)then
     lmax=nlevelmax
  endif
  write(*,*)'time=',t
!   write(*,*)'Working resolution =',2**lmax
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

!   if(TRIM(ordering).eq.'hilbert')then
! 
!      dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
!      do ilevel=1,lmax
!         dx=0.5d0**ilevel
!         if(dx.lt.dmax)exit
!      end do
!      lmin=ilevel
!      bit_length=lmin-1
!      maxdom=2**bit_length
!      imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
!      if(bit_length>0)then
!         imin=int(xxmin*dble(maxdom))
!         imax=imin+1
!         jmin=int(yymin*dble(maxdom))
!         jmax=jmin+1
!         kmin=int(zzmin*dble(maxdom))
!         kmax=kmin+1
!      endif
! 
!      dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
!      ndom=1
!      if(bit_length>0)ndom=8
!      idom(1)=imin; idom(2)=imax
!      idom(3)=imin; idom(4)=imax
!      idom(5)=imin; idom(6)=imax
!      idom(7)=imin; idom(8)=imax
!      jdom(1)=jmin; jdom(2)=jmin
!      jdom(3)=jmax; jdom(4)=jmax
!      jdom(5)=jmin; jdom(6)=jmin
!      jdom(7)=jmax; jdom(8)=jmax
!      kdom(1)=kmin; kdom(2)=kmin
!      kdom(3)=kmin; kdom(4)=kmin
!      kdom(5)=kmax; kdom(6)=kmax
!      kdom(7)=kmax; kdom(8)=kmax
!      
!      do i=1,ndom
!         if(bit_length>0)then
!            call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
!         else
!            order_min=0.0d0
!         endif
!         bounding_min(i)=(order_min(1))*dkey
!         bounding_max(i)=(order_min(1)+1.0D0)*dkey
!      end do
!      
!      cpu_min=0; cpu_max=0
!      do impi=1,ncpu
!         do i=1,ndom
!            if (   bound_key(impi-1).le.bounding_min(i).and.&
!                 & bound_key(impi  ).gt.bounding_min(i))then
!               cpu_min(i)=impi
!            endif
!            if (   bound_key(impi-1).lt.bounding_max(i).and.&
!                 & bound_key(impi  ).ge.bounding_max(i))then
!               cpu_max(i)=impi
!            endif
!         end do
!      end do
!      
!      ncpu_read=0
!      do i=1,ndom
!         do j=cpu_min(i),cpu_max(i)
!            if(.not. cpu_read(j))then
!               ncpu_read=ncpu_read+1
!               cpu_list(ncpu_read)=j
!               cpu_read(j)=.true.
!            endif
!         enddo
!      enddo
!   else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
!   end  if

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
!   nomfich=TRIM(outfich)
!   write(*,*)'Writing file '//TRIM(nomfich)
!  open(unit=20,file=nomfich,form='formatted',access='append')

!   write(*,*) 'ud,ul,ut',ud,ul,ut
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
! ROM: comment the single follwing line for old stuff
     read(10)
!      if(TRIM(ordering).eq.'bisection')then
!         do i=1,5
!            read(10)
!         end do
!      else
        read(10)
!      endif
     read(10)
     read(10)
     read(10)

     ! Open HYDRO file and skip header
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh
     read(11)
     read(11)
     read(11)
     read(11)gamma

     ! Loop over levels
     do ilevel=1,lmax

        ! Geometry
        dx=0.5**ilevel
        dx2=0.5*dx
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(xg(1:ngrida,1:ndim))
           allocate(son(1:ngrida,1:twotondim))
           allocate(var(1:ngrida,1:twotondim,1:nvarh))
           if(.not.allocated(varmax)) allocate(varmax(1:nvarh+8))
           allocate(x  (1:ngrida,1:ndim))
           allocate(rho(1:ngrida))
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
              !   ok_cell= .not.ref(i).and. &
              !        & (x(i,1)+dx2)>=xmin.and.&
              !        & (x(i,2)+dx2)>=ymin.and.&
              !        & (x(i,3)+dx2)>=zmin.and.&
              !        & (x(i,1)-dx2)<=xmax.and.&
              !        & (x(i,2)-dx2)<=ymax.and.&
              !        & (x(i,3)-dx2)<=zmax

                 ok_cell= (x(i,1)+dx2)>=xmin.and.&
                      & (x(i,2)+dx2)>=ymin.and.&
                      & (x(i,3)+dx2)>=zmin.and.&
                      & (x(i,1)-dx2)<=xmax.and.&
                      & (x(i,2)-dx2)<=ymax.and.&
                      & (x(i,3)-dx2)<=zmax

                 if(ok_cell)then
                 
                                   
                    rhoi = var(i,ind,1)*ud

                    bi   = sqrt(0.25d0*((var(i,ind,5)+var(i,ind,8))**2 &
                         &         +(var(i,ind,6)+var(i,ind,9))**2 &
                         &         +(var(i,ind,7)+var(i,ind,10))**2) &
                         &         *4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)

                    bxi = 0.5d0*(var(i,ind,5)+var(i,ind, 8))*sqrt(4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)
                    byi = 0.5d0*(var(i,ind,6)+var(i,ind, 9))*sqrt(4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)
                    bzi = 0.5d0*(var(i,ind,7)+var(i,ind,10))*sqrt(4.0d0*acos(-1.0d0)*ud*(ul/ut)**2)

                    tempi = var(i,ind,17)

                    xi = (x(i,1)-xcenter)*ul*boxlen
                    yi = (x(i,2)-ycenter)*ul*boxlen
                    zi = (x(i,3)-zcenter)*ul*boxlen
                    dxi = dx*ul*boxlen

                    ui = var(i,ind,2)*ul/ut
                    vi = var(i,ind,3)*ul/ut
                    wi = var(i,ind,4)*ul/ut
                    
                    veli = sqrt(ui*ui+vi*vi+wi*wi)

                    icell = icell + 1
                    data_array(icell, 1) = real(ilevel)
                    data_array(icell, 2) = xi
                    data_array(icell, 3) = yi
                    data_array(icell, 4) = zi
                    data_array(icell, 5) = dxi
                    data_array(icell, 6) = rhoi
                    data_array(icell, 7) = veli
                    data_array(icell, 8) = tempi
                    data_array(icell, 9) = bi
                    data_array(icell,10) = ui
                    data_array(icell,11) = vi
                    data_array(icell,12) = yi
                    data_array(icell,13) = bxi
                    data_array(icell,14) = byi
                    data_array(icell,15) = bzi
                                        
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,var,ref,rho,x)
        endif

     end do
     ! End loop over levels

     close(10)
     close(11)

  end do
  ! End loop over cpus
  
  ncells = icell
  write(*,*) 'Read ',ncells,' cells'

end subroutine ramses_data

