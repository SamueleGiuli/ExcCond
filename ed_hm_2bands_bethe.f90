program ed_hm_2bands_bethe
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
  implicit none
  integer                                     :: jo,js,iloop,Le,Nso
  logical                                     :: converged
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,jspin,inso,print_mode
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:),Weiss_(:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable                      :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  complex(8),allocatable,dimension(:)         :: Gtest
  !hamiltonian input:
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(2)                        :: Wband,de
  real(8),dimension(:,:),allocatable          :: Dbands, DbandsLoc
  real(8),dimension(:,:),allocatable          :: Ebands ![Nso][Le]
  real(8),dimension(:),allocatable            :: H0     ![Nso]
  !variables for the model:
  real(8)                                     :: epsbath=1e-4, t1,t2
  real(8)                                     :: mh,wmixing,alpha,excparl,excapar,DelBand,Wuno,Ds,add
  character(len=30)                           :: Params
  character(len=16)                           :: finput
  character(len=32)                           :: hkfile
  logical                                     :: spinsym,mixG0
  !
  real(8),dimension(2)                        :: Eout
  real(8),allocatable                         :: dens(:)
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma5,GammaN
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  complex(8),dimension(4,4)                   :: GammaR0,GammaRx,GammaRy,GammaRz
  real(8),dimension(:),allocatable            :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer                                     :: irank,comm,rank,size2,ierr
  logical                                     :: master

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)

  call cpu_time(t1)
  

  !Parse additional variables && read Input
  call parse_cmd_variable(finput,"FINPUT",default='inputED.in')  
  call parse_input_variable(Params,"Params",finput,default="N",&
       comment="Mh-E0; N-E0; N-Mh; N-Mh-E0; Mh-N-Ex; N-Mh-E0-Ex")
  call parse_input_variable(alpha,"ALPHA",finput,default=1d0,comment="bandwidth ratio W_2 = alpha*W_1=alpha*1.0")
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=1.d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(wuno,"Wuno",finput,default=1.d0)
  call parse_input_variable(delband,"DELBAND",finput,default=0.d0)
  call parse_input_variable(epsbath,"EPSBATH",finput,default=-1.0d-4,comment="If < 1e-4 0 no fit over diagonal part of G")
  call parse_input_variable(excparl,"EXCPARL",finput,default=0.1d0)
  call parse_input_variable(excapar,"EXCAPAR",finput,default=0.1d0)
  !
  !
  call ed_read_input(trim(finput),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(bath_type/="replica")stop "Wrong setup from input file: non replica bath"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if( (excapar.NE.0.d0).AND.(ed_mode=="normal".OR.Nspin/=2) )stop "Wrong setup from input file: chose normal mode or Nspin/=2"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))

  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaN=kron_pauli( pauli_sigma_0, pauli_tau_0)
  !
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron_pauli( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron_pauli( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron_pauli( pauli_sigma_z, pauli_tau_x )
  !
  gammaR0=kron_pauli( pauli_sigma_0, pauli_tau_y )
  gammaRx=kron_pauli( pauli_sigma_x, pauli_tau_y )
  gammaRy=kron_pauli( pauli_sigma_y, pauli_tau_y )
  gammaRz=kron_pauli( pauli_sigma_z, pauli_tau_y )

  
  !Build Hamiltonian structure:
  Wband(1)=1.d0
  Wband(2)=alpha*Wband(1)
  !
  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  !
  do ispin=1,Nspin
     do iorb=1,Norb
        inso = iorb+2*(ispin-1)
        Ebands(inso,:) = linspace(-Wband(iorb), Wband(iorb), Le,mesh=de(iorb))
        Dbands(inso,:) =  dens_bethe(Ebands(inso,:),Wband(iorb))*de(iorb)
     enddo
  enddo
  !
  allocate(DbandsLoc(1,Le))  
  DbandsLoc(1,:) = dens_bethe(Ebands(1,:),Wband(1))*de(1)
  !  
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(H0(Nso))
  Hloc = zero
  H0   = zero
  do js=1,Nspin
     Hloc(js,js,:,:)= Mh*pauli_sigma_z
     do jo=1,Norb
        H0(jo+2*(js-1)) =Hloc(js,js,jo,jo)
     end do
  end do
  !
  
  

  print_mode=3
  if(ed_mode=="nonsu2")print_mode=4

  !Setup solver
  if(bath_type=="replica")then
     select case(trim(Params))
     case default
        stop "Params not in list"

     case("N")
        allocate(lambdasym_vector(1))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,1))
        Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso)) ;lambdasym_vector(1)=0.d0+epsbath

     case("N-Mh")
        allocate(lambdasym_vector(2))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
        Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso)) ;lambdasym_vector(1)=0.0+epsbath
        Hsym_basis(:,:,:,:,2)=j2so(Gamma5(:Nso,:Nso)) ;lambdasym_vector(2)=Mh
        
     case("Mh-E0")
        allocate(lambdasym_vector(2))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5(:Nso,:Nso)) ;lambdasym_vector(1)=0.0d0+epsbath
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0(:Nso,:Nso)) ;lambdasym_vector(2)=-excparl
        
     case("N-E0")
        allocate(lambdasym_vector(2))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
        Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso)) ;lambdasym_vector(1)=0.0d0+epsbath
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0(:Nso,:Nso)) ;lambdasym_vector(2)=-excparl

     case("N-Mh-E0")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso))  ;lambdasym_vector(1)=0.0d0
        Hsym_basis(:,:,:,:,2)=j2so(Gamma5(:Nso,:Nso))  ;lambdasym_vector(2)=Mh
        Hsym_basis(:,:,:,:,3)=j2so(GammaE0(:Nso,:Nso)) ;lambdasym_vector(3)=-excparl


     case("N-Mh-Ex")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso))  ;lambdasym_vector(1)=0.d0
        Hsym_basis(:,:,:,:,2)=j2so(Gamma5(:Nso,:Nso))  ;lambdasym_vector(2)=Mh
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx(:Nso,:Nso)) ;lambdasym_vector(3)=-excapar

     case("N-E0-Ex")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso))  ;lambdasym_vector(1)=0.0d0+epsbath
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0(:Nso,:Nso)) ;lambdasym_vector(2)=-excparl
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx(:Nso,:Nso)) ;lambdasym_vector(3)=-excapar

     case("N-Mh-E0-Ex")
        allocate(lambdasym_vector(4))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,4))
        Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso))  ;lambdasym_vector(1)=0.d0
        Hsym_basis(:,:,:,:,2)=j2so(Gamma5(:Nso,:Nso))  ;lambdasym_vector(2)=Mh
        Hsym_basis(:,:,:,:,3)=j2so(GammaE0(:Nso,:Nso)) ;lambdasym_vector(3)=-excparl
        Hsym_basis(:,:,:,:,4)=j2so(GammaEx(:Nso,:Nso)) ;lambdasym_vector(4)=-excapar

     end select


     call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
     Nb=ed_get_bath_dimension(Hsym_basis)
     allocate(Bath(Nb))
     allocate(Bath_(Nb))
     call ed_init_solver(comm,bath)

  else     
     Nb=ed_get_bath_dimension()
     allocate(Bath(Nb))
     allocate(Bath_(Nb))
     call ed_init_solver(comm,bath)
  endif


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath,Hloc)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_dens(dens)


     !Get GLOC: N.B. here DbandsLoc(1,:) is needed for non diagonal inversion
     call dmft_gloc_matsubara(Ebands,DbandsLoc,H0,Gmats,Smats)
     !Get Weiss to be fitted
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,trim(cg_scheme))

     !Update WeissField:
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=print_mode)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=print_mode)

     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif
     !Fit the new bath, starting from the old bath + the supplied Weiss
     select case(ed_mode)
     case default
        stop "ed_mode!=Normal/Nonsu2"
     case("normal")
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
     case("nonsu2")
        call ed_chi2_fitgf(comm,Weiss,bath)
     end select

     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !
     Gtest=Weiss(1,1,1,1,:)
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)
     call end_loop
  enddo

  !Compute local GFs
  call ed_get_sigma_realaxis(Sreal)
  call dmft_gloc_realaxis(Ebands,DbandsLoc,H0,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=print_mode)

  !Compute Kinetic Energy
  call dmft_kinetic_energy(Ebands,DbandsLoc,H0,Smats)

  call get_Ds(Ds,H0,Ebands,Dbands/de(1),Smats)

  open(8,file="Ds.dat")
  write(8,*) "# Superfluid density"
  write(8,*) Ds
  close(8)
  print*, "Ds:",Ds

  
  call cpu_time(t2)
  print*,"Total Simulation time: ",t2-t1, "ms"
  open(8,file="time.dat")
  write(8,*) "# Simulation time (ms)"
  write(8,*) t2-t1
  close(8)
  
  call finalize_MPI()



contains


  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so


  subroutine get_Ds(Ds,H0,Ebands,Dbands,Smats)
    integer             :: Lw, Lk, Nso
    integer             :: iw, ik, i,j
    real(8)             :: t
    real(8), dimension(:)                  :: H0
    real(8), dimension(:,:)                :: Ebands, Dbands
    complex(8), dimension(:,:,:,:,:)          :: Smats
    real(8), allocatable, dimension(:)     :: Veps, wm
    complex(8), allocatable, dimension(:,:,:) :: Gtot
    complex(8), allocatable, dimension(:,:)   :: Gtmp
    real(8)                                :: Ds, den
    
    Nso=size(Ebands,1)
    Lk=size(Ebands,2)
    Lw=size(Smats,5)
    t=0.25d0*( Ebands(1,Lk)-Ebands(1,1) )
    allocate(Gtmp(Nso,Nso),Veps(Lk),wm(Lw),Gtot(Nso,Nso,Lk))
    Veps = (4.d0*t**2 - Ebands(1,:)**2)/3.d0
    wm = pi/beta*dble(2*arange(1,Lw)-1)
    Ds=0.d0
    
    do iw=1,Lw
       Gtot=0.d0
       do ik=1,Lk
          Gtmp = (xi*wm(iw)*eye(Nso) - diag(Ebands(:,ik)) - diag(H0) -so2j( Smats(:,:,:,:,iw) ) )
          call inv(Gtmp)
          Gtot(:,:,ik) = Gtmp
       end do

       !2 since 2*Lambda_ab
       add = trapz( Dbands(1,:)*Veps*( 2.d0*abs(Gtot(1,2,:))**2 ), Ebands(1,:) )
       Ds = Ds + add
    end do
    add = add*wm(Lw)*2.d0/pi
    ! spin/beta
    Ds = (1-alpha)*Ds*4.d0/beta 

    print*, "Ds :", Ds, " - add:", add
    
      
  end subroutine get_Ds
  
end program ed_hm_2bands_bethe



