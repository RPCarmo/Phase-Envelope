!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!ATOMS - Aplied Thermodynamics and Molecular Simulation!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!UFRJ - Federal University of Rio de Janeiro!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Developed by Rafael Pereira


program main
  
  use GaussianElimination
  use EOS
  
  implicit none
  
  !Continuation Method Variables
  integer :: SpecVar, SpecVar_old, point, flag_crit, maxK_i, flag_error
  real(8) :: S, dS, dSmax, maxTstep, K_CritPoint
  real(8), allocatable :: Var(:), Var_old(:), dfdS(:), sensitivity(:)
  
  !Newton's Method Variables
  real(8) :: diffFrac, diffT, diffP
  real(8), allocatable :: F(:), dF(:), step(:)
  
  !EoS Variables
  real(8) :: Volume(2), amix(2), bmix(2), FugCoef_aux, FugCoef_ref, R = 83.14462175d0!R = cm³.bar/(mol.K)
  real(8), allocatable :: kij(:,:), lij(:,:), Tc(:), Pc(:), acentric(:), b(:), a(:), ac(:)
  
  !Phase Equilibrium Variables
  integer :: comp, phase(2), phase_aux
  real(8) :: P, T
  real(8), allocatable :: Composition(:,:), K(:), z(:)
  
  !Input/Output Variables
  integer :: input_num, output_num, output_num2
  character(LEN=3) :: legend_ELV(2) = ['Vap','Liq']
  
  !Auxiliary Variables
  integer :: i, j, it, maxit
  real(8) :: tol, diff, Gibbs_vap, Gibbs_liq, aux, T_old, maxstep
  
  
  open(NEWUNIT = input_num, FILE = "input.txt", ACTION = "read", STATUS = "old")
  open(NEWUNIT = output_num, FILE = "output.csv", ACTION = "write", STATUS = "unknown")
  open(NEWUNIT = output_num2, FILE = "output.txt", ACTION = "write", STATUS = "unknown")
  
  write(output_num,*) "Incipient Phase Type,Pressure,Temperature,Incipient Phase Composition"
  
  read(input_num,*) comp
  
  allocate(b(comp), ac(comp), a(comp), Tc(comp), Pc(comp), acentric(comp))
  
  allocate(K(comp,9), Composition(comp,10), z(comp), kij(comp,comp), lij(comp,comp))
  
  !defining variables maximum size supposing a maximum number of phases equal to 8
  !Degrees of freedom = (F-1)*C + F + 1
  allocate(dF((8*comp+8+1)*(8*comp+8+1)), step(8*comp+8+1), dfdS(8*comp+8+1), Var(8*comp+8+1), sensitivity(8*comp+8+1), F(8*comp+8+1))
  
  !     Independent Variables     !
  ! (F-1)*C            K          !
  !  (F-1)           beta         !
  !   1           temperature     !
  !   1            pressure       !
  
  !Reading And Calculating Properties********************************************************************************************
  aux = 0.0d0
  do i = 1,comp
      !Reading Global Composition, Critical Temperature, Critical Pressure and Acentric Factor
      read(input_num,*) z(i), Tc(i), Pc(i), acentric(i)
      
      !EoS Parameters Calculation
      b(i) = 0.07780d0*R*Tc(i)/Pc(i) !covolume
      ac(i) = 0.45724*R*R*Tc(i)*Tc(i)/Pc(i)
      aux = aux + z(i)
  enddo
  
  !Binary Interaction Parameter
  kij = 0.d0
  lij = 0.d0

  kij(:,1) = kij(1,:)
  lij(:,1) = lij(1,:)

  close(input_num)
  !******************************************************************************************************************************
  
  
  
  !Initial Settings - Dew point ///////////////////////////////////////////////////////////////////////////////////////////////////
  T = 298.0d0 !Initial Temperature Guess (K)
  P = 1.d0 !Initial Pressure (bar)
  
  do i = 1,comp
    K(i,1) = 1.0d0/dexp(dlog(Pc(i)/P) + 5.373d0*(1.0d0 + Acentric(i))*(1.0d0 - Tc(i)/T)) !Whitson's Approach for Vapor-Liquid Equilibria
    z(i) = z(i)/aux !Normalizing Global Composition
    Composition(i,1) = z(i) !Reference Phase Composition
  enddo
  
  phase(1) = 0 !Reference Phase Index (Vapor)
  phase(2) = 1 !Incipient Phase Index (Liquid)
  !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  !Determining Temperature Initial Guess (Point Near To The Dew Curve)***********************************************************
  T_old = T - 1.0d0
  do while(T_old .ne. T)
      T_old = T
      
      call EoS_param(acentric,Tc,ac,a,T,comp) !Updating Attractive Parameter
      call VdW1fMIX (comp,a,b,kij,lij,z,amix(1),bmix(1)) !Mixing Rule
      call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), 0) !Calculating Vapor Phase Molar Volume
      call EoS_Volume(P, T, bmix(1), amix(1), Volume(2), 1) !Calculating Liquid Phase Molar Volume
      !Calculating Gibbs Energy Of Vapor And Liquid Phases
      Gibbs_vap = 0.0d0
      Gibbs_liq = 0.0d0
      do i = 1,comp
          call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),z,kij(i,:),lij(i,:),i)
          Gibbs_vap = Gibbs_vap + z(i)*FugCoef_ref
          call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_aux,Volume(2),z,kij(i,:),lij(i,:),i)
          Gibbs_liq = Gibbs_liq + z(i)*FugCoef_aux
      enddo
      
      !If Exists A Pure Liquid Or A Liquid-Vapor Equilibrium Nearer To The Bubble Point
      !Than To The Dew Point, The Temperature Is Increased.
      if((Gibbs_liq .lt. Gibbs_vap) .or. (Gibbs_liq .eq. Gibbs_vap .and. (Volume(1)/bmix(1)) .lt. 1.75d0)) then
          T = T + 10.0d0
      endif
      !OBS: The test comparing the ratio between the volume and the covolume to the constant factor 1.75 was proposed by
      !Pedersen and Christensen (Phase Behavior Of Petroleum Reservoir Fluids, 2007 - Chapter 6.5 Phase Identification)
  enddo
  !******************************************************************************************************************************
  
  
  
  !Successive Substitution///////////////////////////////////////////////////////////////////////////////////////////////////////
  tol = 1.0d-5
  diff = 1.0d-6
  step(1) = 1.0d+6
  maxit = 100
  it = 0
  do while(dabs(step(1)) .gt. tol .and. it .lt. maxit)
      it = it + 1
      T_old = T
      
      aux = 0.0d0
      do i = 1,comp
        Composition(i,2) = Composition(i,1)*K(i) !Incipient Phase Composition
        aux = aux + Composition(i,2)
      enddo
      
      F(1) = -1.0d0
      dF(1) = 0.0d0
      
      do i = 1,comp
        Composition(i,2) = Composition(i,2)/aux !Normalizing Composition
        F(1) = F(1) + Composition(i,1)*K(i) !Residual
      enddo
      
      !Numerical Derivative With Respect to Temperature
      T = T_old + diff
      call EoS_param(acentric,Tc,ac,a,T,comp) !Updating Attractive Parameter
      call VdW1fMIX (comp,a,b,kij,lij,Composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
      call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
      call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), 0)
      call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), 1)
      do i = 1,comp
        call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
        call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
        dF(1) = dF(1) + Composition(i,1)*K(i)*(FugCoef_ref - FugCoef_aux)
      enddo
      
      T = T_old - diff
      call EoS_param(acentric,Tc,ac,a,T,comp) !Updating Attractive Parameter
      call VdW1fMIX (comp,a,b,kij,lij,Composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
      call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
      call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), 0)
      call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), 1)
      do i = 1,comp
        call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
        call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
        dF(1) = dF(1) - Composition(i,1)*K(i)*(FugCoef_ref - FugCoef_aux)
      enddo
    
      dF(1) = dF(1)/(2.0d0*diff)
      
      !Temperature Step Calculation
      step(1) = F(1)/dF(1)
      
      !Step Brake
      if(dabs(step(1)) .gt. 0.25d0*T_old) step(1) = 0.25d0*T_old*step(1)/dabs(step(1)) 
      
      !Updating Temperature
      T = T_old - step(1)
      
      
      !Updating K-factors
      call EoS_param(acentric,Tc,ac,a,T,comp)
      call VdW1fMIX (comp,a,b,kij,lij,Composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
      call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
      call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), 0)
      call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), 1)
      do i = 1,comp
        call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
        call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
        K(i) = dexp(FugCoef_ref - FugCoef_aux)
      enddo
  enddo
  
  if(it .eq. maxit .and. dabs(step(1)) .gt. tol) then
      print*, "WARNING: In Successive Substitution Method - Maximum Number of Iterations Reached!"
      print*, "Exiting Program..."
      stop
  endif
  !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  !Continuation Method And Newton's Method Settings******************************************************************************
  S = dlog(P) !Specified Variable Value
  dS = 0.1d0 !Specified Variable Variation
  do i = 1,comp
    Var(i) = dlog(K(i))
  enddo
  Var(comp+1) = dlog(T)
  SpecVar = comp + 2 !Specified Variable Index
  Var(SpecVar) = S !Specified Independent Variable
  do i = 1,(comp+1)
      dfdS(i) = 0.0d0
  enddo
  dfdS(comp+2) = 1.0d0
  tol = 1.0d-6
  diff = 1.0d-6
  maxit = 100 !Maximum Number Of Iteration In Newton's Method
  maxTstep = 5.0d0 !Maximum Temperature Step In Continuation Method
  K_CritPoint = 0.04d0 !K-factor Reference Value Used To Detect Critical Points
  flag_crit = 0 ! if calculating the critical point, flag_crit = 1
  flag_error = 0
  !******************************************************************************************************************************
  
  
  print*, P, T
  !Initializing Continuation Method
  point = 0
  do while(P .ge. 0.5d0 .and. P .lt. 1000.d0 .and. flag_error .eq. 0) !Main Loop
      point = point + 1
      it = 0
      maxstep = 1.d+6
      do while(maxstep .gt. tol .and. it .lt. maxit) !Newton's Method Loop
          it = it + 1
          
          !Calculating Residuals/////////////////////////////////////////////////////////////////////////////////////////////////
          call EoS_param(acentric,Tc,ac,a,T,comp)
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), phase(1))
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), phase(2))
          F(comp+1) = 0.0d0
          do i = 1,comp
            call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
            call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
            !Residual Responsible For The Chemical Equilibrium
            F(i) = Var(i) + (FugCoef_aux - FugCoef_ref)
            !Residual Responsible For Assuring That The Summation Of The Incipient Phase Equals 1
            F(comp+1) = F(comp+1) + Composition(i,2) - Composition(i,1)
          enddo
          !Residual Responsible For Determining The Specified Independent Variable
          F(comp+2) = Var(SpecVar) - S
          !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          
          
          
          !Differentiating The First "C" Residuals With Respect to ln[K(j)]******************************************************
          do i = 1,comp
              do j = 1,comp
                  diffFrac = diff*Composition(j,2)
                  
                  aux = Composition(j,2)
                  !Numerically Differentiating the Fugacity Coefficient of the Incipient Phase
                  Composition(j,2) = aux + diffFrac
                  call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
                  call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), phase(2))
                  call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
                  dF((i-1)*(comp+2)+j) = FugCoef_aux
                  
                  Composition(j,2) = aux - diffFrac
                  call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
                  call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), phase(2))
                  call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
                  dF((i-1)*(comp+2)+j) = dF((i-1)*(comp+2)+j) - FugCoef_aux
                  
                  Composition(j,2) = aux
                  !Derivative of ln[FugacityCoefficient(IncipientPhase,Component i)] With Respect to ln[K(j)]
                  dF((i-1)*(comp+2)+j) = dF((i-1)*(comp+2)+j)*Composition(j,2)/(2.0d0*diffFrac)
                  
                  !Derivative of ln[K(i)] With Respect to ln[K(j)] = Kronecker Delta
                  if(i .eq. j) dF((i-1)*(comp+2)+j) = dF((i-1)*(comp+2)+j) + 1.0d0
              enddo
          enddo
          !**********************************************************************************************************************
          
          
          
          !Differentiating "C+1" Residual With Respect to ln[K(i)]///////////////////////////////////////////////////////////////
          do i = 1,comp
             dF(comp*(comp+2)+i) = Composition(i,2)
          enddo
          !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          
          
          
          !Differentiating The First "C" Residuals With Respect to ln[T]*********************************************************
          diffT = diff*Var(comp+1)
          
          !Numerically Differentiating The ln(FugacityCoefficient) With Respect to ln(T)
          T = dexp(Var(comp+1) + diffT)
          call EoS_param(acentric,Tc,ac,a,T,comp)
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), phase(1))
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), phase(2))
          do i = 1,comp
              call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
              call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
             dF((i-1)*(comp+2)+comp+1) = FugCoef_aux - FugCoef_ref
          enddo
          
          T = dexp(Var(comp+1) - diffT)
          call EoS_param(acentric,Tc,ac,a,T,comp)
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), phase(1))
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), phase(2))
          do i = 1,comp
              call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
              call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
             dF((i-1)*(comp+2)+comp+1) = dF((i-1)*(comp+2)+comp+1) - (FugCoef_aux - FugCoef_ref)
          enddo
          
          do i = 1,comp
             dF((i-1)*(comp+2)+comp+1) = dF((i-1)*(comp+2)+comp+1)/(2.0d0*diffT)
          enddo
          
          T = dexp(Var(comp+1))
          call EoS_param(acentric,Tc,ac,a,T,comp)
          !OBS: The derivative ok ln(K) with respect to ln(T) is null.
          !**********************************************************************************************************************
          
          
          
          !Differentiating The First "C" Residuals With Respect to ln[P]/////////////////////////////////////////////////////////
          diffP = diff*Var(comp+2)
          
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
          call VdW1fMIX (comp,a,b,kij,lij,Composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
          
          !Numerically Differentiating The ln(FugacityCoefficient) With Respect to ln(T)
          P = dexp(Var(comp+2) + diffP)
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), phase(1))
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), phase(2))
          do i = 1,comp
              call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
              call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
             dF((i-1)*(comp+2)+comp+2) = FugCoef_aux - FugCoef_ref
          enddo
          
          P = dexp(Var(comp+2) - diffP)
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), phase(1))
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), phase(2))
          do i = 1,comp
              call fugacity(comp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),Composition(:,1),kij(i,:),lij(i,:),i)
              call fugacity(comp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),Composition(:,2),kij(i,:),lij(i,:),i)
             dF((i-1)*(comp+2)+comp+2) = dF((i-1)*(comp+2)+comp+2) - (FugCoef_aux - FugCoef_ref)
          enddo
          
          do i = 1,comp
             dF((i-1)*(comp+2)+comp+2) = dF((i-1)*(comp+2)+comp+2)/(2.0d0*diffP)
          enddo
          
          P = dexp(Var(comp+2))
          !OBS: The derivative ok ln(K) with respect to ln(P) is null.
          !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          
          
          
          !Derivative of the "C+1" Residual With Respect to ln(T)****************************************************************
          dF(comp*(comp+2)+comp+1) = 0.0d0
          !**********************************************************************************************************************
          
          
          
          !Derivative of the "C+1" Residual With Respect to ln(P)////////////////////////////////////////////////////////////////
          dF(comp*(comp+2)+comp+2) = 0.0d0
          !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          
          
          
          !Derivative of the "C+2" Residual**************************************************************************************
          do i = 1,(comp+2)
              dF((comp+1)*(comp+2)+i) = 0.0d0
          enddo
          dF((comp+1)*(comp+2)+SpecVar) = 1.0d0
          !OBS: The derivative of the specification equation with respect to all independent variables 
          !besides the specified variable are null. Its derivative with respect to the specified variable is 1.
          !**********************************************************************************************************************
          
          
          
          !Solving The System of Equations///////////////////////////////////////////////////////////////////////////////////////
          call GaussElimination(dF,F,step,comp+2)
          !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          
          
          
          !Updating The Independent Variables************************************************************************************
          maxstep = 0.0d0
          do i = 1,(comp+2)
              Var(i) = Var(i) - step(i)
              !Calculating Variable Used As Convergence Criteria
              if(maxstep .lt. dabs(step(i)/Var(i))) maxstep = dabs(step(i)/Var(i))
          enddo
          !**********************************************************************************************************************
          
          
          
          !Calculating The Natural Form Of Independent Variables And Updating Compositions Of The Incipient Phase////////////////
          do i = 1,comp
              K(i) = dexp(Var(i))
              Composition(i,2) = Composition(i,1)*K(i)
          enddo
          T = dexp(Var(comp+1))
          P = dexp(Var(comp+2))
          !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      enddo !End of Newton's Method Loop
      
      write(*,*) "Incipient Phase = ",legend_ELV(phase(2)+1),"    P = ",P,"T = ",T, "SpecVar =",SpecVar

      if(maxstep .gt. tol) then
          flag_error = 1
      else
          do i = 1,(comp+2)
              if(Var(i) .ne. Var(i)) then
                  flag_error = 1
                  exit
              endif
          enddo
      endif

      if(flag_error .eq. 0) then
          if(flag_crit .eq. 2) flag_crit = 0

          write(output_num,*) legend_ELV(phase(2)+1),",",P,",",T,",",(Composition(i,2),",",i=1,comp)
          write(output_num2,*) T, P
          
          !Analyzing Sensitivity Of The Independent Variables************************************************************************
          SpecVar_old = SpecVar
          Var_old = Var
          
          !Sensitivity Vector Calculation
          call GaussElimination(dF,dfdS,sensitivity,comp+2)

          if(flag_crit .eq. 0) then
              !Choosing The New Specified Independent Variable
              do i = 1,(comp+2)
                  if(dabs(sensitivity(i)) .gt. dabs(sensitivity(SpecVar))) SpecVar = i
              enddo
              !OBS: The specified variable is the one with the greatest sensitivity,
              !i.e. the one which its variation makes the system varies more intensely.

              !Updating Specified Variable
              if(SpecVar .ne. SpecVar_old) then
                  dS = dS*sensitivity(SpecVar)
                  do i = 1,(comp+2)
                      if(i .ne. SpecVar) sensitivity(i) = sensitivity(i)/sensitivity(SpecVar)
                  enddo
                  sensitivity(SpecVar) = 1.0d0
                  S = Var(SpecVar)
              endif
              !**************************************************************************************************************************


              !Adjusting Stepsize////////////////////////////////////////////////////////////////////////////////////////////////////////
              dSmax = (dabs(Var(SpecVar))**0.5d0)/10.0d0
              if(dSmax .lt. 0.1d0) dSmax = 0.1d0
              dSmax = dabs(dSmax)*(dabs(dS)/dS)
              dS = dS*4.0d0/it
              if(dabs(dSmax) .lt. dabs(dS)) dS = dSmax
              !Defining Specified Variable Value In The Next Point
              S = S + dS
              !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
              
              
              !Independent Variables Initial Guess For The Next Point********************************************************************
              do i = 1,(comp+2)
                  Var(i) = Var(i) + dS*sensitivity(i)
              enddo
              !**************************************************************************************************************************
              
              
              !Analyzing Temperature Stepsize////////////////////////////////////////////////////////////////////////////////////////////
              T_old = T
              T = dexp(Var(comp+1))
              !Large Temperature Steps Are Not Advisable
              do while(dabs(T-T_old) .gt. maxTstep)
                  dS = dS/2.0d0
                  S = S - dS
                  do i = 1,(comp+2)
                      Var(i) = Var(i) - dS*sensitivity(i)
                  enddo
                  T = dexp(Var(comp+1))
              enddo
              !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                
              
              !Analyzing Proximity to Critical Point*************************************************************************************
              !Seeking The Greatest K-factor
              maxK_i = maxloc(abs(Var(1:comp)),dim=1)
print*, Var(maxK_i)
              !If The ln[K(i)] Stepsize Is Too Big, It Should Be Decreased
              if(abs(Var(maxK_i)) .lt. 0.1d0) then
                  !Analyzing maxK stepsize
                  if(abs(dS*sensitivity(maxK_i)) .gt. K_CritPoint) then
                      S = S - dS
                      do i = 1,(comp+2)
                          Var(i) = Var(i) - dS*sensitivity(i)
                      enddo

                      !Shortening ln(K) Stepsize
                      dS = (dS/abs(dS)) * K_CritPoint/abs(sensitivity(maxK_i))

                      S = S + dS
                      do i = 1,(comp+2)
                          Var(i) = Var(i) + dS*sensitivity(i)
                      enddo
                  endif
                  !OBS: The current point must be near enough to the critical point
                  !so that the algorithm can pass through it without diverging.

                  if(abs(Var(maxK_i)) .lt. K_CritPoint) then
                  !This is gonna be the last point before the critical point
                      S = S - dS
                      do i = 1,(comp+2)
                          Var(i) = Var(i) - dS*sensitivity(i)
                      enddo

                      dS = (K_CritPoint - Var(maxK_i))/sensitivity(maxK_i)

                      S = S + dS
                      do i = 1,(comp+2)
                          Var(i) = Var(i) + dS*sensitivity(i)
                      enddo
                      flag_crit = 1
                  endif
              endif
          else
              !Passing Through The Critical Point
              dS = K_CritPoint*dS/abs(dS)
              S = S + 2.d0*dS
              do i = 1,(comp+2)
                  Var(i) = Var(i) + 2.d0*dS*sensitivity(i)
              enddo

              !Defining Incipient Phase As Vapor - Initializing Bubble Curve
              phase_aux = phase(1)
              phase(1) = phase(2)
              phase(2) = phase_aux
              !OBS: This will cause the definition of the K-factors to change from
              !FugCoef(Vap phase)/FugCoef(Liq phase) to FugCoef(Liq phase)/FugCoef(Vap phase).

              flag_crit = 2
          endif

      elseif(flag_error .eq. 1 .and. flag_crit .eq. 2 .and. K_CritPoint .gt. 0.009d0) then
          K_CritPoint = K_CritPoint - 0.005d0
          Var = Var_old
          S = Var(SpecVar)
          phase_aux = phase(1)
          phase(1) = phase(2)
          phase(2) = phase_aux
          dS = (K_CritPoint - Var(maxK_i))/sensitivity(maxK_i)

          S = S + dS
          do i = 1,(comp+2)
              Var(i) = Var(i) + dS*sensitivity(i)
          enddo

          flag_crit = 1
          flag_error = 0
      elseif(flag_error .eq. 1 .and. flag_crit .eq. 0 .and. abs(dS) .gt. 1.d-6) then
          Var = Var_old
          S = Var(SpecVar)
          dS = dS/4.d0

          S = S + dS
          do i = 1,(comp+2)
              Var(i) = Var(i) + dS*sensitivity(i)
          enddo

          flag_error = 0
      endif

      do i = 1,comp
          K(i) = dexp(Var(i))
          Composition(i,2) = Composition(i,1)*K(i)
      enddo
      T = dexp(Var(comp+1))
      P = dexp(Var(comp+2))
  enddo !End of  the main loop
  
  
  close(output_num)
  close(output_num2)
  
  call system('gnuplot -p plot.gnu')
  
  deallocate(b,ac)
  deallocate(a,Tc)
  deallocate(Pc,acentric)
  deallocate(Var,Composition)
  deallocate(sensitivity,F)
  deallocate(dF,step)
  deallocate(dfdS,K)
  deallocate(z,kij,lij)
  
  
end program


