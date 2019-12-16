    ! ATOMS - Applied Thermodynamics and Molecular Simulation
    !======================================================================================
    !
    ! Nome do arquivo:  TPD.for
    ! Projeto:          Precipitação de parafinas
    ! Conteúdo:         Implementação das subrotinas para análise de estabilidade via TPD
    !
    !======================================================================================
    module TPD

    use EOS
    use LowerUpperDecomposition
    use LiquidModels
    use SolidModels

    implicit none

    contains

    !**************************TPD ANALYSIS****************************************/
    !TPD < 0 --> instability
    !TPD = 0 and MolNum = z (trivial solution)
    !TPD = 0 and MolNum =/= z --> Saturation point (Incipient phase detected)  
    !TPD > 0 and MolNum =/= z --> Shadow region
    !******************************************************************************/

    !==========================================================================================================================================================================
    subroutine TPDliqvap (comp, T, P, a, b, dHf, Tf, dHtr, Ttr, MW, integraldCp1, integraldCp2, BinIntCoef, PhaseNum,&
                        phase, Composition, CompositionTrial, phaseTrial, TPDout, ret, sigma_eos, epsilon_eos)

    !Arguments
    real(8) :: T, P, Composition(comp,comp), CompositionTrial(comp), TPDout, sigma_eos, epsilon_eos
    real(8) :: a(*), b(*), dHf(*), Tf(*), dHtr(*), Ttr(*), MW(*), integraldCp1(*), integraldCp2(*), BinIntCoef(comp, comp)
    integer :: comp, PhaseNum, phase(*), phaseTrial, ret

    !local
    real(8) :: MolNum(comp), d(comp), FugCoef, grad(comp), maxstep
    real(8) :: Hess(comp*comp), Gibbs(2), zaux(comp), step(comp)
    real(8) :: dMol, tol, TPD, TotMolTrial, aux, Vol, amix, bmix, Samix
    integer :: i, j, Trial, Ref, it
    character*6 :: names(2) = (/'Vapor ','Liquid'/)


    !Initializing Variables
    tol = 1.D-10
    dMol = 1.D-6
    it = 1
    ret = 0
    
    !Detecting the reference phase nature: Liquid or vapor
    Ref = phase(PhaseNum)
    if (Ref .eq. 1) then
        Trial = 0
    else 
        Trial = 1
    endif
    phaseTrial = Trial
!    phaseTrial = abs(Ref-1)

    !Composition initial estimates and calculation of TPD's term d(i)
    do i = 1, comp  
      !Trial phase composition initial guess
      MolNum(i) = 1.0d0/abs(Trial*comp-(i-1)-Ref)/comp

      !Feed role on TPD calculation
      call VdW1fMIX (comp, a, b, BinIntCoef, Composition(1,PhaseNum), Samix, amix, bmix, i)

      call fugacity(T, P, a(i), b(i), amix, bmix, Samix, FugCoef, Ref, Vol, sigma_eos, epsilon_eos)

      d(i) = log(Composition(i,PhaseNum)) + FugCoef
    enddo

    !Successive substitution
    do while (maxstep .gt. tol .and. it .le. 20)
      maxstep = 0.0d0
      do i = 1, comp        
        call VdW1fMIX (comp, a, b, BinIntCoef, MolNum, Samix, amix, bmix, i)

        call fugacity(T, P, a(i), b(i), amix, bmix, Samix, FugCoef, Trial, Vol, sigma_eos, epsilon_eos)

        !Calculating new compositions
        step(i) = exp(d(i) - FugCoef)

        !Selecting of the steps squared - Convergence criteria
        maxstep = maxstep + (step(i)-MolNum(i))*(step(i)-MolNum(i))
      enddo

      !Updating compositions
      do i = 1, comp
        MolNum(i) = step(i)
      enddo
      it = it + 1
    enddo

    !Newton Method
    if (it .eq. 21) then
      it = 1
      do while (maxstep.gt.tol .and. it.le.50)
        do i = 1, comp
            
          call VdW1fMIX (comp, a, b, BinIntCoef, MolNum, Samix, amix, bmix, i)

          call fugacity(T, P, a(i), b(i), amix, bmix, Samix, FugCoef, Trial, Vol, sigma_eos, epsilon_eos)

          !Calculating TM gradient vector
          grad(i) = sqrt(MolNum(i))*(log(MolNum(i))+FugCoef-d(i))

          !Calculating TM Hessian matrix
          do j = 1, comp
            !Numerical derivative
            MolNum(j) = MolNum(j) + dMol
            call VdW1fMIX (comp, a, b, BinIntCoef, MolNum, Samix, amix, bmix, i)

            call fugacity(T, P, a(i), b(i), amix, bmix, Samix, FugCoef, Trial, Vol, sigma_eos, epsilon_eos)
            
            Hess(j+(i-1)*comp) = FugCoef

            MolNum(j) = MolNum(j) - 2*dMol
            call VdW1fMIX (comp, a, b, BinIntCoef, MolNum, Samix, amix, bmix, i)

            call fugacity(T, P, a(i), b(i), amix, bmix, Samix, FugCoef, Trial, Vol, sigma_eos, epsilon_eos)
            
            Hess(j+(i-1)*comp) = Hess(j+(i-1)*comp) - FugCoef

            MolNum(j) = MolNum(j) + dMol
            Hess(j+(i-1)*comp) = Hess(j+(i-1)*comp)*(sqrt(MolNum(i)*MolNum(j))/(2.0d0*dMol))

            !Hessian term related to a Kronecker delta
            if(i .eq. j) then
              Hess(j+(i-1)*comp) = Hess(j+(i-1)*comp) + (1.0d0 + grad(i)/(2.0d0*sqrt(MolNum(i))))
            endif
          enddo
        enddo

        !Solving system of equations
        call LUdecomp(Hess,grad,step,comp)

        maxstep = 0.0d0
        do i = 1, comp
            !Updating compositions - Wi
            MolNum(i) = (sqrt(MolNum(i)) - step(i)/2.0d0)*(sqrt(MolNum(i)) - step(i)/2.0d0)

            !Summation of the steps squared - Convergence criteria
            maxstep = maxstep + step(i)*step(i)
        enddo
      enddo 
    endif

    if(it .eq. 51) print*, 'WARNING ON TPDliqvap: Maximum number of iterations reached. Exiting program...'

    !Calculating summation of compositions - WT
    TotMolTrial = 0.0d0
    do i = 1, comp
      TotMolTrial = TotMolTrial + MolNum(i)
    enddo

    !Searching for instability
    if(TotMolTrial .gt. 1.001) then
      print*, 'vapor-liquid stability analysis: Instability detected'
      do j = 1, comp
          CompositionTrial(j) = MolNum(j)
      enddo
      ret = 2
      phaseTrial = Trial

    !Searching for shadow regions
    elseif (TotMolTrial .lt. 0.999) then
      print*, 'vapor-liquid stability analysis: Shadow region detected.'
    endif

    end subroutine TPDliqvap   
    !===========================================================================================================================================================================


    !===========================================================================================================================================================================
    subroutine TPDliqliq (comp, T, P, a, b, dHf, Tf, dHtr, Ttr, MW, integraldCp1, integraldCp2, BinIntCoef, PhaseNum,&
                        phase, Composition, CompositionTrial, ret, sigma_eos, epsilon_eos)

    real(8) :: a(*), b(*), dHf(*), Tf(*), dHtr(*), Ttr(*), MW(*), integraldCp1(*), integraldCp2(*)
    real(8) :: T, P, Composition(comp,comp), CompositionTrial(comp), BinIntCoef(comp,comp), sigma_eos, epsilon_eos
    integer :: comp, PhaseNum, phase(*)

    !The reference phase must be stored in the last line of the Composition matrix.
    real(8) :: MolNum(comp,comp), d(comp), FugCoef(comp), grad(comp,comp), maxstep
    real(8) :: Hess(comp,comp*comp), SAME(comp), aux, TPDaux, step(comp), Vol
    real(8) :: dMol, tol, TPD(comp), TotMolTrial(comp), amix, bmix, Samix
    integer :: i, j, k, it, ret, trivial, liqRef(3), vap, flag    


        ret = 0
        dMol = 1.D-6
        tol = 1.D-6

        !In this routine, the reference phase isn't necessarily the last element of the vector "phase".
        !It's the liquid that already exists in the system.

        !Searching for an already existing liquid phase
        flag = 0
        do i = (PhaseNum-1+1), 1, -1
            if (phase(i).eq.1) then
                liqRef(flag+1) = i 
                flag = flag + 1
            endif
        enddo

        if(flag .eq. 0) liqRef(1) = PhaseNum-1

        !Searching for an already existing vapor phase
        do i = (PhaseNum-1+1), 1, -1 
            if(phase(i).eq.0) then
                vap = i
                exit
            endif
        enddo

        if(i.eq.-1) vap = 1000

        do i = 1, comp
            !Trial phase composition (number of moles) initial guess
            do j = 1, comp
                if (i.eq.j) then
                    MolNum(j,i) = 1
                else
                    MolNum(j,i) = 0.001
                endif
            enddo

            !Feed role on TPD calculation
            call VdW1fMIX (comp, a, b, BinIntCoef, Composition(1,liqRef(1)), Samix, amix, bmix, i)
            call fugacity(T, P, a(i), b(i), amix, bmix, Samix, FugCoef(i), phase(liqRef(1)), Vol, sigma_eos, epsilon_eos)                     

            d(i) = log(Composition(i,liqRef(1))) + FugCoef(i)
        enddo

        !SUCCESSIVE SUBSTITUTION///////////////////////////////////////////////////
        !Successive substitution for all trial phases
        do i = 1, comp
            it = 1
            !Main loop of successive substitution method
            do while (it.le.10 .and. maxstep.gt.tol)
                !Updating compositions (number of moles)///////////////////////////
                maxstep = 0.0d0
                do j = 1, comp
                    call VdW1fMIX (comp, a, b, BinIntCoef, MolNum(1,i), Samix, amix, bmix, j)
                    call fugacity(T, P, a(j), b(j), amix, bmix, Samix, FugCoef(j), 1, Vol, sigma_eos, epsilon_eos)
                    
                    step(j) = exp(d(j) - FugCoef(j))

                    if (dabs(MolNum(j,i)-step(j)).gt.maxstep) then
                        maxstep = dabs(MolNum(j,i)-step(j))
                    else
                        maxstep = maxstep
                    endif
                enddo

                do j = 1, comp
                    MolNum(j,i) = step(j)
                enddo
                !******************************************************************
                it = it + 1
            enddo !while(it <= 10 && maxstep > tol);


            !Modified tangent plane distance (TM)//////////////////////////////////
            TotMolTrial(i) = 0.0d0
            TPD(i) = 1.0d0
            do j = 1, comp
                TotMolTrial(i) = TotMolTrial(i) + MolNum(j,i)
                call VdW1fMIX (comp, a, b, BinIntCoef, MolNum(1,i), Samix, amix, bmix, j)
                call fugacity(T, P, a(j), b(j), amix, bmix, Samix, FugCoef(j), 1, Vol, sigma_eos, epsilon_eos)

                !Actually, it's being calculated the TM.
                TPD(i) = TPD(i) + MolNum(j,i)*(log(MolNum(j,i)) + FugCoef(j) - d(j) - 1.0d0)
            enddo

            TPD(i) = (TPD(i) - (1.0d0 - TotMolTrial(i) + TotMolTrial(i)*log(TotMolTrial(i))))/TotMolTrial(i)
            !printf("\n%d  %g  %5.4lf  %5.4lf\n",it-1,TPD[i],MolNum[i][0],MolNum[i][1]);DEBUG
        enddo
        !printf("\n%5.4lf  %5.4lf\n",Composition[liqRef[0]][0],Composition[liqRef[0]][1]);DEBUG
        !END OF SUCCESSIVE SUBSTITUTION********************************************

        !Calculating the distance of each trial phase from the already existing liquid phase
        if (vap.ne.1000) then
            do i = 1, comp
                SAME(i) = 0.0d0
                do j = 1, comp
                    SAME(i) = SAME(i) + dabs(MolNum(j,i) - Composition(j,(vap+1)))
                enddo
            enddo
        else 
            do i = 1, comp
                SAME(i) = 100.0d0
            enddo
        endif

        do k = 1, flag
            do i = 1, comp
                aux = 0.0d0
                do j = 1, comp
                    aux = aux + dabs(MolNum(j,i) - Composition(j,liqRef(1)))
                enddo

                if (aux.lt.SAME(i)) then
                    SAME(i) = aux
                else
                    SAME(i) = SAME(i)
                endif
            enddo
        enddo

        !Excluding trivial solutions
        trivial = 0
        do i = 1, comp
            if(SAME(i) .lt. 0.01d0) then
                trivial = trivial + 1
            else
                do j = 1, comp
                    MolNum(j,(i-trivial)) = MolNum(j,i)
                enddo
                TPD(i-trivial) = TPD(i)
            endif
        enddo

        if (comp .gt. trivial) then

          !NEWTON METHOD/////////////////////////////////////////////////////////
            if (it.eq.11) then
                tol = 1.D-10
                do i = 1, (comp-trivial)
                    it = 0
                    do while (maxstep.gt.tol .and. it.le.50)
                        it = it + 1
                        do j = 1, comp
                            !Gradient
                            call VdW1fMIX (comp, a, b, BinIntCoef, MolNum(1,i), Samix, amix, bmix, j)
                            call fugacity(T, P, a(j), b(j), amix, bmix, Samix, FugCoef(j), 1, Vol, sigma_eos, epsilon_eos)
                            grad(j,i) = sqrt(MolNum(j,i))*(log(MolNum(j,i))+FugCoef(j)-d(j))

                            !Hessian
                            do k = 1, comp
                                MolNum(k,i) = MolNum(k,i) + dMol
                                call VdW1fMIX (comp, a, b, BinIntCoef, MolNum(1,i), Samix, amix, bmix, j)
                                call fugacity(T, P, a(j), b(j), amix, bmix, Samix, FugCoef(j), 1, Vol, sigma_eos, epsilon_eos)
                                Hess((k+(j-1)*comp),i) = FugCoef(j)

                                MolNum(k,i) = MolNum(k,i) - 2.0d0*dMol
                                call VdW1fMIX (comp, a, b, BinIntCoef, MolNum(1,i), Samix, amix, bmix, j)
                                call fugacity(T, P, a(j), b(j), amix, bmix, Samix, FugCoef(j), 1, Vol, sigma_eos, epsilon_eos) 
                                Hess((k+(j-1)*comp),i) = Hess((k+(j-1)*comp),i) - FugCoef(j)

                                MolNum(k,i) = MolNum(k,i) + dMol
                                Hess((k+(j-1)*comp),i) = Hess((k+(j-1)*comp),i)*(sqrt(MolNum(j,i)*MolNum(k,i))/(2.0d0*dMol))

                                if (j.eq.k) then
                                  Hess((k+(j-1)*comp),i) = Hess((k+(j-1)*comp),i) +&
                                                             (1.0d0 + grad(j,i)/(2.0d0*sqrt(MolNum(j,i))))
                                endif
                            enddo
                        enddo

                        !Calculating the inverse of the Hessian
                        call LUdecomp(Hess(:,i),grad(:,i),step,comp)

                        maxstep = 0.0d0
                        do j = 1, comp
                            !Updating variables
                            MolNum(j,i) = (sqrt(MolNum(j,i)) - step(j)/2.0d0)*(sqrt(MolNum(j,i)) - step(j)/2.0d0)

                            !Selecting maximum step to analyze the convergence of the method
                            maxstep = maxstep + step(j)*step(j)
                        enddo
                    enddo !while(maxstep > tol && it <= 50);

                    if(it.eq.51) print*, 'WARNING ON TPDliqliq: Maximum number of iterations reached. Exiting program...'

                    !Modified tangent plane distance (TM)//////////////////////////
                    TotMolTrial(i) = 0.0d0
                    TPD(i) = 1.0d0
                    do j = 1, comp
                        !Total number of moles////////
                        TotMolTrial(i) = TotMolTrial(i) + MolNum(j,i)

                        call VdW1fMIX (comp, a, b, BinIntCoef, MolNum(1,i), Samix, amix, bmix, j)
                        call fugacity(T, P, a(j), b(j), amix, bmix, Samix, FugCoef(j), 1, Vol, sigma_eos, epsilon_eos)
                                 
                        !Modified tangent plane distance (TM)
                        !At this point, the variable "TPD" is the TM.//////////////
                        TPD(i) = TPD(i) + MolNum(j,i)*(log(MolNum(j,i)) + FugCoef(j) - d(j) - 1.0d0)
                    enddo

                    !TPD calculation
                    TPD(i) = (TPD(i) - (1.0d0 - TotMolTrial(i) + TotMolTrial(i)*log(TotMolTrial(i))))/TotMolTrial(i)

                    !printf("\n%d  %g  %5.4lf  %5.4lf\n",it,TPD[i],MolNum[i][0],MolNum[i][1]);DEBUG
                enddo  

            endif
            !END OF NEWTON METHOD**************************************************



            TPDaux = -1.D-5
            !Searching for instability/////////////////////////////////////////////////
            do i  =  1, (comp-trivial)
              if (TPD(i).lt.TPDaux) then
                  flag = i
                  TPDaux = TPD(i)
              endif
            enddo

            if (TPDaux .lt. -1.D-5) then
              print*, 'liquid-liquid stability analysis: Instability detected!'

              do k = 1, comp
                  CompositionTrial(k) = MolNum(k,flag)
              enddo

              ret = 2
            endif
        endif

    return    

    end subroutine TPDliqliq    
    !==============================================================================================================================================================================


    !==============================================================================================================================================================================
    subroutine TPDliqsol(comp,T,P, a, b, Acentric, Tc, dHf, Tf, dHtr, Ttr, MW, SCN, integraldCp1,integraldCp2,BinIntCoef, PhaseNum,&
                        phase, Composition, CompositionTrial, phaseTrial, TPDout, ret, Uniquac_ri, Uniquac_qi, solidRef,MetCalcLiq,&
                        MetCalcSolid, sigma_eos, epsilon_eos, vki)


      real(8) :: a(*), b(*), Acentric(*), dHf(*), Tf(*), dHtr(*), Ttr(*), MW(*), integraldCp1(*), integraldCp2(*), Tc(*)
      real(8) :: T, P, Composition(comp,comp), CompositionTrial(comp), TPDout, BinIntCoef(comp,comp), Uniquac_ri(*), Uniquac_qi(*)
      integer :: comp, PhaseNum, phase(*), phaseTrial, SCN(*), vki(comp,10)

      real(8) :: MolNum(comp,comp), d(comp), FugCoef(comp), grad(comp,comp), maxstep, FugCoefrefpuroliq(comp)    

      real(8) :: Hess(comp*comp,comp), SAME(comp), aux, TPDaux, step(comp), Vol, TotMolTrialAux, amix, bmix, Samix, sigma_eos, epsilon_eos

      real(8), parameter :: infinite = huge(1.d0)

      real(8) :: dMol, tol, TotMolTrial(20), fugcoefaux, FugCoefrefpuro(comp), gamaaux, SAMEaux, gamaliq(comp), dgamadxk

      integer :: i, j, k, it(comp), ret, trivial, liqRef(3), flag, solidRef, TrialNum, SolidNum, MetCalcLiq, MetCalcSolid

      real(8) :: Pref, Poynting, R2, a4, a5, TermoDeltaCp1, TermoDeltaCp2, R
      real(8), parameter, dimension(3) :: acoef = [-0.05359702d0, 0.00010981d0, 0.00010981d0]      
      real(8), parameter, dimension(3) :: bcoef = [1.37946361d0, 0.05667426d0, 0.95076333d0]         
      real(8), parameter, dimension(3) :: ccoef = [-11.92772402d0, 1.65715239d0, 1.65715239d0]        
      real(8), parameter, dimension(3) :: ecoef = [390.86109841d0, 11.81282730d0, 348.05910146d0]        
      

      ret = 0
      dMol = 1.D-6
      tol = 1.D-6
      MolNum = 0.0d0
      maxstep = 10000d0
      CompositionTrial = 0.0d0
      FugCoefrefpuroliq = 0.0d0
      FugCoefrefpuro = 0.0d0
      gamaliq = 0.0d0
      Pref = 1.0d0
      Poynting = 0.0d0
      
      R2 = 83.14462175d0 !R2 = cm³.bar/(mol.K)
      R = 1.9858775d0 !cal/(mol*K)      
      a4 = 0.3033d0 
      a5 = -4.635d-4
        
      !In this routine, the reference phase isn't necessarily the last element of the vector "phase".
      !It's the liquid that already exists in the system.

      !liqRef(1) = PhaseNum-1

      TrialNum = comp - solidRef + 1
      k = solidRef
      do i = 1, TrialNum    

          !Trial phase composition (number of moles) initial guess
          do j = solidRef, comp

              MolNum(j,i) = 1.d-6
    !            MolNum(i,j) = 0.001                
          enddo

          MolNum(k,i) = 1.0d0
    !        MolNum(i,j) = 1         
          k = k + 1

      enddo

      do i = solidRef, comp    
          !Feed role on TPD calculation

        if (MetCalcLiq.eq.3) then
            call fugacity(T, P, a(i), b(i), a(i), b(i), 2.0d0*a(i), FugCoefrefpuroliq(i), phase(PhaseNum), Vol, sigma_eos, epsilon_eos)                       

            call FloryUnifac (comp, T, P, Composition(1,PhaseNum), gamaliq(i), SCN, vki, i)

            d(i) = log(Composition(i,PhaseNum)) + FugCoefrefpuroliq(i) + gamaliq(i)    !Flory + Unifac

        else
            call VdW1fMIX (comp, a, b, BinIntCoef, Composition(1,PhaseNum), Samix, amix, bmix, i)

            call fugacity(T, P, a(i), b(i), amix, bmix, Samix, FugCoef(i), phase(PhaseNum), Vol, sigma_eos, epsilon_eos)            

            d(i) = log(Composition(i,PhaseNum)) + FugCoef(i) 

        endif


      enddo


      do j = solidRef, comp
          
          call fugacity(T, P, a(j), b(j), a(j), b(j), 2.0d0*a(j), FugCoefrefpuro(j), phase(PhaseNum), Vol, sigma_eos, epsilon_eos)
          
          if (SCN(j).le.20) then
              
              Poynting = (1.0d0/(R2*T))*(P-Pref)*(-(acoef(1)*SCN(j)+bcoef(1))*T + (ccoef(1)*SCN(j)+ecoef(1)))
              
          elseif (SCN(j).gt.20 .and. SCN(j).le.34) then
              
              Poynting = (1.0d0/(R2*T))*(P-Pref)*(-(acoef(2)*SCN(j)+bcoef(2))*T + (ccoef(2)*SCN(j)+ecoef(2)))
              
          elseif (SCN(j).gt.34) then
              
              Poynting = (1.0d0/(R2*T))*(P-Pref)*(-(acoef(3)*SCN(j)+bcoef(3))*T + (ccoef(3)*SCN(j)+ecoef(3)))
              
          endif
          

            TermoDeltaCp1 = (1.0d0/(R*T))*((a4*MW(j)*T+(a5*MW(j)*(T**2)/2.0d0)) -&
                            (a4*MW(j)*Tf(j)+(a5*MW(j)*(Tf(j)**2)/2.0d0)))


            TermoDeltaCp2 = (1.0d0/R)*((a4*MW(j)*log(T)+a5*MW(j)*T) - (a4*MW(j)*log(Tf(j))+a5*MW(j)*Tf(j)))



            FugCoefrefpuro(j) = FugCoefrefpuro(j) + ((-dHf(j)*(1.0d0/T - 1.0d0/Tf(j)) - dHtr(j)*(1.0d0/T - 1.0d0/Ttr(j)))/R -&
                                                     Poynting - TermoDeltaCp1 + TermoDeltaCp2)

          !FugCoefrefpuro(j) = FugCoefrefpuro(j) + (-dHf(j)*(1.0d0/T - 1.0d0/Tf(j)) - dHtr(j)*(1.0d0/T - 1.0d0/Ttr(j)))/R
          
!          FugCoefrefpuro(j) = FugCoefrefpuro(j) + (-dHf(j)*(1.0d0/T - 1.0d0/Tf(j)) - dHtr(j)*(1.0d0/T - 1.0d0/Ttr(j)) + &
!                              MW(j)*(integraldCp1(j) + integraldCp2(j)/T + 0.3033*(log(T) - 1.) + 0.5*(-4.635E-4)*T))/R
          
      enddo

      !SUCCESSIVE SUBSTITUTION///////////////////////////////////////////////////
      !Successive substitution for all trial phases
      do i = 1, TrialNum

          it(i) = 1
          !Main loop of successive substitution method
          do while (it(i).le.10 .and. maxstep.gt.tol)

              !Updating compositions (number of moles)///////////////////////////
              maxstep = 0.0d0
              do j = solidRef, comp

                  if (MetCalcSolid.eq.1) then
                      call Uniquac (comp, T, P, Tc, Acentric, dHf, dHtr, Uniquac_ri, Uniquac_qi, MolNum(1,i), gamaaux, solidRef, j)
                  else
                      gamaaux = 0.0d0
                  endif

                  FugCoef(j) = FugCoefrefpuro(j) + gamaaux

                  step(j) = exp(d(j) - FugCoef(j))

                  if (dabs(MolNum(j,i)-step(j)).gt.maxstep) then
                      maxstep = dabs(MolNum(j,i)-step(j))
    !                if (dabs(MolNum(i,j)-step(j)).gt.maxstep) then
    !                    maxstep = dabs(MolNum(i,j)-step(j))                    
                  endif

              enddo

              do j = solidRef, comp
                  MolNum(j,i) = step(j)
    !                MolNum(i,j) = step(j)                
              enddo
              !******************************************************************
              it(i) = it(i) + 1

          enddo !while(it <= 10 && maxstep > tol);


          !Modified tangent plane distance (TM)//////////////////////////////////
          TotMolTrial(i) = 0.0d0

          do j = solidRef, comp
              TotMolTrial(i) = TotMolTrial(i) + MolNum(j,i)
    !            TotMolTrial(i) = TotMolTrial(i) + MolNum(i,j)                    
          enddo

          maxstep = 10000d0

      enddo
      !printf("\n%5.4lf  %5.4lf\n",Composition[liqRef[0]][0],Composition[liqRef[0]][1]);DEBUG
      !END OF SUCCESSIVE SUBSTITUTION********************************************

      !Calculating the distance of each trial phase from the already existing liquid phase

      do i = 1, TrialNum

          SAME(i) = 1000.0d0

          do j = 1, (PhaseNum-1)

              SAMEaux = 0.0d0
              do k = solidRef, comp                
                  SAMEaux = SAMEaux + dabs(MolNum(k,i) - Composition(k,j))        
              enddo

              if (SAMEaux.lt.SAME(i)) SAME(i) = SAMEaux


          enddo

      enddo


      !Excluding trivial solutions
      trivial = 0

      do i = 1, TrialNum

          if(SAME(i) .lt. 0.01d0) then
              trivial = trivial + 1
          elseif (trivial.ne.0) then
              do j = solidRef, comp
                  MolNum(j,(i-trivial)) = MolNum(j,i)
    !                MolNum((i-trivial),j) = MolNum(i,j)
              enddo
              TotMolTrial(i-trivial) = TotMolTrial(i)
              it(i-trivial) = it(i)

          endif

      enddo

      TrialNum = TrialNum - trivial

      trivial = 0
      if (TrialNum.gt.0) then



          !NEWTON METHOD/////////////////////////////////////////////////////////

          tol = 1.D-7
          do i = 1, TrialNum

              if (it(i).eq.11) then    

                  it(i) = 0

                  do while (maxstep.gt.tol .and. it(i).le.20)

                      it(i) = it(i) + 1
                      do j = 1, (comp-solidRef+1)
    !                    do j = solidRef, comp

                          !Gradient                                      
                          if (MetCalcSolid.eq.1) then
                              call Uniquac (comp, T, P, Tc, Acentric, dHf, dHtr,Uniquac_ri,Uniquac_qi, MolNum(1,i), gamaaux,&
                                            solidRef,(j+solidRef-1))
                          else
                              gamaaux = 0.0d0
                          endif

                          FugCoef(j+solidRef-1) = FugCoefrefpuro(j+solidRef-1) + gamaaux

                          grad(j,i) = sqrt(MolNum(j+solidRef-1,i))*(log(MolNum(j+solidRef-1,i))+FugCoef(j+solidRef-1)&
                                            -d(j+solidRef-1))
    !                        grad(i,j) = sqrt(MolNum(i,j))*(log(MolNum(i,j))+FugCoef(j)-d(j))                        

                          !Hessian
                          do k = 1, (comp-solidRef+1)

                            if (MetCalcSolid.eq.1) then  
                              call DerivadaUniquac (comp, T, P, Tc, Acentric, dHf, dHtr, Uniquac_ri, Uniquac_qi, MolNum(1,i),&
                                                    dgamadxk, solidRef, (j+solidRef-1),(k+solidRef-1))
                            else
                              dgamadxk = 0.0d0  
                            endif

                            Hess((k+(j-1)*(comp-solidRef+1)),i) = dgamadxk*(sqrt(MolNum(j+solidRef-1,i)*&
                                                              MolNum(k+solidRef-1,i)))
                                                              

                            if (j.eq.k) then
                              Hess((k+(j-1)*(comp-solidRef+1)),i) = Hess((k+(j-1)*(comp-solidRef+1)),i) +&
                                                             (1.0d0 + grad(j,i)/(2.0d0*sqrt(MolNum(j+solidRef-1,i))))
                                                           
                            endif

                          enddo

                      enddo

                      !Calculating the inverse of the Hessian
                      call LUdecomp(Hess(:,i),grad(:,i),step,(comp-solidRef+1))

                      maxstep = 0.0d0
                      do j = 1, (comp-solidRef+1)

                          !Updating variables
                          MolNum(j+solidRef-1,i) = (sqrt(MolNum(j+solidRef-1,i)) - step(j)/2.0d0)*(sqrt(MolNum(j+solidRef-1,i))&
                                                   - step(j)/2.0d0)
    !                        MolNum(i,j) = (sqrt(MolNum(i,j)) - step(j)/2.0d0)*(sqrt(MolNum(i,j)) - step(j)/2.0d0)                        

                          !Selecting maximum step to analyze the convergence of the method
                          maxstep = maxstep + step(j)*step(j)

                      enddo

                  enddo    !while(maxstep > tol && it <= 50);

                  if(it(i).eq.21) print*, 'WARNING ON TPDliqsol: Maximum number of iterations reached. Exiting program...'

                  !Modified tangent plane distance (TM)//////////////////////////
                  TotMolTrial(i) = 0.0d0

                  do j = solidRef, comp

                      !Total number of moles////////
                      TotMolTrial(i) = TotMolTrial(i) + MolNum(j,i)
    !                    TotMolTrial(i) = TotMolTrial(i) + MolNum(i,j)                                       

                  enddo

                  !Calculating the distance of each trial phase from the already existing liquid phase

                  SAME(i) = 1000.0d0

                  do j = 1, (PhaseNum-1)

                      SAMEaux = 0.0d0
                      do k = solidRef, comp                
                          SAMEaux = SAMEaux + dabs(MolNum(k,i) - Composition(k,j))        
                      enddo

                      if (SAMEaux.lt.SAME(i)) SAME(i) = SAMEaux

                  enddo

                  !Excluding trivial solutions

                  if(SAME(i) .lt. 0.01d0) then
                      trivial = trivial + 1
                  elseif (trivial.ne.0) then
                      do j = solidRef, comp
                          MolNum(j,(i-trivial)) = MolNum(j,i)
          !                MolNum((i-trivial),j) = MolNum(i,j)  
                      enddo
                          TotMolTrial(i-trivial) = TotMolTrial(i)
                  endif


              else

                  do j = solidRef, comp
                      MolNum(j,(i-trivial)) = MolNum(j,i)
      !                MolNum((i-trivial),j) = MolNum(i,j)  
                  enddo                
                      TotMolTrial(i-trivial) = TotMolTrial(i)
              endif    

          enddo  

          !END OF NEWTON METHOD**************************************************

          TrialNum = TrialNum - trivial


          if ((TrialNum).gt.0) then
              phaseTrial = 2
              flag = 10000
              TotMolTrialAux  = 1.001d0
              !TotMolTrialAux  = 1.1d0              
              !Searching for instability/////////////////////////////////////////////////
              do i  =  1, (TrialNum)
                  if (TotMolTrial(i).gt.TotMolTrialAux) then
                      flag = i
                      TotMolTrialAux = TotMolTrial(i)
                  endif
              enddo

              if (flag.ne.10000) then

                  print*, 'solid-liquid stability analysis: Instability detected!'

                  aux = 0.0d0
                  do k = solidRef, comp
                      aux = aux + MolNum(k,flag)
                  enddo

                  do k = solidRef, comp
                      CompositionTrial(k) = MolNum(k,flag)/aux
      !                    CompositionTrial(k) = MolNum(flag,k)                    
                  enddo

                  ret = 2

              endif
              !**********************************************************************

    !            if (ret.eq.0) then
    !
    !                do i = 1, (TrialNum)
    !
    !                    if (TotMolTrial(i).gt.0.9999d0 .and. TotMolTrial(i).lt.1.00001d0) then
    !
    !                    !Searching for saturation points///////////////////////////////
    !                        ret = 2
    !
    !                        print*, 'solid-liquid stability analysis: Saturation point detected.'
    !                        exit
    !
    !                    endif
    !
    !                enddo
    !
    !            endif

              do i = solidRef, comp

                  if (CompositionTrial(i).ne.CompositionTrial(i) .or. CompositionTrial(i).eq.Infinite) then
                      ret = -1
                      exit
                  endif

              enddo

          endif

      endif



    return    

    end subroutine TPDliqsol    

    
!==============================================================================================================================================================================    
   
    subroutine TPDmultisolid(comp, T, P, a, b, Acentric, Tc, dHf, Tf, dHtr, Ttr, MW, SCN, integraldCp1,integraldCp2,BinIntCoef,&
                              PNA, PhaseNum, phase, Composition, CompositionTrial, phaseTrial, TPDout, ret, sigma_eos, epsilon_eos,&
                              solidRef)


    real(8) :: a(*), b(*), Acentric(*), dHf(*), Tf(*), dHtr(*), Ttr(*), MW(*), integraldCp1(*), integraldCp2(*), Tc(*)
    real(8) :: T, P, Composition(comp,comp), CompositionTrial(comp), TPDout, BinIntCoef(comp,comp), sigma_eos, epsilon_eos
    integer :: comp, PhaseNum, phase(*), phaseTrial, PNA(*), SCN(*), solidRef

    real(8), parameter :: R = 1.9858775d0
    real(8) :: TPDprecipitate(comp), TPDlower, FugCoef, trial(comp,comp), Vol, amix, bmix, Samix
    integer :: i, j, k, solidIndex, solid(comp), precipitate(comp), flagNoSolid, ret, flag
    
    solidIndex = 1
    ret = 0
    flagNoSolid = 100000
    solid = 0.0d0
    flag = flagNoSolid
    trial = 0.0d0
    precipitate = 0
    TPDprecipitate = 0.0d0
    
    !Initializing variables
    do i = solidRef, comp
        precipitate(i) = flagNoSolid
        TPDprecipitate(i) = real(flagNoSolid)
    enddo    
    TPDlower = real(flagNoSolid)    
    
    
    !Finding already precipitated compounds
    solidIndex = 0
    do j = 1, PhaseNum    
        if (phase(j) .eq. 2) then
           do i = comp, solidRef, -1
               if ((Composition(i,j) - 1.0d0) .gt. -0.999d0) then
                   solidIndex = solidIndex + 1
                   solid(solidIndex) = i                   
                   exit
               endif
           enddo
        endif
    enddo    
    
    
    !Determining which are the next compounds to precipitate
    k = 0
    do i = comp, solidRef, -1
        
        flag = flagNoSolid
               
        if (PNA(i).eq.1) then
            
            do j = 1, solidIndex
                if (solid(j).eq.i) then
                    flag = 1
                    exit
                endif            
            enddo
            
            if (flag.ne.1) then 
                k = k + 1
                precipitate(k) = i
                do j = solidRef, comp
                    trial(j,k) = 0.0d0
                enddo
                trial(i,k) = 1.0d0
                
            endif
            
        endif
        
        if (k.eq.comp) exit
        
    enddo
    
    do i = 1, comp
        CompositionTrial(i) = 0.0d0
    enddo
    
    
    !stability analysis
    flag = flagNoSolid
    if (k.ne.0) then
    
        do i = 1, k
            !Calculates the natural logarithm of pure solid fugacity coefficient of the component "precipitate(i)"
            call fugacity(T, P, a(precipitate(i)), b(precipitate(i)), a(precipitate(i)), b(precipitate(i)), 2.0d0*a(precipitate(i)), FugCoef, 2, Vol, sigma_eos, epsilon_eos)

            FugCoef = FugCoef + (-dHf(precipitate(i))*(1.0d0/T - 1.0d0/Tf(precipitate(i))) - dHtr(precipitate(i))*(1.0d0/T - 1.0d0/Ttr(precipitate(i))))/R 
            
            TPDprecipitate(i) = FugCoef + log(P)            
                        
            !Calculates the natural logarithm of reference phase fugacity coefficient of the component "precipitate(i)"
            call VdW1fMIX (comp, a, b, BinIntCoef, Composition(1,PhaseNum), Samix, amix, bmix, precipitate(i))

            call fugacity(T, P, a(precipitate(i)), b(precipitate(i)), amix, bmix, Samix, FugCoef, phase(PhaseNum), Vol, sigma_eos, epsilon_eos)
            
            TPDprecipitate(i) = TPDprecipitate(i) - (FugCoef + log(Composition(precipitate(i),PhaseNum)*P))           
            
            
            if (TPDprecipitate(i).lt.TPDlower .and. TPDprecipitate(i).lt.0.0d0) then            
                TPDlower = TPDprecipitate(i)
                flag = precipitate(i)
                CompositionTrial(flag) = 1.0d0
                exit
            endif
            
            if (TPDprecipitate(i).ne.TPDprecipitate(i)) then
                ret = -1 
                exit
            endif
            
        enddo
        
        !CompositionTrial(flag) = 1.0d0        
    endif              
    
    !Searching for negative TPDs
    if (TPDlower.lt.0.0d0 .and. ret.ne.-1) then
        TPDout = TPDlower
        phaseTrial = 2                
        ret = 2
        print*, 'multisolid stability analysis: Instability detected!'
    endif
    
    end subroutine TPDmultisolid      
    
!==============================================================================================================================================================================    
    
    end module TPD      
