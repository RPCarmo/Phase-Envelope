

module SuccessiveSubst

    use EOS

implicit none

    subroutine  SuccessiveSubstitution(Ncomp, Nphase, Nphase_inc, T, P, acentric, Tc, ac, composition, beta, K, a, b, kij lij)
      !Arguments
      real(8) :: composition(:,:), K(:,:), beta(:), a(:), b(:), kij(:,:), lij(:,:)
      integer :: Ncomp
      
      !Local
      real(8), parameter :: tol = 1.d-5, diff = 1.d-6
      integer, parameter :: maxit = 100
      real(8) :: T_old, aux, F(:), dF(:), step(:), bmix(2), amix(2), Volume(2), FugCoef_ref, FugCoef_aux
      integer :: it, i
      
      step(1) = 1.0d+6
      it = 0
      do while(dabs(step(1)) .gt. tol .and. it .lt. maxit)
          it = it + 1
          T_old = T
          
          do j = 2, Nphase
            composition(:,j) = composition(:,1)*K(:,j)
            F(j-1) = sum(composition(:,j)) - 1.d0 !Residuals
          enddo
          
          !Numerical Derivative With Respect to Temperature
          T = T_old + diff
          call EoS_param(acentric,Tc,ac,a,T,Ncomp) !Updating Attractive Parameter
          call VdW1fMIX (Ncomp,a,b,kij,lij,composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
          call VdW1fMIX (Ncomp,a,b,kij,lij,composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), 0)
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), 1)
          do i = 1,Ncomp
            call fugacity(Ncomp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),composition(:,1),kij(i,:),lij(i,:),i)
            call fugacity(Ncomp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),composition(:,2),kij(i,:),lij(i,:),i)
            dF(1) = dF(1) + composition(i,1)*K(i)*(FugCoef_ref - FugCoef_aux)
          enddo
          
          T = T_old - diff
          call EoS_param(acentric,Tc,ac,a,T,Ncomp) !Updating Attractive Parameter
          call VdW1fMIX (Ncomp,a,b,kij,lij,composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
          call VdW1fMIX (Ncomp,a,b,kij,lij,composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), 0)
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), 1)
          do i = 1,Ncomp
            call fugacity(Ncomp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),composition(:,1),kij(i,:),lij(i,:),i)
            call fugacity(Ncomp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),composition(:,2),kij(i,:),lij(i,:),i)
            dF(1) = dF(1) - composition(i,1)*K(i)*(FugCoef_ref - FugCoef_aux)
          enddo
        
          dF(1) = dF(1)/(2.0d0*diff)
          
          !Temperature Step Calculation
          step(1) = F(1)/dF(1)
          
          !Step Brake
          if(dabs(step(1)) .gt. 0.25d0*T_old) step(1) = 0.25d0*T_old*step(1)/dabs(step(1)) 
          
          !Updating Temperature
          T = T_old - step(1)
          
          
          !Updating K-factors
          call EoS_param(acentric,Tc,ac,a,T,Ncomp)
          call VdW1fMIX (Ncomp,a,b,kij,lij,composition(:,1),amix(1),bmix(1)) !Mixing Rule - Reference Phase
          call VdW1fMIX (Ncomp,a,b,kij,lij,composition(:,2),amix(2),bmix(2)) !Mixing Rule - Incipient Phase
          call EoS_Volume(P, T, bmix(1), amix(1), Volume(1), 0)
          call EoS_Volume(P, T, bmix(2), amix(2), Volume(2), 1)
          do i = 1,Ncomp
            call fugacity(Ncomp,T,P,a,b,amix(1),bmix(1),FugCoef_ref,Volume(1),composition(:,1),kij(i,:),lij(i,:),i)
            call fugacity(Ncomp,T,P,a,b,amix(2),bmix(2),FugCoef_aux,Volume(2),composition(:,2),kij(i,:),lij(i,:),i)
            K(i) = dexp(FugCoef_ref - FugCoef_aux)
          enddo
      enddo
      
      if(it .eq. maxit .and. dabs(step(1)) .gt. tol) then
          print*, "WARNING: In Successive Substitution Method - Maximum Number of Iterations Reached!"
          print*, "Exiting Program..."
          stop
      endif
     
    end subroutine

end module
