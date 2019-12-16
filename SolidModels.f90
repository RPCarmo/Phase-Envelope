! ATOMS - Applied Thermodynamics and Molecular Simulation
!====================================================================
!
! Nome do arquivo:  Uniquac.for
! Projeto:          Precipitação de parafinas
! Conteúdo:         Implementação de modelos para fase sólida: UNIQUAC
!
!====================================================================

module SolidModels

  implicit none

  contains

  subroutine Uniquac (comp, T, P, Tc, Acentric, dHf, dHtr, Uniquac_ri, Uniquac_qi, MoleFrac, gama, solidRef, index)
     
    !Arguments
    real(8) :: T, P, Tc(*), Acentric(*), dHf(*), dHtr(*), Uniquac_ri(*), Uniquac_qi(*), MoleFrac(comp), gama
    integer :: solidRef, index, comp
    
    !Local    
    real(8) :: Uniquac_dHvap, Uniquac_dH0, Uniquac_dH1, Uniquac_dH2, Uniquac_dHsub(comp), MoleFracAux(comp), sumMole 
    real(8) :: R, sumphi, sumteta, x, phi(comp), teta(comp), lambda(comp,comp),  lj(comp)
    real(8) :: termo1, termo2, termo3, termo4, termo5, termo6, termo7, sumlj, Tau(comp,comp), somatorio, den, sumTauij, Z
    integer :: i, j, k
    real(8) :: Param
    
    !Param = 1.3730468750000000d-4 !Bim 0
    !Param = -2.4511718749999998d-3 !Bim 3
    !Param = -2.0710937500000005d-3 !Bim 5
    !Param = -1.2937304687499976d-2 !Bim 9
    !Param = -4.7722656250000000d-3 !Bim 13
    Param = 0.0d0
    
    
    !Initializing variables
    lambda = 0.0d0
    Tau = 0.0d0
    R = 1.9858775d0 !cal/(mol*K)
    Z = 10.0d0

    sumMole = 0.0d0
    do i = solidRef,comp
      sumMole  = sumMole + MoleFrac(i)
    enddo

    sumphi = 0.0d0
    sumteta = 0.0d0
    sumlj = 0.0d0
    do i = solidRef,comp
      MoleFracAux(i) = MoleFrac(i)/sumMole
      sumphi = sumphi + MoleFracAux(i)*Uniquac_ri(i)
      sumteta = sumteta + MoleFracAux(i)*Uniquac_qi(i)  
      lj(i) = (Z/2.0d0)*(Uniquac_ri(i)-Uniquac_qi(i))-(Uniquac_ri(i)-1.0d0)    
      sumlj = sumlj + MoleFracAux(i)*lj(i)
    enddo
    
    do i = solidRef, comp
      !Calculating vaporization enthalpy   
      x = 1.d0 - T/Tc(i)
      Uniquac_dH0 = 5.2804d0*x**0.3333d0 + 12.8650d0*x**0.8333d0 + 1.1710d0*x**1.2083d0 - 13.1160d0*x + 0.4858d0*x**2 - 1.0880d0*x**3  
      Uniquac_dH1 = 0.80022d0*x**0.3333d0 + 273.23d0*x**0.8333d0 + 465.08d0*x**1.2083d0 - 638.51d0*x - 145.12d0*x**2 + 74.049d0*x**3
      Uniquac_dH2 = 7.2543d0*x**0.3333d0 - 346.45d0*x**0.8333d0 - 610.48d0*x**1.2083d0 + 839.89d0*x + 160.05d0*x**2 - 50.711d0*x**3
      Uniquac_dHvap = Uniquac_dH0 + Acentric(i)*Uniquac_dH1 + Acentric(i)*Acentric(i)*Uniquac_dH2
      Uniquac_dHvap = Uniquac_dHvap*R*Tc(i)
      
      !Calculating sublimation enthalpy
      Uniquac_dHsub(i) = Uniquac_dHvap + dHf(i) + dHtr(i)
      
      !Calculating interaction parameter (lambda)    
      do j = i, comp                    
        lambda(i,j) = (-2.0d0/Z)*(Uniquac_dHsub(i) - R*T)
        lambda(j,i) = lambda(i,j)
      enddo
        
      phi(i) = (MoleFracAux(i)*Uniquac_ri(i))/sumphi
      teta(i) = (MoleFracAux(i)*Uniquac_qi(i))/sumteta
    enddo

    do i = solidRef, comp    
      do j = solidRef, comp                    
        Tau(i,j) = exp(-(lambda(i,j)*(1.0d0 + Param)-lambda(i,i))/(Uniquac_qi(i)*R*T))
      enddo 
    enddo
    
    !Calculating first term of UNIQUAC equation
    termo1 = 0.0d0
    termo1 = log(Uniquac_ri(index)/sumphi)
    
    !Calculating second term of UNIQUAC equation
    termo2 = 0.0d0
    termo2 = (Z/2)*Uniquac_qi(index)*log((Uniquac_qi(index)*sumphi)/(Uniquac_ri(index)*sumteta))
    
    !Calculating third term of UNIQUAC equation
    termo3 = 0.0d0
    termo3 =  lj(index)
    
    !Calculating fourth term of UNIQUAC equation
    termo4 = 0.0d0
    termo4 = (Uniquac_ri(index)/sumphi)*sumlj
    
    !Calculating fifth term of UNIQUAC equation
    termo5 = 0.0d0
    sumTauij = 0.0d0
    do j = solidRef,comp
      sumTauij  = sumTauij + teta(j)*Tau(index,j)      
    enddo
    termo5 = Uniquac_qi(index)*log(sumTauij)
    
    !Calculating sixth term of UNIQUAC equation
    termo6 = 0.0d0
    termo6 = Uniquac_qi(index)
    
    !Calculating seventh term of UNIQUAC equation
    termo7 = 0.0d0
    somatorio = 0.0d0
    do j = solidRef,comp
      den = 0.0d0
      do k = solidRef, comp                    
        den = den + teta(k)*Tau(j,k)    
      enddo     
      somatorio = somatorio + (teta(j)*Tau(j,index))/den  
    enddo
    termo7 = Uniquac_qi(index)*somatorio
    
    !Calculation of activity coefficient
    gama = termo1 + termo2 + termo3 - termo4 - termo5 + termo6 - termo7
  
    
  end subroutine

!===================================================================================================================================
  
  subroutine DerivadaUniquac (comp, T, P, Tc, Acentric, dHf, dHtr, Uniquac_ri, Uniquac_qi, MoleFrac, dlngamadxk, solidRef,&
                              index, indexk)
     
    !Arguments
    real(8) :: T, P, Tc(*), Acentric(*), dHf(*), dHtr(*), Uniquac_ri(*), Uniquac_qi(*), MoleFrac(comp), dlngamadxk
    integer :: solidRef, index, indexk, comp
    
    !Local    
    real(8) :: Uniquac_dHvap, Uniquac_dH0, Uniquac_dH1, Uniquac_dH2, Uniquac_dHsub(comp), MoleFracAux(comp), sumMole,termo7parcela2 
    real(8) :: R, sumphi, sumteta, x, phi(comp), teta(comp), lambda(comp,comp),  lj(comp), sumtetaTji, sumtetaTmk, sumtetaTmj
    real(8) :: termo1, termo2, termo3, termo4, termo5, termo6, termo7, sumlj, Tau(comp,comp), somatorio, den, sumTauij, Z
    integer :: i, j, k, m, n
    real(8) :: Param
    
    !Param = 1.3730468750000000d-4 !Bim 0
    !Param = -2.4511718749999998d-3 !Bim 3
    !Param = -2.0710937500000005d-3 !Bim 5
    !Param = -1.2937304687499976d-2 !Bim 9
    !Param = -4.7722656250000000d-3 !Bim 13
    Param = 0.0d0
    
    !Initializing variables
    lambda = 0.0d0
    Tau = 0.0d0
    R = 1.9858775d0 !cal/(mol*K)
    Z = 10.0d0

    sumMole = 0.0d0
    do i = solidRef,comp
      sumMole  = sumMole + MoleFrac(i)
    enddo

    sumphi = 0.0d0
    sumteta = 0.0d0
    sumlj = 0.0d0
    do i = solidRef,comp
      MoleFracAux(i) = MoleFrac(i)
      sumphi = sumphi + MoleFracAux(i)*Uniquac_ri(i)
      sumteta = sumteta + MoleFracAux(i)*Uniquac_qi(i)  
      lj(i) = (Z/2.0d0)*(Uniquac_ri(i)-Uniquac_qi(i))-(Uniquac_ri(i)-1.0d0)    
      sumlj = sumlj + MoleFracAux(i)*lj(i)
    enddo
    
    do i = solidRef, comp
      !Calculating vaporization enthalpy   
      x = 1.d0 - T/Tc(i)
      Uniquac_dH0 = 5.2804d0*x**0.3333d0 + 12.8650d0*x**0.8333d0 + 1.1710d0*x**1.2083d0 - 13.1160d0*x + 0.4858d0*x**2 - 1.0880d0*x**3  
      Uniquac_dH1 = 0.80022d0*x**0.3333d0 + 273.23d0*x**0.8333d0 + 465.08d0*x**1.2083d0 - 638.51d0*x - 145.12d0*x**2 + 74.049d0*x**3
      Uniquac_dH2 = 7.2543d0*x**0.3333d0 - 346.45d0*x**0.8333d0 - 610.48d0*x**1.2083d0 + 839.89d0*x + 160.05d0*x**2 - 50.711d0*x**3
      Uniquac_dHvap = Uniquac_dH0 + Acentric(i)*Uniquac_dH1 + Acentric(i)*Acentric(i)*Uniquac_dH2
      Uniquac_dHvap = Uniquac_dHvap*R*Tc(i)
      
      !Calculating sublimation enthalpy
      Uniquac_dHsub(i) = Uniquac_dHvap + dHf(i) + dHtr(i)
      
      !Calculating interaction parameter (lambda)    
      do j = i, comp                    
        lambda(i,j) = (-2.0d0/Z)*(Uniquac_dHsub(i) - R*T)
        lambda(j,i) = lambda(i,j)
      enddo
        
      phi(i) = (MoleFracAux(i)*Uniquac_ri(i))/sumphi
      teta(i) = (MoleFracAux(i)*Uniquac_qi(i))/sumteta
    enddo

    do i = solidRef, comp    
      do j = solidRef, comp                    
        Tau(i,j) = exp(-(lambda(i,j)*(1.0d0 + Param)-lambda(i,i))/(Uniquac_qi(i)*R*T))
      enddo 
    enddo
           
    !Calculation of dln(gama)/dxi   
    termo1 = 0.0d0
    termo2 = 0.0d0
    termo3 = 0.0d0
    termo4 = 0.0d0
    termo5 = 0.0d0
    termo6 = 0.0d0
    termo7 = 0.0d0    
    
    !termo1
    termo1 = -Uniquac_ri(indexk)/sumphi
    
    !termo2
    termo2 = (Z/2.0d0)*Uniquac_qi(index)*((Uniquac_ri(indexk)/sumphi) - (Uniquac_qi(indexk)/sumteta))
    
    !termo3
    termo3 = 0.0d0
    
    !termo4
    termo4 = Uniquac_ri(index)*((lj(indexk)*sumphi-Uniquac_ri(indexk)*sumlj)/(sumphi**2))
    
    !Termo5
    sumtetaTji = 0.0d0
    do j = solidRef, comp        
        sumtetaTji = sumtetaTji + MoleFracAux(j)*Uniquac_qi(j)*Tau(index,j)        
    enddo
    
    termo5 = Uniquac_qi(index)*((Uniquac_qi(indexk)*Tau(index,indexk)/sumtetaTji) - (Uniquac_qi(indexk)/sumteta))    
    
    !termo6
    termo6 = 0.0d0
    
    !termo7
    sumtetaTmk = 0.0d0
    do m = solidRef, comp
        sumtetaTmk = sumtetaTmk + MoleFracAux(m)*Uniquac_qi(m)*Tau(indexk,m)        
    enddo    
    
    termo7parcela2 = 0.0d0
    do j = solidRef, comp
      
        sumtetaTmj = 0.0d0
        do m = solidRef, comp
            sumtetaTmj = sumtetaTmj + MoleFracAux(m)*Uniquac_qi(m)*Tau(j,m)            
        enddo        
        
        termo7parcela2 = termo7parcela2 + MoleFracAux(j)*Uniquac_qi(j)*Tau(j,index)*Uniquac_qi(indexk)*Tau(j,indexk)/(sumtetaTmj**2)        
        
    enddo
               
    termo7 = Uniquac_qi(index)*((Uniquac_qi(indexk)*Tau(indexk,index))/(sumtetaTmk) - termo7parcela2)    
    
    dlngamadxk = termo1 + termo2 + termo3 - termo4 - termo5 + termo6 - termo7
    
  end subroutine  

!===================================================================================================================================

  subroutine DerivadaTempUniquac (comp, T, P, Tc, Acentric, dHf, dHtr, Uniquac_ri, Uniquac_qi, MoleFrac, dlngamaidT, solidRef)
     
    !Arguments
    real(8) :: T, P, Tc(*), Acentric(*), dHf(*), dHtr(*), Uniquac_ri(*), Uniquac_qi(*), MoleFrac(comp), dlngamaidT(comp)
    integer :: solidRef, index, comp
    
    !Local    
    real(8) :: Uniquac_dHvap, Uniquac_dH0, Uniquac_dH1, Uniquac_dH2, Uniquac_dHsub(comp), MoleFracAux(comp), sumMole,termo7parcela2 
    real(8) :: R, sumphi, sumteta, x, phi(comp), teta(comp), lambda(comp,comp),  lj(comp), sumtetaTji, sumtetaTmk, sumtetaTmj
    real(8) :: termo1, termo2, sumlj, Tau(comp,comp), somatorio, sumTauij, Z
    real(8) :: d_dH0, d_dH1, d_dH2, d_dHvap, d_dHsub(comp), d_Tau(comp,comp), d_lambda(comp,comp), sumTauji, d_sumTauji, sumTaukj
    real(8) :: d_sumTaukj
    integer :: i, j, k, m, n
    real(8) :: Param
    
    !Param = 1.3730468750000000d-4 !Bim 0
    !Param = -2.4511718749999998d-3 !Bim 3
    !Param = -2.0710937500000005d-3 !Bim 5
    !Param = -1.2937304687499976d-2 !Bim 9
    !Param = -4.7722656250000000d-3 !Bim 13
    Param = 0.0d0
    
    !Initializing variables
    lambda = 0.0d0
    Tau = 0.0d0
    R = 1.9858775d0 !cal/(mol*K)
    Z = 10.0d0

    sumMole = 0.0d0
    do i = solidRef,comp
      sumMole  = sumMole + MoleFrac(i)
    enddo

    sumphi = 0.0d0
    sumteta = 0.0d0
    sumlj = 0.0d0
    do i = solidRef,comp
      MoleFracAux(i) = MoleFrac(i)/sumMole
      sumphi = sumphi + MoleFracAux(i)*Uniquac_ri(i)
      sumteta = sumteta + MoleFracAux(i)*Uniquac_qi(i)  
      lj(i) = (Z/2.0d0)*(Uniquac_ri(i)-Uniquac_qi(i))-(Uniquac_ri(i)-1.0d0)    
      sumlj = sumlj + MoleFracAux(i)*lj(i)
    enddo
    
    do i = solidRef, comp
      !Calculating vaporization enthalpy   
      x = 1.d0 - T/Tc(i)
      Uniquac_dH0 = 5.2804d0*x**0.3333d0 + 12.8650d0*x**0.8333d0 + 1.1710d0*x**1.2083d0 - 13.1160d0*x + 0.4858d0*x**2 - 1.0880d0*x**3  
      Uniquac_dH1 = 0.80022d0*x**0.3333d0 + 273.23d0*x**0.8333d0 + 465.08d0*x**1.2083d0 - 638.51d0*x - 145.12d0*x**2 + 74.049d0*x**3
      Uniquac_dH2 = 7.2543d0*x**0.3333d0 - 346.45d0*x**0.8333d0 - 610.48d0*x**1.2083d0 + 839.89d0*x + 160.05d0*x**2 - 50.711d0*x**3
      Uniquac_dHvap = Uniquac_dH0 + Acentric(i)*Uniquac_dH1 + Acentric(i)*Acentric(i)*Uniquac_dH2
      Uniquac_dHvap = Uniquac_dHvap*R*Tc(i)
      
      !Calculating sublimation enthalpy
      Uniquac_dHsub(i) = Uniquac_dHvap + dHf(i) + dHtr(i)
      
      !Calculating interaction parameter (lambda)    
      do j = i, comp                    
        lambda(i,j) = (-2.0d0/Z)*(Uniquac_dHsub(i) - R*T)
        lambda(j,i) = lambda(i,j)
      enddo
        
      phi(i) = (MoleFracAux(i)*Uniquac_ri(i))/sumphi
      teta(i) = (MoleFracAux(i)*Uniquac_qi(i))/sumteta
      
      
      !Derivada do deltaH de sublimação com a temperatura
      d_dH0 = -1.7600d0*(x**(-0.6667))/Tc(i) - 10.7204d0*(x**(-0.1667d0))/Tc(i) - 1.4149d0*(x**(0.2083))/Tc(i) + 13.116d0/Tc(i) -&
               0.9716d0*x/Tc(i) + 3.264d0*(x**2)/Tc(i)      
      d_dH1 = -0.2667d0*(x**(-0.6667))/Tc(i) - 227.68d0*(x**(-0.1667d0))/Tc(i) - 561.96d0*(x**(0.2083))/Tc(i) + 638.51d0/Tc(i) +&
               290.24d0*x/Tc(i) - 222.147d0*(x**2)/Tc(i)            
      d_dH2 = -2.4179d0*(x**(-0.6667))/Tc(i) + 288.70d0*(x**(-0.1667d0))/Tc(i) + 737.64d0*(x**(0.2083))/Tc(i) - 839.89/Tc(i) -&
               320.10d0*x/Tc(i) + 152.133d0*(x**2)/Tc(i)
               
      d_dHvap = d_dH0 + Acentric(i)*d_dH1 + Acentric(i)*Acentric(i)*d_dH2
      d_dHvap = d_dHvap*R*Tc(i)         
      d_dHsub(i) = d_dHvap
      
      !Calculating derivative of interaction parameter (lambda)    
      do j = i, comp                    
        d_lambda(i,j) = (-2.0d0/Z)*(d_dHsub(i) - R*T)
        d_lambda(j,i) = d_lambda(i,j)
      enddo      
           
    enddo   
    
    do i = solidRef, comp    
      do j = solidRef, comp                    
        Tau(i,j) = exp(-(lambda(i,j)*(1.0d0 + Param)-lambda(i,i))/(Uniquac_qi(i)*R*T))
      enddo 
    enddo

    do i = solidRef, comp    
      do j = solidRef, comp                    
        d_Tau(i,j) = (exp(-(lambda(i,j)*(1.0d0 + Param)-lambda(i,i))/(Uniquac_qi(i)*R*T)))*(((d_lambda(i,i)-d_lambda(i,j)*&
                            (1.0d0 + Param))*Uniquac_qi(i)*R*T -&
                            (lambda(i,i)-lambda(i,j)*(1.0d0 + Param))*Uniquac_qi(i)*R)/((Uniquac_qi(i)*R*T)**2))
      enddo 
    enddo    
     
    !Calculation of dln(gamai)/dT   
    termo1 = 0.0d0
    termo2 = 0.0d0
                
    do i = solidRef, comp
        
        index = i
        
        sumTauji = 0.0d0
        do j = solidRef,comp
          sumTauji  = sumTauji + teta(j)*Tau(index,j)      
        enddo        
        
        d_sumTauji = 0.0d0
        do j = solidRef,comp
          d_sumTauji  = d_sumTauji + teta(j)*d_Tau(index,j)      
        enddo         
        
        !Termo1
        termo1 = sumTauji/(d_sumTauji)
        
        !termo2
        somatorio = 0.0d0
        do j = solidRef,comp
          sumTaukj = 0.0d0
          do k = solidRef, comp                    
            sumTaukj = sumTaukj + teta(k)*Tau(j,k)    
          enddo
          
          d_sumTaukj = 0.0d0
          do k = solidRef, comp                    
            d_sumTaukj = d_sumTaukj + teta(k)*d_Tau(j,k)    
          enddo          
          somatorio = somatorio + (teta(j)*d_Tau(j,index)*sumTaukj - teta(j)*Tau(j,index)*d_sumTaukj)/(sumTaukj**2)
        enddo
        
        termo2 = somatorio                
        
        dlngamaidT(index) = -Uniquac_qi(index)*(termo1+termo2)
        
    enddo
        
  end subroutine  

!===================================================================================================================================  
    
end module


