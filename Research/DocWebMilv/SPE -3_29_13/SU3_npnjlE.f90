!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! - EDITED - !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************
!	added subroutines: DNSQE, QK15I 
!
!
!
!
!
!
!
!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************
		
!!!!!!!!!!!!!!!!!PROGRAMA NON LOCAL NJL CON POLYAKOV EN SU(3)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Basado en el articulo ''Nonlocal SU(3) chiral quark models at finite temperature:
!The role of the Polyakov loop'' Gustavo A. Contrera, Daniel G�mez Dumm, Norberto N. Scoccola

!PROGRAMA A TEMPERATURA Y POTENCIAL QUIMICO FINITOS
      
!Recordar cambiar la libreria : imsl.lib imsls_err.lib IMSLMPISTUB.LIB 

!BEGIN MAIN***************************************************************************************
	Implicit none
	Integer itmax,l,nreg,N,Nmax,k,i,j,Mmax,m 
	Integer title1,title2,title3,title4,poly,INTERV
	Parameter (N=3)
	Parameter (Nmax=70) !Para el potencial quimico
	Parameter (Mmax=199) !Para la Temperatura
	Complex*16 rg,w2
	Real*8 N_c,pi,pi2,Lambda,G_s,H 
	Real*8 FNORM, mu, muinf, musup, dmu,F(N)
	Real*8 X(N),sigui,sigsi,phi3i,XGUESS(N),erro
	Real*8 Tinf,Tsup,dT,T,xT,dxT, xm,dxm, xms,dxms,ms,mc
	Real*8 sigu,sigs,phi3,dxmu,xmu,rho,hbarc,Om_pot,Denso
	Real*8 BOUND,ERRABS,ERREL,ERREST,C_ss,qq_s,qq_u,M_s,M_u
	Real*8 Q_uu,Q_ss,cond_u,cond_s,Omega,Pot_U,Denu,Dens,rho_B
	Real*8 sigumax,sigsmax,p_0,chi_u,chi_s
	Real*8 Susc_u,Susc_s,qq_ushift,qq_sshift
	
	
	Integer IOPT, LWA, NPRINT, INFO !DNSQE
	Parameter (LWA = 180) !DNSQE
	Real*8 TOL, D1MACH, FVEC(N), WA(LWA)  !DNSQE
	
	Real*8 A, B, ABSERR, RESABS, RESASC !QK15I
	Parameter (A = 0)!QK15I
	Parameter (B = 1)!QK15I


	EXTERNAL FCN,DNEQNF,C_ss,C_uu
	
	EXTERNAL JAC !DNSQE

      	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
	Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/fcnsu3/sigu,sigs,phi3
	Common/polyakov/poly
	Common/param_DQDAGI/BOUND,ERRABS,ERREL,ERREST,INTERV
	Common/maxvalues/sigumax,sigsmax 
	Common/vacuo/p_0

      
	Data hbarc/197.327d0/  !MeV.fm
	Data dxm/1.d-7/  !.9d-5
	Data dxms/1.d-5/ !.9d-5
	Data dxmu/1.d-6/ !1.d-5
	Data dxT/1.d-4/  !.9d-4
	Data rho/1.3d6/
	Data mu/0.0d0/


	pi=4.0D0*datan(1.0D0)   !	pi=3.141592654D0

	pi2=pi**2
      
	title1=1  !Nombre de los archivos de salida
	title2=1
	title3=1
	title4=1
	
      !poly=1 NO incluye Polyakov
	poly=1
	nreg=1 !marca los sets de parametros


	If(poly.eq.1)then
	  N_c=3.d0
	else
	  N_c=1.d0
	endif

!+++++++++PARAMETROS PARA LAS INTEGRALES SEMI-INFINITAS(DQDAGI)+++
!Set limits of integration
      BOUND  = 0.0
      INTERV = 1

!Set error tolerances for the integral
      ERRABS =0.0d0 !1.d-5
      ERREL =1.d-6      !1.d-8
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
         !Datos para el SET I
!      If (nreg.eq.1) then
!	Lambda=842.9517109432D0        !MeV
!      mc=5.0d0              !MeV 
!	ms=119.2474094840D0              !MeV 
!      G_s=13.3450291392D0/Lambda**2   !1/MeV^2!
!	H=-273.7477267621D0/Lambda**5   !1/MeV^5
!	sigumax=366.469274690640d0	       
!	sigsmax=486.191144251437d0
!	Endif

	         !Datos para el SET I
!      If (nreg.eq.1) then
!	Lambda=709.D0        !MeV
!      mc=8.5d0              !MeV 
!	ms=223.D0              !MeV 
!      G_s=10.99D0/Lambda**2   !1/MeV^2!
!	H=-295.3D0/Lambda**5   !1/MeV^5
!	Endif
	!Datos para el SET COMPARATIVO
	If (nreg.eq.1) then
	Lambda=705.996432415474D0        !MeV
	mc=6.19711753633038d0              !MeV 
	ms=140.7D0              !MeV 
	G_s=15.0398354983559D0/Lambda**2   !1/MeV^2!
	H=-337.711951782219D0/Lambda**5   !1/MeV^5
	p_0=-1789388888.53626d0

      endif



	
      If(nreg.eq.1.and.poly.eq.1)then
      	sigui=450.d0   !366.9d0  !chutes iniciais p/sist de eq. n-lineares set I
      	sigsi=600.d0      !424.d0  !542.553d0
	phi3i=0.d0
      else
	sigui=350.d0  !chutes iniciais p/sist de eq. n-lineares set I
	sigsi=500.d0
	phi3i=60.d0
      endif

!**************************************************!        
! FIJO LA TEMPERATURA Y VARIO EL POTENCIAL QUIMICO !
!**************************************************!

      muinf=100.d0    
      musup=350.d0       !132.55d0         
      dmu=(musup-muinf)/Nmax
  
      xmu=muinf
!      Do l=0,Nmax
!**************************************************!        
! FIJO EL POTENCIAL QUIMICO Y VARIO LA TEMPERATURA !
!**************************************************!

      Tinf=50.d0   
      Tsup=200.d0        
      dT=(Tsup-Tinf)/Mmax

      xT=Tinf

!      Do j=0,Mmax

!	Do 50 m=1,2
!	xT=T+(m-2)*dxT

      
!	Do 30 i=1,2
!	xmu=mu+(i-2)*dxmu

      Do 80 k=1,2
	xm=mc+(k-1)*dxm

      Do 30 i=1,2
	xms=ms+(i-1)*dxms
	
!	Do 50 l=1,2
!	xms=ms  !+(i-2)*dxm


 
!*****************************************************************************

!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************
!	
! subroutine DNSQE does not use: erro, xguess or Itmax 
!
	
		
!Set values of initial guess for the non-linear system
   !   If(T.eq.30.d0)then
      
! 	If(xT.eq.Tinf.and.mu.eq.muinf)then
!     If(T.eq.Tinf)then
!
!     If(xmu.eq.muinf.or.xmu.eq.musup.and.xT.eq.Tinf)then
!	xguess(1)=sigui
!     	xguess(2)=sigsi
!    	xguess(3)=phi3i
!     else
!	xguess(1)=sigu
!	xguess(2)=sigs
!      	xguess(3)=phi3
!     endif

!	Itmax=1000      !(maximum number of iterations for the non-linear system)

!	erro=1.d-8    !1.d-10
	
!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************
!		
!	CALL DNEQNF (FCN,erro, N,itmax, xguess, x,fnorm) 
!			-hopefully replaced by DNSQE	
!
!Set values for DNSQE variables:
	IOPT = 2
        NPRINT = 0

!      SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
!      UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
!      THIS IS THE RECOMMENDED SETTING.
 
        TOL = SQRT(D1MACH(4))

 CALL DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA) 
!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************	
	
!      call FCN (X, F, N)
 !     call Densidad(x,Denu,Dens,rho_B)

!      CALL DQDAGI (Funs_bar1,BOUND,INTERV,ERRABS1,ERREL1,xint1,ERREST)

!	ERRABS2 =0.0d0  ! 1.d-3
!      ERREL2 =1.d-5   !1.d-8
!      CALL DQDAGI (Funs_bar2,BOUND,INTERV,ERRABS1,ERREL2,xint2,ERREST)
       
!      ERRABSCON=0.0D0
!	ERRELCON=1.D-5    !1.d-8 

!	CALL DQDAGI (Fun_qq,BOUND,INTERV,ERRABSCON,ERRELCON,xpair,ERREST)

	!Set error tolerances for the integral
!      ERRAB =0.0d0     !1.d-5
!      ERREL = 1.d-5   !1.d-8

!	CALL DQDAGI(Deriv_mu,BOUND,INTERV,ERRAB,ERREL,xdensi,ERREST)

    !**********************************************
    
!50    continue



!********************************************************************************
!--------------------------------------------------------------------------------
!	QK15I(F, BOUND, INTERV, A, B, RESULT, ABSERR, RESABS, RESASC) 
!********************************************************************************	
			
     ! CALL DQDAGI (C_ss,BOUND,INTERV,ERRABS,ERREL,qq_s,ERREST)     	     
 CALL QK15I(C_ss,BOUND,INTERV,A,B,qq_s,ABSERR,RESABS,RESASC) 

      If (i.eq.1)qq_sshift=qq_s
30    continue

      !CALL DQDAGI (C_uu,BOUND,INTERV,ERRABS,ERREL,qq_u,ERREST)
 CALL QK15I(C_uu,BOUND,INTERV,A,B,qq_u,ABSERR,RESABS,RESASC) 
           
	If (k.eq.1)qq_ushift=qq_u
80    continue
      
!	Call Densidad(x,den)


!	If (i.eq.1)xdenmenos=xdensi
!      If (i.eq.1)Om_pot=Omega(sigu,sigs,phi3)
!30    continue
!      If (i.eq.1)Denso=den
!30    continue


!50    If (m.eq.1)Om_pot=Omega(sigu,sigs,phi3) !,xT)



!	conden=-(4.d0*N_c*xT/pi2)*xpair
!     cond_u=(-Q_uu(x(1)))**(1.0D0/3.0D0)
!	cond_s= (-Q_ss(x(2)))**(1.0D0/3.0D0) 

	Q_uu=-4.d0*xT*N_c*qq_u/pi2
	Q_ss=-4.d0*xT*N_c*qq_s/pi2
   

	If (Q_uu.lt. 0.0D0) then
		cond_u=(-Q_uu)**(1.0D0/3.0D0)
	else
		cond_u=-(abs(Q_uu)**(1.0D0/3.0D0))
	endif

	If (Q_ss.lt. 0.0D0) then
		cond_s=(-Q_ss)**(1.0D0/3.0D0)
	else
		cond_s=-(abs(Q_ss)**(1.0D0/3.0D0))
	endif

!****************************************************
      !EXPRESION PARA EL CALOR ESPECIFICO
!	dif=Inpo_temp-Po_temp(phi3,s_bar1,s_bar2) !,xT)

!	Ca=xT*(dif/dxT)

!	write(*,*)'primera derivada'
!      write(*,*)Po_temp(phi3,s_bar1,s_bar2,xT)
!	write(*,*)'primera derivada'
!	Pause
!****************************************************
      !EXPRESION PARA LA SUSCEPTIBILIDAD QUIRAL

	Susc_u=-(4.d0*N_c*xT/pi2)*((qq_ushift-qq_u)/dxm)
      Susc_s=-(4.d0*N_c*xT/pi2)*((qq_sshift-qq_s)/dxms)

!****************************************************
      !EXPRESION PARA LA SUSCEPTIBILIDAD BARIONICA	

!	Denso=(4.d0*N_c*xT/pi2)*(Om_pot-Omega(sigu,sigs,phi3))/dxmu

!      Susbarion=(den-denso)/dxmu

!*****************************************************
      !SUSCEPTIBILIDAD ASOCIADA AL LOOP DE POLYAKOV

!!	JISPHI3=(2.D0/3.D0*xT**2)*DSIN(phi3/xT)
!	trace=(1.d0/3.d0)*(1.d0+2.d0*dcos(phi3/xT))
   
!*****************SALIDA******************************
      Open(3,file='sigmas_mu50.dat')
	Open(4,file='suscep_mu50.dat')

8       Format(3x,7D14.6)   !
  

      M_u=dreal(xm+x(1)*rg(w2))
      M_s=dreal(xms+x(2)*rg(w2))


	
	If (title1.eq.1)then
!	write(3,*)'T,sigu,sigs,<uu>,<ss>,Omega'
      write(3,*)'mu,T,sigu,sigs,cond_u,cond_s, Granpot'
	write(3,*)''
	endif
	title1=2

	If (title2.eq.1)then
!	write(3,*)'T,sigu,sigs,<uu>,<ss>,Omega'
      write(4,*)'mu,T,chi_u,chi_s, Granpot'
	write(4,*)''
	endif
	title2=2


      

	write(3,8)xmu,xT,x(1),x(2),cond_u,cond_s,Omega(x(1),x(2),x(3))
      write(4,8)xmu,xT,Susc_u,Susc_s,Omega(x(1),x(2),x(3))


!	write(3,8)xT,x(1),x(2),cond_u,cond_s,Omega(sigu,sigs,phi3)
   
      write(*,*)'T=',xT,'xmu=',xmu
	write(*,*)''
	write(*,*)'condensados       ','<uu>=',cond_u,'<ss>=', cond_s
	write(*,*)''
	write(*,*)'M_u',M_u,'M_s',M_s	 
      write(*,*)''
	write(*,*)'Granpotencial=',Omega(x(1),x(2),x(3)) !,'Densidad',Denso
	write(*,*)''
	write(*,*)'Pot_U=',Pot_U(x(3),xT)   !,'denu',denu,'dens',dens
	write(*,*)''
      write(*,*)'sigu=',x(1),'sigs=',x(2),'phi3=',x(3)
	write(*,*)'suscept','x_u=',Susc_u,'x_s=',Susc_s
	write(*,*)''
 !     pause
!	xT=xT+dT			
!      Enddo

!      xmu=xmu+dmu      
!      Enddo
!      Pause
      Stop
 !     pause
	end
!***********************************END MAIN*************************************
!regulador g interacciones no locales
!****************************************************************     
      Function rg(w2)
      Implicit none
	Integer nreg
      Complex*16 rg,w2
      Real*8 pi2,N_c,xT,mu,G_s,H,sigu,sigs,phi3,Lambda

      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
	Common/fcnsu3/sigu,sigs,phi3
      
	If(nreg.eq.1)then
	rg = cdexp(-w2/Lambda**2)
	Endif

	Return
	 End
!***************************************************************
!Evaluacion de las integrales en un intervalo semi-infinito 
!***************************************************************
!INTEGRAL S PARA LOS QUARKS u y d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Function S_u(p)
	  Implicit none 
        Integer j,nreg,n3,poly,n3_inic
	  Complex*16 rg,w2,M 
        Real*8 Sum,dSum, p,xmu,mu,sigu,sigs,phi3,suma_ant
        Real*8 Lambda,freq,xm,xms
        Real*8 pi2,N_c,xT,S_u,G_s,H
	  Real*8 Color,dColor
      
      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

!	write(*,*)'sigu',sigu

	  
! 	p=500.d0 
	 Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

        Do n3=n3_inic,1

	  
	  Sum=0.d0
       Do j = 0,3000
	suma_ant=sum
 
		freq=((2*j+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara

      w2=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))**2+(1.,0.)*p**2
	
	If(nreg.eq.1)then
!		write(*,*)'sigu',sigu
		M=xm+sigu*rg(w2)
      endif

       dSum=dreal((rg(w2)*M)/(w2+M**2))
       Sum=Sum+dSum
	if(sum .eq. suma_ant) exit
       Enddo
!      write(*,*)''
!	write(*,*)'num',dreal(rg(w2)*M),'denom',dreal((w2+M**2)),'n',j
!	write(*,*)'M',M,'rg(w2)',rg(w2)
!	write(*,*)'xm',xm,'sigu',sigu
!	write(*,*)'Sum',Sum,'suma_ant',suma_ant,'n',j
!      pause
!	write(*,*)'w2',w2,'freq',freq,'rg',rg(w2),'M',M
!	pause

!      do n=0,3000
!      suma_ant=suma
!      wn=(2*n+1)*pi*temp
!      w2=(wn+i*pot)**2+p**2
!      denom=dreal((w2+sigmag_compleja(mu,sigu,w2)**2))
!      numerad=dreal(g_compleja(w2)*sigmag_compleja(mu,sigu,w2))
!      suma=suma+(numerad/denom)
!      if (suma_ant .eq. suma) exit
!      end do


	 dColor=Sum
	 Color=Color+dColor
	 Enddo      
	  
	    S_u= p**2*Color
	        
       RETURN           
       END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INTEGRAL S PARA EL QUARK s !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Function S_s(p)
	  Implicit none 
        Integer j,nreg,n3,poly,n3_inic
	  Complex*16 rg,w2,M 
        Real*8 Sum,dSum, p, xmu,mu,sigu,sigs,phi3,suma_ant
        Real*8 Lambda,freq,xm,xms
        Real*8 pi2, N_c,xT,S_s ,G_s,H
	  Real*8 Color,dColor
      
      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3
	  
! 	p=500.d0 
	 Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

        Do n3=n3_inic,1
	       
	  
	  Sum=0.d0
       Do j = 0,3000
	suma_ant=sum
 
		freq=((2*j+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara

      w2=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))**2+(1.,0.)*p**2
	
	If(nreg.eq.1)then
	
		M=xms+sigs*rg(w2)
      endif

       dSum=dreal(rg(w2)*M/(w2+M**2))
       Sum=Sum+dSum
	if(sum .eq. suma_ant) exit
       Enddo
!	write(*,*)''
!	write(*,*)'sigu=',sigu,'sigs=',sigs
!	write(*,*)''
!	write(*,*)'num',dreal(rg(w2)*M),'denom',dreal((w2+M**2))
!	write(*,*)'M',M,'rg(w2)',rg(w2)
!	write(*,*)'xms',xms,'sigs',sigs
!	write(*,*)'Sum',Sum,'n',j

!	pause


	 dColor=Sum
	 Color=Color+dColor
	 Enddo      
        
	    S_s= p**2*Color
	        
       RETURN           
       END
!**************************************************************************************
!INTEGRAL CONDENSADO PARA LOS QUARKS u y d
        FUNCTION C_uu(p)
	  Implicit none
        Integer i,nreg,n3,n3_inic,poly
	  Complex*16 w2,rg,e1,e2,M,n1,n2
        Real*8 pi2,N_c,xT,mu,G_s,H,Lambda,xmu,Ep,suma_ant
        Real*8 p,sigu,sigs,phi3,xm,xms,freq,C_uu
	  Real*8 Color,dColor,Sum,dSum,termo

      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

!   	xT=220.d0
!	p=200.d0
	  Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	Do n3=n3_inic,1

	 	   
        Sum=0.d0
       Do i = 0,3000
	suma_ant=sum
  
	   freq=((2*i+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara
	       
	w2=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))**2+(1.,0.)*p**2
	

	If(nreg.eq.1)then
	M=(xm+sigu*rg(w2))

	endif  
	     		       
          dSum=dreal(M/(w2+M**2)-xm/(w2+xm**2))
          Sum=Sum+dSum
	if(sum .eq. suma_ant) exit
	Enddo
      
	!Termino de regularizacion!!!!!!!!!!!!!!!!!!!!!!!      
!	Ep=dsqrt(p**2+xm**2)
      

!	e1=cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
!	e2=cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT)
!	n1=1.d0+e1
!	n2=1.d0+e2
!      termo=(xm/(2.d0*xT*Ep))*dreal((e1/n1)+(e2/n2))

!	

!	write(*,*)'Sum',Sum,'term',-termo,'n',i
!      pause    
		dColor=Sum  !-termo	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    Color=Color+dColor

      Enddo
     
	 C_uu =p**2*Color

   
	   RETURN	   
         END
!**************************************************************************************
!INTEGRAL CONDENSADO PARA EL QUARK s
        FUNCTION C_ss(p)
	  Implicit none
        Integer i,nreg,n3,n3_inic,poly
	  Complex*16 w2,rg,e1,e2,M,n1,n2
        Real*8 pi2,N_c,xT,mu,G_s,H,Lambda,xmu,Ep,suma_ant
        Real*8 p,sigu,sigs,phi3,xm,xms,freq,C_ss
	  Real*8 Color,dColor,Sum,dSum,termo

      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

	   	   
	  Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	Do n3=n3_inic,1

	 	   
        Sum=0.d0
       Do i = 0,3000
	suma_ant=sum
  
	   freq=((2*i+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara
	       
	w2=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))**2+(1.,0.)*p**2
	

	If(nreg.eq.1)then
	M=(xms+sigs*rg(w2))

	endif  
	     		       
          dSum=dreal(M/(w2+M**2)-xms/(w2+xms**2))
          Sum=Sum+dSum
	if(sum .eq. suma_ant) exit
	Enddo
      
	!Termino de regularizacion!!!!!!!!!!!!!!!!!!!!!!!      
!	Ep=dsqrt(p**2+xms**2)
      

!	e1=cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
!	e2=cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT)
!	n1=1.d0+e1
!	n2=1.d0+e2
!      termo=(xms/(2.d0*xT*Ep))*dreal((e1/n1)+(e2/n2))
!          
		dColor=Sum  !-termo	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    Color=Color+dColor

      Enddo
     
	 C_ss =p**2*Color

   
	   RETURN	   
         END
!*************************************************************************
!CONTRIBUCION DE LOS QUARKS u Y d AL GRAN POTENCIAL
!*********************************
        Function Fun_Omu(p)
	  Implicit none
	  Integer k,nreg,n3,poly,n3_inic
        Complex*16 rg,w2, M,logui,n1,n2, loga
        Real*8 Suma, dSuma,sigu,sigs,phi3,p,Ep,suma_ant
        Real*8 pi2,N_c,xT,mu,xmu,xm,xms,G_s,H,Lambda,freq
        Real*8 Fun_Omu,Color,dColor

	Intrinsic cdlog, cdexp,dsqrt

	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      

	Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	Do n3=n3_inic,1

        Suma=0.d0

       Do k = 0,3000
	suma_ant=suma
		freq=((2*k+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara

      w2=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))**2+(1.,0.)*p**2

		
	If(nreg.eq.1)then
		 M=(xm+sigu*rg(w2))
	endif
      
	logui=cdlog((w2+M**2)/(w2+xm**2))
 
       dSuma=dreal(logui)
       
	 Suma=Suma+dSuma
      
	if(suma .eq. suma_ant) exit

	 Enddo
	  
      Ep=dsqrt(p**2+xm**2)

      n1=(1.d0+cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT))
	
      n2=1.d0+cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
	
       
	loga=cdlog(n1*n2)
 	    dColor=Suma
	    Color=Color+dColor+dreal(loga)/2.0d0
	 Enddo
	 
	
      Fun_Omu= p**2*Color

	 
       RETURN           
       END
!********************************************************************************
!CONTRIBUCION DEL QUARK s al GRANPOTENCIAL
!*********************************
        Function Fun_Oms(p)
	  Implicit none
	  Integer k,nreg,n3,poly,n3_inic
        Complex*16 rg,w2, M,logui,n1,n2, loga
        Real*8 Suma, dSuma,sigu,sigs,phi3,p,Ep,suma_ant
        Real*8 pi2,N_c,xT,mu,xmu,xm,xms,G_s,H,Lambda,freq
        Real*8 Fun_Oms,Color,dColor


	Intrinsic cdlog, cdexp,dsqrt

	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      

	Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	Do n3=n3_inic,1

        Suma=0.d0

       Do k = 0,3000
	suma_ant=suma
		freq=((2*k+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara

      w2=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))**2+(1.,0.)*p**2

		
	If(nreg.eq.1)then
		 M=(xms+sigs*rg(w2))
	endif
      
	logui=cdlog((w2+M**2)/(w2+xms**2))
 
       dSuma=dreal(logui)
       
	 Suma=Suma+dSuma
      
	if(suma .eq. suma_ant) exit

	 Enddo
	  
      Ep=dsqrt(p**2+xms**2)

      n1=(1.d0+cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT))
	
      n2=1.d0+cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
	
       
	loga=cdlog(n1*n2)
 	    dColor=Suma
	    Color=Color+dColor+dreal(loga)/2.0d0
	 Enddo
	 
	
      Fun_Oms= p**2*Color

	 
       RETURN           
       END
!********************************************************************************
!POTENCIAL OMEGA TOTAL: u,d,s!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Function Omega(su,ss,fi3)
	Implicit none
	Integer poly, nreg, INTERV
	Real*8 su,ss,fi3,sigu,sigs,phi3
	Real*8 BOUND,ERRABS,ERREL,ERREST
	Real*8 xm,xms,xmu,pi2,N_c,xT,mu,G_s,H,Lambda,Omega
	Real*8 Fun_Omu, Fun_Oms, xomu,xoms,ud_cont,s_cont,uds,S_s,S_u
	Real*8 int_u,int_s,Ese_u,Ese_s,ter,Pot_U,U_poly

			
	External DQDAGI, Fun_Omu, Fun_Oms,S_s,S_u

	Common/param_DQDAGI/BOUND,ERRABS,ERREL,ERREST,INTERV 
	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
	Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

	sigu=su
	sigs=ss
	phi3=fi3
	
!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************
		
      ! CALL DQDAGI (Fun_Omu,BOUND,INTERV,ERRABS,ERREL,xomu,ERREST)
      ! CALL DQDAGI (Fun_Oms,BOUND,INTERV,ERRABS,ERREL,xoms,ERREST)
      ! CALL DQDAGI (S_u,BOUND,INTERV,ERRABS,ERREL,int_u,ERREST)
      ! CALL DQDAGI (S_s,BOUND,INTERV,ERRABS,ERREL,int_s,ERREST)
      	      
 CALL QK15I(Fun_Omu,BOUND,INTERV,A,B,xomu,ABSERR,RESABS,RESASC)
 CALL QK15I(Fun_Oms,BOUND,INTERV,A,B,xoms,ABSERR,RESABS,RESASC)      
 CALL QK15I(S_u,BOUND,INTERV,A,B,int_u,ABSERR,RESABS,RESASC)
 CALL QK15I(S_s,BOUND,INTERV,A,B,int_s,ABSERR,RESABS,RESASC)
      	
	Ese_u=(-8.d0*N_c*xT*int_u)/pi2
	Ese_s=(-8.d0*N_c*xT*int_s)/pi2


	ud_cont=2.d0*(sigu*Ese_u+G_s*Ese_u**2/2.d0)
	s_cont=sigs*Ese_s+G_s*Ese_s**2/2.d0
	uds=H*Ese_u**2*Ese_s/2.d0
	ter=(ud_cont+s_cont+uds)/2.d0                !+U_poly(phi3,xT)

      Omega=(-2.d0*N_c*xT/pi2)*(2.d0*xomu+xoms)-ter+Pot_U(phi3,xT)
!      write(*,*)'xoms',2.d0*xoms
!      write(*,*)'Omega',Omega
!	pause

      return
	end
!*******************************************************************************
!INCLUSI�N DEL POLYAKOV
!********************************************************************************
!Potencial efectivo logaritmico debido a la inclusion del Poliakov
!*****************************************************       
	 Function Pot_U(xphi3,temp)
	Implicit none
	Integer nreg, poly
	 Real*8 pi2,N_c,xT,mu,G_s,H,Lambda 
       Real*8 sigu,sigs,phi3,xphi3,temp,xmu,xm,xms
	 Real*8 a_0,a_1,a_2,b_3,a_temp,b_temp,T_0,Traphi,loga
	 Real*8 Pot_U 

	Intrinsic dcos,dlog

      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

	 phi3=xphi3
	 xT=temp

	T_0=187.d0  !temperatura critica de deconfinamiento

	a_0=3.51d0
	a_1=-2.47d0
	a_2=15.2d0
	b_3=-1.75d0
   
	a_temp=a_0+a_1*(T_0/xT)+a_2*(T_0/xT)**2
	b_temp=b_3*(T_0/xT)**3
     	

	If(poly.eq.1)then
	Traphi=0.0d0
	loga=0.0d0
	else
	Traphi=(1.0d0/3.0d0)*(2.0d0*dcos(xphi3/xT)+1.0d0)
	loga=dlog(1.0d0-6.0d0*Traphi**2+8.0d0*Traphi**3-3.0d0*Traphi**4)
	endif

	
      Pot_U=((-1.0d0/2.0d0)*a_temp*Traphi**2+b_temp*loga)*xT**4


	RETURN	   

         END
!****************************************************************
!Potencial efectivo polinomico debido a la inclusion del Poliakov
!****************************************************************       
	 Function U_poly(xphi3,temp)
	Implicit none
	Integer nreg, poly
	 Real*8 pi2,N_c,xT,mu,G_s,H,Lambda 
       Real*8 sigu,sigs,phi3,xphi3,temp,xmu,xm,xms
	 Real*8 a_0,a_1,a_2,a_3,b_temp,T_0,Traphi,loga
	 Real*8 U_poly,b_3,b_4,bet2

	Intrinsic dcos,dlog

      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

	 phi3=xphi3
	 xT=temp

	T_0=187.d0  !temperatura critica de deconfinamiento (2+1)

	a_0=6.75d0
	a_1=-1.95d0
	a_2=2.625d0
	a_3=-7.44d0
	b_3=0.75d0
	b_4=7.5d0
   
	b_temp=a_0+a_1*(T_0/xT)+a_2*(T_0/xT)**2+a_3*(T_0/xT)**3
     	
	If(poly.eq.1)then
	Traphi=0.0d0
	else
	Traphi=(1.0d0/3.0d0)*(2.0d0*dcos(xphi3/xT)+1.0d0)
	endif

      bet2=-b_temp*Traphi**2/2.d0
	
      U_poly=(bet2-b_3*Traphi**3/3.d0+b_4*Traphi**4/4.d0)*xT**4

	RETURN	   

         END
!*******************************************************************************
!primera DERIVADA DEL GRANPOT RESPECTO DE PHI3 (QUARK u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Function Der1_Om_u(p)
      Implicit none
      Integer k,nreg,n3,n3_inic,poly
      Complex*16 rg,w2,M,p0,e1,e2,n1,n2
      Real*8 Sum, dSum,p,sigu,sigs,phi3,sum_ant,Color,dColor
      Real*8 pi2,N_c,xT,mu,G_s,H,Lambda 
	Real*8 Ep,xmu,xm,xms,Der1_Om_u,termo,freq
	Complex*16 dT1dphi3,dT2dphi3,Contrib1,Contrib2,paren

	Intrinsic cdlog, cdexp,dsqrt,dexp,cdsqrt
      
	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      

	Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	Do n3=n3_inic,1

        Sum=0.d0

       Do k = 0,4000
	sum_ant=sum

		freq=((2*k+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara
      
	p0=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))
      w2=p0**2+(1.,0.)*p**2

	If(nreg.eq.1)then
      M=(xm+sigu*rg(w2))
      paren=-rg(w2)*(sigu/(Lambda**2))

      Contrib1=2.d0*((M/(w2+M**2))*paren)
	Contrib2=(xm**2-M**2)/((w2+M**2)*(w2+xm**2))
	dT1dphi3=2.d0*n3*p0*(Contrib1+Contrib2)  !*0.0d0

	Endif

	!**************************************************************
	
       dSum=dreal(dT1dphi3)   !*0.d0
       
	 Sum=Sum+dSum
      
	if(sum .eq. sum_ant) exit
      
	 Enddo
	  
      Ep=dsqrt(p**2+xm**2)
!	write(*,*)Ey
!	pause

      e1=cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT)
	e2=cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
	n1=1.d0+e1
	n2=1.d0+e2
	
	dT2dphi3=(0.,1.)*n3*((e1/n1)-(e2/n2))
     
	termo=dreal(dT2dphi3)
    
	    dColor=Sum   
    
	    Color=Color+dColor+termo/(xT*2.d0)
	 Enddo
	
      Der1_Om_u = p**2*Color
	 

       RETURN           
       END
!*******************************************************************************
!primera DERIVADA DEL GRANPOT RESPECTO DE PHI3 (QUARK s)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Function Der1_Om_s(p)
      Implicit none
      Integer k,nreg,n3,n3_inic,poly
      Complex*16 rg,w2,M,p0,e1,e2,n1,n2
      Real*8 Sum, dSum,p,sigu,sigs,phi3,sum_ant,Color,dColor
      Real*8 pi2,N_c,xT,mu,G_s,H,Lambda 
	Real*8 Ep,xmu,xm,xms,Der1_Om_s,termo,freq
	Complex*16 dT1dphi3,dT2dphi3,Contrib1,Contrib2,paren

	Intrinsic cdlog, cdexp,dsqrt,dexp,cdsqrt
      
	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      

	Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	Do n3=n3_inic,1

        Sum=0.d0

       Do k = 0,4000
	sum_ant=sum

		freq=((2*k+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara
      
	p0=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))
      w2=p0**2+(1.,0.)*p**2

	If(nreg.eq.1)then
      M=(xms+sigs*rg(w2))
      paren=-rg(w2)*(sigs/(Lambda**2))

      Contrib1=2.d0*((M/(w2+M**2))*paren)
	Contrib2=(xms**2-M**2)/((w2+M**2)*(w2+xms**2))
	dT1dphi3=2.d0*n3*p0*(Contrib1+Contrib2)

	Endif

	!**************************************************************
	
       dSum=dreal(dT1dphi3)
       
	 Sum=Sum+dSum
      
	if(sum .eq. sum_ant) exit
      
	 Enddo
	  
      Ep=dsqrt(p**2+xms**2)
!	write(*,*)Ey
!	pause

      e1=cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT)
	e2=cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
	n1=1.d0+e1
	n2=1.d0+e2
	
	dT2dphi3=(0.,1.)*n3*((e1/n1)-(e2/n2))
     
	termo=dreal(dT2dphi3)
    
	    dColor=Sum   
    
	    Color=Color+dColor+termo/(xT*2.d0)
	 Enddo
	
      Der1_Om_s = p**2*Color
	 

       RETURN           
       END

!*******************************************************************************
!DERIVADA TOTAL DE OMEGA RESPECTO DE PHI3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Function dOmega_dphi3(su,ss,fi3)
	Implicit none
	Integer nreg, poly,INTERV,pot1
	Real*8 su,ss,fi3,sigu,sigs,phi3,dOmega_dphi3
	Real*8 pi2,N_c,xT,mu,G_s,H,Lambda,xm,xms,xmu
	Real*8 BOUND,ERRABS,ERREL,ERREST
	Real*8 Der1_Om_u,Der1_Om_s,x1u_dphi3,x1s_dphi3
	Real*8 dertrace,trace,paren,Valor,T_0,a_0,a_1,a_2,a_3
	Real*8 a_temp,b_temp,denom,Potenderiv,dUdphi3,S_s,S_u
	Real*8 dUpolyfi3,b_3,b_4,bet2


	External DQDAGI,Der1_Om_u,Der1_Om_s

      Common/param_DQDAGI/BOUND,ERRABS,ERREL,ERREST,INTERV 
      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      sigu=su
	sigs=ss
	phi3=fi3

      Pot1=1

      Valor=(-2.d0*N_c*xT/pi2)

      trace=(1.d0/3.d0)*(1.d0+2.d0*dcos(phi3/xT))
      dertrace=(-2.d0/(3.d0*xT))*dsin(phi3/xT)


	If(Pot1.eq.1)then
	T_0=187.d0  !temperatura critica de deconfinamiento

	a_0=3.51d0
	a_1=-2.47d0
	a_2=15.2d0
	b_3=-1.75d0
   
	a_temp=a_0+a_1*(T_0/xT)+a_2*(T_0/xT)**2
	b_temp=b_3*(T_0/xT)**3

      denom=(1.d0-6.d0*trace**2+8.d0*trace**3-3.d0*trace**4)

	paren=-a_temp-b_temp*(12.d0-24.d0*trace+12.d0*trace**2)/denom

	dUdphi3=xT**4*dertrace*trace*paren

	else

	T_0=187.d0  !temperatura critica de deconfinamiento (2+1)

	a_0=6.75d0
	a_1=-1.95d0
	a_2=2.625d0
	a_3=-7.44d0
	b_3=0.75d0
	b_4=7.5d0
   
	b_temp=a_0+a_1*(T_0/xT)+a_2*(T_0/xT)**2+a_3*(T_0/xT)**3
	bet2=-b_temp*trace**2
	
      dUdphi3=(-b_3*trace**3+b_4*trace**4)*dertrace*xT**4
	endif
!	write(*,*)dUdphi3
!	pause
	
	

	! CALL DQDAGI (Der1_Om_u,BOUND,INTERV,ERRABS,ERREL,x1u_dphi3,ERREST)
	! CALL DQDAGI (Der1_Om_s,BOUND,INTERV,ERRABS,ERREL,x1s_dphi3,ERREST)
	
 CALL QK15I(Der1_Om_u,BOUND,INTERV,A,B,x1u_dphi3,ABSERR,RESABS,RESASC)	
 CALL QK15I(Der1_Om_s,BOUND,INTERV,A,B,x1s_dphi3,ABSERR,RESABS,RESASC)	
	
	
	
	
!	write(*,*)int_u,int_s,x1u_dphi3,x1s_dphi3
!	pause
!      	write(*,*)'integral', Valor*x1u_dphi3
!	pause

	Potenderiv=(Valor)*(2.d0*x1u_dphi3+x1s_dphi3)

	dOmega_dphi3=Potenderiv+dUdphi3

!	write(*,*)dOmega_dphi3/2.d0
!	pause

	return
	end

!********************************************************************************
!Derivada analitica del Granpotencial de u respecto de mu
      Function Derivu_mu(p)
      Implicit none
      Integer k,nreg,n3,n3_inic,poly
      Complex*16 rg,w2,M,p0,e1,e2,n1,n2
      Real*8 Sum, dSum,p,sigu,sigs,phi3,sum_ant,Color,dColor
      Real*8 pi2,N_c,xT,mu,G_s,H,Lambda 
	Real*8 Ep,xmu,xm,xms,Derivu_mu,termo,freq
	Complex*16 dT1dmu,dT2dmu,Contrib1,Contrib2,paren

	Intrinsic cdlog, cdexp,dsqrt,dexp,cdsqrt
      
	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      

	Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	 Do n3=n3_inic,1	  
        Sum=0.d0

       Do k = 0,4000
	sum_ant=sum

		freq=((2*k+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara
      
	p0=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))
      w2=p0**2+(1.,0.)*p**2

	If(nreg.eq.1)then

      M=(xm+sigu*rg(w2))
       paren=-rg(w2)*(sigu/(Lambda**2))

      Contrib1=2.d0*((M/(w2+M**2))*paren)
	Contrib2=(xm**2-M**2)/((w2+M**2)*(w2+xm**2))
	dT1dmu=(2.d0*(0.,1.)*p0)*(Contrib1+Contrib2)

	Endif
	!**************************************************************
	
       dSum=dreal(dT1dmu)
       
	 Sum=Sum+dSum
      
	if(sum .eq. sum_ant) exit
      
	 Enddo
	  
      Ep=dsqrt(p**2+xm**2)
!	write(*,*)Ey
!	pause

      e1=cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT)
	e2=cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
	n1=1.d0+e1
	n2=1.d0+e2
	
	dT2dmu=((e1/n1)-(e2/n2))
!	write(*,*)dT2dphi3
!	pause
     
	termo=dreal(dT2dmu)
    
	    dColor=Sum   
    
	    Color=Color+dColor-termo/(xT*2.d0)
	 Enddo
	
      Derivu_mu = p**2*Color
	 

       RETURN           
       END     
!********************************************************************************
!Derivada analitica del Granpotencial de s respecto de mu
      Function Derivs_mu(p)
      Implicit none
      Integer k,nreg,n3,n3_inic,poly
      Complex*16 rg,w2,M,p0,e1,e2,n1,n2
      Real*8 Sum, dSum,p,sigu,sigs,phi3,sum_ant,Color,dColor
      Real*8 pi2,N_c,xT,mu,G_s,H,Lambda 
	Real*8 Ep,xmu,xm,xms,Derivs_mu,termo,freq
	Complex*16 dT1dmu,dT2dmu,Contrib1,Contrib2,paren

	Intrinsic cdlog, cdexp,dsqrt,dexp,cdsqrt
      
	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/activia/nreg
      Common/variablem/xm,xms
	Common/variablemu/xmu
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      

	Color=0.d0
       If(poly.eq.1) then
               n3_inic=1
       else
               n3_inic=-1
       endif

	 Do n3=n3_inic,1	  
        Sum=0.d0

       Do k = 0,4000
	sum_ant=sum

		freq=((2*k+1)*dsqrt(pi2)*xT)!Frecuencias de Matsubara
      
	p0=(freq*(1.,0.)+xmu*(0.,1.)+n3*phi3*(1.,0.))
      w2=p0**2+(1.,0.)*p**2

	If(nreg.eq.1)then

      M=(xms+sigs*rg(w2))
       paren=-rg(w2)*(sigs/(Lambda**2))

      Contrib1=2.d0*((M/(w2+M**2))*paren)
	Contrib2=(xms**2-M**2)/((w2+M**2)*(w2+xms**2))
	dT1dmu=(2.d0*(0.,1.)*p0)*(Contrib1+Contrib2)

	Endif

	!**************************************************************
	
       dSum=dreal(dT1dmu)
       
	 Sum=Sum+dSum
      
	if(sum .eq. sum_ant) exit
      
	 Enddo
	  
      Ep=dsqrt(p**2+xm**2)
!	write(*,*)Ey
!	pause

      e1=cdexp(-(Ep+xmu-n3*phi3*(0.,1.))/xT)
	e2=cdexp(-(Ep-xmu+n3*phi3*(0.,1.))/xT)
	n1=1.d0+e1
	n2=1.d0+e2
	
	dT2dmu=((e1/n1)-(e2/n2))
!	write(*,*)dT2dphi3
!	pause
     
	termo=dreal(dT2dmu)
    
	    dColor=Sum   
    
	    Color=Color+dColor-termo/(xT*2.d0)

	 Enddo
	
      Derivs_mu = p**2*Color
	

       RETURN           
       END     
!*********************************************************************************
!Calculo de la densidad    
      Subroutine Densidad(x,Denu,Dens,rho_B)
	Implicit none
	Integer INTERV,N,poly
	Parameter (N=3)
	Real*8 ERRABS,ERREL,ERREST,X(3)
	Real*8 sigu,sigs,phi3,Derivu_mu,Denu,Dens,rho_B
	Real*8 pi2,N_c,xT,mu,G_s,H,Lambda,xdensiu,xdensis
	Real*8 BOUND,Valor

	
	EXTERNAL DQDAGI, Derivu_mu, Derivs_mu

	Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
      Common/param_DQDAGI/BOUND,ERRABS,ERREL,ERREST,INTERV 
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3

      
	Valor=(-2.d0*N_c*xT/pi2)
	
      ! CALL DQDAGI(Derivu_mu,BOUND,INTERV,ERRABS,ERREL,xdensiu,ERREST)
      ! CALL DQDAGI(Derivs_mu,BOUND,INTERV,ERRABS,ERREL,xdensis,ERREST)
     
 CALL QK15I(Derivu_mu,BOUND,INTERV,A,B,xdensiu,ABSERR,RESABS,RESASC)
 CALL QK15I(Derivs_mu,BOUND,INTERV,A,B,xdensis,ABSERR,RESABS,RESASC)
      	

	If(poly.eq.1)then

	Denu=0.0d0
	Dens=0.0d0
	else	  
	Denu=-Valor*xdensiu
	Dens=-Valor*xdensis

	endif

	rho_B=(1.d0/3.d0)*(2.d0*Denu+Dens)


	RETURN
	END
!********************************************************************************
!Resolucion del sistema de ecuaciones no lineales
        SUBROUTINE FCN (X, F, N)
        Implicit none
        INTEGER N, INTERV,nreg,poly
        Real*8 X(3),F(N),sigu,sigs,phi3,campu,camps,sigumax,sigsmax 
	  Real*8 BOUND,ERRABS,ERREL,ERREST,pi2,N_c,xT,mu,G_s,H,Lambda
	  Real*8 int_u,int_s,S_u,S_s,Ese_u,Ese_s,xmu,dOmega_dphi3

      EXTERNAL DQDAGI,S_u,S_s  
   
    
      Common/dat1/pi2,N_c,xT,mu,G_s,H,Lambda
	Common/param_DQDAGI/BOUND,ERRABS,ERREL,ERREST,INTERV 
	Common/polyakov/poly
	Common/fcnsu3/sigu,sigs,phi3
	Common/variablemu/xmu
	Common/maxvalues/sigumax,sigsmax 

!      call Map(x,campu,camps)


	sigu=x(1)   !campu     !x(1)   !
	sigs=x(2)    !camps      !x(2)      !
	phi3=x(3)


      ! CALL DQDAGI (S_u,BOUND,INTERV,ERRABS,ERREL,int_u,ERREST)
      ! CALL DQDAGI (S_s,BOUND,INTERV,ERRABS,ERREL,int_s,ERREST)
      
 CALL QK15I(S_u,BOUND,INTERV,A,B,int_u,ABSERR,RESABS,RESASC)
 CALL QK15I(S_s,BOUND,INTERV,A,B,int_s,ABSERR,RESABS,RESASC)


!      pause
!	write(*,*)''
!	write(*,*)'antes','sigu=',sigu,'sigs=',sigs,'phi3=',phi3
!	Write(*,*)'mu=',xmu,'T=',xt
!      write(*,*)''

	Ese_u=(-8.d0*N_c*xT*int_u)/pi2
	Ese_s=(-8.d0*N_c*xT*int_s)/pi2

!      write(*,*)'S_u=',Ese_u,'S_s=',Ese_s
!      pause

!	write(*,*)Ese_u,N_c
!	pause

	  
      If(poly.eq.1)then

	x(3)=0.0d0

	F(1)=sigu+G_s*Ese_u+H*Ese_u*Ese_s/2.d0
	F(2)=sigs+G_s*Ese_s+H*Ese_u**2/2.d0
	F(3)=0.0d0
	else
	F(1)=sigu/G_s+Ese_u+H*Ese_u*Ese_s/(2.d0*G_s)
!	F(1)=sigu+G_s*Ese_u+H*Ese_u*Ese_s/2.d0
	F(2)=sigs+G_s*Ese_s+H*Ese_u**2/2.d0
	F(3)=dOmega_dphi3(sigu,sigs,phi3)
      endif
!      write(*,*)''
!	write(*,*)'f1=',f(1),'f2=',f(2),'f3=',f(3)
!      write(*,*)''
!	write(*,*)'Ese_u=',Ese_u,'Ese_=',Ese_s
!	write(*,*)''
!	write(*,*)'despues','sigu=',sigu,'sigs=',sigs,'phi3=',phi3
!	write(*,*)''
!	pause

       RETURN           
       END
!********************************************************************************
      Subroutine Map(x,campu,camps)
	Real*8 x(3),campu,camps,sigumax,sigsmax
		
	Common/maxvalues/sigumax,sigsmax 

      campu=sigumax*(dsin(x(1)))**2
	camps=sigsmax*(dsin(x(2)))**2

	return
	end
!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************
	
!*DECK DNSQE
      SUBROUTINE DNSQE (FCN, JAC, IOPT, N, X, FVEC, TOL, NPRINT, INFO, &
     &   WA, LWA)
! C***BEGIN PROLOGUE  DNSQE
! C***PURPOSE  An easy-to-use code to find a zero of a system of N
! C            nonlinear functions in N variables by a modification of
! C            the Powell hybrid method.
! C***LIBRARY   SLATEC
! C***CATEGORY  F2A
! C***TYPE      DOUBLE PRECISION (SNSQE-S, DNSQE-D)
! C***KEYWORDS  EASY-TO-USE, NONLINEAR SQUARE SYSTEM,
! C             POWELL HYBRID METHOD, ZEROS
! C***AUTHOR  Hiebert, K. L. (SNLA)
! C***DESCRIPTION
! C
! C 1. Purpose.
! C
! C       The purpose of DNSQE is to find a zero of a system of N
! C       nonlinear functions in N variables by a modification of the
! C       Powell hybrid method.  This is done by using the more general
! C       nonlinear equation solver DNSQ.  The user must provide a
! C       subroutine which calculates the functions.  The user has the
! C       option of either to provide a subroutine which calculates the
! C       Jacobian or to let the code calculate it by a forward-difference
! C       approximation.  This code is the combination of the MINPACK
! C       codes (Argonne) HYBRD1 and HYBRJ1.
! C
! C 2. Subroutine and Type Statements.
! C
! C       SUBROUTINE DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,
! C      *                  WA,LWA)
! C       INTEGER IOPT,N,NPRINT,INFO,LWA
! C       DOUBLE PRECISION TOL
! C       DOUBLE PRECISION X(N),FVEC(N),WA(LWA)
! C       EXTERNAL FCN,JAC
! C
! C 3. Parameters.
! C
! C       Parameters designated as input parameters must be specified on
! C       entry to DNSQE and are not changed on exit, while parameters
! C       designated as output parameters need not be specified on entry
! C       and are set to appropriate values on exit from DNSQE.
! C
! C       FCN is the name of the user-supplied subroutine which calculates
! C         the functions.  FCN must be declared in an external statement
! C         in the user calling program, and should be written as follows.
! C
! C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
! C         INTEGER N,IFLAG
! C         DOUBLE PRECISION X(N),FVEC(N)
! C         ----------
! C         Calculate the functions at X and
! C         return this vector in FVEC.
! C         ----------
! C         RETURN
! C         END
! C
! C         The value of IFLAG should not be changed by FCN unless the
! C         user wants to terminate execution of DNSQE.  In this case set
! C         IFLAG to a negative integer.
! C
! C       JAC is the name of the user-supplied subroutine which calculates
! C         the Jacobian.  If IOPT=1, then JAC must be declared in an
! C         external statement in the user calling program, and should be
! C         written as follows.
! C
! C         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
! C         INTEGER N,LDFJAC,IFLAG
! C         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
! C         ----------
! C         Calculate the Jacobian at X and return this
! C         matrix in FJAC.  FVEC contains the function
! C         values at X and should not be altered.
! C         ----------
! C         RETURN
! C         END
! C
! C         The value of IFLAG should not be changed by JAC unless the
! C         user wants to terminate execution of DNSQE. In this case set
! C         IFLAG to a negative integer.
! C
! C         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
! C
! C       IOPT is an input variable which specifies how the Jacobian will
! C         be calculated.  If IOPT=1, then the user must supply the
! C         Jacobian through the subroutine JAC.  If IOPT=2, then the
! C         code will approximate the Jacobian by forward-differencing.
! C
! C       N is a positive integer input variable set to the number of
! C         functions and variables.
! C
! C       X is an array of length N.  On input X must contain an initial
! C         estimate of the solution vector.  On output X contains the
! C         final estimate of the solution vector.
! C
! C       FVEC is an output array of length N which contains the functions
! C         evaluated at the output X.
! C
! C       TOL is a nonnegative input variable.  Termination occurs when
! C         the algorithm estimates that the relative error between X and
! C         the solution is at most TOL.  Section 4 contains more details
! C         about TOL.
! C
! C       NPRINT is an integer input variable that enables controlled
! C         printing of iterates if it is positive.  In this case, FCN is
! C         called with IFLAG = 0 at the beginning of the first iteration
! C         and every NPRINT iterations thereafter and immediately prior
! C         to return, with X and FVEC available for printing. Appropriate
! C         print statements must be added to FCN(see example).  If NPRINT
! C         is not positive, no special calls of FCN with IFLAG = 0 are
! C         made.
! C
! C       INFO is an integer output variable.  If the user has terminated
! C         execution, INFO is set to the (negative) value of IFLAG.  See
! C         description of FCN and JAC. Otherwise, INFO is set as follows.
! C
! C         INFO = 0  Improper input parameters.
! C
! C         INFO = 1  Algorithm estimates that the relative error between
! C                   X and the solution is at most TOL.
! C
! C         INFO = 2  Number of calls to FCN has reached or exceeded
! C                   100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2.
! C
! C         INFO = 3  TOL is too small.  No further improvement in the
! C                   approximate solution X is possible.
! C
! C         INFO = 4  Iteration is not making good progress.
! C
! C         Sections 4 and 5 contain more details about INFO.
! C
! C       WA is a work array of length LWA.
! C
! C       LWA is a positive integer input variable not less than
! C         (3*N**2+13*N))/2.
! C
! C 4. Successful Completion.
! C
! C       The accuracy of DNSQE is controlled by the convergence parameter
! C       TOL.  This parameter is used in a test which makes a comparison
! C       between the approximation X and a solution XSOL.  DNSQE
! C       terminates when the test is satisfied.  If TOL is less than the
! C       machine precision (as defined by the  function D1MACH(4)), then
! C       DNSQE only attempts to satisfy the test defined by the machine
! C       precision.  Further progress is not usually possible.  Unless
! C       high precision solutions are required, the recommended value
! C       for TOL is the square root of the machine precision.
! C
! C       The test assumes that the functions are reasonably well behaved,
! C       and, if the Jacobian is supplied by the user, that the functions
! C       and the Jacobian are coded consistently. If these conditions are
! C       not satisfied, then DNSQE may incorrectly indicate convergence.
! C       The coding of the Jacobian can be checked by the subroutine
! C       DCKDER.  If the Jacobian is coded correctly or IOPT=2, then
! C       the validity of the answer can be checked, for example, by
! C       rerunning DNSQE with a tighter tolerance.
! C
! C       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
! C         vector Z, then this test attempts to guarantee that
! C
! C               DENORM(X-XSOL) .LE. TOL*DENORM(XSOL).
! C
! C         If this condition is satisfied with TOL = 10**(-K), then the
! C         larger components of X have K significant decimal digits and
! C         INFO is set to 1.  There is a danger that the smaller
! C         components of X may have large relative errors, but the fast
! C         rate of convergence of DNSQE usually avoids this possibility.
! C
! C 5. Unsuccessful Completion.
! C
! C       Unsuccessful termination of DNSQE can be due to improper input
! C       parameters, arithmetic interrupts, an excessive number of
! C       function evaluations, errors in the functions, or lack of good
! C       progress.
! C
! C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1, or
! C         IOPT .GT. 2, or N .LE. 0, or TOL .LT. 0.E0, or
! C         LWA .LT. (3*N**2+13*N)/2.
! C
! C       Arithmetic Interrupts.  If these interrupts occur in the FCN
! C         subroutine during an early stage of the computation, they may
! C         be caused by an unacceptable choice of X by DNSQE.  In this
! C         case, it may be possible to remedy the situation by not
! C         evaluating the functions here, but instead setting the
! C         components of FVEC to numbers that exceed those in the initial
! C         FVEC.
! C
! C       Excessive Number of Function Evaluations.  If the number of
! C         calls to FCN reaches 100*(N+1) for IOPT=1 or 200*(N+1) for
! C         IOPT=2, then this indicates that the routine is converging
! C         very slowly as measured by the progress of FVEC, and INFO is
! C         set to 2.  This situation should be unusual because, as
! C         indicated below, lack of good progress is usually diagnosed
! C         earlier by DNSQE, causing termination with INFO = 4.
! C
! C       Errors In the Functions.  When IOPT=2, the choice of step length
! C         in the forward-difference approximation to the Jacobian
! C         assumes that the relative errors in the functions are of the
! C         order of the machine precision.  If this is not the case,
! C         DNSQE may fail (usually with INFO = 4).  The user should
! C         then either use DNSQ and set the step length or use IOPT=1
! C         and supply the Jacobian.
! C
! C       Lack of Good Progress.  DNSQE searches for a zero of the system
! C         by minimizing the sum of the squares of the functions.  In so
! C         doing, it can become trapped in a region where the minimum
! C         does not correspond to a zero of the system and, in this
! C         situation, the iteration eventually fails to make good
! C         progress.  In particular, this will happen if the system does
! C         not have a zero.  If the system has a zero, rerunning DNSQE
! C         from a different starting point may be helpful.
! C
! C 6. Characteristics of The Algorithm.
! C
! C       DNSQE is a modification of the Powell Hybrid method.  Two of
! C       its main characteristics involve the choice of the correction as
! C       a convex combination of the Newton and scaled gradient
! C       directions, and the updating of the Jacobian by the rank-1
! C       method of Broyden.  The choice of the correction guarantees
! C       (under reasonable conditions) global convergence for starting
! C       points far from the solution and a fast rate of convergence.
! C       The Jacobian is calculated at the starting point by either the
! C       user-supplied subroutine or a forward-difference approximation,
! C       but it is not recalculated until the rank-1 method fails to
! C       produce satisfactory progress.
! C
! C       Timing.  The time required by DNSQE to solve a given problem
! C         depends on N, the behavior of the functions, the accuracy
! C         requested, and the starting point.  The number of arithmetic
! C         operations needed by DNSQE is about 11.5*(N**2) to process
! C         each evaluation of the functions (call to FCN) and 1.3*(N**3)
! C         to process each evaluation of the Jacobian (call to JAC,
! C         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
! C         the timing of DNSQE will be strongly influenced by the time
! C         spent in FCN and JAC.
! C
! C       Storage.  DNSQE requires (3*N**2 + 17*N)/2 single precision
! C         storage locations, in addition to the storage required by the
! C         program.  There are no internally declared storage arrays.
! C
! C *Long Description:
! C
! C 7. Example.
! C
! C       The problem is to determine the values of X(1), X(2), ..., X(9),
! C       which solve the system of tridiagonal equations
! C
! C       (3-2*X(1))*X(1)           -2*X(2)                   = -1
! C               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
! C                                   -X(8) + (3-2*X(9))*X(9) = -1
! C
! C       **********
! C
! C       PROGRAM TEST
! C C
! C C     DRIVER FOR DNSQE EXAMPLE.
! C C
! C       INTEGER J,N,IOPT,NPRINT,INFO,LWA,NWRITE
! C       DOUBLE PRECISION TOL,FNORM
! C       DOUBLE PRECISION X(9),FVEC(9),WA(180)
! C       DOUBLE PRECISION DENORM,D1MACH
! C       EXTERNAL FCN
! C       DATA NWRITE /6/
! C C
! C       IOPT = 2
! C       N = 9
! C C
! C C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
! C C
! C       DO 10 J = 1, 9
! C          X(J) = -1.E0
! C    10    CONTINUE
! C
! C       LWA = 180
! C       NPRINT = 0
! C C
! C C     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
! C C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
! C C     THIS IS THE RECOMMENDED SETTING.
! C C
! C       TOL = SQRT(D1MACH(4))
! C C
! C       CALL DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
! C       FNORM = DENORM(N,FVEC)
! C       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
! C       STOP
! C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
! C      *        5X,' EXIT PARAMETER',16X,I10 //
! C      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
! C       END
! C       SUBROUTINE FCN(N,X,FVEC,IFLAG)
! C       INTEGER N,IFLAG
! C       DOUBLE PRECISION X(N),FVEC(N)
! C       INTEGER K
! C       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
! C       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
! C C
! C       DO 10 K = 1, N
! C          TEMP = (THREE - TWO*X(K))*X(K)
! C          TEMP1 = ZERO
! C          IF (K .NE. 1) TEMP1 = X(K-1)
! C          TEMP2 = ZERO
! C          IF (K .NE. N) TEMP2 = X(K+1)
! C          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
! C    10    CONTINUE
! C       RETURN
! C       END
! C
! C       RESULTS OBTAINED WITH DIFFERENT COMPILERS OR MACHINES
! C       MAY BE SLIGHTLY DIFFERENT.
! C
! C       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07
! C
! C       EXIT PARAMETER                         1
! C
! C       FINAL APPROXIMATE SOLUTION
! C
! C       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
! C       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
! C       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
! C
! C***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
! C                 tions. In Numerical Methods for Nonlinear Algebraic
! C                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
! C                 1988.
! C***ROUTINES CALLED  DNSQ, XERMSG
! C***REVISION HISTORY  (YYMMDD)
! C   800301  DATE WRITTEN
! C   890531  Changed all specific intrinsics to generic.  (WRB)
! C   890831  Modified array declarations.  (WRB)
! C   890831  REVISION DATE from Version 3.2
! C   891214  Prologue converted to Version 4.0 format.  (BAB)
! C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
! C   920501  Reformatted the REFERENCES section.  (WRB)
! C***END PROLOGUE  DNSQE
      INTEGER INDEX, INFO, IOPT, J, LR, LWA, MAXFEV, ML, MODE, MU, N,&
     &     NFEV, NJEV, NPRINT !deleted a 1 at the beggining of this line
      DOUBLE PRECISION EPSFCN, FACTOR, FVEC(*), ONE, TOL, WA(*),&
     &     X(*), XTOL, ZERO
      EXTERNAL FCN, JAC
      SAVE FACTOR, ONE, ZERO
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
! C     BEGIN BLOCK PERMITTING ...EXITS TO 20
! C***FIRST EXECUTABLE STATEMENT  DNSQE
         ! INFO = 0
! C
! C        CHECK THE INPUT PARAMETERS FOR ERRORS.
! C
! C     ...EXIT
         IF (IOPT .LT. 1 .OR. IOPT .GT. 2 .OR. N .LE. 0 &
     &       .OR. TOL .LT. ZERO .OR. LWA .LT. (3*N**2 + 13*N)/2)
     2      GO TO 20
! C
! C        CALL DNSQ.
! C
         MAXFEV = 100*(N + 1)
         IF (IOPT .EQ. 2) MAXFEV = 2*MAXFEV
         XTOL = TOL
         ML = N - 1
         MU = N - 1
         EPSFCN = ZERO
         MODE = 2
         DO 10 J = 1, N
            WA(J) = ONE
   10    CONTINUE
         LR = (N*(N + 1))/2
         INDEX = 6*N + LR
         CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,WA(INDEX+1),N,XTOL,MAXFEV,ML,&
          &        MU,EPSFCN,WA(1),MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
          &       WA(6*N+1),LR,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),
          &        WA(5*N+1))
         IF (INFO .EQ. 5) INFO = 4
   20 CONTINUE
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DNSQE', &
     &   'INVALID INPUT PARAMETER.', 2, 1)
      RETURN
! C
! C     LAST CARD OF SUBROUTINE DNSQE.
! C
      END  
!********************************************************************************
!--------------------------------------------------------------------------------
!********************************************************************************
!*DECK QK15I
      SUBROUTINE QK15I (F, BOUN, INF, A, B, RESULT, ABSERR, RESABS,&
     &   RESASC)
! C***BEGIN PROLOGUE  QK15I
! C***PURPOSE  The original (infinite integration range is mapped
! C            onto the interval (0,1) and (A,B) is a part of (0,1).
! C            it is the purpose to compute
! C            I = Integral of transformed integrand over (A,B),
! C            J = Integral of ABS(Transformed Integrand) over (A,B).
! C***LIBRARY   SLATEC (QUADPACK)
! C***CATEGORY  H2A3A2, H2A4A2
! C***TYPE      SINGLE PRECISION (QK15I-S, DQK15I-D)
! C***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
! C***AUTHOR  Piessens, Robert
! C             Applied Mathematics and Programming Division
! C             K. U. Leuven
! C           de Doncker, Elise
! C             Applied Mathematics and Programming Division
! C             K. U. Leuven
! C***DESCRIPTION
! C
! C           Integration Rule
! C           Standard Fortran subroutine
! C           Real version
! C
! C           PARAMETERS
! C            ON ENTRY
! C              F      - Real
! C                       Function subprogram defining the integrand
! C                       FUNCTION F(X). The actual name for F needs to be
! C                       Declared E X T E R N A L in the calling program.
! C
! C              BOUN   - Real
! C                       Finite bound of original integration
! C                       Range (SET TO ZERO IF INF = +2)
! C
! C              INF    - Integer
! C                       If INF = -1, the original interval is
! C                                   (-INFINITY,BOUND),
! C                       If INF = +1, the original interval is
! C                                   (BOUND,+INFINITY),
! C                       If INF = +2, the original interval is
! C                                   (-INFINITY,+INFINITY) AND
! C                       The integral is computed as the sum of two
! C                       integrals, one over (-INFINITY,0) and one over
! C                       (0,+INFINITY).
! C
! C              A      - Real
! C                       Lower limit for integration over subrange
! C                       of (0,1)
! C
! C              B      - Real
! C                       Upper limit for integration over subrange
! C                       of (0,1)
! C
! C            ON RETURN
! C              RESULT - Real
! C                       Approximation to the integral I
! C                       Result is computed by applying the 15-POINT
! C                       KRONROD RULE(RESK) obtained by optimal addition
! C                       of abscissae to the 7-POINT GAUSS RULE(RESG).
! C
! C              ABSERR - Real
! C                       Estimate of the modulus of the absolute error,
! C                       WHICH SHOULD EQUAL or EXCEED ABS(I-RESULT)
! C
! C              RESABS - Real
! C                       Approximation to the integral J
! C
! C              RESASC - Real
! C                       Approximation to the integral of
! C                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) over (A,B)
! C
! C***REFERENCES  (NONE)
! C***ROUTINES CALLED  R1MACH
! C***REVISION HISTORY  (YYMMDD)
! C   800101  DATE WRITTEN
! C   890531  Changed all specific intrinsics to generic.  (WRB)
! C   890531  REVISION DATE from Version 3.2
! C   891214  Prologue converted to Version 4.0 format.  (BAB)
! C***END PROLOGUE  QK15I
! C
      REAL A,ABSC,ABSC1,ABSC2,ABSERR,B,BOUN,CENTR,&
     &  DINF,R1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,
     &  FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,TABSC1,TABSC2,
     &  UFLOW,WG,WGK,XGK
      INTEGER INF,J
      EXTERNAL F
! C
      DIMENSION FV1(7),FV2(7),XGK(8),WGK(8),WG(8)
! C
! C           THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
! C           (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
! C           THEIR CORRESPONDING WEIGHTS ARE GIVEN.
! C
! C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
! C                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
! C                    GAUSS RULE
! C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
! C                    ADDED TO THE 7-POINT GAUSS RULE
! C
! C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
! C
! C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING
! C                    TO THE ABSCISSAE XGK(2), XGK(4), ...
! C                    WG(1), WG(3), ... ARE SET TO ZERO.
! C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),
     1  XGK(8)/
     2     0.9914553711208126E+00,     0.9491079123427585E+00,
     3     0.8648644233597691E+00,     0.7415311855993944E+00,
     4     0.5860872354676911E+00,     0.4058451513773972E+00,
     5     0.2077849550078985E+00,     0.0000000000000000E+00/
! C
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),
     1  WGK(8)/
     2     0.2293532201052922E-01,     0.6309209262997855E-01,
     3     0.1047900103222502E+00,     0.1406532597155259E+00,
     4     0.1690047266392679E+00,     0.1903505780647854E+00,
     5     0.2044329400752989E+00,     0.2094821410847278E+00/
! C
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.0000000000000000E+00,     0.1294849661688697E+00,
     2     0.0000000000000000E+00,     0.2797053914892767E+00,
     3     0.0000000000000000E+00,     0.3818300505051189E+00,
     4     0.0000000000000000E+00,     0.4179591836734694E+00/
! C
! C
! C           LIST OF MAJOR VARIABLES
! C           -----------------------
! C
! C           CENTR  - MID POINT OF THE INTERVAL
! C           HLGTH  - HALF-LENGTH OF THE INTERVAL
! C           ABSC*  - ABSCISSA
! C           TABSC* - TRANSFORMED ABSCISSA
! C           FVAL*  - FUNCTION VALUE
! C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
! C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
! C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
! C                    INTEGRAND OVER (A,B), I.E. TO I/(B-A)
! C
! C           MACHINE DEPENDENT CONSTANTS
! C           ---------------------------
! C
! C           EPMACH IS THE LARGEST RELATIVE SPACING.
! C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
! C
! C***FIRST EXECUTABLE STATEMENT  QK15I
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
      DINF = MIN(1,INF)
! C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      TABSC1 = BOUN+DINF*(0.1E+01-CENTR)/CENTR
      FVAL1 = F(TABSC1)
      IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
      FC = (FVAL1/CENTR)/CENTR
! C
! C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
! C           THE INTEGRAL, AND ESTIMATE THE ERROR.
! C
      RESG = WG(8)*FC
      RESK = WGK(8)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,7
        ABSC = HLGTH*XGK(J)
        ABSC1 = CENTR-ABSC
        ABSC2 = CENTR+ABSC
        TABSC1 = BOUN+DINF*(0.1E+01-ABSC1)/ABSC1
        TABSC2 = BOUN+DINF*(0.1E+01-ABSC2)/ABSC2
        FVAL1 = F(TABSC1)
        FVAL2 = F(TABSC2)
        IF(INF.EQ.2) FVAL1 = FVAL1+F(-TABSC1)
        IF(INF.EQ.2) FVAL2 = FVAL2+F(-TABSC2)
        FVAL1 = (FVAL1/ABSC1)/ABSC1
        FVAL2 = (FVAL2/ABSC2)/ABSC2
        FV1(J) = FVAL1
        FV2(J) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(J)*FSUM
        RESABS = RESABS+WGK(J)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(8)*ABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESASC = RESASC*HLGTH
      RESABS = RESABS*HLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.E0) ABSERR = RESASC*
     1 MIN(0.1E+01,(0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1 ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
      
