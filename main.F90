!************************************************************


Fixed bug in First release!!!!!!!!
Dut dut dut release 2!!:wq
Okay final release 2
Ding!
!*
!*  MAIN
!*
!************************************************************
	
program main
      implicit none
      include 'runhydro.h'
      
!************************************************************
!*
!*  Global Variables

      real, dimension(numr,numz,numphi) :: pot, rho
      common /poisson/ pot, rho

      real, dimension(numr,numz,numphi) :: psi, enth

!*
!************************************************************      
!*
!*   Local variables
      real :: w, phi_a, phi_b, h_a, h_b, psi_a,psi_b,phi_c, rho_c, h_0, c_0, rho_norm, h_0sq, h_max
      integer :: i,j,k,count
      real :: delta_c,delta_h, c_prev, h_prev,re
!* 
!************************************************************    
      print*, "The polytropic index = ", np

!Guess the initial density
      call guessrho
      
!Find rotational potential	
      do i=1,numr
        do j=1,numz
          do k=1,numphi
            w=1.0/ax*(i-1.5)
            psi(i,j,k)=-w**2/2
            if (j==1) then
              print*,"w=",w, "psi=",psi(i,j,k), i,j,k    
            endif
          enddo
        enddo
      enddo  

      
      call print2d(psi,psi.txt)
!      open(unit=10,file='psix.txt')
!        do i=1,numr  
!          write(10,*) psi(i,2,1) 
!        enddo
!      close(10)
!      print*,"Cross section file psix.txt printed"

     stop 

!Find potential and normalize     
      call poisson_solve
      Re=ax*1.0/(numr)
      pot=pot/Re**2
      
!      open(unit=10,file='solx.txt')
!        do i=1,numr  
!          write(10,*) pot(i,2,1) 
!        enddo
!      close(10)
!      print*,"Cross section file solx.txt printed"
      
!      open(unit=10,file='soly.txt')
!        do i=1,numz  
!          write(10,*) pot(2,i,1) 
!        enddo
!      close(10)
!      print*,"Cross section file soly.txt printed"      
      
      
      
      
!Find the constants h_0 and C_0      
      phi_a=pot(ax,ay,1) 
      phi_b=pot(bx,by,1)
      psi_a=psi(ax,ay,1)
      psi_b=psi(bx,by,1)
      
!      print*,"phi_a",phi_a,"phi_b",phi_b, "phi_a-phi_b",phi_a-phi_b
!      print*,"psi_a",psi_a,"psi_b",psi_b, "psi_a-psi_b",psi_a-psi_b
!      print*, "(phi_a-phi_b)/(psi_a-psi_b)",(phi_a-phi_b)/(psi_a-psi_b)
      
!      if (phi_a.eq.phi_b) then 
!        c_0=phi_b
!        h_0=((c_0-h_a-phi_a)/psi_a)**0.5
!      else
	
	
!        h_0= (-1.0*(phi_a-phi_b)/(psi_a-psi_b))**(0.5)
        h_0= ((-1.0*(phi_a-phi_b)/(psi_a-psi_b)))**(0.5)
!	h_0=0.0
        c_0=phi_a+h_0**2*psi_a
        
        print*,"phi_a",phi_a
        
        
!      endif
	
	
!      c_0=phi_b
!      h_0=((c_0-h_a-phi_a)/psi_a)**0.5
	
      print*,"C_00",C_0, "h_00",h_0
      
      
!Get enthalpy      
      do i=1,numr
        do j=1,numz
          do k=1,numphi
            enth(i,j,k)=  C_0 - pot(i,j,k) - h_0**2 * psi(i,j,k)
          enddo
        enddo
      enddo  
      h_max=maxval(enth)
   
   
!Find the new normalized density      
      phi_c=pot(1,1,1)
      rho_c=rho(1,1,1)
!      print*,"phi_c",phi_c,"rho_c",rho_c
      
      K= (C_0 -phi_c)/(np+1)/rho_c**(1.0/np)
      
      do i=1,numr
        do j=1,numz
          do k=1,numphi
         !   rho(i,j,k)=(enth(i,j,k)/(np+1.0)/K)**np
	    if (enth(i,j,k).gt.0) then 
	      rho(i,j,k)=(enth(i,j,k)/h_max)**np
	    else
          !  if ((rho(i,j,k).lt.0).or.(i.gt.ax).or.(j.gt.by)) then
              rho(i,j,k)=0.0
            endif      
          enddo
        enddo
      enddo        
      
      rho_norm=rho(1,1,1)
!      rho=rho/rho_norm
      
      

      
      
!     open(unit=10,file='rho1.txt')
!      do j=1,numz
!        do i=1,numr  
!          write(10,*) i,j,rho(i,j,1) 
!        enddo
!        write(10,*)
!!      enddo
 !     close(10)
 !     print*,"First iteration rho1.txt printed"
      
      
      
      
      
      
!!!!!!!Iterate till Convergence!!!!!!!
      delta_c=1.0
      delta_h=1.0
      count=0
      
      do while ((delta_c .gt. 1d-4).and.(delta_h.gt.1d-4))
        count=count+1
        
        
        !Find rotational potential	
        !Generally changes with rho not with const omega
        
        
        !Poisson solve for density      
        call poisson_solve
        pot=pot/Re**2

!      open(unit=10,file='pot.txt')
!      do j=1,numz
 !       do i=1,numr  
  !        write(10,*) i,j,pot(i,j,1) 
   !     enddo
    !    write(10,*)
     ! enddo
      !close(10)
      !print*,"Intermediate potential pot.txt printed" 
        
!Find the constants h_0 and C_0      
      phi_a=pot(ax,ay,1) 
      phi_b=pot(bx,by,1)
      psi_a=psi(ax,ay,1)
      psi_b=psi(bx,by,1)
      
!      print*,"phi_a",phi_a,"phi_b",phi_b, "phi_a-phi_b",phi_a-phi_b
!      print*,"psi_a",psi_a,"psi_b",psi_b, "psi_a-psi_b",psi_a-psi_b
 !     print*, "(phi_a-phi_b)/(psi_a-phi_b)",(phi_a-phi_b)/(psi_a-phi_b)
      
      
!      if (phi_a.eq.phi_b) then 
!        c_0=phi_b
!        h_0=((c_0-h_a-phi_a)/psi_a)**0.5
!      else
        h_0= ((-1.0*(phi_a-phi_b)/(psi_a-psi_b)))**(0.5)
!        h_0=0.0 !Enabled for non rotation
        c_0=phi_a+h_0**2*psi_a
        
        print*,"phi_a",phi_a
        
!      endif
      
!      c_0=phi_b
!      h_0=((c_0-h_a-phi_a)/psi_a)**0.5
	
!      print*,"C_0",C_0, "h_0",h_0
      
      
      
        !Get enthalpy      
        do i=1,numr
          do j=1,numz
            do k=1,numphi
              enth(i,j,k)=  C_0 - pot(i,j,k) - h_0**2 * psi(i,j,k)
            enddo
          enddo
        enddo  
        
        h_max=maxval(enth)
   
        !Find the new normalized density      
        phi_c=pot(1,1,1)
        rho_c=rho(1,1,1)
!        print*,"phi_c",phi_c,"rho_c",rho_c
      
        K= (C_0 -phi_c)/(np+1)/rho_c**(1.0/np)
      
      do i=1,numr
        do j=1,numz
          do k=1,numphi
         !   rho(i,j,k)=(enth(i,j,k)/(np+1.0)/K)**np
	    if (enth(i,j,k).gt.0) then 
	      rho(i,j,k)=(enth(i,j,k)/h_max)**np
	    else
          !  if ((rho(i,j,k).lt.0).or.(i.gt.ax).or.(j.gt.by)) then
              rho(i,j,k)=0.0
            endif      
          enddo
        enddo
      enddo          
      
        rho_norm=rho(1,1,1)
      
!        rho=rho/rho_norm
      
        delta_c=abs(c_prev-c_0)
        delta_h=abs(h_prev-h_0)
        
        c_prev=c_0
        h_prev=h_0
        
        print*, "Iteration number = ",count
        print*,"C_0 = ",c_0, "h_0 = ",h_0
        print*,"delta_c = ",delta_c, "delta_h = ",delta_h
        
        
        
      enddo
      
      
      open(unit=10,file='res.txt')
      do j=1,numz
        do i=1,numr  
          write(10,*) i,j,rho(i,j,1) 
        enddo
        write(10,*)
      enddo
      close(10)
      print*,"First iteration density res.txt printed"
 
      
      open(unit=10,file='ii.txt')
        do i=1,numr  
          write(10,*) rho(1,i,1) 
        enddo
      close(10)
      print*,"Cross section file ii.txt printed"
      
      
      stop
      end program main

