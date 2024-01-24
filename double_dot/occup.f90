module occup

implicit none

public :: dirac, fermi  ! Make the function public

contains

	!ALWAYS DEFINE THE FUNCTIONS AT FIRST
	function dirac(a,b)
	integer a, b, dirac
	if(a.eq.b) then
		dirac=1
	else
		dirac=0
	endif
	end function

	function fermi(x,kT)
	real(8) x, kT, fermi

	fermi=1.d0/(dexp(x/kT)+1.d0)

	end function
	
	subroutine normal_distribution_example(mu, sigma, normal_values)

	  real(kind=8), intent(in) :: mu
	  real(kind=8), intent(in) :: sigma
	  integer, parameter :: num_values = 2
	  real(kind=8), dimension(num_values) :: random_values
	  real(kind=8), intent(out) :: normal_values
	  
	  ! Generate uniform random numbers
	  call random_number(random_values)

	  ! Transform uniform random numbers into normal random numbers
	  normal_values = mu + sigma * sqrt(-2.0d0 * log(random_values(1))) * cos(2.0d0 * acos(-1.0d0) * random_values(2))
	  !normal_values(2) = mu + sigma * sqrt(-2.0d0 * log(random_values(1))) * sin(2.0d0 * acos(-1.0d0) * random_values(2))

	end subroutine normal_distribution_example
	
	subroutine normal_noise(img_size, noise)
	  
	  integer(kind=8), intent(in) :: img_size
	  real(kind=8), intent(out) :: noise(img_size,img_size)
	  real(kind=8) :: n, random, sigma
	  integer(kind=8) :: i,j, mu
	  	   
	  !mu = 0
	  !sigma = 0.25
	  do i = 1, img_size
	  	do j = 1, img_size
	  		call random_number(random)
	  		n = (2.0*random)-1.0
	  		noise(i,j) = n
	  	enddo
	  enddo 
	
	end subroutine normal_noise
	
	subroutine ocupaciones(d, Ic)

	real(kind=8), intent(in) :: d(4)
	real(kind=8), intent(out) :: Ic(30,30)
	integer, parameter :: ne=29
	complex(8) im,Qx(9,9)
	real(8) pi,u,vb,eL,eR,tc,ep,ed,kT
	real(8) en(9),mu(2),v(2,2),Tm(9,9,2,2,-1:1),tL,tR
	real(8) WL(9,9),WR(9,9),Wx(9,9),WI(9,9),vx(9),p(9), Ic1
	real(8) a1,a2,emin,emax,m1,m2,nt
	integer ie1,ie2,i,j,nc(9)


	pi=acos(-1.d0)
	im=(0.d0,1.d0)



	kT=1

	u =30.*kT
	tc= 5.*kT

	nc(1)=0
	nc(2)=1
	nc(3)=1
	nc(4)=1
	nc(5)=1
	nc(6)=2
	nc(7)=2
	nc(8)=2
	nc(9)=2

	emin=-40*kT
	emax=+20*kT

	vb=15*kT

	mu(1)=+.5*vb
	mu(2)=-.5*vb
	
	open(1,file='q02.dat')

	do ie1=0,ne
		! write(*,*) ie1
		m1=emin+(emax-emin)*ie1/dble(ne)

		do ie2=0,ne
			
			m2=emin+(emax-emin)*ie2/dble(ne)
				
			eL=d(1)*m1+d(2)*m2
			eR=d(3)*m1+d(4)*m2
			!eL=0.83*m1+0.32*m2
			!eR=0.11*m1+1.19*m2

			ep=.5*(eL+eR)
			ed=eL-eR

			a1=1.-ed/sqrt(ed**2+tc**2)
			a2=1.+ed/sqrt(ed**2+tc**2)

			!Tunnel amplitudes

			Tm=0

			tL=.1*sqrt(kT/(8*pi))
			tR=.1*sqrt(kT/(8*pi))

			v(1,1)=+tL/sqrt(2.)*sqrt(a1)
			v(1,2)=-tL/sqrt(2.)*sqrt(a2)
			v(2,1)=+tR/sqrt(2.)*sqrt(a2)
			v(2,2)=+tR/sqrt(2.)*sqrt(a1)

			do i=1,2
				Tm(2,1,i,1,+1)=+v(i,1)
				Tm(1,2,i,1,-1)=+v(i,1)
				Tm(3,1,i,2,+1)=+v(i,1)
				Tm(1,3,i,2,-1)=+v(i,1)	
				Tm(4,1,i,1,+1)=+v(i,2)
				Tm(1,4,i,1,-1)=+v(i,2)
				Tm(5,1,i,2,+1)=+v(i,2)
				Tm(1,5,i,2,-1)=+v(i,2)

				Tm(6,2,i,1,+1)=-v(i,2)
				Tm(2,6,i,1,-1)=-v(i,2)
				Tm(7,2,i,2,+1)=-v(i,2)*dirac(i,2)
				Tm(2,7,i,2,-1)=-v(i,2)*dirac(i,2)
				Tm(8,2,i,2,+1)=-v(i,2)*dirac(i,1)
				Tm(2,8,i,2,-1)=-v(i,2)*dirac(i,1)

				Tm(7,3,i,1,+1)=-v(i,2)*dirac(i,1)
				Tm(3,7,i,1,-1)=-v(i,2)*dirac(i,1)
				Tm(8,3,i,1,+1)=-v(i,2)*dirac(i,2)
				Tm(3,8,i,1,-1)=-v(i,2)*dirac(i,2)
				Tm(9,3,i,2,+1)=-v(i,2)
				Tm(3,9,i,2,-1)=-v(i,2)

				Tm(6,4,i,1,+1)=+v(i,1)
				Tm(4,6,i,1,-1)=+v(i,1)
				Tm(7,4,i,2,+1)=+v(i,1)*dirac(i,2)
				Tm(4,7,i,2,-1)=+v(i,1)*dirac(i,2)
				Tm(8,4,i,2,+1)=+v(i,1)*dirac(i,1)
				Tm(4,8,i,2,-1)=+v(i,1)*dirac(i,1)

				Tm(7,5,i,1,+1)=+v(i,1)*dirac(i,1)
				Tm(5,7,i,1,-1)=+v(i,1)*dirac(i,1)
				Tm(8,5,i,1,+1)=+v(i,1)*dirac(i,2)
				Tm(5,8,i,1,-1)=+v(i,1)*dirac(i,2)
				Tm(9,5,i,2,+1)=+v(i,1)
				Tm(5,9,i,2,-1)=+v(i,1)	
			enddo

			! Energies

			en(1)=0.
			en(2)=ep-.5*sqrt(ed**2+tc**2)
			en(3)=ep-.5*sqrt(ed**2+tc**2)
			en(4)=ep+.5*sqrt(ed**2+tc**2)
			en(5)=ep+.5*sqrt(ed**2+tc**2)
			en(6)=2*ep+u
			en(7)=2*ep+u
			en(8)=2*ep+u
			en(9)=2*ep+u

			call kernel(en,kT,mu,1,Tm,WL)
			call kernel(en,kT,mu,2,Tm,WR)

			do i=1,9
				do j=1,9
					Wx(i,j)=WL(i,j)+WR(i,j)-WL(i,i)-WR(i,i)
					WI(i,j)=(nc(i)-nc(j))*WL(i,j)
				enddo
				vx(i)=-WL(i,i)-WR(i,i)
			enddo

			Qx=Wx
			do i=1,9
				Qx(i,i)=im*1.e-5
			enddo
			call invert(9,Qx)		
				
			p=matmul(dble(Qx),vx)
				
			Ic1=0
			!nt=0
				
			do i=1,9			
				do j=1,9
					Ic1=Ic1+WI(i,j)*p(j)			
				enddo
				
				!nt=nt+nc(i)*p(i)
					
			enddo
				
			Ic(ie1+1,ie2+1)=Ic1 
			!write(1,*) m1, m2, Ic1
				
		enddo
			
		!write(1,*)
			
	enddo
	write(1,*) Ic

	end subroutine ocupaciones


	subroutine kernel(en,kT,mu,r,T,W)

		implicit none

		real(8) W(9,9),kT
		real(8) pi,E(9,9),en(9),mu(2)
		real(8) T(9,9,2,2,-1:1),Q
		integer a0(-1:1),a1(-1:1),a2(-1:1)
		integer i1m,i1p,i0,i2
		integer D0,p1,p2,eta1,r,s1
		integer i,j

		pi=acos(-1.d0)
		E=0

		do i=1,9
			do j=1,9
				E(i,j)=en(i)-en(j)
			enddo
		enddo

		W=0

		do i0=1,9
		do i2=1,9

			a0(-1)=i0
			a0(+1)=i0
			a2(-1)=i2
			a2(+1)=i2

			do p1=-1,1,2
			do p2=-1,1,2

			!do r1=1,2
			do eta1=-1,1,2

			do i1m=1,9
			do i1p=1,9
			
			a1(-1)=i1m
			a1(+1)=i1p
			
			D0=dirac(a2(-p2),a1(-p2))*dirac(a1(-p1),a0(-p1))

			if(D0.ne.0) then

				Q=0
				do s1=1,2
					Q=Q+T(a2(p2),a1(p2),r,s1,-eta1*p2)*T(a1(p1),a0(p1),r,s1,eta1*p1)
				enddo

				if(Q.ne.0) W(i2,i0)=W(i2,i0)-D0*Q*p2*p1*pi*fermi(p1*(E(a1(1),a1(-1))-eta1*mu(r)),kT)			

			endif

			enddo
			enddo
			
			enddo
			!enddo
			
			enddo
			enddo
			
		enddo
		enddo

	end subroutine

	subroutine mult(N,a1,a2,v)

	integer N
	real(8) a1(N,N),a2(N),v(N)
	integer i,j

	v=0
	do i=1,N
		do j=1,N
			v(i)=v(i)+a1(i,j)*a2(j)
		enddo
	enddo

	end subroutine


	subroutine invert(N,A)

	integer INFO, LDA, LWORK, N
	integer IPIV(N)
	complex(8) A(N,N), WORK(N)

	LDA=N
	LWORK=N
	      	
	CALL ZGETRF(N,N,A,LDA,IPIV,INFO)
	CALL ZGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
	    	    
	end subroutine

end module occup
