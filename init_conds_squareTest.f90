!===============================================================================
!
!*** 				INITIAL CONDITIONS			     ***
!
!-------------------------------------------------------------------------------
!Subroutine for initial conditions functions for density, velocity, and pressure
!& as a function of x.
!-------------------------------------------------------------------------------
subroutine init_conds(x,y,z,rho,ux,uy,uz,p)			
	use constants								!Tells program to use global constants
	implicit none

		Real*8, intent(in) :: x,y,z
		Real*8, intent(out) :: rho,ux,uy,uz,p
		Real*8 :: delx,dely,delz					!Ranges in x, y and z for the maximum range limit

!	Cells format: 	cells(i,j,k,1:14) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz/)
!	uy=0.d0
	uz=0.d0


!Sets values of pressure, x-velocity, and density according to starting conditions for 2D Sod Shock
	If ((x>=0.25d0).and.(x<=0.75d0).and.(y>=0.25d0).and.(y<=0.75d0)) then
		p=2.5d0
		rho=4.0d0
	else
		p=2.5d0
		rho=1.0d0
	endif
	
	ux=0.5d0
	uy=0.5d0

end subroutine init_conds
!===============================================================================
