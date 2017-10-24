!==============================================================================
!# Subroutine to calculate gradients in all **REAL** cells. All variables of W
!# change in each axis, where W:
!# W=(/rho,ux,uy,uz,p/)
!#
!# Thus we Have
!#	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!#						& grad_d_x,grad_d_y,grad_d_z, &
!#						& grad_ux_x,grad_uy_x,grad_uz_x, &
!#						& grad_ux_y,grad_uy_y,grad_uz_y, &
!#						& grad_ux_z,grad_uy_z,grad_uz_z, &
!#						& grad_p_x,grad_p_y,grad_p_z /)
!==============================================================================
subroutine get_gradients(time)
	use constants
	implicit none

		Real*8, intent(in) :: time
		Integer :: i,j,k						!Indices of cell being analysed
		Integer :: idum,jdum,kdum					!Dummy i,j,k indiced for inner loop
		Integer :: iL,iR,jL,jR,kL,kR					!Left and right values of i,j,k
		Character (len=1) :: xyz


!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)		

! Scroll over all **REAL** cells:
	i=2
	j=2
	k=1!2

	Do while (k==1)
		j=1
		Do while(j<=1)
			i=2
			Do while (i<=NcX-1)

!===============================================================================
!#####				X-Axis Analysis				   #####
!===============================================================================
!		Set axis selector to x-axis
				xyz="x"

!		Set dummy left and right indices
!		(x-axis):
				iL=i-1
				iR=i+1
				jL=j
				jR=j
				kL=k
				kR=k

				call gradients(xyz,i,j,k,iL,jL,kL,iR,jR,kR,time)

!===============================================================================
!#####				Y-axis Analysis				   #####
!===============================================================================
!!		Set axis selector to y-axis
!				xyz="y"
!
!!		Set dummy left and right indices
!!		(y-axis):
!				iL=i
!				iR=i
!				jL=j-1
!				jR=j+1
!				kL=k
!				kR=k
!
!				call gradients(xyz,i,j,k,iL,jL,kL,iR,jR,kR,time)
!
!===============================================================================
!#####				Z-Axis Analysis				   #####
!===============================================================================
!!		Set axis selector to z-axis
!				xyz="z"
!
!!		Set dummy left and right indices
!!		(z-axis):
!				iL=i
!				iR=i
!				jL=j
!				jR=j
!				kL=k-1
!				kR=k+1
!
!				call gradients(xyz,i,j,k,iL,jL,kL,iR,jR,kR,time)


!===============================================================================
!#####				Apply Slope Limiter			   #####
!===============================================================================

				Call slope_limiter(i,j,k,time)

				i=i+1
			enddo
			j=j+1
		enddo
		k=k+1
	enddo


end subroutine get_gradients

subroutine gradients(xyz,i,j,k,iL,jL,kL,iR,jR,kR,time)
	use constants
	implicit none

		Integer,intent(in) :: i,j,k					!Indices of cell being analysed
		Integer,intent(in) :: iL,iR,jL,jR,kL,kR				!Left and right values of i,j,k
		Real*8, intent(in) :: time
		Character(len=1), intent(in):: xyz
	
		Real*8 :: xM,yM,zM,dM,pM,uxM,uyM,uzM				!Quantities for analysed cell
		Real*8 :: xmin,ymin,zmin,dmin,pmin,uxmin,uymin,uzmin		!Quantities for cells to the left of analysed cell
		Real*8 :: xplus,yplus,zplus,dplus,pplus,uxplus,uyplus,uzplus	!Quantities for cell to right of analysed cell
		Real*8 :: grad_d_dum,grad_ux_dum,grad_uy_dum,grad_uz_dum,grad_p_dum	!Gradients of pressure, velocity and densit. Dummy values that are to be passed into the slope limiter
		Real*8 :: grad_d,grad_ux,grad_uy,grad_uz,grad_p			!FINAL Gradients of pressure, velocity and density
		Real*8 :: Dq							!Distance between left and right cell for any given axis 


!#	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!#						& grad_d_x,grad_d_y,grad_d_z, &
!#						& grad_ux_x,grad_uy_x,grad_uz_x, &
!#						& grad_ux_y,grad_uy_y,grad_uz_y, &
!#						& grad_ux_z,grad_uy_z,grad_uz_z, &
!#						& grad_p_x,grad_p_y,grad_p_z /)

!	Set values for analysed cell:
	dM=cells(i,j,k,rhon)

	uxM=cells(i,j,k,uxn)
	uyM=cells(i,j,k,uyn)
	uzM=cells(i,j,k,uzn)

	pM=cells(i,j,k,pn)

	xM=cells(i,j,k,xn)
	yM=cells(i,j,k,yn)
	zM=cells(i,j,k,zn)


!	Set values for cells left of analysed cell:
	dmin=cells(iL,jL,kL,rhon)

	uxmin=cells(iL,jL,kL,uxn)
	uymin=cells(iL,jL,kL,uyn)
	uzmin=cells(iL,jL,kL,uzn)

	pmin=cells(iL,jL,kL,pn)

	xmin=cells(iL,jL,kL,xn)
	ymin=cells(iL,jL,kL,yn)
	zmin=cells(iL,jL,kL,zn)

!	Set values for cells Right of analysed cell:
	dplus=cells(iR,jR,kR,rhon)

	uxplus=cells(iR,jR,kR,uxn)
	uyplus=cells(iR,jR,kR,uyn)
	uzplus=cells(iR,jR,kR,uzn)

	pplus=cells(iR,jR,kR,pn)

	xplus=cells(iR,jR,kR,xn)
	yplus=cells(iR,jR,kR,yn)
	zplus=cells(iR,jR,kR,zn)


!===============================================================================
!#####				Gradient Analysis			   #####
!===============================================================================

!-------------------------------------------------------------------------------
!#####				x-Axis					   #####
!-------------------------------------------------------------------------------
	If (xyz=="x") then
!	Calculate distance between  left and right cell:
		Dq = xplus - xmin

!	Calculate Cell wide gradients:
		grad_d = (dplus - dMin)/Dq
		grad_ux = (uxplus - uxMin)/Dq
		grad_uy = (uyplus - uyMin)/Dq
		grad_uz = (uzplus - uzMin)/Dq
		grad_p = (pplus - pMin)/Dq

!#	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!#						& grad_d_x,grad_d_y,grad_d_z, &
!#						& grad_ux_x,grad_uy_x,grad_uz_x, &
!#						& grad_ux_y,grad_uy_y,grad_uz_y, &
!#						& grad_ux_z,grad_uy_z,grad_uz_z, &
!#						& grad_p_x,grad_p_y,grad_p_z /)

!	Update cells:
		!grad_d_x:
		cells(i,j,k,gdxn) = grad_d
		!grad_ux_x:	
		cells(i,j,k,guxxn) = grad_ux
		!grad_uy_x:
		cells(i,j,k,guyxn) = grad_uy
		!grad_uz_x:	
		cells(i,j,k,guzxn) = 0.d0!grad_uz
		!grad_p_x
		cells(i,j,k,gpxn) = grad_p

	else if (xyz=="y") then
!-------------------------------------------------------------------------------
!#####				y-Axis					   #####
!-------------------------------------------------------------------------------
!	Calculate distance between  left and right cell:
		Dq = Abs(yplus - ymin)

!	Calculate Cell wide gradients:
		grad_d = (dplus - dMin)/Dq
		grad_ux = (uxplus - uxMin)/Dq
		grad_uy = (uyplus - uyMin)/Dq
		grad_uz = (uzplus - uzMin)/Dq
		grad_p = (pplus - pMin)/Dq


!#	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!#						& grad_d_x,grad_d_y,grad_d_z, &
!#						& grad_ux_x,grad_uy_x,grad_uz_x, &
!#						& grad_ux_y,grad_uy_y,grad_uz_y, &
!#						& grad_ux_z,grad_uy_z,grad_uz_z, &
!#						& grad_p_x,grad_p_y,grad_p_z /)

!	Update cells:
		!grad_d_y:
		cells(i,j,k,gdyn) = grad_d
		!grad_ux_y:	
		cells(i,j,k,guxyn) = grad_ux
		!grad_uy_y:	
		cells(i,j,k,guyyn) = grad_uy
		!grad_uz_y:	
		cells(i,j,k,guzyn) = 0.d0!grad_uz
		!grad_p_y
		cells(i,j,k,gpyn) = grad_p

	else if (xyz=="z") then
!-------------------------------------------------------------------------------
!#####				z-Axis					   #####
!-------------------------------------------------------------------------------
!	Calculate distance between  left and right cell:
		Dq = Abs(zplus-zmin)

!	Calculate Cell wide gradients:
		grad_d = (dplus - dMin)/Dq
		grad_ux = (uxplus - uxMin)/Dq
		grad_uy = (uyplus - uyMin)/Dq
		grad_uz = (uzplus - uzMin)/Dq
		grad_p = (pplus - pMin)/Dq


!#	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!#						& grad_d_x,grad_d_y,grad_d_z, &
!#						& grad_ux_x,grad_uy_x,grad_uz_x, &
!#						& grad_ux_y,grad_uy_y,grad_uz_y, &
!#						& grad_ux_z,grad_uy_z,grad_uz_z, &
!#						& grad_p_x,grad_p_y,grad_p_z /)

!	Update cells:
		!grad_d_z:
		cells(i,j,k,gdzn) = grad_d
		!grad_ux_z:	
		cells(i,j,k,guxzn) = grad_ux
		!grad_uy_z:	
		cells(i,j,k,guyzn) = grad_uy
		!grad_uz_z:	
		cells(i,j,k,guzzn) = 0.d0!grad_uz
		!grad_p_z
		cells(i,j,k,gpzn) = grad_p

!-------------------------------------------------------------------------------
!#####				Axis Indicator Failure			   #####
!-------------------------------------------------------------------------------
	else

		Print*,"Axis indicator failure [@Gradients]! Program will Abort!"
		STOP
	endif

end subroutine gradients
