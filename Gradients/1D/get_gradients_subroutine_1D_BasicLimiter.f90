!********************************************************************************
!* This file is part of:
!* BAFFLE (https://github.com/ATHannington/BAFFLE)
!* BAFFLE: BAsic Fortran FLuids Engine
!*
!* BAFFLE:
!* Copyright (C) 2017 Andrew Hannington (ath4@st-andrews.ac.uk)
!*
!* This software is a derivative of work by Bert Vandenbroucke
!* (bert.vandenbroucke@gmail.com)
!* Find more work by Bert Vandenbroucke at: (https://github.com/bwvdnbro)
!*
!* BAFFLE is a free software: you can redistribute it and/or modify it under the 
!* terms of the GNU Affero General Public License
!* as published by the Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* BAFFLE is distributed in the hope that it will 
!* be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!* GNU Affero General Public License for more details.
!*
!* You should have received a copy of the GNU Affero General Public License
!* along with BAFFLE. If not, see
!* <http://www.gnu.org/licenses/>.
!********************************************************************************

!********************************************************************************
!* @file get_gradients_subroutine_1D_BasicLimiter.f90
!*
!* @Gets gradients of primitive variables of cells for second order functionality.
!* Basic slope limiter algorithm version!
!* 1D Version: Fortran 90 version.
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!********************************************************************************


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
	j=1!2
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
		Real*8 :: xM,yM,zM,dM,pM,uxM,uyM,uzM				!Quantities for cell being analysed
		Real*8 :: xmin,ymin,zmin,dmin,pmin,uxmin,uymin,uzmin		!Quantities for cells to the left of analysed cell
		Real*8 :: xplus,yplus,zplus,dplus,pplus,uxplus,uyplus,uzplus	!Quantities for cell to right of analysed cell
		Real*8 :: grad_d_L,grad_ux_L,grad_uy_L,grad_uz_L,grad_p_L	!Gradients of pressure, velocity and density, for LEFT cell
		Real*8 :: grad_d_M,grad_ux_M,grad_uy_M,grad_uz_M,grad_p_M	!Gradients of pressure, velocity and density,for MIDDLE cell, i.e. the cell being analysed
		Real*8 :: grad_d_R,grad_ux_R,grad_uy_R,grad_uz_R,grad_p_R	!Gradients of pressure, velocity and density, for RIGHT cell
		Real*8 :: grad_d,grad_ux,grad_uy,grad_uz,grad_p			!FINAL Gradients of pressure, velocity and density
		Real*8 :: Dqmin,Dqplus						!Left and right distance to cell for any given axis 


!#	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!#						& grad_d_x,grad_d_y,grad_d_z, &
!#						& grad_ux_x,grad_uy_x,grad_uz_x, &
!#						& grad_ux_y,grad_uy_y,grad_uz_y, &
!#						& grad_ux_z,grad_uy_z,grad_uz_z, &
!#						& grad_p_x,grad_p_y,grad_p_z /)


!	Set values for cell being analysed:
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
!	Calculate distance between cell analysed and left and right cell:
		Dqmin = Abs(xM - xmin)
		Dqplus = Abs(xPlus - xM)

!	Calculate Left gradients:
		grad_d_L = (dM - dMin)/DqMin
		grad_ux_L = (uxM - uxMin)/DqMin
		grad_uy_L = (uyM - uyMin)/DqMin
		grad_uz_L = (uzM - uzMin)/DqMin
		grad_p_L =(pM - pMin)/DqMin

!	Calculate Right gradients:
		grad_d_R = (dPlus - dM)/Dqplus
		grad_ux_R = (uxPlus - uxM)/Dqplus
		grad_uy_R = (uyPlus - uyM)/Dqplus
		grad_uz_R = (uzPlus - uzM)/Dqplus
		grad_p_R = (pPlus - pM)/Dqplus

!-------------------------------------------------------------------------------
!	Gradient management. This step compares the two gradients and takes the smaller of the two
!	as the final gradient. This *should* help prevent over-estimation of variables to non-physical
!	values. This is the least invasive method, as all other methods interfere with the hydrodynamics.
!	
!	Implementinmg slope_limiter slope limiter. Returns "smallest" of absolute values if signs are same,
!	and zero if signs are opposite
!-------------------------------------------------------------------------------
!	Density Gradient:
		call slope_limiter(i,j,k,grad_d_L,grad_d_R,grad_d)
!	Ux Gradient:
		call slope_limiter(i,j,k,grad_ux_L,grad_ux_R,grad_ux)
!	Uy Gradient:
		call slope_limiter(i,j,k,grad_uy_L,grad_uy_R,grad_uy)
!	Uz Gradient:
		call slope_limiter(i,j,k,grad_uz_L,grad_uz_R,grad_uz)
!	Pressure Gradient:
		call slope_limiter(i,j,k,grad_p_L,grad_p_R,grad_p)
!-------------------------------------------------------------------------------



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
Print*,
Print*,i,j,k
Print*,grad_p

	else if (xyz=="y") then
!-------------------------------------------------------------------------------
!#####				y-Axis					   #####
!-------------------------------------------------------------------------------
!	Calculate distance between cell analysed and left and right cell:
		Dqmin = yM - ymin
		Dqplus = yPlus - yM

!	Calculate Left gradients:
		grad_d_L = (dM - dMin)/DqMin
		grad_ux_L = (uxM - uxMin)/DqMin
		grad_uy_L = (uyM - uyMin)/DqMin
		grad_uz_L = (uzM - uzMin)/DqMin
		grad_p_L = (pM - pMin)/DqMin

!	Calculate Right gradients:
		grad_d_R = (dPlus - dM)/Dqplus
		grad_ux_R = (uxPlus - uxM)/Dqplus
		grad_uy_R = (uyPlus - uyM)/Dqplus
		grad_uz_R = (uzPlus - uzM)/Dqplus
		grad_p_R = (pPlus - pM)/Dqplus

!-------------------------------------------------------------------------------
!	Gradient management. This step compares the two gradients and takes the smaller of the two
!	as the final gradient. This *should* help prevent over-estimation of variables to non-physical
!	values. This is the least invasive method, as all other methods interfere with the hydrodynamics.
!	
!	Implementinmg slope_limiter slope limiter. Returns "smallest" of absolute values if signs are same,
!	and zero if signs are opposite
!-------------------------------------------------------------------------------
!	Density Gradient:
		call slope_limiter(i,j,k,grad_d_L,grad_d_R,grad_d)
!	Ux Gradient:
		call slope_limiter(i,j,k,grad_ux_L,grad_ux_R,grad_ux)
!	Uy Gradient:
		call slope_limiter(i,j,k,grad_uy_L,grad_uy_R,grad_uy)
!	Uz Gradient:
		call slope_limiter(i,j,k,grad_uz_L,grad_uz_R,grad_uz)
!	Pressure Gradient:
		call slope_limiter(i,j,k,grad_p_L,grad_p_R,grad_p)
!-------------------------------------------------------------------------------

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
!	Calculate distance between cell analysed and left and right cell:
		Dqmin = zM - zmin
		Dqplus = zPlus - zM

!	Calculate Left gradients:
		grad_d_L = (dM - dMin)/DqMin
		grad_ux_L = (uxM - uxMin)/DqMin
		grad_uy_L = (uyM - uyMin)/DqMin
		grad_uz_L = (uzM - uzMin)/DqMin
		grad_p_L = (pM - pMin)/DqMin

!	Calculate Right gradients:
		grad_d_R = (dPlus - dM)/Dqplus
		grad_ux_R = (uxPlus - uxM)/Dqplus
		grad_uy_R = (uyPlus - uyM)/Dqplus
		grad_uz_R = (uzPlus - uzM)/Dqplus
		grad_p_R = (pPlus - pM)/Dqplus

!-------------------------------------------------------------------------------
!	Gradient management. This step compares the two gradients	 and takes the smaller of the two
!	as the final gradient. This *should* help prevent over-estimation of variables to non-physical
!	values. This is the least invasive method, as all other methods interfere with the hydrodynamics.
!	
!	Implementinmg slope_limiter slope limiter. Returns "smallest" of absolute values if signs are same,
!	and zero if signs are opposite
!-------------------------------------------------------------------------------
!	Density Gradient:
		call slope_limiter(i,j,k,grad_d_L,grad_d_R,grad_d)
!	Ux Gradient:
		call slope_limiter(i,j,k,grad_ux_L,grad_ux_R,grad_ux)
!	Uy Gradient:
		call slope_limiter(i,j,k,grad_uy_L,grad_uy_R,grad_uy)
!	Uz Gradient:
		call slope_limiter(i,j,k,grad_uz_L,grad_uz_R,grad_uz)
!	Pressure Gradient:
		call slope_limiter(i,j,k,grad_p_L,grad_p_R,grad_p)
!-------------------------------------------------------------------------------

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
		cells(i,j,k,guzzn) = grad_uz
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

!-------------------------------------------------------------------------------
! Implementinmg slope_limiter slope limiter. Returns "smallest" of absolute values if signs are same,
! and zero if signs are opposite
!-------------------------------------------------------------------------------
subroutine slope_limiter(i,j,k,grad_L,grad_R,grad)
	use constants
	implicit none

		Integer,intent(in) :: i,j,k
		Real*8, intent(in) :: grad_L, grad_R					
		Real*8, intent(out) :: grad

	If ((i==2).or.(j==2)) then !.or.(k==2)) then
		grad=grad_R
	else if ((i==NcX-1).or.(j==NcY-1)) then !.or.(k==NcZ-1)) then
		grad=grad_L
	else if ((grad_L<0.d0).and.(grad_R>0.d0)) then
		grad=0.d0
	else if ((grad_R<0.d0).and.(grad_L>0.d0)) then
		grad=0.d0
	else if (abs(grad_L)<abs(grad_R)) then
		grad= grad_L
	else if (abs(grad_L)>abs(grad_R)) then
		grad = grad_R
	else if ((grad_L).eq.(grad_R)) then
		grad=grad_L
	else
		Print*,"Slope Limiter Failure! Program will Abort!"
		STOP
	endif

end subroutine slope_limiter
