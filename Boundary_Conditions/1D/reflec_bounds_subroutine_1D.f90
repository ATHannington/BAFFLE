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
!* 
!* Find more work by Bert Vandenbroucke at: (https://github.com/bwvdnbro)
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
!* @file reflec_bounds_subroutine_1D.f90
!*
!* @Reflective Boundary Conditions Subroutine in 1D: Fortran 90 version.
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!********************************************************************************

!-------------------------------------------------------------------------------
!Compact subroutine for reflective boundary conditions. Reverses the direction of
!velocity and momentum
!-------------------------------------------------------------------------------
subroutine bounds							

	use constants
	implicit none
		integer :: i,j,k,d
		Real*8 :: m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
			& grad_d_x,grad_d_y,grad_d_z, &
			& grad_ux_x,grad_uy_x,grad_uz_x, &
			& grad_ux_y,grad_uy_y,grad_uz_y, &
			& grad_ux_z,grad_uy_z,grad_uz_z, &
			& grad_p_x,grad_p_y,grad_p_z

!	Cells format: 	cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)	

	d=1									!Counter to move over variables within a cell for two ghost cells
	i=1
	j=1
	k=1

!==============================================================================
!*****			X-Axis Face Ghost Cells:			  ******
!==============================================================================

!	Face normal to: -ve x-axis
	i=1
	j=1
	k=1
	Do while(k<=1)
		j=1
		Do while(j<=1)!NcY)

!	Set values of this face equal to those +2 cell in x-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i+1,j,k,mn)
			qx=-1.d0*cells(i+1,j,k,qxn)
			qy=+1.d0*cells(i+1,j,k,qyn)
			qz=+1.d0*cells(i+1,j,k,qzn)
			E = cells(i+1,j,k,En)
			rho = cells(i+1,j,k,rhon)
			ux = -1.d0*cells(i+1,j,k,uxn)
			uy = +1.d0*cells(i+1,j,k,uyn)
			uz = +1.d0*cells(i+1,j,k,uzn)
			p = cells(i+1,j,k,pn)
			V = cells(i+1,j,k,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients to zero for reflective boundary:
			grad_d_x= 0.d0
			grad_d_y=0.d0
			grad_d_z=0.d0
			grad_ux_x=0.d0
			grad_uy_x=0.d0
			grad_uz_x=0.d0
			grad_ux_y=0.d0
			grad_uy_y=0.d0
			grad_uz_y=0.d0
			grad_ux_z=0.d0
			grad_uy_z=0.d0
			grad_uz_z=0.d0
			grad_p_x=0.d0
			grad_p_y=0.d0
			grad_p_z=0.d0

			cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
						& grad_d_x,grad_d_y,grad_d_z, &
						& grad_ux_x,grad_uy_x,grad_uz_x, &
						& grad_ux_y,grad_uy_y,grad_uz_y, &
						& grad_ux_z,grad_uy_z,grad_uz_z, &
						& grad_p_x,grad_p_y,grad_p_z /)
			j=j+1
		enddo
		k=k+1
	enddo

!	Face normal to: +ve x-axis
	i=NcX
	j=1
	k=1
	Do while(k<=1)
		j=1
		Do while(j<=1)!NcY)

!	Set values of this face equal to those NcX-1 cell in x-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i-1,j,k,mn)
			qx=-1.d0*cells(i-1,j,k,qxn)
			qy=+1.d0*cells(i-1,j,k,qyn)
			qz=+1.d0*cells(i-1,j,k,qzn)
			E = cells(i-1,j,k,En)
			rho = cells(i-1,j,k,rhon)
			ux = -1.d0*cells(i-1,j,k,uxn)
			uy = +1.d0*cells(i-1,j,k,uyn)
			uz = +1.d0*cells(i-1,j,k,uzn)
			p = cells(i-1,j,k,pn)
			V = cells(i-1,j,k,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients to zero for reflective boundary:
			grad_d_x= 0.d0
			grad_d_y=0.d0
			grad_d_z=0.d0
			grad_ux_x=0.d0
			grad_uy_x=0.d0
			grad_uz_x=0.d0
			grad_ux_y=0.d0
			grad_uy_y=0.d0
			grad_uz_y=0.d0
			grad_ux_z=0.d0
			grad_uy_z=0.d0
			grad_uz_z=0.d0
			grad_p_x=0.d0
			grad_p_y=0.d0
			grad_p_z=0.d0

			cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
						& grad_d_x,grad_d_y,grad_d_z, &
						& grad_ux_x,grad_uy_x,grad_uz_x, &
						& grad_ux_y,grad_uy_y,grad_uz_y, &
						& grad_ux_z,grad_uy_z,grad_uz_z, &
						& grad_p_x,grad_p_y,grad_p_z /)
			j=j+1
		enddo
		k=k+1
	enddo

!==============================================================================
!*****			y-Axis Face Ghost Cells:			  ******
!==============================================================================
!
!!	Face normal to: -ve y-axis
!	i=1
!	j=1
!	k=1
!	Do while(k<=1)
!		i=1
!		Do while(i<=NcX)
!
!!	Set values of this face equal to those +2 cell in y-direction. Equal and opposite for x-vel. and x-mom.: 
!			m=cells(i,j+1,k,mn)
!			qx=+1.d0*cells(i,j+1,k,qxn)
!			qy=-1.d0*cells(i,j+1,k,qyn)
!			qz=+1.d0*cells(i,j+1,k,qzn)
!			E = cells(i,j+1,k,En)
!			rho = cells(i,j+1,k,rhon)
!			ux = +1.d0*cells(i,j+1,k,uxn)
!			uy = -1.d0*cells(i,j+1,k,uyn)
!			uz = +1.d0*cells(i,j+1,k,uzn)
!			p = cells(i,j+1,k,pn)
!			V = cells(i,j+1,k,vn)
!
!!	Set position of faces equal to those already assigned in Cell Setup:
!			Midx=cells(i,j,k,xn)
!			Midy=cells(i,j,k,yn)
!			Midz=cells(i,j,k,zn)
!
!!	Set gradients to zero for reflective boundary:
!			grad_d_x= 0.d0
!			grad_d_y=0.d0
!			grad_d_z=0.d0
!			grad_ux_x=0.d0
!			grad_uy_x=0.d0
!			grad_uz_x=0.d0
!			grad_ux_y=0.d0
!			grad_uy_y=0.d0
!			grad_uz_y=0.d0
!			grad_ux_z=0.d0
!			grad_uy_z=0.d0
!			grad_uz_z=0.d0
!			grad_p_x=0.d0
!			grad_p_y=0.d0
!			grad_p_z=0.d0
!
!			cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)
!			i=i+1
!		enddo
!		k=k+1
!	enddo
!
!!	Face normal to: +ve y-axis
!	i=1
!	j=NcY
!	k=1
!	Do while(i<=1)
!		i=1
!		Do while(i<=NcX)
!
!!	Set values of this face equal to those NcY-1 cell in y-direction. Equal and opposite for x-vel. and x-mom.: 
!			m=cells(i,j-1,k,mn)
!			qx=+1.d0*cells(i,j-1,k,qxn)
!			qy=-1.d0*cells(i,j-1,k,qyn)
!			qz=+1.d0*cells(i,j-1,k,qzn)
!			E = cells(i,j-1,k,En)
!			rho = cells(i,j-1,k,rhon)
!			ux = +1.d0*cells(i,j-1,k,uxn)
!			uy = -1.d0*cells(i,j-1,k,uyn)
!			uz = +1.d0*cells(i,j-1,k,uzn)
!			p = cells(i,j-1,k,pn)
!			V = cells(i,j-1,k,Vn)
!
!!	Set position of faces equal to those already assigned in Cell Setup:
!			Midx=cells(i,j,k,xn)
!			Midy=cells(i,j,k,yn)
!			Midz=cells(i,j,k,zn)
!!
!!	Set gradients to zero for reflective boundary:
!			grad_d_x= 0.d0
!			grad_d_y=0.d0
!			grad_d_z=0.d0
!			grad_ux_x=0.d0
!			grad_uy_x=0.d0
!			grad_uz_x=0.d0
!			grad_ux_y=0.d0
!			grad_uy_y=0.d0
!			grad_uz_y=0.d0
!			grad_ux_z=0.d0
!			grad_uy_z=0.d0
!			grad_uz_z=0.d0
!			grad_p_x=0.d0
!			grad_p_y=0.d0
!			grad_p_z=0.d0
!
!			cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)
!			i=i+1
!		enddo
!		k=k+1
!	enddo
!
!
!==============================================================================
!*****			z-Axis Face Ghost Cells:			  ******
!==============================================================================
!
!!	Face normal to: -ve z-axis
!	i=1
!	j=1
!	k=1
!	Do while(j<=NcY)
!		i=1
!		Do while(i<=NcX)
!
!!	Set values of this face equal to those +2 cell in z-direction. Equal and opposite for x-vel. and x-mom.: 
!			m=cells(i,j,k+1,mn)
!			qx=+1.d0*cells(i,j,k+1,qxn)
!			qy=+1.d0*cells(i,j,k+1,qyn)
!			qz=-1.d0*cells(i,j,k+1,qzn)
!			E = cells(i,j,k+1,En)
!			rho = cells(i,j,k+1,rhon)
!			ux = +1.d0*cells(i,j,k+1,uxn)
!			uy = +1.d0*cells(i,j,k+1,uyn)
!			uz = -1.d0*cells(i,j,k+1,uzn)
!			p = cells(i,j,k+1,pn)
!			V = cells(i,j,k+1,Vn)
!
!!	Set position of faces equal to those already assigned in Cell Setup:
!			Midx=cells(i,j,k,xn)
!			Midy=cells(i,j,k,yn)
!			Midz=cells(i,j,k,zn)
!
!!	Set gradients to zero for reflective boundary:
!			grad_d_x= 0.d0
!			grad_d_y=0.d0
!			grad_d_z=0.d0
!			grad_ux_x=0.d0
!			grad_uy_x=0.d0
!			grad_uz_x=0.d0
!			grad_ux_y=0.d0
!			grad_uy_y=0.d0
!			grad_uz_y=0.d0
!			grad_ux_z=0.d0
!			grad_uy_z=0.d0
!			grad_uz_z=0.d0
!			grad_p_x=0.d0
!			grad_p_y=0.d0
!			grad_p_z=0.d0
!
!			cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)
!
!			i=i+1
!		enddo
!		j=j+1
!	enddo
!
!!	Face normal to: +ve z-axis
!	i=1
!	j=1dot_product(vector_a, vector_b)
!	k=NcZ
!	Do while(j<=NcY)
!		i=1
!		Do while(i<=NcX)
!
!!	Set values of this face equal to those NcZ-1 cell in z-direction. Equal and opposite for x-vel. and x-mom.: 
!			m=cells(i,j,k-1,mn)
!			qx=+1.d0*cells(i,j,k-1,qxn)
!			qy=+1.d0*cells(i,j,k-1,qyn)
!			qz=-1.d0*cells(i,j,k-1,qzn)
!			E = cells(i,j,k-1,En)
!			rho = cells(i,j,k-1,rhon)
!			ux = +1.d0*cells(i,j,k-1,uxn)
!			uy = +1.d0*cells(i,j,k-1,uyn)
!			uz = -1.d0*cells(i,j,k-1,uzn)
!			p = cells(i,j,k-1,pn)
!			V = cells(i,j,k-1,Vn)
!
!!	Set position of faces equal to those already assigned in Cell Setup:
!			Midx=cells(i,j,k,xn)
!			Midy=cells(i,j,k,yn)
!			Midz=cells(i,j,k,zn)
!
!!	Set gradients to zero for reflective boundary:
!			grad_d_x= 0.d0
!			grad_d_y=0.d0
!			grad_d_z=0.d0
!			grad_ux_x=0.d0
!			grad_uy_x=0.d0
!			grad_uz_x=0.d0
!			grad_ux_y=0.d0
!			grad_uy_y=0.d0
!			grad_uz_y=0.d0
!			grad_ux_z=0.d0
!			grad_uy_z=0.d0
!			grad_uz_z=0.d0
!			grad_p_x=0.d0
!			grad_p_y=0.d0
!			grad_p_z=0.d0
!
!			cells(i,j,k,1:29) = (/m,qx,qy,qz,E,rho,ux,uy,uz,p,V,Midx,Midy,Midz, &
!						& grad_d_x,grad_d_y,grad_d_z, &
!						& grad_ux_x,grad_uy_x,grad_uz_x, &
!						& grad_ux_y,grad_uy_y,grad_uz_y, &
!						& grad_ux_z,grad_uy_z,grad_uz_z, &
!						& grad_p_x,grad_p_y,grad_p_z /)
!
!			i=i+1
!		enddo
!		j=j+1
!	enddo

end subroutine bounds
