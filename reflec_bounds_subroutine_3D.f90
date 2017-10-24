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
	Do while(k<=NcZ)
		j=1
		Do while(j<=NcY)

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
	Do while(k<=NcZ)
		j=1
		Do while(j<=NcY)

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

!	Face normal to: -ve y-axis
	i=1
	j=1
	k=1
	Do while(k<=NcZ)
		i=1
		Do while(i<=NcX)

!	Set values of this face equal to those +2 cell in y-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,j+1,k,mn)
			qx=+1.d0*cells(i,j+1,k,qxn)
			qy=-1.d0*cells(i,j+1,k,qyn)
			qz=+1.d0*cells(i,j+1,k,qzn)
			E = cells(i,j+1,k,En)
			rho = cells(i,j+1,k,rhon)
			ux = +1.d0*cells(i,j+1,k,uxn)
			uy = -1.d0*cells(i,j+1,k,uyn)
			uz = +1.d0*cells(i,j+1,k,uzn)
			p = cells(i,j+1,k,pn)
			V = cells(i,j+1,k,vn)

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
			i=i+1
		enddo
		k=k+1
	enddo

!	Face normal to: +ve y-axis
	i=1
	j=NcY
	k=1
	Do while(k<=NcZ)
		i=1
		Do while(i<=NcX)

!	Set values of this face equal to those NcY-1 cell in y-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,j-1,k,mn)
			qx=+1.d0*cells(i,j-1,k,qxn)
			qy=-1.d0*cells(i,j-1,k,qyn)
			qz=+1.d0*cells(i,j-1,k,qzn)
			E = cells(i,j-1,k,En)
			rho = cells(i,j-1,k,rhon)
			ux = +1.d0*cells(i,j-1,k,uxn)
			uy = -1.d0*cells(i,j-1,k,uyn)
			uz = +1.d0*cells(i,j-1,k,uzn)
			p = cells(i,j-1,k,pn)
			V = cells(i,j-1,k,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)
!
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
			i=i+1
		enddo
		k=k+1
	enddo


!==============================================================================
!*****			z-Axis Face Ghost Cells:			  ******
!==============================================================================

!	Face normal to: -ve z-axis
	i=1
	j=1
	k=1
	Do while(j<=NcY)
		i=1
		Do while(i<=NcX)

!	Set values of this face equal to those +2 cell in z-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,j,k+1,mn)
			qx=+1.d0*cells(i,j,k+1,qxn)
			qy=+1.d0*cells(i,j,k+1,qyn)
			qz=-1.d0*cells(i,j,k+1,qzn)
			E = cells(i,j,k+1,En)
			rho = cells(i,j,k+1,rhon)
			ux = +1.d0*cells(i,j,k+1,uxn)
			uy = +1.d0*cells(i,j,k+1,uyn)
			uz = -1.d0*cells(i,j,k+1,uzn)
			p = cells(i,j,k+1,pn)
			V = cells(i,j,k+1,Vn)

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

			i=i+1
		enddo
		j=j+1
	enddo

!	Face normal to: +ve z-axis
	i=1
	j=1
	k=NcZ
	Do while(j<=NcY)
		i=1
		Do while(i<=NcX)

!	Set values of this face equal to those NcZ-1 cell in z-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,j,k-1,mn)
			qx=+1.d0*cells(i,j,k-1,qxn)
			qy=+1.d0*cells(i,j,k-1,qyn)
			qz=-1.d0*cells(i,j,k-1,qzn)
			E = cells(i,j,k-1,En)
			rho = cells(i,j,k-1,rhon)
			ux = +1.d0*cells(i,j,k-1,uxn)
			uy = +1.d0*cells(i,j,k-1,uyn)
			uz = -1.d0*cells(i,j,k-1,uzn)
			p = cells(i,j,k-1,pn)
			V = cells(i,j,k-1,Vn)

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

			i=i+1
		enddo
		j=j+1
	enddo

end subroutine bounds
