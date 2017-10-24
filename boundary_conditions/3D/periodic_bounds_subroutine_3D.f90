
!-------------------------------------------------------------------------------
!Compact subroutine for periodic boundary conditions. Reverses the direction of
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

!	Set values of this face equal to those NcX-1 cell in x-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(NcX-1,j,k,mn)
			qx=+1.d0*cells(NcX-1,j,k,qxn)
			qy=+1.d0*cells(NcX-1,j,k,qyn)
			qz=+1.d0*cells(NcX-1,j,k,qzn)
			E = cells(NcX-1,j,k,En)
			rho = cells(NcX-1,j,k,rhon)
			ux = +1.d0*cells(NcX-1,j,k,uxn)
			uy = +1.d0*cells(NcX-1,j,k,uyn)
			uz = +1.d0*cells(NcX-1,j,k,uzn)
			p = cells(NcX-1,j,k,pn)
			V = cells(NcX-1,j,k,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients equal to those in NcX-1 cell in x-direction for periodic boundary:
			grad_d_x= cells(NcX-1,j,k,gdxn)
			grad_d_y=cells(NcX-1,j,k,gdyn)
			grad_d_z=cells(NcX-1,j,k,gdzn)
			grad_ux_x=cells(NcX-1,j,k,guxxn)
			grad_uy_x=cells(NcX-1,j,k,guyxn)
			grad_uz_x=cells(NcX-1,j,k,guzxn)
			grad_ux_y=cells(NcX-1,j,k,guxyn)
			grad_uy_y=cells(NcX-1,j,k,guyyn)
			grad_uz_y=cells(NcX-1,j,k,guzyn)
			grad_ux_z=cells(NcX-1,j,k,guxzn)
			grad_uy_z=cells(NcX-1,j,k,guyzn)
			grad_uz_z=cells(NcX-1,j,k,guzzn)
			grad_p_x=cells(NcX-1,j,k,gpxn)
			grad_p_y=cells(NcX-1,j,k,gpyn)
			grad_p_z=cells(NcX-1,j,k,gpzn)

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

!	Set values of this face equal to those +2 cell in x-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(2,j,k,mn)
			qx=+1.d0*cells(2,j,k,qxn)
			qy=+1.d0*cells(2,j,k,qyn)
			qz=+1.d0*cells(2,j,k,qzn)
			E = cells(2,j,k,En)
			rho = cells(2,j,k,rhon)
			ux = +1.d0*cells(2,j,k,uxn)
			uy = +1.d0*cells(2,j,k,uyn)
			uz = +1.d0*cells(2,j,k,uzn)
			p = cells(2,j,k,pn)
			V = cells(2,j,k,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients equal to those in +2 cell in x-direction for periodic boundary:

			grad_d_x= cells(2,j,k,gdxn)
			grad_d_y=cells(2,j,k,gdyn)
			grad_d_z=cells(2,j,k,gdzn)
			grad_ux_x=cells(2,j,k,guxxn)
			grad_uy_x=cells(2,j,k,guyxn)
			grad_uz_x=cells(2,j,k,guzxn)
			grad_ux_y=cells(2,j,k,guxyn)
			grad_uy_y=cells(2,j,k,guyyn)
			grad_uz_y=cells(2,j,k,guzyn)
			grad_ux_z=cells(2,j,k,guxzn)
			grad_uy_z=cells(2,j,k,guyzn)
			grad_uz_z=cells(2,j,k,guzzn)
			grad_p_x=cells(2,j,k,gpxn)
			grad_p_y=cells(2,j,k,gpyn)
			grad_p_z=cells(2,j,k,gpzn)

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

!	Set values of this face equal to those NcY-1 cell in y-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,NcY-1,k,mn)
			qx=+1.d0*cells(i,NcY-1,k,qxn)
			qy=+1.d0*cells(i,NcY-1,k,qyn)
			qz=+1.d0*cells(i,NcY-1,k,qzn)
			E = cells(i,NcY-1,k,En)
			rho = cells(i,NcY-1,k,rhon)
			ux = +1.d0*cells(i,NcY-1,k,uxn)
			uy = +1.d0*cells(i,NcY-1,k,uyn)
			uz = +1.d0*cells(i,NcY-1,k,uzn)
			p = cells(i,NcY-1,k,pn)
			V = cells(i,NcY-1,k,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients equal to those in NcY-1 cell in y-direction for periodic boundary:

			grad_d_x= cells(i,NcY-1,k,gdxn)
			grad_d_y=cells(i,NcY-1,k,gdyn)
			grad_d_z=cells(i,NcY-1,k,gdzn)
			grad_ux_x=cells(i,NcY-1,k,guxxn)
			grad_uy_x=cells(i,NcY-1,k,guyxn)
			grad_uz_x=cells(i,NcY-1,k,guzxn)
			grad_ux_y=cells(i,NcY-1,k,guxyn)
			grad_uy_y=cells(i,NcY-1,k,guyyn)
			grad_uz_y=cells(i,NcY-1,k,guzyn)
			grad_ux_z=cells(i,NcY-1,k,guxzn)
			grad_uy_z=cells(i,NcY-1,k,guyzn)
			grad_uz_z=cells(i,NcY-1,k,guzzn)
			grad_p_x=cells(i,NcY-1,k,gpxn)
			grad_p_y=cells(i,NcY-1,k,gpyn)
			grad_p_z=cells(i,NcY-1,k,gpzn)

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

!	Set values of this face equal to those +2 cell in y-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,2,k,mn)
			qx=+1.d0*cells(i,2,k,qxn)
			qy=+1.d0*cells(i,2,k,qyn)
			qz=+1.d0*cells(i,2,k,qzn)
			E = cells(i,2,k,En)
			rho = cells(i,2,k,rhon)
			ux = +1.d0*cells(i,2,k,uxn)
			uy = +1.d0*cells(i,2,k,uyn)
			uz = +1.d0*cells(i,2,k,uzn)
			p = cells(i,2,k,pn)
			V = cells(i,2,k,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients equal to those in +2 cell in y-direction for periodic boundary:

			grad_d_x= cells(i,2,k,gdxn)
			grad_d_y=cells(i,2,k,gdyn)
			grad_d_z=cells(i,2,k,gdzn)
			grad_ux_x=cells(i,2,k,guxxn)
			grad_uy_x=cells(i,2,k,guyxn)
			grad_uz_x=cells(i,2,k,guzxn)
			grad_ux_y=cells(i,2,k,guxyn)
			grad_uy_y=cells(i,2,k,guyyn)
			grad_uz_y=cells(i,2,k,guzyn)
			grad_ux_z=cells(i,2,k,guxzn)
			grad_uy_z=cells(i,2,k,guyzn)
			grad_uz_z=cells(i,2,k,guzzn)
			grad_p_x=cells(i,2,k,gpxn)
			grad_p_y=cells(i,2,k,gpyn)
			grad_p_z=cells(i,2,k,gpzn)

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

!	Set values of this face equal to those NcZ-1 cell in z-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,j,NcZ-1,mn)
			qx=+1.d0*cells(i,j,NcZ-1,qxn)
			qy=+1.d0*cells(i,j,NcZ-1,qyn)
			qz=+1.d0*cells(i,j,NcZ-1,qzn)
			E = cells(i,j,NcZ-1,En)
			rho = cells(i,j,NcZ-1,rhon)
			ux = +1.d0*cells(i,j,NcZ-1,uxn)
			uy = +1.d0*cells(i,j,NcZ-1,uyn)
			uz = +1.d0*cells(i,j,NcZ-1,uzn)
			p = cells(i,j,NcZ-1,pn)
			V = cells(i,j,NcZ-1,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients equal to those in NcZ-1 cell in z-direction for periodic boundary:

			grad_d_x= cells(i,j,NcZ-1,gdxn)
			grad_d_y=cells(i,j,NcZ-1,gdyn)
			grad_d_z=cells(i,j,NcZ-1,gdzn)
			grad_ux_x=cells(i,j,NcZ-1,guxxn)
			grad_uy_x=cells(i,j,NcZ-1,guyxn)
			grad_uz_x=cells(i,j,NcZ-1,guzxn)
			grad_ux_y=cells(i,j,NcZ-1,guxyn)
			grad_uy_y=cells(i,j,NcZ-1,guyyn)
			grad_uz_y=cells(i,j,NcZ-1,guzyn)
			grad_ux_z=cells(i,j,NcZ-1,guxzn)
			grad_uy_z=cells(i,j,NcZ-1,guyzn)
			grad_uz_z=cells(i,j,NcZ-1,guzzn)
			grad_p_x=cells(i,j,NcZ-1,gpxn)
			grad_p_y=cells(i,j,NcZ-1,gpyn)
			grad_p_z=cells(i,j,NcZ-1,gpzn)


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

!	Set values of this face equal to those +2 cell in z-direction. Equal and opposite for x-vel. and x-mom.: 
			m=cells(i,j,2,mn)
			qx=+1.d0*cells(i,j,2,qxn)
			qy=+1.d0*cells(i,j,2,qyn)
			qz=+1.d0*cells(i,j,2,qzn)
			E = cells(i,j,2,En)
			rho = cells(i,j,2,rhon)
			ux = +1.d0*cells(i,j,2,uxn)
			uy = +1.d0*cells(i,j,2,uyn)
			uz = +1.d0*cells(i,j,2,uzn)
			p = cells(i,j,2,pn)
			V = cells(i,j,2,Vn)

!	Set position of faces equal to those already assigned in Cell Setup:
			Midx=cells(i,j,k,xn)
			Midy=cells(i,j,k,yn)
			Midz=cells(i,j,k,zn)

!	Set gradients equal to those in +2 cell in z-direction for periodic boundary:
			grad_d_x= cells(i,j,2,gdxn)
			grad_d_y=cells(i,j,2,gdyn)
			grad_d_z=cells(i,j,2,gdzn)
			grad_ux_x=cells(i,j,2,guxxn)
			grad_uy_x=cells(i,j,2,guyxn)
			grad_uz_x=cells(i,j,2,guzxn)
			grad_ux_y=cells(i,j,2,guxyn)
			grad_uy_y=cells(i,j,2,guyyn)
			grad_uz_y=cells(i,j,2,guzyn)
			grad_ux_z=cells(i,j,2,guxzn)
			grad_uy_z=cells(i,j,2,guyzn)
			grad_uz_z=cells(i,j,2,guzzn)
			grad_p_x=cells(i,j,2,gpxn)
			grad_p_y=cells(i,j,2,gpyn)
			grad_p_z=cells(i,j,2,gpzn)

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
