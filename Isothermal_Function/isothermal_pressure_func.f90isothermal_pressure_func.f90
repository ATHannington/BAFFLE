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
!* @file isothermal_pressure_func.f90
!*
!* @Isothermal, two-level temperature system function: Fortran 90 version.
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!********************************************************************************

!===============================================================================
!
!*** 			ISOTHERMAL PRESSURE FUNCTION			     ***
!
!-------------------------------------------------------------------------------
!Function to calculate Isotherml pressure from density.
!Uses Specific heat ratio (gamma), and Global Isothermal Sound Speed (CsIso) 
!written in constants module.
!
!-------------------------------------------------------------------------------
Real*8 function isotherm_p (rho)

	Use constants
	Implicit None
		Real*8, intent(in) :: rho					!Density
		Real*8 :: THigh,TLow						!Two Phase temperature
		Real*8 :: TRatio						!Scaling Ratio between minimum and maximum TEMPERATURE
		Real*8 :: T_Scale						!Scaling factor of pressure BY TEMPERATURE

	TRatio = RhoMaxSet/RhoMinSet

	THigh = TRatio
	TLow = 1.0d0

	If (TwoPhaseBool .eqv. .TRUE.) then
		If (rho.lt.RhoCrit) then
			T_Scale = THigh
		else if (rho.ge.RhoCrit) then
			T_Scale = TLow
		endif
	else
		T_Scale = 1.d0
	endif



!	Isothermal Sound Speed:
!	CsIso^2 = U*(Gamma-1)

!	Internal Energy:
!	U = (k * T)/(gamma-1)
!	U = CsIso^2 / (Gamma-1)

!	Pressure:
!	P = rho * kb * T / mu
!	P = rho*(gam-1)*U
!	P = rho * Cs^2
	isotherm_p = rho*(CsIso**2)*T_Scale

end function isotherm_p
