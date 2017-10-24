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
