!********************************************************************************
!* This file is adapted from part of python_finite_volume_solver
!* (https://github.com/bwvdnbro/python_finite_volume_solver), and is part of 
!* BAFFLE (https://github.com/ATHannington/BAFFLE)
!*
!* BAFFLE:
!* Copyright (C) 2017 Andrew Hannington (ath4@st-andrews.ac.uk)
!*
!* python_finite_volume_solver:
!* Copyright (C) 2017 Kenneth Wood (kw25@st-andrews.ac.uk)
!*                    Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
!*
!* python_finite_volume_solver and BAFFLE are free software: you can redistribute 
!* them and/or modify them under the terms of the GNU Affero General Public License
!* as published by the Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* python_finite_volume_solver and BAFFLE are distributed in the hope that they will 
!* be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!* GNU Affero General Public License for more details.
!*
!* You should have received a copy of the GNU Affero General Public License
!* along with python_finite_volume_solver and BAFFLE. If not, see
!* <http://www.gnu.org/licenses/>.
!********************************************************************************

!********************************************************************************
!* @file riemann.f90
!*
!* @brief Standalone Riemann solver library: Fortran 90 version.
!*
!* This Riemann solver is the original exact Riemann solver as it is presented in
!* Toro, E., Riemann Solvers and Numerical Methods for Fluid Dynamics, 3rd
!* edition (Springer, 2009), chapter 4.
!*
!* It was carefully copied (by hand) from that book to actual Fortran 77 code by
!* Kenneth Wood, and slightly modified to be used as a Riemann solver routine in
!* a 1-3D finite volume code. It was then adapted to Fortran 90, and combined with
!* Bert Vandenbroucke's python riemann solver, by Andrew Hannington.
!*
!* Additional changes (including this documentation) were made by Bert
!* Vandenbroucke to make the Python, C++ and Fortran version of the Riemann
!* solver more homogeneous. Subsequent changes (including edits to this documentation)
!* were made by Andrew Hannington, including adding vacuum capability from Bert
!* Vandenbroucke's own python Riemann solver.
!*
!* @author Andrew Hannington (ath4@st-andrews.ac.uk)
!* @author Kenneth Wood (kw25@st-andrews.ac.uk)
!* @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
!********************************************************************************

!********************************************************************************
!* @brief Solve the Riemann problem with the given left and right state.
!*
!* @pCRam d1 Left state density.
!* @pCRam u1 Left state fluid velocity.
!* @pCRam p1 Left state pressure.
!* @pCRam d2 Right state density.
!* @pCRam u2 Right state fluid velocity.
!* @pCRam p2 Right state pressure.
!* @pCRam gg Adiabatic index of the gas.
!* @pCRam ds VCRiable to store the resulting density in.
!* @pCRam us VCRiable to store the resulting fluid velocity in.
!* @pCRam ps VCRiable to store the resulting pressure in.
!* @pCRam uflag Flag vCRiable used to signal if the solution was sampled from the
!* left state (-1) or the right state (1). A vacuum state (0) is currently not
!* supported by this version of the Riemann solver.
!* @pCRam x_over_t Point (in velocity space) where we want to sample the
!* solution. In a finite volume scheme you usually want to set this value to 0.
!********************************************************************************
module gamma_module
	REAL*8 :: G1, G2, G3, G4, G5, G6, G7, G8
end module gamma_module

subroutine initialise_gammas
	USE CONSTANTS
	use gamma_module
     	implicit none

      G1 = (gam - 1.0d0)/(2.0d0*gam)
      G2 = (gam + 1.0d0)/(2.0d0*gam)
      G3 = 2.0d0*gam/(gam - 1.0d0)
      G4 = 2.0d0/(gam - 1.0d0)
      G5 = 2.0d0/(gam + 1.0d0)
      G6 = (gam - 1.0d0)/(gam + 1.0d0)
      G7 = (gam - 1.0d0)/2.0d0
      G8 = gam - 1.0d0

end subroutine initialise_gammas
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************

subroutine riemann(d1, u1, p1, d2, u2, p2, gg, ds,us, ps, uflag, x_over_t,i,j,k,xyz,time)

	USE CONSTANTS
	use gamma_module
     	implicit none

	Integer,intent(in) :: i,j,k
	Real*8, intent(in) :: time
	Character(len=1),intent(in) :: xyz					!ChCRacter passed in from cells to Riemann that tells this subroutine whether we CRe looking at x, y, or z direction.

!     DeclCRation of vCRiables:
      integer uflag
      Real*8 ::  AGAMMA, DL, UL, PL, CL, DR, UR, PR, CR, &
          DIAPH, DOMLEN, DS, DX, PM, MPA, PS, S, &
          TIMEOUT, UM, US, XPOS

      real*8 gg, d1, u1, p1, d2, u2, p2
      real*8 x_over_t

      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR

!     Compute gamma related constants

      AGAMMA = gg
      DL = d1
      UL = u1
      PL = p1
      DR = d2
      UR = u2
      PR = p2

      MPA=1.0d0  ! Not sure about this????

!     Compute sound speeds

	If (IsoBool .eqv. .FALSE.) then
		CL = DSQRT(AGAMMA*PL/DL)
		CR = DSQRT(AGAMMA*PR/DR)
	else if (IsoBool .eqv. .TRUE.) then
		CL = DSQRT(PL/DL)
		CR = DSQRT(PR/DR)
	ENDIF

!     The pressure positivity condition is tested for

      IF(G4*(CL+CR).LE.(UR-UL))THEN
!		Print*,
!		Print*,"***!*!*!*!*!***"
!		Print*, "Vacuum generated by data [@RIEMANN SOLVER]."
!		Print*,"Whilst solving in ",xyz," direction."
!		Print*,"i=",i,"j=",j,"k=",k
!		Print*,"Time=",time
!		Print*,
!		Print*,"Attempting to Solve!"
		call solve_vacuum(dL, uL, PL, CL, dR, uR, PR, CR, x_over_t, ds, us, Ps, uflag)
	!	     solve_vacuum(rhoL, uL, PL, CL, rhoR, uR, PR, CR, x_over_t, rhosol, usol, Psol, flag)

	else

!     Exact solution for pressure and velocity in star region is found
		CALL STCRPU(PM, UM, MPA)



	      if(pm.ne.pm) then
		Print*,"NAN DETECTED! [@RIEMANN SOLVER]. PROGRAM WILL ABORT!"
		Print*,"Whilst solving in ",xyz," direction."
		Print*,"i=",i,"j=",j,"k=",k
		Print*,"Time=",time
		print*,"PM=",pm,"PL=",pl,"DL=",dl,"uL=",ul,"PR=",pr,"DR=",dr,"uR=",ur
		stop
	      endif
	      if(UM .gt.x_over_t) then
		uflag=-1 ! "left state"
	      else
		uflag=1 ! "right state
	      endif


	      CALL SAMPLE(PM, UM, x_over_t, DS, US, PS)

	endif

      RETURN
      END

!********************************************************************************
!* @brief Get the pressure and velocity in the middle state.
!*
!* @pCRam p VCRiable to store the resulting pressure in.
!* @pCRam u VCRiable to store the resulting fluid velocity in.
!* @pCRam mpa Mach number? Not used.
!********************************************************************************
      SUBROUTINE STCRPU(P, U, MPA)

	USE CONSTANTS
	use gamma_module
     	implicit none

!     Purpose: to compute the solution for pressure and velocity in stCR region

      INTEGER I, NRITER

      REAL*8 DL, UL, PL, CL, DR, UR, PR, CR, &
         CHANGE, FL, FLD, FR, FRD, P, POLD, PSTCRT, &
         TOLPRE, U, UDIFF, MPA

      COMMON/STATES/ DL, UL, PL, CL, DR, UR, PR, CR
      DATA TOLPRE, NRITER/1.0d-09, 20/

!     Guessed value PSTCRT is computed

      CALL GUESSP(PSTCRT)

      POLD = PSTCRT
      UDIFF = UR - UL
!      WRITE(6,*)'Iteration number Change '

      DO 10 I = 1, NRITER

         CALL PREFUN(FL, FLD, POLD, DL, PL, CL)
         CALL PREFUN(FR, FRD, POLD, DR, PR, CR)
         P = POLD - (FL + FR + UDIFF)/(FLD + FRD)

         CHANGE = 2.0d0*DABS((P - POLD)/(P + POLD))
!         WRITE(6, 30)I, CHANGE
         IF(CHANGE.LE.TOLPRE)GOTO 20
         IF(P.LT.0.0)P = TOLPRE
         POLD = P

 10   CONTINUE

!      WRITE(6,*)'Divergence in Newton-Raphson iteration'

 20   CONTINUE

!     Compute velocity in stCR region

      U = 0.5d0*(UL + UR + FR - FL)

!      WRITE(6,*)'Pressure	Velocity'
!      WRITE(6,40)P/MPA, U

 30   FORMAT(5X, I5,15X, F12.7)
 40   FORMAT(2(F14.6, 5X))

      RETURN
      END

!********************************************************************************
!* @brief Get an initial guess for the pressure in the middle state, as a
!* stCRting point for the Newton-Raphson iteration.
!*
!* @pCRam pm VCRiable to store the resulting pressure guess in.
!********************************************************************************
      SUBROUTINE GUESSP(PM)

!     Purpose: to provide a guess value for pressure PM in the StCR Region.
!     The choice is made according to adaptive Riemann solver using the PVRS,
!     TRRS, and TSRS approximate Riemann solvers

	USE CONSTANTS
	use gamma_module
     	implicit none

      REAL*8 DL, UL, PL, CL, DR, UR, PR, CR, &
          AGAMMA, &
          CUP, GEL, GER, PM, PMAX, PMIN, PPV, PQ, &
          PTL, PTR, QMAX, QUSER, UM


      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR

      QUSER = 2.0d0

!     Compute guess pressure for PVRS Riemann solver

      CUP = 0.25d0*(DL + DR)*(CL + CR)
      PPV = 0.5d0*(PL + PR) + 0.5d0*(UL - UR)*CUP
      PPV = AMAX1(0.0d0, PPV)
      PMIN = AMIN1(PL, PR)
      PMAX = AMAX1(PL, PR)
      QMAX = PMAX/PMIN

      IF(QMAX.LE.QUSER.AND. &
       (PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN

!        Select PVRS Riemann solver

         PM = PPV

      ELSE
         IF(PPV.LT.PMIN)THEN

!           Select Two-RCRefaction Riemann solver

            PQ = (PL/PR)**G1
            UM = (PQ*UL/CL + UR/CR + &
                 G4*(PQ - 1.0d0))/(PQ/CL + 1.0d0/CR)
            PTL = 1.0d0 + G7*(UL - UM)/CL
            PTR = 1.0d0 + G7*(UM - UR)/CR
            PM = 0.5d0*(PL*PTL**G3 + PR*PTR**G3)
         ELSE

!           Select Two=Shock Riemann solver with PVRS as estimate

            GEL = DSQRT((G5/DL)/(G6*PL + PPV))
            GER = DSQRT((G5/DR)/(G6*PR + PPV))
            PM = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
         ENDIF
      ENDIF

      RETURN
      END

!********************************************************************************
!* @brief Evaluate the pressure functions for the given pressure guess and the
!* given left or right state vCRiables.
!*
!* @pCRam f VCRiable to store the pressure function value in.
!* @pCRam fd VCRiable to store the value of the derivative of the presure
!* function in.
!* @pCRam p Current pressure guess.
!* @pCRam dk Left or right state density.
!* @pCRam pk Left or right state pressure.
!* @pCRam ck Left or right state sound speed.
!********************************************************************************
      SUBROUTINE PREFUN(F, FD, P, DK, PK, CK)

!     Purpose: to evaluate the pressure functions FL and FR in the exact Riemann
!     solverPM

	USE CONSTANTS
	use gamma_module
     	implicit none

      REAL*8 AK, BK, CK, DK, F, FD, P, PK, PRAT, QRT, &
          AGAMMA

      COMMON /GAMMAS/ AGAMMA

      IF(P.LE.PK) THEN

!        RCRefaction wave

         PRAT = P/PK
         F = G4*CK*(PRAT**G1 - 1.0d0)
         FD = (1.0d0/(DK*CK))*PRAT**(-G2)

      ELSE

!        Shock wave
         AK = G5/DK
         BK = G6*PK
         QRT = DSQRT(AK/(BK + P))
         F = (P - PK)*QRT
         FD = (1.0d0 - 0.5d0*(P - PK)/(BK + P))*QRT

      ENDIF

      RETURN
      END

!********************************************************************************
!* @brief Sample the Riemann solution with the given middle state pressure and
!* fluid velocity at the given point in velocity space.
!*
!* @pCRam pm Middle state pressure.
!* @pCRam um Middle state fluid velocity.
!* @pCRam s Sampling point in velocity space.
!* @pCRam d VCRiable to store the resulting density in.
!* @pCRam u VCRiable to store the resulting fluid velocity in.
!* @pCRam p VCRiable to store the resulting pressure in.
!********************************************************************************
      SUBROUTINE SAMPLE(PM, UM, S, D, U, P)

!     Purpose to sample the solution throughout the wave pattern. Pressure PM
!     and velocity UM in the StCR Region CRe known. Sampling is performed in
!     terms of the 'speed' S = X/T. Sampled values CRe D, U, P

!     Input vCRiables: PM, UM, S, /GAMMAS. /STATES/
!     Output vCRiables: D, U, P
	USE CONSTANTS
	use gamma_module
     	implicit none

      REAL*8 DL, UL, PL, CL, DR, UR, PR, CR, &
          AGAMMA, &
          C, CML, CMR, D, P, PM, PML, PMR, S, &
          SHL, SHR, SL, SR, STL, STR, U, UM

      COMMON /GAMMAS/ AGAMMA
      COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR

      IF(S.LE.UM)THEN

!     Sampling point lies to the left of the contact discontinuity

         IF(PM.LE.PL)THEN

!          Left rCRefaction

            SHL = UL - CL

            IF(S.LE.SHL)THEN

!              Sampled point is left data state

               D = DL
               U = UL
               P = PL
            ELSE
               CML = CL*(PM/PL)**G1
               STL = UM - CML

               IF(S.GT.STL)THEN

!                 Sampled point is StCR Left state

                  D = DL*(PM/PL)**(1.0d0/AGAMMA)
                  U = UM
                  P = PM
               ELSE

!                 Sampled point is inside left fan

                  U = G5*(CL + G7*UL + S)
                  C = G5*(CL + G7*(UL - S))
                  D = DL*(C/CL)**G4
                  P = PL*(C/CL)**G3
               ENDIF
            ENDIF
!         ENDIF
         ELSE

!            Left shock

            PML = PM/PL
            SL = UL - CL*DSQRT(G2*PML + G1)

            IF(S.LE.SL)THEN

!              Sampled point is left data state

               D = DL
               U = UL
               P = PL

             ELSE

!              Sampled point is StCR Left state

               D = DL*(PML+G6)/(PML*G6 + 1.0d0)
               U = UM
               P = PM
             ENDIF
           ENDIF
          ELSE

!    Sampling point lies to the right of the contact discontinuity

         IF(PM.GT.PR)THEN

!           Right shock

            PMR = PM/PR
            SR = UR + CR*DSQRT(G2*PMR + G1)

            IF(S.GE.SR)THEN

!              Sampled point is right data state

               D = DR
               U = UR
               P = PR
            ELSE

!              Sampled point is StCR Right state

               D = DR*(PMR + G6)/(PMR*G6 + 1.0d0)
               U = UM
               P = PM
            ENDIF
         ELSE


!           Right rCRefaction

            SHR = UR + CR

            IF(S.GE.SHR)THEN

!              Sampled point is right data state

               D = DR
               U = UR
               P = PR
            ELSE
               CMR = CR*(PM/PR)**G1
               STR = UM + CMR

               IF(S.LE.STR)THEN

!                 Sampled point is StCR Right state

                  D = DR*(PM/PR)**(1.0d0/AGAMMA)
                  U = UM
                  P = PM
               ELSE

!                 Sampled point is inside left fan

                  U=G5*(-CR + G7*UR + S)
                  C = G5*(CR - G7*(UR-S))
                  D = DR*(C/CR)**G4
                  P = PR*(C/CR)**G3
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
END

!  ##############################################################################
  ! @brief Vacuum Riemann solver.
!
  ! This solver is called when one or both states have a zero density, or when
  ! the vacuum generation condition is satisfied (meaning vacuum is generated
  ! in the middle state, although strictly speaking there is no "middle"
  ! state if vacuum is involved).
!
  ! @pCRam rhoL Density of the left state.
  ! @pCRam uL Velocity of the left state.
  ! @pCRam PL Pressure of the left state.
  ! @pCRam CL Soundspeed of the left state.
  ! @pCRam rhoR Density of the right state.
  ! @pCRam uR Velocity of the right state.
  ! @pCRam PR Pressure of the right state.
  ! @pCRam CR Soundspeed of the right state.
  ! @pCRam x_over_t Point in velocity space where we want to sample the solution.
  ! @return Density, velocity and pressure solution, and a flag indicating
  ! wether the left state (-1), the right state (1), or a vacuum state (0) was
  ! sampled.
!  ##############################################################################
subroutine solve_vacuum(rhoL, uL, PL, CL, rhoR, uR, PR, CR, x_over_t, rhosol, usol, Psol, flag)

	USE CONSTANTS
	use gamma_module
     	implicit none

		Real*8,intent(in) :: rhoL,uL,PL,CL
		Real*8,intent(in) :: rhoR,uR,PR,CR
		Real*8,intent(in) :: x_over_t 
	   	Real*8,intent(out) :: rhosol, usol, Psol
		Integer,intent(out) :: flag
		Real*8 :: SR,SL,Base

! if both states are vacuum, the solution is also vacuum
	If ((rhoL == 0.d0) .and. (rhoR == 0.d0)) then
     		RhoSol= 0.d0
		uSol= 0.d0
		PSol= 0.d0
		flag=0

	else if (rhoR == 0.d0) then 
		! vacuum right state
		call sample_right_vacuum(rhoL, uL, PL, CL, x_over_t, rhosol, usol, Psol, flag)
	else if (rhoL == 0.d0) then
		! vacuum left state
		call sample_left_vacuum(rhoR, uR, PR, CR, x_over_t, rhosol, usol, Psol, flag)
	else
		! vacuum "middle" state
		call sample_vacuum_generation(rhoL, uL, PL, CL, rhoR, uR, PR, CR, x_over_t, rhosol, usol, Psol, flag)
	endif

end subroutine solve_vacuum

!  ##############################################################################
!  ! @brief Sample the vacuum Riemann problem if the right state is a vacuum.
  ! @pCRam rhoL Density of the left state.
  ! @pCRam uL Velocity of the left state.
  ! @pCRam PL Pressure of the left state.
  ! @pCRam CL Soundspeed of the left state.
  ! @pCRam x_over_t Point in velocity space where we want to sample the solution.
  ! @return Density, velocity and pressure solution, and a flag indicating
  ! wether the left state (-1), the right state (1), or a vacuum state (0) was
  ! sampled.
 ! ##############################################################################
subroutine sample_right_vacuum(rhoL, uL, PL, CL, x_over_t, rhosol, usol, Psol, flag)
	USE CONSTANTS
	use gamma_module
     	implicit none

		Real*8,intent(in) :: rhoL,uL,PL,CL,x_over_t 
	   	Real*8,intent(out) :: rhosol, usol, Psol
		Integer,intent(out) :: flag
		Real*8 :: SL,Base

	If ((uL - CL).lt.x_over_t) then

	! vacuum regime
	!get the vacuum rCRefaction wave speed
		SL = uL + G4 * CL

      		If (SL .gt. x_over_t) then
	! rCRefaction wave regime
	! vCRiable used twice below
			base = G5 + G6* (uL - x_over_t) / CL

			rhosol = rhoL * base**G4
			usol = G5 * (CL + (G7 * uL) + x_over_t)
			Psol = PL * base**G3
			flag = -1
      		else
			! vacuum
			rhosol = 0.d0
			usol = 0.d0
			Psol = 0.d0
			flag = 0
		endif
	else
		! left state regime
     		rhosol = rhoL
      		usol = uL
     		Psol = PL
      		flag = -1
	endif

end subroutine sample_right_vacuum
!##############################################################################
  ! @brief Sample the vacuum Riemann problem if the left state is a vacuum.
!
  ! @pCRam rhoR Density of the right state.
  ! @pCRam uR Velocity of the right state.
  ! @pCRam PR Pressure of the right state.
  ! @pCRam CR Soundspeed of the right state.
  ! @pCRam x_over_t Point in velocity space where we want to sample the solution.
  ! @return Density, velocity and pressure solution, and a flag indicating
  ! wether the left state (-1), the right state (1), or a vacuum state (0) was
  ! sampled.
!  ##############################################################################
subroutine sample_left_vacuum(rhoR, uR, PR, CR, x_over_t, rhosol, usol, Psol, flag)
	USE CONSTANTS
	use gamma_module
     	implicit none

		Real*8,intent(in) :: rhoR,uR,PR,CR,x_over_t 
	   	Real*8,intent(out) :: rhosol, usol, Psol
		Integer,intent(out) :: flag
		Real*8 :: SR,Base

	If (x_over_t .lt. (uR + CR)) then
	! vacuum regime
	! get the vacuum rarefaction wave speed
		SR = uR - G4 * CR

		if (SR .lt. x_over_t) then
		! rarefaction wave regime
			! vCRiable used twice below
			base = G5 - G6* (uR - x_over_t) / CR
			rhosol = rhoR * base**G4
			usol = G5 * (-1.d0*CR + G4 * uR + x_over_t)
			Psol = PR * base**G3
			flag = 1
		else
		! vacuum
			rhosol = 0.d0
			usol = 0.d0
			Psol = 0.d0
			flag = 0
		endif

	else
	! right state regime
	      rhosol = rhoR
	      usol = uR
	      Psol = PR
	      flag = 1
	endif

end subroutine sample_left_vacuum

!  ##############################################################################
  ! @brief Sample the vacuum Riemann problem in the case vacuum is generated in
  ! between the left and right state.
!
  ! @pCRam rhoL Density of the left state.
  ! @pCRam uL Velocity of the left state.
  ! @pCRam PL Pressure of the left state.
  ! @pCRam CL Soundspeed of the left state.
  ! @pCRam rhoR Density of the right state.
  ! @pCRam uR Velocity of the right state.
  ! @pCRam PR Pressure of the right state.
  ! @pCRam CR Soundspeed of the right state.
  ! @pCRam x_over_t Point in velocity space where we want to sample the solution.
  ! @return Density, velocity and pressure solution, and a flag indicating
  ! wether the left state (-1), the right state (1), or a vacuum state (0) was
  ! sampled.
!  ##############################################################################
subroutine sample_vacuum_generation(rhoL, uL, PL, CL, rhoR, uR, PR, CR, x_over_t, rhosol, usol, Psol, flag)

	USE CONSTANTS
	use gamma_module
     	implicit none

		Real*8,intent(in) :: rhoL,uL,PL,CL
		Real*8,intent(in) :: rhoR,uR,PR,CR
		Real*8,intent(in) :: x_over_t 
	   	Real*8,intent(out) :: rhosol, usol, Psol
		Integer,intent(out) :: flag
		Real*8 :: SR,SL,Base

! get the speeds of the left and right rCRefaction waves
    SR = uR - G4 * CR
    SL = uL + G4 * CL

	If ((SR .gt. x_over_t) .and. (SL .lt. x_over_t)) then
		! vacuum
	      rhosol = 0.d0
	      usol = 0.d0
	      Psol = 0.d0
	      flag = 0

	else if (SL .lt. x_over_t) then
	! right state
		If (x_over_t .lt. (uR + CR)) then
		! right rCRefaction wave regime
		! vCRiable used twice below
			base = G5 - G6* (uR - x_over_t) / CR
			rhosol = rhoR * base**G4
			usol = G5 * (-1.d0*CR + G4 * uR + x_over_t)
			Psol = PR * base**G3
   		else
		! right state regime
			rhosol = rhoR
			usol = uR
			Psol = PR
			flag = 1
		endif
	else
	! left state
		If (x_over_t .gt. (uL - CL)) then
		! left rCRefaction wave regime
		! vCRiable used twice below
			base = G5 + G6* (uL - x_over_t) / CL
			rhosol = rhoL * base**G4
			usol = G5 * (CL + G4 * uL + x_over_t)
			Psol = PL * base**G3
        	else
			! left state regime
			rhosol = rhoL
			usol = uL
			Psol = PL
			flag = -1
		endif
	endif

END subroutine sample_vacuum_generation
