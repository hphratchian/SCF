      module scf_mod
!
!     This module supports the program SCF.
!
!     -H. P. Hratchian, 2020.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
!
!
!     Module Procedures
!
      CONTAINS
!
!
      subroutine formFock(nBasis,density,ERIs,coulomb)
!
!     This subroutine forms a Coulomb matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::coulomb
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempCoulomb
!
      call ERIs%print(IOut,' In formCoulomb: ERIs=',blankAtTop=.True.)
!
!     Do the work...
!
      tempCoulomb = float(0)
      do iSigma    = 1,nBasis
        do iLambda = 1,nBasis
          do iNu   = 1,nBasis
            do iMu = 1,nBasis
              write(IOut,'(/,1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3,2x,I3,2x,I3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])),  &
                float(density%getVal([iLambda,iSigma])),  &
                MQC_Variable_getArrayPosition(ERIs,[iMu,iNu,iLambda,iSigma])
              tempCoulomb(iMu,iNu) = tempCoulomb(iMu,iNu) +  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      coulomb = tempCoulomb
!
      return
      end subroutine formFock
!
!
      subroutine formCoulomb(nBasis,density,ERIs,coulomb,initialize)
!
!     This subroutine forms a Coulomb matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::coulomb
      logical,optional::initialize
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempCoulomb
      logical::init
!
      call ERIs%print(IOut,' In formCoulomb: ERIs=',blankAtTop=.True.)
!
!     Do initial set-up work.
!
      init = .true.
      if(Present(initialize)) init = initialize
!
!     Work through the integral loops to build Coulomb contributions.
!
      if(init) tempCoulomb = float(0)
      do iSigma    = 1,nBasis
        do iLambda = 1,nBasis
          do iNu   = 1,nBasis
            do iMu = 1,nBasis
              write(IOut,'(/,1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3,2x,I3,2x,I3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])),  &
                float(density%getVal([iLambda,iSigma])),  &
                MQC_Variable_getArrayPosition(ERIs,[iMu,iNu,iLambda,iSigma])
              tempCoulomb(iMu,iNu) = tempCoulomb(iMu,iNu) +  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      coulomb = tempCoulomb
!
      return
      end subroutine formCoulomb
!
!
      subroutine formExchange(nBasis,density,ERIs,exchange,initialize)
!
!     This subroutine forms an Exchange matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::exchange
      logical,optional::initialize
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempExchange
      logical::init
!
!     Do initial set-up work.
!
      init = .true.
      if(Present(initialize)) init = initialize
!
!     Work through the integral loops to build Exchange contributions.
!
      if(init) tempExchange = float(0)
      do iMu = 1,nBasis
        do iNu = 1,nBasis
          do iLambda = 1,nBasis
            do iSigma = 1,nBasis
              write(IOut,'(1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iSigma,iLambda,iNu])),  &
                float(density%getVal([iLambda,iSigma]))
              tempExchange(iMu,iNu) = tempExchange(iMu,iNu) -  &
                float(ERIs%getVal([iMu,iSigma,iLambda,iNu])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      exchange = tempExchange
!
      return
      end subroutine formExchange
!
!
      end module scf_mod
