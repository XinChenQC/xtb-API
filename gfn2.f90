program test
   call test_gfn2_scc
   call test_gfn2_api
end program test

subroutine test_gfn2_scc
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_basisset
   use tbdef_param
   use tbdef_data
   use tbdef_pcem
   use setparam
   use aoparam
   use xbasis
   use scf_module
   use scc_core

   implicit none
   real(wp),parameter :: thr = 1.0e-7_wp
   real(wp),parameter :: thr2 = 1.0e-5_wp
   integer(8), parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))
   real(wp),parameter :: et = 300.0_wp
   integer, parameter :: maxiter = 50
   integer, parameter :: prlevel = 2
   logical, parameter :: lgrad = .true.
   logical, parameter :: restart = .false.
   real(wp),parameter :: acc = 1.0_wp

   type(tb_molecule)     :: mol
   type(scc_results)     :: res
   type(tb_basisset)     :: basis
   type(tb_wavefunction) :: wfn
   type(scc_parameter)   :: param
   type(tb_pcem)         :: pcem
   real(wp) :: etot,egap
   real(wp), allocatable :: g(:,:)
   real(wp) :: globpar(25)
   logical  :: okpar,okbas
   gfn_method = 2

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   mol%chrg = 0.0_wp
   call mol%set_nuclear_charge
   call mol%update

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )
   call use_parameterset('.param_gfn2.xtb',globpar,okpar)
   call set_gfn2_parameter(param,globpar,mol%n,mol%at)
   call xbasis0(mol%n,mol%at,basis)
   call xbasis_gfn2(mol%n,mol%at,basis,okbas)
   call xbasis_cao2sao(mol%n,mol%at,basis)


   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   call iniqshell(mol%n,mol%at,mol%z,basis%nshell,wfn%q,wfn%qsh,gfn_method)

   g = 0.0_wp

   call scf(output_unit,mol,wfn,basis,param,pcem, &
      &   egap,et,maxiter,prlevel,restart,lgrad,acc,etot,g,res)


   call mol%deallocate
   call wfn%deallocate
   call basis%deallocate

end subroutine test_gfn2_scc


subroutine test_gfn2_api
   use iso_fortran_env, wp => real64, istdout => output_unit

!   use assertion

   use tbdef_options
   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_pcem

   use tb_calculators

   implicit none

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,6,6,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[0.00000000000000_wp, 0.00000000000000_wp,-1.79755622305860_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 0.95338756106749_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 3.22281255790261_wp, &
      &-0.96412815539807_wp,-1.66991895015711_wp,-2.53624948351102_wp, &
      &-0.96412815539807_wp, 1.66991895015711_wp,-2.53624948351102_wp, &
      & 1.92825631079613_wp, 0.00000000000000_wp,-2.53624948351102_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 5.23010455462158_wp], shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true.,&
      &  solvent = "none")

   type(tb_molecule)    :: mol
   type(tb_wavefunction):: wfn
   type(tb_environment) :: env
   type(gfn_parameter)  :: gfn
   type(tb_pcem)        :: pcem

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call env%setup

   call mol%allocate(nat)
   mol%at  = at
   mol%xyz = xyz
   call mol%set_nuclear_charge
   call mol%update

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call gfn2_calculation &
      (istdout,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

   print *,gradient(2,3),'0.000000'
   print *,gradient(3,1),'-0.746499'
   print *,gradient(1,4),'-0.284337'
   print *,gradient(3,7),'0.997505'

end subroutine test_gfn2_api
