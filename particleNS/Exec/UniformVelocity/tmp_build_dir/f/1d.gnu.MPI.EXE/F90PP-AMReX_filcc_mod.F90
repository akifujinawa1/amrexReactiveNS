











module amrex_filcc_module

  use amrex_fort_module, only : amrex_real, amrex_spacedim
  use amrex_bc_types_module
  use amrex_constants_module

  implicit none

  interface amrex_filcc
     module procedure amrex_filcc_1
     module procedure amrex_filcc_n
  end interface amrex_filcc

  private
  public :: amrex_filcc, amrex_fab_filcc, amrex_hoextraptocc
  public :: amrex_filccn

  public :: filccn

contains

  subroutine amrex_filcc_n(q,qlo,qhi,domlo,domhi,dx,xlo,bclo,bchi)
    integer, intent(in) :: qlo(4), qhi(4)
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bclo(amrex_spacedim,*), bchi(amrex_spacedim,*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qlo(4):qhi(4))
    integer :: i, bc(amrex_spacedim,2)
    do i = qlo(4), qhi(4)
       bc(:,1) = bclo(:,i)
       bc(:,2) = bchi(:,i)
       call amrex_filccn(qlo(1:3), qhi(1:3), q(:,:,:,i), qlo(1:3), qhi(1:3), 1, &
            domlo, domhi, dx, xlo, bc);
    end do
  end subroutine amrex_filcc_n


  subroutine amrex_filcc_1(q,qlo1,qhi1,domlo,domhi,dx,xlo,bc)
    integer, intent(in) :: qlo1,qhi1,domlo(amrex_spacedim),domhi(amrex_spacedim)
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim, 2)
    real(amrex_real), intent(inout) :: q(qlo1:qhi1)
    integer :: q_lo(3), q_hi(3)
    q_lo = [qlo1,0,0]
    q_hi = [qhi1,0,0]
    call amrex_filccn(q_lo, q_hi, q, q_lo, q_hi, 1, domlo, domhi, dx, xlo, bc);
  end subroutine amrex_filcc_1


  subroutine amrex_fab_filcc (q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc) &
       bind(c, name='amrex_fab_filcc')

    implicit none

    integer, intent(in) :: qlo(3), qhi(3), nq
    integer, dimension(amrex_spacedim), intent(in) :: domlo, domhi
    real(amrex_real), intent(in) :: dx(amrex_spacedim), xlo(amrex_spacedim)
    integer, intent(in) :: bc(amrex_spacedim,2,nq)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nq)

    integer :: lo(3), hi(3)

    lo = qlo
    hi = qhi

    call amrex_filccn(lo, hi, q, qlo, qhi, nq, domlo, domhi, dx, xlo, bc)

  end subroutine amrex_fab_filcc

  subroutine filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)
    implicit none
    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)
    call amrex_filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)
  end subroutine filccn

  subroutine amrex_filccn(lo, hi, q, q_lo, q_hi, ncomp, domlo, domhi, dx, xlo, bc)

    implicit none

    integer,          intent(in   ) :: lo(3), hi(3)
    integer,          intent(in   ) :: q_lo(3), q_hi(3)
    integer,          intent(in   ) :: ncomp
    integer,          intent(in   ) :: domlo(amrex_spacedim), domhi(amrex_spacedim)
    real(amrex_real), intent(in   ) :: xlo(amrex_spacedim), dx(amrex_spacedim)
    real(amrex_real), intent(inout) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),ncomp)
    integer,          intent(in   ) :: bc(amrex_spacedim,2,ncomp)

    integer :: is, ie, ilo, ihi, imin, imax
    integer :: i, j, k, n

    is = max(q_lo(1), domlo(1))
    ie = min(q_hi(1), domhi(1))
    ilo = domlo(1)
    ihi = domhi(1)



    do n = 1, ncomp

       if (lo(1) < ilo) then
          imin = lo(1)
          imax = ilo-1

          if (bc(1,1,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,1,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

                      if (i < ilo - 1) then
                         q(i,j,k,n) = q(ilo,j,k,n)
                      else if (i == ilo - 1) then
                         if (ilo+2 <= ie) then
                            q(i,j,k,n) = eighth * (15*q(ilo,j,k,n) - 10*q(ilo+1,j,k,n) + 3*q(ilo+2,j,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(ilo,j,k,n) - q(ilo+1,j,k,n))
                         end if
                      end if

                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. amrex_bc_hoextrapcc) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = (ilo-i)*(q(ilo,j,k,n) - q(ilo+1,j,k,n)) + q(ilo,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,1,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ilo+(ilo-i)-1,j,k,n)
                   end do
                end do
             end do

          end if

       end if

       if (hi(1) > ihi) then
          imin = ihi+1
          imax = hi(1)

          if (bc(1,2,n) .eq. amrex_bc_ext_dir) then

             ! Do nothing.

          else if (bc(1,2,n) .eq. amrex_bc_foextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,2,n) .eq. amrex_bc_hoextrap) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax

                      if (i > ihi + 1) then
                         q(i,j,k,n) = q(ihi,j,k,n)
                      else if (i == ihi + 1) then
                         if (ihi-2 >= is) then
                            q(i,j,k,n) = eighth * (15*q(ihi,j,k,n) - 10*q(ihi-1,j,k,n) + 3*q(ihi-2,j,k,n))
                         else
                            q(i,j,k,n) = half * (3*q(ihi,j,k,n) - q(ihi-1,j,k,n))
                         end if
                      end if

                   end do
                end do
             end do

          else if (bc(1,2,n) .eq. amrex_bc_hoextrapcc) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = (i-ihi)*(q(ihi,j,k,n) - q(ihi-1,j,k,n)) + q(ihi,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,2,n) .eq. amrex_bc_reflect_even) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do

          else if (bc(1,2,n) .eq. amrex_bc_reflect_odd) then

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = imin, imax
                      q(i,j,k,n) = -q(ihi-(i-ihi)+1,j,k,n)
                   end do
                end do
             end do

          end if

       end if










    end do

  end subroutine amrex_filccn


  subroutine amrex_hoextraptocc (q, qlo, qhi, domlo, domhi, dx, xlo) &
       bind(c,name='amrex_hoextraptocc')
    integer, intent(in) :: qlo(3), qhi(3), domlo(*), domhi(*)
    real(amrex_real), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    real(amrex_real), intent(in) :: dx(*), xlo(*)

  end subroutine amrex_hoextraptocc




end module amrex_filcc_module
