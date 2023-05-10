
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the state
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag        <=  integer tag array
! ::: tag_lo,hi   => index extent of tag array
! ::: state       => state array
! ::: state_lo,hi => index extent of state array
! ::: set         => integer value to tag cell for refinement
! ::: clear       => integer value to untag cell
! ::: lo,hi       => work region we are allowed to change
! ::: dx          => cell size
! ::: problo      => phys loc of lower left corner of prob domain
! ::: time        => problem evolution time
! ::: level       => refinement level of this array
! ::: -----------------------------------------------------------

! // will edit this file to tag high error cells based on density, momentum, and energy

subroutine state_error(tag,tag_lo,tag_hi, &
                       state,state_lo,state_hi, &
                       set,clear,&
                       lo,hi,&
                       dx,problo,time,level) bind(C, name="state_error")

  use tagging_params_module, only : phierr, phigrad, max_phierr_lev, max_phigrad_lev
  implicit none
  
  integer          :: lo(3),hi(3)
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  double precision :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3), &
                            1:4)
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  double precision :: problo(3),dx(3),time
  integer          :: level,set,clear

  double precision :: ax, ay, az, left, right, center
  integer          :: i, j, k, h, dim

  if (state_lo(2) .eq. state_hi(2)) then
     dim = 1
  else if (state_lo(3) .eq. state_hi(3)) then
     dim = 2
  else
     dim = 3
  end if

  ! Tag on regions of high phi gradient
  if (level .lt. max_phigrad_lev) then
     do          k = lo(3), hi(3)
        do       j = lo(2), hi(2)
           do    i = lo(1), hi(1)
                  h = 5   ! Tag on refinement criterion based on density gradient, hence, h=1  - 2023W2

                  ! The gradients are calculated via first order differences at the boundary, while for all other points
                  ! in the domain, the left and right gradient is calculated, and the maximum is taken for a given cell.
                  ! Depending on the dimensionality of the simulation, the irrelevant quantities (dy, dz if 1-D) are 
                  ! manually set to zero. -2023W2

                  if (dim .eq. 1) then    

                     left   = state(i-1,j,k,h)/state(i-1,j,k,1)
                     center = state(i,j,k,h)/state(i,j,k,1)
                     right  = state(i+1,j,k,h)/state(i+1,j,k,1)

                     if (i .eq. lo(1)) then
                        ax = abs(center-right)
                     else if (i .eq. hi(1)) then
                        ax = abs(left-center)
                     else
                        ax = abs(left-center)
                        ax = max(ax, abs(center-right))
                     end if

                     ay = 0.d0
                     az = 0.d0

                  end if

                  if (dim .eq. 2) then  
                     if (i .eq. lo(1)) then
                        ax = abs(state(i,j,k,h)-state(i+1,j,k,h))
                     else if (i .eq. hi(1)) then
                        ax = abs(state(i-1,j,k,h)-state(i,j,k,h))
                     else
                        ax = abs(state(i-1,j,k,h)-state(i,j,k,h))
                        ax = max(ax, abs(state(i,j,k,h)-state(i+1,j,k,h)))
                     end if

                     if (j .eq. lo(2)) then
                        ay = abs(state(i,j,k,h)-state(i,j+1,k,h))
                     else if (j .eq. hi(2)) then
                        ay = abs(state(i,j-1,k,h)-state(i,j,k,h))
                     else 
                        ay = abs(state(i,j-1,k,h)-state(i,j,k,h))
                        ay = max(ay, abs(state(i,j,k,h)-state(i,j+1,k,h)))
                     end if

                     az = 0.d0

                  end if

                  if (dim .eq. 3) then

                     if (i .eq. lo(1)) then
                        ax = abs(state(i,j,k,h)-state(i+1,j,k,h))
                     else if (i .eq. hi(1)) then
                        ax = abs(state(i-1,j,k,h)-state(i,j,k,h))
                     else
                        ax = abs(state(i-1,j,k,h)-state(i,j,k,h))
                        ax = max(ax, abs(state(i,j,k,h)-state(i+1,j,k,h)))
                     end if
                     
                     if (j .eq. lo(2)) then
                        ay = abs(state(i,j,k,h)-state(i,j+1,k,h))
                     else if (j .eq. hi(2)) then
                        ay = abs(state(i,j-1,k,h)-state(i,j,k,h))
                     else 
                        ay = abs(state(i,j-1,k,h)-state(i,j,k,h))
                        ay = max(ay, abs(state(i,j,k,h)-state(i,j+1,k,h)))
                     end if

                     if (k .eq. lo(3)) then
                        az = abs(state(i,j,k,h)-state(i,j,k+1,h))
                     else if (k .eq. hi(3)) then
                        az = abs(state(i,j,k-1,h)-state(i,j,k,h))
                     else 
                        az = abs(state(i,j,k-1,h)-state(i,j,k,h))
                        az = max(az, abs(state(i,j,k,h)-state(i,j,k+1,h)))
                     end if

                  end if

                  if (max(ax,ay,az) .ge. phigrad(level)) then
                     tag(i,j,k) = set
                  end if

                  !write(*,*) "phigrad is ", phigrad(level), " maxgrad = ", max(ax,ay,az), "tag is ", tag(i,j,k)

            end do
         end do
      end do
   end if
  
end subroutine state_error

