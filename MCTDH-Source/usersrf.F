C#######################################################################
C    UDEFSRF    User-DEFinition-of-SuRFace
C
C Multi-dimensional surfaces defined by the user.
C Do NOT put 1D functions here! Implement them on user1d.F .
C Compare with funcsrf.F
C
C ifile=203
C
C#######################################################################

      subroutine udefsrf(hopilab,label,ifile,pesopts,zpesopts,hopipar,
     +           hoprpar,maxhpar,pesmin,pesmax,lcutpes,ndof,ierr)

      implicit none

#include "global.inc"
#include "channels.inc"

      integer maxhpar,ndof,ilbl,jlbl,ll
      integer zpesopts(2*maxhpar,*)
      integer hopilab,ifile,ierr,hopipar(*)
      real*8  hoprpar(*),pesmin,pesmax
      character*(*) label,pesopts(*)
      logical lcutpes

      hopilab=0
      ifile=0
      ierr=0

C-----------------------------------------------------------------------
C Determine label length
C-----------------------------------------------------------------------
         ilbl=index(label,'{')-1
         jlbl=index(label,' ')-1
         if(ilbl.le.0) ilbl=99
         if(jlbl.le.0) jlbl=99
         ll = min(ilbl,jlbl)

C-----------------------------------------------------------------------
C FOR EXAMPLES SEE defsrf AND surfinfo ON opfuncs/funcsrf.F
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C MYSURF surface.
C-----------------------------------------------------------------------
      if (label(1:ll) .eq. 'mysrf') then
         hopilab=1

C-----------------------------------------------------------------------
C EXAMPLE: YOURSRF surface.
C NOTE: You also have to modify the subroutines uvpoint and usersurfinfo
C       which are on this file (opfuncs/usersrf.F) below.
C-----------------------------------------------------------------------
C     elseif (label(1:ll) .eq. 'yoursrf') then
C        hopilab=2
       else if (label(1:ll) .eq. 'vhsticking') then
        hopilab=2
       else if (label(1:ll) .eq. 'vhcollinear') then
        hopilab=3
       else if (label(1:ll) .eq. 'vhcartesian') then
        hopilab=4

C-----------------------------------------------------------------------
C end of if loop
C-----------------------------------------------------------------------
      endif

C-----------------------------------------------------------------------
C if label is assigned, assign ifile number.
C-----------------------------------------------------------------------
      if (hopilab .ne. 0) ifile=203

      return
      end


C##########################################################################
C    UVPOINT   User-V-POINT
C
C This routine calculates the value of a potential at a grid point.
C The point on the MCTDH grid is contained in gpoint. The modelabel
C flags allow a check so that gpoint(f) is passed as the correct degree
C of freedom to the potential calculation routine.
C
C To add new surfaces, add the new ifunc number, and then any code
C needed to manipulate the grid point and call the potential surface
C routine. Note that mysrf and DEF_mysrf already exist on source/surfaces,
C the Fortran text of mysrf, however, has to be adapted to the needs.
C##########################################################################

      subroutine uvpoint (v,gpoint,ndof1,ifunc,hopipar,hoprpar)

      implicit none

#include "global.inc"

      integer ifunc,ndof1,hopipar(*)
      real*8  gpoint(ndof1),hoprpar(*)
      real*8  v, r1, r2, theta, phi
      real*8  rho, zc, zh, xh, yh

C Decide which surface to call, accroding to ifunc

      if (ifunc .eq. 1) then
         r1   =gpoint(1)
         r2   =gpoint(2)     ! You may have to edit these lines.
         theta=gpoint(3)     ! You may have to edit these lines.
         phi  =gpoint(4)     ! You may have to edit these lines.
         call mysrf(r1,r2,theta,phi,v)
C ----   The routine mysrf is on source/surfaces/mysrf.f
C ----   One has to edit mysrf.f as well as to adapt usersurfinfo below.

C      elseif(ifunc .eq. 2) then
C ----   add the code for calling another surface here.
        else if(ifunc .eq. 2) then
           zh  = gpoint(1)
           zc  = gpoint(2)
           rho = gpoint(3)
           call vhsticking( rho, 0.0d0, zh, zc, v )
        else if(ifunc .eq. 3) then
           call vhsticking( 0.d0, 0.d0, gpoint(1),gpoint(2) ,v)
        else if(ifunc .eq. 4) then
           zh  = gpoint(1)
           zc  = gpoint(2)
           xh = gpoint(3)
           yh = gpoint(4)
           call vhsticking( xh, yh, zh, zc, v )

      else
         routine='uvpoint (opfuncs/usersrf.F)'
         write(message,'(a,i5)') 'Surface not found :',ifunc
         call errormsg
      endif

      return
      end


C#######################################################################
C  User Surface Info. This routine is called from surfinfo.
C  Description of the user implemented surface and their coordinates.
C#######################################################################

      subroutine usersurfinfo(hopilab,hopipar,hoprpar,mlabel,nn,ilog)

      implicit none

#include "global.inc"

      integer   hopilab,ilog,nn,i,jj
      integer   hopipar(*)
      real*8    hoprpar(*)
      character*16 mlabel(*)

      if (hopilab .eq. 1) then
         write(ilog,'(a)')
     +   'MYSRF (change this text accordingly).'
         write(ilog,'(a)') 'Mysrf Coordinates '
         write(ilog,'(2a)') 'r1    : ',mlabel(1)
         write(ilog,'(2a)') 'r2    : ',mlabel(2)
         write(ilog,'(2a)') 'theta : ',mlabel(3)
         write(ilog,'(2a)') 'phi   : ',mlabel(4)
         jj=5   ! This is dimension+1

C --- The following is just an example ! See also surfinfo on funcsrf.F
C     elseif (hopilab .eq. 2) then
C        write(ilog,'(a)') '(Tell what surface it is).'
C        write(ilog,'(2a)') '1. DOF : ',mlabel(1)
C        write(ilog,'(2a)') '2. DOF : ',mlabel(2)
C        jj=3   ! This is dimension+1
       else if (hopilab .eq. 2) then
        write(ilog,'(a)') 
     + 'VHSTICKING (Jackson PES for Gr-H sticking, 3D).'
        write(ilog,'(2a)') '1. zH  : ',mlabel(1)
        write(ilog,'(2a)') '2. zC  : ',mlabel(2)
        write(ilog,'(2a)') '3. rho : ',mlabel(3)
        jj=4   ! This is dimension+1

       else if (hopilab .eq. 3) then
        write(ilog,'(a)') 
     + 'VHCOLLINEAR (Jackson PES for Gr-H sticking, 2D).'
        write(ilog,'(2a)') '1. zH  : ',mlabel(1)
        write(ilog,'(2a)') '2. zC  : ',mlabel(2)
        jj=3   ! This is dimension+1

       else if (hopilab .eq. 4) then
        write(ilog,'(a)')
     + 'VHCARTESIAN (Jackson PES for Gr-H sticking, 3D).'
        write(ilog,'(2a)') '1. zH  : ',mlabel(1)
        write(ilog,'(2a)') '2. zC  : ',mlabel(2)
        write(ilog,'(2a)') '3. xH  : ',mlabel(3)
        write(ilog,'(2a)') '4. yH  : ',mlabel(4)
        jj=5   ! This is dimension+1

      else
         routine='usersurfinfo (opfuncs/usersrf.F)'
         write(message,'(a,i5)') 'Surface not found :', hopilab
         call errormsg
      endif


      do i =jj,nn
         write(ilog,'(2a)') 'Dummy Variable     : ',mlabel(i)
      enddo

      return
      end

C#######################################################################

