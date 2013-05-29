
      subroutine vhsticking( xh, yh, zh, zc ,vv )
      real*8   xh, yh, zh, zc, vv
      write(6,'(a,/a)') '###########################################',
     +  'v:H2 needed. lsth is not linked. Run compile with -i option.'
      write(2,'(a,/a)') '###########################################',
     +  'v:H2 needed. lsth is not linked. Run compile with -i option.'
      stop
      end
