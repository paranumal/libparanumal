      program interfaceTest
      implicit none

      integer setuph, n, ierr
      character basis*32

      call paranumalsetupaidecreate(setuph,ierr)

      call paranumalsetupaidesetarg('BASIS'//char(0),
     $  'NODAL'//char(0),setuph,ierr)
      call paranumalsetupaidegetarg(basis,n,'BASIS'//char(0),setuph,
     $  ierr)

      write(6,*) 'Basis is:',basis(1:n)

      call paranumalsetupaidedestroy(setuph,ierr)

      end
