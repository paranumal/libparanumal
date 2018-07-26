      program interfaceTest
      implicit none

      integer setuph, n, ierr
      character basis*32

      call holmessetupaidecreate(setuph,ierr)

      call holmessetupaidesetarg('BASIS'//char(0),
     $  'NODAL'//char(0),setuph,ierr)
      call holmessetupaidegetarg(basis,n,'BASIS'//char(0),setuph,ierr)

      write(6,*) 'Basis is:',basis(1:n)

      call holmessetupaidedestroy(setuph,ierr)

      end
