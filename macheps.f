      program macheps
      implicitnone
      doubleprecision dble_epsi
      doubleprecision dble_eps

      doubleprecision dble_dble
      doubleprecision dble_one / 1d0 /

      integer*8 i8_dble
      integer*8 i8_one

      equivalence (dble_dble,i8_dble)
      equivalence (dble_one,i8_one)

      real*4 real_eps
      real*4 real_epsi

      real*4 real_real
      real*4 real_one / 1.0 /

      integer*4 i4_real
      integer*4 i4_one

      equivalence (real_real,i4_real)
      equivalence (real_one,i4_one)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      dble_eps = 1d0
      dble_dble = 1.5d0
      do while (dble_dble .ne. dble_one)
        dble_eps = dble_eps * 0.5d0
        dble_dble = 1d0+(0.5d0*dble_eps)
      enddo

      dble_epsi = 1d0
      dble_dble = 1.5d0
      do while (i8_dble .ne. i8_one)
        dble_epsi = dble_epsi * 0.5d0
        dble_dble = 1d0+(0.5d0*dble_epsi)
      enddo

      real_eps = 1.0
      real_real = 1.5
      do while (real_real .ne. real_one)
        real_eps = real_eps * 0.5
        real_real = 1.0+(0.5*real_eps)
      enddo

      real_epsi = 1.0
      real_real = 1.5
      do while (i4_real .ne. i4_one)
        real_epsi = real_epsi * 0.5
        real_real = 1.0+(0.5*real_epsi)
      enddo

      print'(1x,e24.16,a)',dble_eps ,'~Double epsilon (DOUBLE compare)'
     &                    ,dble_epsi,'~Double epsilon (INT*8 compare)'
     &                    ,real_eps ,'~Single epsilon (REAL compare)'
     &                    ,real_epsi,'~Single epsilon (INT*4 compare)'

      dble_eps = 1d0
      do while ((1.0d0+(0.5d0*dble_eps)) .ne. dble_one)
        dble_eps = dble_eps * 0.5d0
      enddo

      real_eps = 1.0
      do while ((1.0+(0.5*real_eps)) .ne. real_one)
        real_eps = real_eps * 0.5
      enddo

      print'(1x,e24.16,a)',dble_eps ,'~Double epsilon (Expr. compare)'
     &                    ,real_eps ,'~Single epsilon (Expr. compare)'

      end
