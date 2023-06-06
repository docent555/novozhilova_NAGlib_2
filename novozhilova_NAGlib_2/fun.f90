module fun
   use, intrinsic :: iso_c_binding
   use ifcore

   integer(c_int) ne, nt, nz, freq_out, neqp, lenwrk, l, method, neqf, lwork, liwork, nrd, it, ifailp, inharm
   real(c_double) zex, dz, tend, dtr(2), q(3), icu(2), th(2), a(2), dcir(2), r(2), &
      f0(3), dt, pitch, f10, f20, f30, p10, p20, p30, ftol, ptol, hstart, zstart, xout, nharm, &
      gamma, ukv, betta, betta2, betta_z, betta_z2, betta_perp, betta_perp2, gmp, w_op, c, e, m, b(2), w_n(2), ia(2), norm, nfac
   complex(c_double_complex) fp(2)
   logical(c_bool) errass, wc, fok, lensm, btod, iatoi

   common xout, it

   integer(c_int) breaknum(3)
   real(c_double) phitmp0(3), phitmp1(3)
   complex(c_double_complex) fc, fcomp(3)

   integer(c_int), allocatable, target :: idxre(:, :), idxim(:, :), iwork(:), idxp(:, :)
   complex(c_double_complex), allocatable, target :: u(:), mean(:)
   real(c_double), allocatable, target :: tax(:), zax(:), eta(:, :), etag(:, :), w(:, :), f(:, :), p(:, :), &
                                          phi(:, :), phios(:, :), wos(:, :), &
                                          pgot(:), ppgot(:), pmax(:), thres(:), workp(:), work(:), dp2dz(:, :), mean2(:, :), &
                                          cl1(:), lhs1(:), rhs1(:), cl2(:), lhs2(:), rhs2(:)

   include 'z.inc'
   include 're.inc'
   include 'im.inc'

   complex(c_double_complex), parameter :: ic = (0.0d0, 1.0d0)
   real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)

   private freq_out, zex, tend, dtr, q, i, th, a, dcir, r, f0, pitch

contains

   subroutine init()
      implicit none

      integer(c_int) ii

      call read_param()

      nharm = dble(inharm)
      gamma = 1.0 + ukv/511.0
      betta = dsqrt(1.0d0 - 1.0d0/(gamma*gamma))
      betta2 = 1.0d0 - 1.0d0/(gamma*gamma)
      betta_z = betta/dsqrt(pitch*pitch + 1.0d0)
      betta_z2 = betta2/(pitch*pitch + 1.0d0)
      betta_perp2 = betta2 - betta_z2
      gmp = 0.048715056967419
      w_op = 2*pi*w_op*1e9
      c = 29979245800.0d0
      e = 4.803e-10
      m = 9.1093837015e-28

      if (lensm .eq. .true.) zex = betta_perp2/2.0d0/betta_z*w_op*zex/nharm/c
      print *, 'Zex = ', zex

      if (btod .eq. .true.) then
         w_n(1) = e*B(1)/(m*c)*10000.0d0
         dtr(1) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(1))/(gamma*w_op))
         w_n(2) = e*B(2)/(m*c)*10000.0d0
         dtr(2) = (2.0/betta_perp2)*(1.0 - (2.0*w_n(2))/(gamma*w_op))
      end if            
      
      call norma(norm)

      if (iatoi .eq. .true.) then         
         nfac = factorial(inharm)         
         icu(1) = 2.35/10000*IA(1)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(1)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
         icu(2) = 2.35/10000*IA(2)*(nharm**(inharm + 1)/2.0**(inharm - 1)/nfac)**2 &
                  *(Q(2)*gmp*betta**(2*inharm - 4)/gamma/betta_z)/norm
      end if      
      
      print *, 'dtr1 = ', dtr(1), 'dtr2 = ', dtr(2)
      print *, 'icu1 = ', icu(1), 'icu2 = ', icu(2)      

      nt = tend/dt + 1
      nz = zex/dz + 1

      neqp = 4*ne
      lenwrk = 32*neqp
      zstart = 0.0d0
      errass = .false.
      hstart = 0.0d0
      method = 2
      ifailp = 1

      neqf = 6
      nrd = 6
      lwork = 8*neqf + 5*nrd + 21
      liwork = nrd + 21

      call allocate_arrays()

      do l = 1, neqp
         thres(l) = 1.0d-8
      end do

      f(1, 1) = f10
      f(2, 1) = p10
      f(3, 1) = f20
      f(4, 1) = p20
      f(5, 1) = f30
      f(6, 1) = p30

      do ii = 1, nt
         tax(ii) = (ii - 1)*dt
      end do

      do ii = 1, nz
         zax(ii) = (ii - 1)*dz
      end do

      call calc_u(u, zex, nz, zax)

      do ii = 1, 2
         idxre(ii, :) = (/2*(ii - 1)*ne + 1:(2*ii - 1)*ne/)
         idxim(ii, :) = (/(2*ii - 1)*ne + 1:2*ii*ne/)
      end do

      idxp(1, :) = (/1:ne/)
      idxp(2, :) = (/ne + 1:2*ne/)

      do i = 1, ne
         p(i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(2*ne + i, 1) = dreal(cdexp(ic*(i - 1)/dble(ne)*2*pi))
         p(3*ne + i, 1) = dimag(cdexp(ic*(i - 1)/dble(ne)*2*pi))
      end do

   end subroutine init
   
   subroutine norma(norm)
      use, intrinsic :: iso_c_binding, only: c_double, c_double_complex, c_int
      import, only:rea, ima, betta_perp2, betta_z, inharm, c, w_op!, omega
      implicit none

      real(c_double) :: dzz = 0.0280211, norm
      complex(c_double_complex) :: u(663)

      dzz = betta_perp2/2.0d0/betta_z*w_op*dzz/inharm/c
      u = dcmplx(rea, ima)

      norm = sum(cdabs(u(:))*cdabs(u(:)))*dzz

      print *, 'N = ', norm

   end subroutine norma

   recursive function factorial(p) result(l)
      integer, intent(in) :: p
      integer l
      if (p == 1) then
         l = 1
      else
         l = p*factorial(p - 1)
      end if
   end function

   function squval(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) squval
      real(c_double) z, re, im, dz, z1
      integer(c_int) l

      dz = 0.280211
      z = zz/zex*185.5d0
      l = z/dz

      if (l .eq. 0) then
         l = 2
      elseif (l .ge. 662) then
         l = 662
      else
         if ((z - za(l)) .gt. 0.5*dz) l = l + 1
      end if

      z1 = za(l)
      z = z - 8.5d0

      re = rea(l - 1) + ((-rea(l - 1) + rea(l))*(dz + z - z1))/dz + ((rea(l - 1)/2.0d0 - rea(l) + &
                                                                      rea(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      im = ima(l - 1) + ((-ima(l - 1) + ima(l))*(dz + z - z1))/dz + ((ima(l - 1)/2.0d0 - ima(l) + &
                                                                      ima(l + 1)/2.0d0)*(z - z1)*(dz + z - z1))/dz/dz
      !!!NE RABOTAET
      !re = ((rea(l - 1) - 2*rea(l) + rea(l + 1))*z**2)/(2.*dz**2) &
      !     + (z*(dz*(-rea(l - 1) + rea(l + 1)) - 2*(rea(l - 1) - 2*rea(l) + rea(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*rea(l) + dz*(rea(l - 1) - rea(l + 1))*z1 + (rea(l - 1) - 2*rea(l) + rea(l + 1))*z1**2)/(2.*dz**2)
      !im = ((ima(l - 1) - 2*ima(l) + ima(l + 1))*z**2)/(2.*dz**2) + &
      !     (z*(dz*(-ima(l - 1) + ima(l + 1)) - 2*(ima(l - 1) - 2*ima(l) + ima(l + 1))*z1))/(2.*dz**2) + &
      !     -(2*dz**2*ima(l) + dz*(ima(l - 1) - ima(l + 1))*z1 + (ima(l - 1) - 2*ima(l) + ima(l + 1))*z1**2)/(2.*dz**2)

      squval = dcmplx(re, im)

   end function squval

   function uval(zz)

      implicit none

      real(c_double), intent(in) :: zz

      complex(c_double_complex) uval
      real(c_double) z, re, im, d
      integer(c_int) l

      z = zz/zex*185.5 - 8.5
      l = (z + 8.5)/0.28021 + 1
      d = z - za(l)

      !print *, z, l, d

      if (d .gt. 0.0 .and. l /= 663) then
         re = (rea(l)*za(l + 1) - rea(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (rea(l + 1) - rea(l))/(za(l + 1) - za(l))*z
         im = (ima(l)*za(l + 1) - ima(l + 1)*za(l))/(za(l + 1) - za(l)) + &
              (ima(l + 1) - ima(l))/(za(l + 1) - za(l))*z
      else if (d .lt. 0.0 .and. l /= 1) then
         re = (rea(l - 1)*za(l) - rea(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (rea(l) - rea(l - 1))/(za(l) - za(l - 1))*z
         im = (ima(l - 1)*za(l) - ima(l)*za(l - 1))/(za(l) - za(l - 1)) + &
              (ima(l) - ima(l - 1))/(za(l) - za(l - 1))*z
      else
         re = rea(l)
         im = ima(l)
      end if

      uval = dcmplx(re, im)

   end function uval

   subroutine allocate_arrays()
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int) err_alloc

      allocate (f(6, nt), p(neqp, nz), u(nz), tax(nt), zax(nz), mean(nz), eta(2, nt), etag(2, nt), w(3, nt), &
                idxre(2, ne), idxim(2, ne), pgot(neqp), ppgot(neqp), pmax(neqp), thres(neqp), workp(lenwrk), work(lwork), &
                iwork(liwork), wos(3, nt), phi(3, nt), phios(3, nt), dp2dz(2*ne, nz), idxp(2, ne), mean2(2, nz - 1), &
                cl1(nt), lhs1(nt), rhs1(nt), cl2(nt), lhs2(nt), rhs2(nt), stat=err_alloc)

      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if
   end subroutine allocate_arrays

   subroutine deallocate_arrays()
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int) err_dealloc

      deallocate (f, p, u, tax, zax, mean, eta, etag, w, stat=err_dealloc)

      if (err_dealloc /= 0) then
         print *, "deallocation error"
         pause
         stop
      end if
   end subroutine deallocate_arrays

   subroutine read_param() bind(c, name='read_param')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2

      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2

      open (unit=1, file='input_fortran.in', status='old', err=101)
      read (unit=1, nml=param, err=102)
      close (unit=1)

      q(1) = q1
      q(2) = q2
      q(3) = q3
      icu(1) = i1
      icu(2) = i2
      th(1) = th1
      th(2) = th2
      a(1) = a1
      a(2) = a2
      dtr(1) = dtr1
      dtr(2) = dtr2
      dcir(1) = dcir1
      dcir(2) = dcir2
      r(1) = r1
      r(2) = r2
      b(1) = b1
      b(2) = b2
      ia(1) = ia1
      ia(2) = ia2

      write (*, nml=param)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine read_param

   subroutine write_param() bind(c, name='write_param')
      use, intrinsic :: iso_c_binding
      import
      implicit none

      namelist /param/ ne, tend, zex, q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, &
         dcir1, dcir2, r1, r2, f10, f20, f30, p10, p20, p30, dt, dz, pitch, ftol, ptol, wc, fok, inharm, ukv, &
         w_op, lensm, btod, b1, b2, iatoi, ia1, ia2

      real(c_double) q1, q2, q3, i1, i2, th1, th2, a1, a2, dtr1, dtr2, dcir1, dcir2, r1, r2, b1, b2, ia1, ia2
      logical wc

      open (unit=1, file='input_fortran.in', status='old', err=101)
      !print *, 'OK1'
      read (unit=1, nml=param, err=102)
      close (unit=1)

      f10 = f(1, nt)
      f20 = f(3, nt)
      f30 = f(5, nt)
      p10 = mod(f(2, nt), 2*pi)
      p20 = mod(f(4, nt), 2*pi)
      p30 = mod(f(6, nt), 2*pi)

      open (unit=1, file='input_fortran_next.in', err=101)
      !print *, 'OK2'
      write (1, nml=param)
      close (unit=1)

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of reading file "input_fortran.in"'; pause; stop
   end subroutine write_param

   subroutine write_results()

      implicit none

      integer i, j

      if (wc .eq. .true.) then
         w(:, 1) = 0
         do i = 2, nt
            do j = 1, 3
               !w(j, i - 1) = dimag(log(f(2*j - 1, i)*cdexp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*cdexp(ic*f(2*j, i - 1)))))/dt
               w(j, i) = (f(2*j, i) - f(2*j, i - 1))/dt
            end do
         end do
         print *, 'Frequency calculated from phase. ( WC = ', wc, ')'
      elseif (wc .eq. .false.) then
         call freq()
         print *, 'Frequency calculated from RHS. ( WC = ', wc, ')'
      end if

      phi(:, 1) = 0; 
      do i = 2, nt
         do j = 1, 3
            phi(j, i) = phi(j, i - 1) + dimag(log(f(2*j - 1, i)*cdexp(ic*f(2*j, i))/(f(2*j - 1, i - 1)*cdexp(ic*f(2*j, i - 1)))))
         end do
      end do

      breaknum(:) = 0
      fcomp(1) = f(2*1 - 1, 1)*cdexp(ic*f(2*1, 1))
      fcomp(2) = f(2*2 - 1, 1)*cdexp(ic*f(2*2, 1))
      fcomp(3) = f(2*3 - 1, 1)*cdexp(ic*f(2*3, 1))
      phitmp0(:) = datan2(dimag(fcomp(:)), dreal(fcomp(:)))
      !phitmp0(:) = datan2(dimag(f(:, 1)), dreal(f(:, 1)))
      phios(:, 1) = phitmp0(:)
      do i = 2, nt
         do j = 1, 3
            fc = f(2*j - 1, i)*cdexp(ic*f(2*j, i))
            phitmp1(j) = datan2(dimag(fc), dreal(fc))
            if ((phitmp1(j) - phitmp0(j)) .gt. pi) breaknum(j) = breaknum(j) - 1
            if ((phitmp1(j) - phitmp0(j)) .lt. -pi) breaknum(j) = breaknum(j) + 1
            phios(j, i) = phitmp1(j) + 2.*pi*breaknum(j)
            !phios(j, i) = phitmp1(j)
            phitmp0(j) = phitmp1(j)
         end do
      end do

      do i = 1, nt - 1
         do j = 1, 3
            wos(j, i) = (phios(j, i + 1) - phios(j, i))/dt
         end do
      end do

      write (*, '(/)')

      !pause

      open (3, file='cl1.dat')
      do i = 1, nt
         write (3, '(5e17.8)') tax(i), cl1(i), lhs1(i), rhs1(i), abs(cl1(i)/lhs1(i))
      end do
      close (3)

      open (3, file='cl2.dat')
      do i = 1, nt
         write (3, '(5e17.8)') tax(i), cl2(i), lhs2(i), rhs2(i), abs(cl2(i)/lhs2(i))
      end do
      close (3)

      open (1, file='F.dat')
      do i = 1, nt
         !write (1, '(4e17.8)') tax(i), dabs(f(1, i)), dabs(f(3, i)), dabs(f(5, i))
         write (1, '(4e17.8)') tax(i), f(1, i), f(3, i), f(5, i)
      end do
      close (1)

      open (13, file='FCMPLX.dat')
      do i = 1, nt
         fcomp(1) = f(2*1 - 1, i)*cdexp(ic*f(2*1, i))
         fcomp(2) = f(2*2 - 1, i)*cdexp(ic*f(2*2, i))
         fcomp(3) = f(2*3 - 1, i)*cdexp(ic*f(2*3, i))
         write (13, '(7e17.8)') tax(i), dreal(fcomp(1)), dimag(fcomp(1)), dreal(fcomp(2)), dimag(fcomp(2)), &
            dreal(fcomp(3)), dimag(fcomp(3))
      end do
      close (13)

      open (2, file='E.dat')
      do i = 1, nt
         write (2, '(5e17.8)') tax(i), eta(1, i), etag(1, i), eta(2, i), etag(2, i)
      end do
      close (2)

      open (3, file='W.dat')
      do i = 1, nt
         write (3, '(4e17.8)') tax(i), w(1, i), w(2, i), w(3, i)
      end do
      close (3)

      open (1, file='P.dat')
      do i = 1, nt
         !write (1, '(4e17.8)') tax(i), phi(1, i), phi(2, i), phi(3, i)
         write (1, '(4e17.8)') tax(i), f(2, i), f(4, i), f(6, i)
      end do
      close (1)

      open (1, file='POS.dat')
      do i = 1, nt
         write (1, '(4e17.8)') tax(i), phios(1, i), phios(2, i), phios(3, i)
      end do
      close (1)

      open (3, file='WOS.dat')
      do i = 1, nt - 1
         write (3, '(4e17.8)') tax(i + 1), wos(1, i), wos(2, i), wos(3, i)
      end do
      close (3)

      call write_param()

      return

101   print *, 'error of file open.'
      pause
      stop
102   print *, 'error of file reading.'
      pause
      stop
103   print *, 'error of file writing.'
      pause
      stop
   end subroutine write_results

   subroutine ode4f()
      import
      implicit none

      external d02nbz, d02nby

      integer(c_int) i, j, ifailf
      real(c_double) :: h, t, zwant, zgot
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27

      integer(c_int), parameter :: neq = 6, neqmax = neq, nrw = 50 + 4*neqmax, ninf = 23, nwkjac = neqmax*(neqmax + 1), &
                                   maxord = 5, ny2dim = maxord + 1, maxstp = 0, mxhnil = 5
      real(c_double), parameter :: h0 = 0.0d0, hmax = 10.0d0, hmin = 1.0d-10
      integer(c_int) itask, itol, itrace, inform(ninf)
      real(c_double) atolnag(neqmax), const(6), rtolnag(neqmax), rwork(nrw), wkjac(nwkjac), y(neqmax), ydot(neqmax), &
         ysave(neqmax, ny2dim), tout, atol(neqmax), rtol(neqmax), tcrit, pex(neqp)
      logical(c_bool), parameter ::    petzld = .false.

      !solve eq. at t=0
      call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifailp)

      fp(1) = f(1, 1)*cdexp(ic*f(2, 1))
      fp(2) = f(3, 1)*cdexp(ic*f(4, 1))

      do i = 1, nz - 1
         !zwant = i*dz
         zwant = zax(i + 1)
         call d02pcf(dpdz, zwant, zgot, pgot, ppgot, pmax, workp, ifailp)

         if (ifailp .ne. 0) then
            write (*, *)
            write (*, 99998) 'exit d02pcf with ifail = ', ifailp, '  and t = ', t
            pause
            stop
         end if

         p(:, i + 1) = pgot
      end do
99998 format(1x, a, i2, a, d12.5)

      eta(:, 1) = eff(p(:, nz))
      etag(:, 1) = pitch**2/(pitch**2 + 1.0d0)*eta(:, 1)

      !open (1, file='test2.dat')
      !do i = 1, ne
      !    write(1,'(i,f17.8)') i, p(i,nz)
      !end do
      !close (1)
      !stop

      call d02nvf(neqmax, ny2dim, maxord, 'newton', petzld, const, tcrit, hmin, hmax, &
                  h0, maxstp, mxhnil, 'average-l2', rwork, ifailf)
      call d02nsf(neq, neqmax, 'a', nwkjac, rwork, ifailf)

      t = 0.0d0
      tout = tend
      itask = 1
      itol = 1
      rtol(1) = ftol
      atol(1) = 0.001*ftol
      do i = 1, 6
         const(i) = 0.0d0
      end do
      tcrit = tout
      y(:) = f(:, 1)
      itrace = -1
      xout = dt
      it = 2
      ifailf = 1

      call d02nbf(neq, neqmax, t, tout, y, ydot, rwork, rtol, atol, itol, inform, dfdt, ysave, &
                  ny2dim, jac, wkjac, nwkjac, monitr, itask, itrace, ifailf)

      if (ifailf .ne. 0) then
         write (*, *)
         write (*, 99998) 'exit d02nbf with ifail = ', ifailf, '  and t = ', t
         pause
         stop
      end if

      !call calcpex(y, pex, cl1(it), lhs1(it), rhs1(it), cl2(it), lhs2(it), rhs2(it))
      !eta(:, nt) = eff(p(:, nz))
      !etag(:, nt) = pitch**2/(pitch**2 + 1)*eta(:, nt)
      !do j = 1, neqf
      !   f(j, nt) = y(j)
      !end do

   end subroutine ode4f

   function eff(pex) result(eta)
      use, intrinsic :: iso_c_binding, only: c_double, c_int
      import, only:ne, idxre, idxim

      implicit none

      integer(c_int) i
      real(c_double) eta(2)
      real(c_double), intent(in) :: pex(:)

      do i = 1, 2
         eta(i) = 1 - sum(cdabs(dcmplx(pex(idxre(i, :)), pex(idxim(i, :))))**2)/ne
      end do

      !open(1, file='test1.dat')
      !do i=1,ne
      !   write(1,'(i,f17.8)') i, pex(i)
      !enddo
      !close(1)
      !stop

   end function eff

   subroutine dpdz(z, p, prhs)
      import :: ne, zex, f, ic, dtr
      implicit none

      real(c_double) z, p(*), prhs(*)

      integer(c_int) i, reidx(ne), imidx(ne)
      complex(c_double_complex) :: u
      complex(c_double_complex) s(ne), ptmp(ne)

      if (fok .eq. .false.) then
         u = dexp(-3*((z - zex/2)/(zex/2))**2)
      else
         u = squval(z)
      end if

      do i = 1, 2
         ptmp = dcmplx(p(idxre(i, :)), p(idxim(i, :)))

         !s = ic*(fp(i)*u*dconjg(ptmp) - (dtr(i) + cdabs(ptmp)**2 - 1)*ptmp)
         s = ic*(fp(i)*u*dconjg(ptmp)**(inharm - 1) - (dtr(i) + cdabs(ptmp)**2 - 1)*ptmp)

         prhs(idxre(i, :)) = dreal(s)
         prhs(idxim(i, :)) = dimag(s)
      end do
   end subroutine dpdz

   complex(c_double_complex) function xi(p, num)
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
      import, only:ne, nz, mean, u, dz, idxre, idxim, inharm

      implicit none

      integer(c_int) i, num
      real(c_double) p(:, :)

      do i = 1, nz
         mean(i) = sum(dcmplx(p(idxre(num, :), i), p(idxim(num, :), i))**inharm, 1)/ne
      end do

      mean = dconjg(u)*mean

      !xi1 = (0.5d0*(mean(1) + mean(2)) + sum(mean(2:nz - 1)))*dz
      xi = (0.5d0*(mean(1) + mean(nz)) + sum(mean(2:(nz - 1))))*dz

   end function

   subroutine dfdt(neqf, t, f, s, ires) ! nag
      implicit none

      integer(c_int) :: ii, jj, iter_num = 1, time_num = 1, neqf, ires
      real(c_double) t, f(neqf), s(neqf), &
         x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2, zwant, zgot, e1, e2
      complex(c_double_complex) x1, x2
      logical bp

      call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifailp)

      fp(1) = f(1)*cdexp(ic*f(2))
      fp(2) = f(3)*cdexp(ic*f(4))

      do jj = 1, nz - 1
         zwant = zax(jj + 1)
         call d02pcf(dpdz, zwant, zgot, pgot, ppgot, pmax, workp, ifailp)

         if (ifailp .ne. 0) then
            write (*, *)
            write (*, 99998) 'exit d02pcf with ifail = ', ifailp, '  and t = ', t
            pause
            stop
         end if
99998    format(1x, a, i2, a, d12.5)

         p(:, jj + 1) = pgot
      end do

      !if (t>10) then
      !open (1, file='test.dat')
      !do jj = 1, nz
      !    write(1,'(3f17.8)') zax(jj), p(100,jj), p(105,jj)
      !end do
      !close (1)
      !stop
      !endif

      !x1 = xi(p(1:2*ne, :), 1)
      !x2 = xi(p(2*ne + 1:4*ne, :), 1)
      x1 = xi(p, 1)
      x2 = xi(p, 2)

      x1r = dreal(x1)
      x1i = dimag(x1)
      x2r = dreal(x2)
      x2i = dimag(x2)

      f1 = f(1)
      phi1 = f(2)
      f2 = f(3)
      phi2 = f(4)
      f3 = f(5)
      phi3 = f(6)

      q31 = q(3)/q(1)
      i1 = icu(1)
      r1 = r(1)
      th1 = th(1)
      dcir1 = dcir(1)
      cos1 = dcos(phi1)
      sin1 = dsin(phi1)

      q32 = q(3)/q(2)
      i2 = icu(2)
      r2 = r(2)
      th2 = th(2)
      dcir2 = dcir(2)
      cos2 = dcos(phi2)
      sin2 = dsin(phi2)

      q3 = q(3)
      a1 = a(1)
      a2 = a(2)

      s(1) = (-nharm*f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*nharm*f3*dcos(phi3 - phi1 - th1))*q31
      s(2) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*nharm*(f3/f1)*dsin(phi3 - phi1 - th1))*q31

      s(3) = (-nharm*f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*nharm*f3*dcos(phi3 - phi2 - th2))*q32
      s(4) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*nharm*(f3/f2)*dsin(phi3 - phi2 - th2))*q32

      s(5) = -f3 + a1*f1*dcos(phi1 - phi3) + a2*f2*dcos(phi2 - phi3)
      s(6) = a1*f1/f3*dsin(phi1 - phi3) + a2*f2/f3*dsin(phi2 - phi3)

   end subroutine dfdt

   subroutine calc_u(u, zex, nz, zax)
      import
      implicit none

      integer(c_int), intent(in) :: nz
      real(c_double), intent(in) :: zex, zax(nz)
      complex(c_double_complex), intent(out) :: u(:)

      integer(c_int) i

      if (fok .eq. .false.) then
         do i = 1, nz
            u(i) = dexp(-3*((zax(i) - zex/2)/(zex/2))**2)
         end do
      else
         do i = 1, nz
            u(i) = squval(zax(i))
         end do
      end if

      !open(1, file = 'test.dat')
      !do i = 1,nz
      !   write(1, '(i,2f12.6)') i, dreal(u(i)), dimag(u(i))
      !end do
      !close(1)
      !stop

   end subroutine

   subroutine jac(neq, t, y, h, d, pp)
      implicit none
!     .. scalar arguments ..
      double precision d, h, t
      integer neq
!     .. array arguments ..
      double precision pp(neq, neq), y(neq)
!     .. local scalars ..
      double precision hxd
!     .. executable statements ..

      real(c_double) x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q333, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2, ph311, ph322, phi13, phi23, &
         sin311, sin322, cos311, cos322, cos13, sin13, cos23, sin23, phi131, phi232
      complex(c_double_complex) x1, x2

      hxd = h*d

      !x1 = xi(p(1:2*ne, :), 1)
      !x2 = xi(p(2*ne + 1:4*ne, :), 1)
      x1 = xi(p, 1)
      x2 = xi(p, 2)

      x1r = dreal(x1)
      x1i = dimag(x1)
      x2r = dreal(x2)
      x2i = dimag(x2)

      f1 = y(1)
      phi1 = y(2)
      f2 = y(3)
      phi2 = y(4)
      f3 = y(5)
      phi3 = y(6)

      q31 = q(3)/q(1)
      i1 = icu(1)
      r1 = r(1)
      th1 = th(1)
      dcir1 = dcir(1)
      cos1 = dcos(phi1)
      sin1 = dsin(phi1)

      q32 = q(3)/q(2)
      i2 = icu(2)
      r2 = r(2)
      th2 = th(2)
      dcir2 = dcir(2)
      cos2 = dcos(phi2)
      sin2 = dsin(phi2)

      ph311 = phi3 - phi1 - th1
      ph322 = phi3 - phi2 - th2
      phi13 = phi1 - phi3
      phi23 = phi2 - phi3
      phi131 = phi1 - phi3 + th1
      phi232 = phi2 - phi3 + th2

      cos311 = dcos(ph311)
      sin311 = dsin(ph311)
      cos322 = dcos(ph322)
      sin322 = dsin(ph322)
      cos13 = dcos(phi13)
      sin13 = dsin(phi13)
      cos23 = dcos(phi23)
      sin23 = dsin(phi23)

      !pp(1, 1) = 1.0d0 + hxd*q31*nharm
      !pp(1, 2) = -hxd*q31*(i1*(x1i*sin1 + x1r*cos1) + 2.0d0*nharm*r1*f3*sin311)
      !!pp(1,3) = 0
      !!pp(1,4) = 0
      !pp(1, 5) = -hxd*2.0d0*nharm*r1*q31*cos311
      !pp(1, 6) = hxd*2.0d0*nharm*r1*q31*f3*sin311
      !
      !pp(2, 1) = hxd*q31/(f1**2)*(i1*(x1r*cos1 + x1i*sin1) + 2.0d0*nharm*r1*f3*sin311)
      !pp(2, 2) = 1.0d0 - hxd*q31/f1*(i1*(-x1r*sin1 + x1i*cos1) - 2.0d0*nharm*r1*f3*cos311)
      !!pp(2,3) = 0
      !!pp(2,4) = 0
      !pp(2, 5) = -hxd*2.0d0*nharm*r1*q31/f1*sin311
      !pp(2, 6) = -hxd*2.0d0*nharm*r1*q31*f3/f1*cos311
      !
      !!pp(3,1) = 0
      !!pp(3,2) = 0
      !pp(3, 3) = 1.0d0 + hxd*q32*nharm
      !pp(3, 4) = -hxd*q32*(i2*(x2i*sin2 + x2r*cos2) + 2.0d0*nharm*r2*f3*sin322)
      !pp(3, 5) = -hxd*2.0d0*nharm*r2*q32*cos322
      !pp(3, 6) = hxd*2.0d0*nharm*r2*q32*f3*sin322
      !
      !!pp(4,1) = 0
      !!pp(4,2) = 0
      !pp(4, 3) = hxd*q32/(f2**2)*(i2*(x2r*cos2 + x2i*sin2) + 2.0d0*nharm*r2*f3*sin322)
      !pp(4, 4) = 1.0d0 - hxd*q32/f2*(i2*(-x2r*sin2 + x2i*cos2) - 2.0d0*nharm*r2*f3*cos322)
      !pp(4, 5) = -hxd*2.0d0*nharm*r2*q32/f2*sin322
      !pp(4, 6) = -hxd*2.0d0*nharm*r2*q32*f3/f2*cos322
      !
      !pp(5, 1) = -hxd*a1*cos13
      !pp(5, 2) = hxd*a1*f1*sin13
      !pp(5, 3) = -hxd*a2*cos23
      !pp(5, 4) = hxd*a2*f2*sin23
      !pp(5, 5) = 1.0d0 + hxd
      !pp(5, 6) = -hxd*(a1*f1*sin13 + a2*f2*sin23)
      !
      !pp(6, 1) = -hxd*a1/f3*sin13
      !pp(6, 2) = -hxd*a1*f1/f3*cos13
      !pp(6, 3) = -hxd*a2/f3*sin23
      !pp(6, 4) = -hxd*a2*f2/f3*cos23
      !pp(6, 5) = hxd/(f3**2)*(a1*f1*sin13 + a2*f2*sin23)
      !pp(6, 6) = 1.0d0 + hxd/f3*(a1*f1*cos13 + a2*f2*cos23)

      pp(1, 1) = hxd*nharm*q31 + 1.0D0
      pp(1, 2) = -hxd*q31*(i1*(x1r*cos1 + x1i*sin1) - f3*nharm*r1*sin(phi131)*2.0D0)
      !pp(1, 3) = 0.0D0
      !pp(1, 4) = 0.0D0
      pp(1, 5) = hxd*nharm*q31*r1*cos(phi131)*(-2.0D0)
      pp(1, 6) = f3*hxd*nharm*q31*r1*sin(phi131)*(-2.0D0)

      pp(2, 1) = hxd*q31*(1.0D0/f1**2*i1*(x1r*cos1 + x1i*sin1) - 1.0D0/f1**2*f3*nharm*r1*sin(phi131)*2.0D0)
      pp(2, 2) = -hxd*q31*((i1*(x1i*cos1 - x1r*sin1))/f1 - (f3*nharm*r1*cos(phi131)*2.0D0)/f1) + 1.0D0
      !pp(2, 3) = 0.0D0
      !pp(2, 4) = 0.0D0
      pp(2, 5) = (hxd*nharm*q31*r1*sin(phi131)*2.0D0)/f1
      pp(2, 6) = (f3*hxd*nharm*q31*r1*cos(phi131)*(-2.0D0))/f1

      !pp(3, 1) = 0.0D0
      !pp(3, 2) = 0.0D0
      pp(3, 3) = hxd*nharm*q32 + 1.0D0
      pp(3, 4) = -hxd*q32*(i2*(x2r*cos2 + x2i*sin2) - f3*nharm*r2*sin(phi232)*2.0D0)
      pp(3, 5) = hxd*nharm*q32*r2*cos(phi232)*(-2.0D0)
      pp(3, 6) = f3*hxd*nharm*q32*r2*sin(phi232)*(-2.0D0)

      !pp(4, 1) = 0.0D0
      !pp(4, 2) = 0.0D0
      pp(4, 3) = hxd*q32*(1.0D0/f2**2*i2*(x2r*cos2 + x2i*sin2) - 1.0D0/f2**2*f3*nharm*r2*sin(phi232)*2.0D0)
      pp(4, 4) = -hxd*q32*((i2*(x2i*cos2 - x2r*sin2))/f2 - (f3*nharm*r2*cos(phi232)*2.0D0)/f2) + 1.0D0
      pp(4, 5) = (hxd*nharm*q32*r2*sin(phi232)*2.0D0)/f2
      pp(4, 6) = (f3*hxd*nharm*q32*r2*cos(phi232)*(-2.0D0))/f2

      pp(5, 1) = -a1*hxd*cos13
      pp(5, 2) = a1*f1*hxd*sin13
      pp(5, 3) = -a2*hxd*cos23
      pp(5, 4) = a2*f2*hxd*sin23
      pp(5, 5) = hxd + 1.0D0
      pp(5, 6) = -hxd*(a1*f1*sin13 + a2*f2*sin23)

      pp(6, 1) = -(a1*hxd*sin13)/f3
      pp(6, 2) = -(a1*f1*hxd*cos13)/f3
      pp(6, 3) = -(a2*hxd*sin23)/f3
      pp(6, 4) = -(a2*f2*hxd*cos23)/f3
      pp(6, 5) = hxd*(a1*f1*1.0D0/f3**2*sin13 + a2*f2*1.0D0/f3**2*sin23)
      pp(6, 6) = hxd*((a1*f1*cos13)/f3 + (a2*f2*cos23)/f3) + 1.0D0

      return
   end

   subroutine monitr(n, nmax, t, hlast, h, y, ydot, ysave, r, acor, imon, inln, hmin, hmxi, nqu)
      implicit none
      !..parameters..
      integer nout, it, j
      parameter(nout=6)
      integer ny2dim
      parameter(ny2dim=6)
      !..scalar arguments..
      double precision h, hlast, hmin, hmxi, t
      integer imon, inln, n, nmax, nqu
      !..array arguments..
      double precision acor(nmax, 2), r(n), y(n), ydot(n), ysave(nmax, *)
      !..scalars in common..
      double precision xout
      !..local scalars..
      integer i, ifail
      logical(4) pressed
      character(1) key
      integer(c_int), parameter :: esc = 27

      real(c_double) pex(neqp)

      !..external subroutines..
      external d02xkf
      !..common blocks..
      common xout, it
      !..executable statements..
      if (imon .ne. 1) return
20    if (.not. (t - hlast .lt. xout .and. xout .le. t)) return
      ifail = 1
      !c1 interpolation
      call d02xkf(xout, r, n, ysave, nmax, ny2dim, acor(1, 2), n, t, nqu, hlast, h, ifail)

      if (ifail .ne. 0) then
         imon = -2
      else
         !write (nout, 99999) xout, (r(i), i=1, n)
         call calcpex(r, pex, cl1(it), lhs1(it), rhs1(it), cl2(it), lhs2(it), rhs2(it))
         !call zs(r, pex, cl1(it), lhs1(it), rhs1(it), cl2(it), lhs2(it), rhs2(it))
         eta(:, it) = eff(pex)
         etag(:, it) = pitch**2/(pitch**2 + 1)*eta(:, it)
         write (*, '(a,f8.3,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a,f5.3,a,f5.3,a,\,a)') 't =', xout, &
            '  |f1| =', r(1), '  |f2| =', r(3), '  |f3| =', r(5), '  e1 =', eta(1, it), '  e2 =', eta(2, it), &
            '  w1 = ', ydot(2), '  w2 = ', ydot(4), '  w3 = ', ydot(6), '  c1 =', dabs(cl1(it)/rhs1(it))*100, '%  c2 =', dabs(cl2(it)/rhs2(it)*100), '%', char(13)
         do j = 1, n
            f(j, it) = r(j)
         end do
         xout = xout + dt
         nt = it
         it = it + 1
         pressed = peekcharqq()
         if (pressed) then
            key = getcharqq()
            if (ichar(key) .eq. esc) then
               write (*, '(/,a)') 'quit?'
               key = getcharqq()
               if (ichar(key) .eq. 121 .or. ichar(key) .eq. 89) then
                  nt = it - 1
                  call write_results()
                  !imon = -2
                  !return
               end if
            end if
         end if
         if (xout .lt. tend) go to 20
      end if

      return

99999 format(1x, f8.3, 6(f13.5, 2x))
   end subroutine monitr

   subroutine calcpex(fpex, pex, c1, lhs1, rhs1, c2, lhs2, rhs2)

      implicit none

      real(8), intent(out) :: c1, lhs1, rhs1, c2, lhs2, rhs2

      integer(c_int) i, idx(ne)
      real(c_double) :: zwant, zgot, pex(neqp), fpex(6), p2ex_mean(2), p20_mean(2)
      complex(8) ptmp(ne)

      call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifailp)

      fp(1) = fpex(1)*cdexp(ic*fpex(2))
      fp(2) = fpex(3)*cdexp(ic*fpex(4))

      call d02pcf(dpdz, zex, zgot, pex, ppgot, pmax, workp, ifailp)

      do i = 1, 2
         idx = idxp(i, :)

         ptmp(:) = dcmplx(p(idxre(i, :), 1), p(idxim(i, :), 1))
         p20_mean(i) = sum(cdabs(ptmp(:)*cdabs(ptmp(:))))/ne

         ptmp(:) = dcmplx(pex(idxre(i, :)), pex(idxim(i, :)))
         p2ex_mean(i) = sum(cdabs(ptmp(:)*cdabs(ptmp(:))))/ne
      end do

      lhs1 = 4*fpex(1)**2 - 8*r(1)*fpex(5)*fpex(1)*dcos(th(1) - fpex(6) + fpex(2))
      rhs1 = -icu(1)*(p2ex_mean(1) - p20_mean(1))
      c1 = lhs1 - rhs1

      lhs2 = 4*fpex(3)**2 - 8*r(2)*fpex(5)*fpex(3)*dcos(th(2) - fpex(6) + fpex(4))
      rhs2 = -icu(2)*(p2ex_mean(2) - p20_mean(2))
      c2 = lhs2 - rhs2

   end subroutine calcpex

   subroutine zs(f, pex, c1, lhs1, rhs1, c2, lhs2, rhs2)

      implicit none

      real(8), intent(in) :: f(neqf)
      real(8), intent(inout) :: c1, lhs1, rhs1, c2, lhs2, rhs2, pex(neqp)

      integer(4) j, i, idx(ne)
      real(c_double) :: zwant, zgot
      complex(c_double_complex) ptmp(2, ne), ptmp_old(2, ne)

      fp(1) = f(1)*cdexp(ic*f(2))
      fp(2) = f(3)*cdexp(ic*f(4))

      call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifailp)

      do i = 1, 2
         ptmp_old(i, :) = dcmplx(p(idxre(i, :), 1), p(idxim(i, :), 1))
      end do

      do j = 1, nz - 1
         zwant = zax(j + 1)
         call d02pcf(dpdz, zwant, zgot, pgot, ppgot, pmax, workp, ifailp)

         if (ifailp .ne. 0) then
            write (*, *)
            write (*, 99998) 'exit d02pcf with ifail = ', ifailp, '  and z = ', zgot
            pause
            stop
         end if
99998    format(1x, a, i2, a, d12.5)

         p(:, j + 1) = pgot

         do i = 1, 2
            idx = idxp(i, :)
            ptmp(i, :) = dcmplx(p(idxre(i, :), j + 1), p(idxim(i, :), j + 1))
            dp2dz(idxp(i, :), j) = (cdabs(ptmp(i, :))**2 - cdabs(ptmp_old(i, :))**2)/dz
            mean2(i, j) = sum(dp2dz(idx, j))/ne
            ptmp_old(i, :) = ptmp(i, :)
         end do
      end do

      pex(:) = pgot

      lhs1 = 4*f(1)**2 - 8*r(1)*f(5)*f(1)*dcos(th(1) - f(6) + f(2))
      rhs1 = -icu(1)*(0.5d0*(mean2(1, 1) + mean2(1, nz - 1)) + sum(mean2(1, 2:(nz - 2))))*dz
      c1 = lhs1 - rhs1

      lhs2 = 4*f(3)**2 - 8*r(2)*f(5)*f(3)*dcos(th(2) - f(6) + f(4))
      rhs2 = -icu(2)*(0.5d0*(mean2(2, 1) + mean2(2, nz - 1)) + sum(mean2(2, 2:(nz - 2))))*dz
      c2 = lhs2 - rhs2

      !open (1, file='test.dat')
      !do j = 1, nz - 1
      !   write (1, '(2f17.8)') zax(j), mean2(2, j)
      !end do
      !close (1)
      !stop

   end subroutine zs

   subroutine freq()

      implicit none

      integer(c_int) :: ii, iparf, aiparp(1), itp, j
      real(c_double) t, z, zwant, zgot, &
         x1r, x1i, q31, i1, r1, th1, dcir1, cos1, sin1, &
         x2r, x2i, q32, i2, r2, th2, dcir2, cos2, sin2, q3, &
         f1, f2, f3, phi1, phi2, phi3, a1, a2
      complex(c_double_complex) x1, x2

      do ii = 2, nt

         fp(1) = f(1, ii)*exp(ic*f(2, ii))
         fp(2) = f(3, ii)*exp(ic*f(4, ii))

         call d02pvf(neqp, zstart, p(:, 1), zex, ptol, thres, method, 'usual task', errass, hstart, workp, lenwrk, ifailp)

         do j = 1, nz - 1
            zwant = zax(j + 1)
            call d02pcf(dpdz, zwant, zgot, pgot, ppgot, pmax, workp, ifailp)

            if (ifailp .ne. 0) then
               write (*, *)
               write (*, 99998) 'exit d02pcf with ifail = ', ifailp, '  and z = ', zgot
               pause
               stop
            end if
99998       format(1x, a, i2, a, d12.5)

            p(:, j + 1) = pgot
         end do

         !rparp = 0.0
         !iparp = 0
         !itolp = 0
         !!rtolp = 1.0d-7
         !rtolp = ptol
         !atolp = rtolp
         !ioutp = neqp
         !z = zax(1)
         !xoutp = z
         !itp = 0
         !yp = p(:, 1)
         !iworkp(:) = 0
         !workp(:) = 0.0d0
         !iworkp(5) = neqp
         !
         !artolp(1) = rtolp
         !aatolp(1) = atolp
         !arparp(1) = rparp
         !aiparp(1) = iparp
         !
         !call dopri5_p(neqp, dpdz, z, yp, zex, artolp, aatolp, itolp, soloutp, ioutp, &
         !              workp, lworkp, iworkp, liworkp, arparp, aiparp, ididp)
         !
         !p(:, nz) = yp(:)

         !open (1, file='test.dat')
         !do ii = 1, nz
         !    write (1,'(5f17.8)') zax(ii), p(1,ii), p(ne+1,ii), p(10,ii), p(ne+10,ii)
         !end do
         !close (1)
         !stop

         x1 = xi(p(1:2*ne, :), 1)
         x2 = xi(p(2*ne + 1:4*ne, :), 1)

         x1r = dreal(x1)
         x1i = dimag(x1)
         x2r = dreal(x2)
         x2i = dimag(x2)

         f1 = f(1, ii)
         phi1 = f(2, ii)
         f2 = f(3, ii)
         phi2 = f(4, ii)
         f3 = f(5, ii)
         phi3 = f(6, ii)

         q31 = q(3)/q(1)
         i1 = icu(1)
         r1 = r(1)
         th1 = th(1)
         dcir1 = dcir(1)
         cos1 = dcos(phi1)
         sin1 = dsin(phi1)

         q32 = q(3)/q(2)
         i2 = icu(2)
         r2 = r(2)
         th2 = th(2)
         dcir2 = dcir(2)
         cos2 = dcos(phi2)
         sin2 = dsin(phi2)

         q3 = q(3)
         a1 = a(1)
         a2 = a(2)

         !s(1) = (-nharm*f1 + i1*(-x1i*cos1 + x1r*sin1) + 2*r1*nharm*f3*dcos(phi3 - phi1 - th1))*q31
         w(1, ii) = -2*dcir1*q3 + (i1/f1*(x1r*cos1 + x1i*sin1) + 2*r1*nharm*(f3/f1)*dsin(phi3 - phi1 - th1))*q31

         !s(3) = (-nharm*f2 + i2*(-x2i*cos2 + x2r*sin2) + 2*r2*nharm*f3*dcos(phi3 - phi2 - th2))*q32
         w(2, ii) = -2*dcir2*q3 + (i2/f2*(x2r*cos2 + x2i*sin2) + 2*r2*nharm*(f3/f2)*dsin(phi3 - phi2 - th2))*q32

         !s(5) = -f3 + a1*f1*dcos(phi1 - phi3) + a2*f2*dcos(phi2 - phi3)
         w(3, ii) = a1*f1/f3*dsin(phi1 - phi3) + a2*f2/f3*dsin(phi2 - phi3)
      end do
   end subroutine freq

end module fun
