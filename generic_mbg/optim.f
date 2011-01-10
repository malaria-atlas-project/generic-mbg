      subroutine splev(t,n,c,k,x,y,m,ier)

cf2py intent(hide) n, m
cf2py intent(out) y, ier 
      
c  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
c  a spline s(x) of degree k, given in its b-spline representation.
c
c  calling sequence:
c     call splev(t,n,c,k,x,y,m,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    x    : array,length m, which contains the points where s(x) must
c           be evaluated.
c    m    : integer, giving the number of points where s(x) must be
c           evaluated.
c
c  output parameter:
c    y    : array,length m, giving the value of s(x) at the different
c           points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl.
c
c  references :
c    de boor c  : on calculating with b-splines, j. approximation theory
c                 6 (1972) 50-62.
c    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
c                 applics 10 (1972) 134-149.
c    dierckx p. : curve and surface fitting with splines, monographs on
c                 numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,m,ier
c  ..array arguments..
      real t(n),c(n),x(m),y(m)
c  ..local scalars..
      integer i,j,k1,l,ll,l1,nk1
      real arg,sp,tb,te
c  ..local array..
      real h(6)
c  ..
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(m-1) 100,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 100
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
      l = k1
      l1 = l+1
c  main loop for the different points.
      do 80 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
  40    if(arg.lt.t(l1) .or. l.eq.nk1) go to 50
        l = l1
        l1 = l+1
        go to 40
c  evaluate the non-zero b-splines at arg.
  50    call fpbspl(t,n,k,arg,l,h)
c  find the value of s(x) at x=arg.
        sp = 0.
        ll = l-k1
        do 60 j=1,k1
          ll = ll+1
          sp = sp+c(ll)*h(j)
  60    continue
        y(i) = sp
  80  continue
 100  return
      end
      
      
      subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
c  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
c  the derivative of order nu of a spline s(x) of degree k,given in
c  its b-spline representation.
c
c  calling sequence:
c     call splder(t,n,c,k,nu,x,y,m,wrk,ier)
c
c  input parameters:
c    t    : array,length n, which contains the position of the knots.
c    n    : integer, giving the total number of knots of s(x).
c    c    : array,length n, which contains the b-spline coefficients.
c    k    : integer, giving the degree of s(x).
c    nu   : integer, specifying the order of the derivative. 0<=nu<=k
c    x    : array,length m, which contains the points where the deriv-
c           ative of s(x) must be evaluated.
c    m    : integer, giving the number of points where the derivative
c           of s(x) must be evaluated
c    wrk  : real array of dimension n. used as working space.
c
c  output parameters:
c    y    : array,length m, giving the value of the derivative of s(x)
c           at the different points.
c    ier  : error flag
c      ier = 0 : normal return
c      ier =10 : invalid input data (see restrictions)
c
c  restrictions:
c    0 <= nu <= k
c    m >= 1
c    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
c
c  other subroutines required: fpbspl
c
c  references :
c    de boor c : on calculating with b-splines, j. approximation theory
c                6 (1972) 50-62.
c    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
c                applics 10 (1972) 134-149.
c   dierckx p. : curve and surface fitting with splines, monographs on
c                numerical analysis, oxford university press, 1993.
c
c  author :
c    p.dierckx
c    dept. computer science, k.u.leuven
c    celestijnenlaan 200a, b-3001 heverlee, belgium.
c    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
c
c  latest update : march 1987
c
c  ..scalar arguments..
      integer n,k,nu,m,ier
c  ..array arguments..
      real t(n),c(n),x(m),y(m),wrk(n)
c  ..local scalars..
      integer i,j,kk,k1,k2,l,ll,l1,l2,nk1,nk2,nn
      real ak,arg,fac,sp,tb,te
c  ..local arrays ..
      real h(6)
c  before starting computations a data check is made. if the input data
c  are invalid control is immediately repassed to the calling program.
      ier = 10
      if(nu.lt.0 .or. nu.gt.k) go to 200
      if(m-1) 200,30,10
  10  do 20 i=2,m
        if(x(i).lt.x(i-1)) go to 200
  20  continue
  30  ier = 0
c  fetch tb and te, the boundaries of the approximation interval.
      k1 = k+1
      nk1 = n-k1
      tb = t(k1)
      te = t(nk1+1)
c  the derivative of order nu of a spline of degree k is a spline of
c  degree k-nu,the b-spline coefficients wrk(i) of which can be found
c  using the recurrence scheme of de boor.
      l = 1
      kk = k
      nn = n
      do 40 i=1,nk1
         wrk(i) = c(i)
  40  continue
      if(nu.eq.0) go to 100
      nk2 = nk1
      do 60 j=1,nu
         ak = kk
         nk2 = nk2-1
         l1 = l
         do 50 i=1,nk2
            l1 = l1+1
            l2 = l1+kk
            fac = t(l2)-t(l1)
            if(fac.le.0.) go to 50
            wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
  50     continue
         l = l+1
         kk = kk-1
  60  continue
      if(kk.ne.0) go to 100
c  if nu=k the derivative is a piecewise constant function
      j = 1
      do 90 i=1,m
         arg = x(i)
  70     if(arg.lt.t(l+1) .or. l.eq.nk1) go to 80
         l = l+1
         j = j+1
         go to 70
  80     y(i) = wrk(j)
  90  continue
      go to 200
 100  l = k1
      l1 = l+1
      k2 = k1-nu
c  main loop for the different points.
      do 180 i=1,m
c  fetch a new x-value arg.
        arg = x(i)
        if(arg.lt.tb) arg = tb
        if(arg.gt.te) arg = te
c  search for knot interval t(l) <= arg < t(l+1)
 140    if(arg.lt.t(l1) .or. l.eq.nk1) go to 150
        l = l1
        l1 = l+1
        go to 140
c  evaluate the non-zero b-splines of degree k-nu at arg.
 150    call fpbspl(t,n,kk,arg,l,h)
c  find the value of the derivative at x=arg.
        sp = 0.
        ll = l-k1
        do 160 j=1,k2
          ll = ll+1
          sp = sp+wrk(ll)*h(j)
 160    continue
        y(i) = sp
 180  continue
 200  return
      end
      
      subroutine lmdif(t,c,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
     *                 factor,nprint,info,nfev,fjac,ldfjac,
     *                 ipvt,qtf,wa1,wa2,wa3,wa4,norms,nt,nn)
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
      integer ipvt(n)
      double precision ftol,xtol,gtol,epsfcn,factor
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),
     *                 wa1(n),wa2(n),wa3(n),wa4(m)
      double precision t(nt), c(nt), norms(nn)
      integer nt, nn
c     **********
c
c     subroutine lmdif
c
c     the purpose of lmdif is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the levenberg-marquardt algorithm. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
c                        factor,nprint,info,nfev,fjac,
c                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
c
c     where
c
c     Modified 17 May, 2010 by A. Patil to take the t, c, and norms arguments.
c     t and c are the locations and coefficients of spline knots.
c     norms is a vector of standard normal random variables.
c 
cf2py intent(hide) nt, nn
c
c       m is a positive integer input variable set to the number
c         of functions.
c
cf2py intent(hide) m
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
cf2py intent(hide) n
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
cf2py intent(inplace) x
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
cf2py intent(out) fvec
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
cf2py intent(optional) ftol=0
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
cf2py intent(optional) xtol=0
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
cf2py intent(optional) gtol=1
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn is at least
c         maxfev by the end of an iteration.
c
cf2py intent(optional) maxfev=1000
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
cf2py intent(in) epsfcn
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.). 100. is a generally recommended value.
c
cf2py intent(optional) factor=100
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x and fvec available
c         for printing. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
cf2py intent(optional) nprint=0
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn has reached or
c                   exceeded maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
cf2py intent(out) input
c
c       nfev is an integer output variable set to the number of
c         calls to fcn.
c
cf2py intent(out) nfev
c
c       fjac is an output m by n array. the upper n by n submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
cf2py intent(out) fjac
c
c       ldfjac is a positive integer input variable not less than m
c         which specifies the leading dimension of the array fjac.
c
cf2py intent(out) ldfjac
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
cf2py intent(out) ipvt
c
c       qtf is an output array of length n which contains
c         the first n elements of the vector (q transpose)*fvec.
c
cf2py intent(out) qtf
c
c       wa1, wa2, and wa3 are work arrays of length n.
c
c       wa4 is a work array of length m.
c
cf2py intent(hide) wa1,wa2,wa3,wa4
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
c
c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,
     *                 one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,
     *                 sum,temp,temp1,temp2,xnorm,zero
      double precision dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
c
c     epsmch is the machine precision.
c
      mode = 1

      epsmch = dpmpar(1)
c
      info = 0
      iflag = 0
      nfev = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m
     *    .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero
     *    .or. maxfev .le. 0 .or. factor .le. zero) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
      iflag = 1
      call fcn(m,n,x,fvec,iflag)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(m,fvec)
c
c     initialize levenberg-marquardt parameter and iteration counter.
c
      par = zero
      iter = 1
c
c     beginning of the outer loop.
c
   30 continue
c
c        calculate the jacobian matrix.
c
         iflag = 2
         call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
         nfev = nfev + n
         if (iflag .lt. 0) go to 300
c
c        if requested, call fcn to enable printing of iterates.
c
         if (nprint .le. 0) go to 40
         iflag = 0
         if (mod(iter-1,nprint) .eq. 0) call fcn(m,n,x,fvec,iflag)
         if (iflag .lt. 0) go to 300
   40    continue
c
c        compute the qr factorization of the jacobian.
c
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 80
         if (mode .eq. 2) go to 60
         do 50 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   50       continue
   60    continue
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
         do 70 j = 1, n
            wa3(j) = diag(j)*x(j)
   70       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   80    continue
c
c        form (q transpose)*fvec and store the first n components in
c        qtf.
c
         do 90 i = 1, m
            wa4(i) = fvec(i)
   90       continue
         do 130 j = 1, n
            if (fjac(j,j) .eq. zero) go to 120
            sum = zero
            do 100 i = j, m
               sum = sum + fjac(i,j)*wa4(i)
  100          continue
            temp = -sum/fjac(j,j)
            do 110 i = j, m
               wa4(i) = wa4(i) + fjac(i,j)*temp
  110          continue
  120       continue
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
  130       continue
c
c        compute the norm of the scaled gradient.
c
         gnorm = zero
         if (fnorm .eq. zero) go to 170
         do 160 j = 1, n
            l = ipvt(j)
            if (wa2(l) .eq. zero) go to 150
            sum = zero
            do 140 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  140          continue
            gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
  150       continue
  160       continue
  170    continue
c
c        test for convergence of the gradient norm.
c
         if (gnorm .le. gtol) info = 4
         if (info .ne. 0) go to 300
c
c        rescale if necessary.
c
         if (mode .eq. 2) go to 190
         do 180 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  180       continue
  190    continue
c
c        beginning of the inner loop.
c
  200    continue
c
c           determine the levenberg-marquardt parameter.
c
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,
     *                 wa3,wa4)
c
c           store the direction p and x + p. calculate the norm of p.
c
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = enorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
c
c           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
            call fcn(m,n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(m,wa4)
c
c           compute the scaled actual reduction.
c
            actred = -one
            if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction and
c           the scaled directional derivative.
c
            do 230 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 220 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
            ratio = zero
            if (prered .ne. zero) ratio = actred/prered
c
c           update the step bound.
c
            if (ratio .gt. p25) go to 240
               if (actred .ge. zero) temp = p5
               if (actred .lt. zero)
     *            temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par .ne. zero .and. ratio .lt. p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
c
c           test for successful iteration.
c
            if (ratio .lt. p0001) go to 290
c
c           successful iteration. update x, fvec, and their norms.
c
            do 270 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  270          continue
            do 280 i = 1, m
               fvec(i) = wa4(i)
  280          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  290       continue
c
c           tests for convergence.
c
            if (dabs(actred) .le. ftol .and. prered .le. ftol
     *          .and. p5*ratio .le. one) info = 1
            if (delta .le. xtol*xnorm) info = 2
            if (dabs(actred) .le. ftol .and. prered .le. ftol
     *          .and. p5*ratio .le. one .and. info .eq. 2) info = 3
            if (info .ne. 0) go to 300
c
c           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 5
            if (dabs(actred) .le. epsmch .and. prered .le. epsmch
     *          .and. p5*ratio .le. one) info = 6
            if (delta .le. epsmch*xnorm) info = 7
            if (gnorm .le. epsmch) info = 8
            if (info .ne. 0) go to 300
c
c           end of the inner loop. repeat if iteration unsuccessful.
c
            if (ratio .lt. p0001) go to 200
c
c        end of the outer loop.
c
         go to 30
  300 continue
c
c     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(m,n,x,fvec,iflag)
      return
c
c     last card of subroutine lmdif.
c
      end
      

      SUBROUTINE postparams(mu,C,nmu,mul,Vl,ind,nind,d)
c Finds the approximate likelihood mean and variance
c at indices 'ind'.
cf2py intent(hide) nmu, nt, nind
cf2py intent(inplace) mul, Cl
cf2py threadsafe
      INTEGER nmu, nt, i, j, nind, ind(nind), d, k, 
      DOUBLE PRECISION mu(nmu), C(nmu,nmu), mul(nmu), Vl(nmu)
      DOUBLE PRECISION t(nt,nmu), k(nt,nmu)
      
      do k=1,nind
          
      end do
      
      RETURN
      END
      
!       SUBROUTINE postparams_s(mu,C,t,k,nmu,nt,mul,Vl,ind,nind,d)
! c Finds the approximate likelihood mean and variance
! c at indices 'ind'.
! cf2py intent(hide) nmu, nt, nind
! cf2py intent(inplace) mul, Cl
! cf2py threadsafe
!       INTEGER nmu, nt, i, j, nind, ind(nind), d, k, 
!       DOUBLE PRECISION mu(nmu), C(nmu,nmu), mul(nmu), Vl(nmu)
!       DOUBLE PRECISION t(nt,nmu), k(nt,nmu)
!       
!       do k=1,nind
!           
!       
!       RETURN
!       END