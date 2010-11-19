! Copyright (C) 2009 Anand Patil
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


      SUBROUTINE mtcimt(m,s,nm,out)
cf2py intent(copy) m
cf2py intent(out) out
cf2py intent(hide) nm
      EXTERNAL DTRSM
! DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
      DOUBLE PRECISION m(nm), s(nm,nm), out
      INTEGER i, nm
      EXTERNAL dtrsm

      ! pm.flib.dtrsm_wrap(s,m_,side='l',transa='n',uplo='l')
      CALL dtrsm('L','L','N','N',nm,1,1.0D0,s,nm,m,nm)
      out = 0.0D0
      do i=1,nm
         out = out + m(i)*m(i)
      end do

      RETURN
      END


      SUBROUTINE impw(mr, cr, ml, vl, mp, cp, out, n)
cf2py intent(out) out
cf2py intent(hide) n
cf2py intent(copy) cr, cp
      DOUBLE PRECISION mr(n), cr(n,n), ml(n), vl(n), mp(n), cp(n,n), out
      DOUBLE PRECISION sp,sr
      INTEGER n, i, info
      EXTERNAL DPOTRF

      call DPOTRF( 'L', n, cr, n, info )
      call DPOTRF( 'L', n, cp, n, info )
      call mtcimt(mp,cp,n,sp)
      call mtcimt(mr,cr,n,sr)
      
      out = sp-sr
            
      do i=1,n
         out = out - ml(i)*ml(i)/vl(i)
         out = out + 2.*(dlog(cp(i,i))-dlog(cr(i,i)))
      end do
      
      out = out * 0.5D0

      RETURN
      END

      SUBROUTINE meanl(ml, vl, mp, vp, n, out)
cf2py intent(out) out
cf2py intent(hide) n
      DOUBLE PRECISION ml(n), vl(n), mp(n), vp(n), out
      DOUBLE PRECISION a, b, one, mpi, vpi, mli, vli
      INTEGER n, i

      one = 1.0D0
      out = 0.0D0

      do i=1,n
         mli = ml(i)
         vli = vl(i)
         mpi = mp(i)
         vpi = vp(i)

         a = one/vpi+one/vli
         b = mpi/vpi+mli/vli
         
         out = out + DLOG(a*vpi)
         out = out-b*b/a+mli*mli/vli+mpi*mpi/vpi
      end do

      out = -0.5D0*out

      RETURN
      END

      SUBROUTINE obsc(mp, cp, ml, vl, i, n)
cf2py intent(inplace) mp, cp
cf2py intent(hide) n
      
      DOUBLE PRECISION mp(n), cp(n,n), cps(n), cpo(n)
      DOUBLE PRECISION dm, ml, vl, d     
      INTEGER n,i,j,k
      
      i=i+1
      
      d = cp(i,i)+vl
      dm = ml-mp(i)
      
      
      do j=1,n
         cpo(j) = cp(i,j)
         cps(j) = cp(i,j)/d
         mp(j) = mp(j)+dm*cps(j)
      end do

      do j=1,n
         do k=1,n
            cp(j,k) = cp(j,k)-cpo(j)*cps(k)
         end do
      end do


      RETURN
      END



      SUBROUTINE linit(l, norms,n, ml, vl)
cf2py intent(out) ml, vl
cf2py intent(hide) n
cf2py intent(inplace) l
      DOUBLE PRECISION l(n), norms(n), mp, vp, ml, vl
      INTEGER n, i
      DOUBLE PRECISION ws, wm, dn, no, no2
      
      wm = l(1)
      ws = 0.0D0
      ml = 0.0D0
      vl = 0.0D0
      mp = 0.0D0
      vp = 0.0D0
      dn = n

      do i=2,n
         wm = max(l(i), wm)
      end do
      do i=1,n
         no = norms(i)
         no2 = no*no
         l(i) = dexp(l(i)-wm)
         ws = ws + l(i)
         ml = ml + l(i)*no
         vl = vl + l(i)*no2
         mp = mp + no
         vp = vp + no2
      end do
      
      ml = ml / ws
      vl = vl / ws - ml*ml

      mp = mp / n
      vp = vp / n - mp*mp
      
      vl = vl*vp/(vp-vl)
      ml = (ml*(vl+vp)-mp*vl)/vp

      RETURN
      END

      SUBROUTINE meshmatch(tempvals,mesh,cache,cacheval,nm,nc,nd)
cf2py intent(inplace) tempvals
cf2py intent(hide) nm,nd,nc      
      INTEGER nm,nc,nd,j,k,l
      DOUBLE PRECISION mesh(nm,nd), cache(nc,nd)
      DOUBLE PRECISION tempvals(nm), cacheval(nc)
      LOGICAL match
      
      do j=1,nm
          if (.TRUE.) then
              do k=1,nc
                  match=.TRUE.
                  do l=1,nd
                      if (DABS(mesh(j,l)-cache(k,l)).GE.1.0D-10) then
                          match=.FALSE.
                          go to 101
                       end if
101                end do
                   if (match) then
                       tempvals(j) = cacheval(k)
                       go to 102
                   end if
102           end do
          else
            continue
          end if
      end do

      RETURN
      END

      SUBROUTINE bufster(arr,ci,barr,ni,nj,nci)
cf2py intent(out) barr
cf2py intent(hide) ni, nj, nci
cf2py intent(in) arr, ci
      INTEGER ci(nci,2)
      LOGICAL arr(ni,nj)
      INTEGER ni,nj,nci
      INTEGER i,j,k,barr(ni,nj), this
      
      do j=1,nj
          do i=1,ni
            this = 0
            do k=1,nci
                if (arr(i+ci(k,1),j+ci(k,2)+1)) then
                    this = 1
                end if
            end do
            barr(i,j) = this
        end do
      end do    
      
      RETURN
      END


      SUBROUTINE subset_eq(m1,m2,start,stopp,n1,n2,nd)
cf2py intent(hide) n1,n2,nd
cf2py intent(out) start,stopp
      INTEGER n1,n2,start,stopp,i,j
      DOUBLE PRECISION m1(n1,nd), m2(n2,nd)
      LOGICAL all

!      If the second is bigger than the first, no match is possible.
      if (n2.GT.n1) then
        start = -1
        stopp = -1
        RETURN

      else if (n2.EQ.n1) then
        do j=1,nd
          do i=1,n1
!           If they're the same size, no mismatch is tolerated.
            if (DABS(m1(i,j)-m2(i,j)).GE.1.0D-10) then
              start = -1
              stopp = -1
              RETURN
            end if
          end do
        end do
        start = 0
        stopp = n2
        RETURN

      else
!       Look for a match between any in m1 and the first row in m2.
        do i=1,(n1-n2+1)
          all=.TRUE.
          do j=1,nd
            if (DABS(m1(i,j)-m2(1,j)).GE.1.0D-10) then
              all=.FALSE.
            end if
          end do
!         If a match is found, check the next n2 rows of m1 against m2.
          if (all) then
!             print*,'First row matches at ',i
            do j=1,nd
              all = .TRUE.
!             Check the next n2 rows. If any mismatch is found, drop out of the loop.
              do k=1,n2-1
                if (DABS(m1(i+k,j)-m2(k+1,j)).GE.1.0D-10) then
                  all = .FALSE.
!                   print *,'Mismatch at ',k
                  goto 1
                end if
              end do
!             If all the next n2 rows match, return.
    1         if (all) then
!                 print *,'All match, returning ',i
                start=i-1
                stopp=i+n2-1
                RETURN
              end if
            end do
          end if
        end do
        start = -1
        stopp = -1
        RETURN
      end if
      
!       print*, 'Control reached EOF'
      
      END

      SUBROUTINE multiinc(x,ind,nxi,nxk)
cf2py intent(hide) nxi,nxk
cf2py intent(inplace) x
cf2py threadsafe
      INTEGER x(nxi,nxk), ind(nxi), nxi
      INTEGER nxk, i, k
      
      do i=1,nxi
          k = max(1, min(ind(i)+1, nxk))
          x(i,k) = x(i,k) + 1
      end do

      RETURN
      END
      
      SUBROUTINE qextract(x,n,q,out,bin,nxi,nxk,nq)
cf2py intent(hide) nxi,nxk,nq
cf2py intent(out) out
cf2py threadsafe
      INTEGER x(nxi,nxk), nxi, i, k, l
      INTEGER nxk, nq, n 
      DOUBLE PRECISION q(nq), bin(nxk), out(nq, nxi)
      DOUBLE PRECISION cusum, next

      do i=1,nxi
          cusum = 0.0D0
          l = 0      
!           print *,i,nxi,cusum,l
!           print *,
          do k=1,nq
              out(k,i) = 0.0D0
              next = q(k)*n
              do while (cusum.LT.next)
!                   print *,l,k,cusum,next
                  l = l + 1
                  cusum = cusum + x(i,l)
              end do
!               print *,k,next,cusum,l,n
!               print *,
              out(k,i) = bin(l)
          end do
      end do

      RETURN
      END
     

      SUBROUTINE iinvlogit(C,nx,cmin,cmax)

cf2py intent(inplace) C
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny
cf2py threadsafe

      DOUBLE PRECISION C(nx)
      INTEGER nx, i, cmin, cmax

      if (cmax.EQ.-1) then
          cmax = nx
      end if


        do i=cmin+1,cmax
            C(i) = 1.0D0 / (1.0D0 + dexp(-C(i)))
        end do


      RETURN
      END   


      SUBROUTINE isinvlogit(C,a1,a2,nx,cmin,cmax)

cf2py intent(inplace) C
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny
cf2py threadsafe

      DOUBLE PRECISION C(nx), a1, a2
      INTEGER nx, i, cmin, cmax

      if (cmax.EQ.-1) then
          cmax = nx
      end if


        do i=cmin+1,cmax
            if (C(i).GT.0.0D0) then
                if (a1.GT.0.0D0) then
                    C(i) = (dexp(a1*C(i))-1.0D0)/a1
                else if (a1.LT.0.0D0) then
                    C(i) = -dlog(1.0D0-a1*C(i))/a1
                end if

            else if (C(i).LT.0.0D0) then
                if (a2.GT.0.0D0) then
                    C(i) = -(dexp(-a2*C(i))-1.0D0)/a2
                else if (a2.LT.0.0D0) then
                    C(i) = dlog(1.0D0+a2*C(i))/a2
                end if

            end if
        end do
        
        CALL iinvlogit(C,nx,cmin,cmax)


      RETURN
      END   


      SUBROUTINE iaaxpy(a,x,y,n,cmin,cmax)

cf2py intent(inplace) y
cf2py intent(in) x,a
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) n
cf2py threadsafe

      DOUBLE PRECISION y(n)
      DOUBLE PRECISION x(n),a
      INTEGER n, cmin, cmax, i
!      EXTERNAL DAXPY

      if (cmax.EQ.-1) then
          cmax = n
      end if

      do i=cmin+1,cmax
          y(i)=a*x(i)+y(i)
      end do
      !CALL DAXPY(cmax-cmin,a,x(cmin+1),1,y(cmin+1),1)


      RETURN
      END


      SUBROUTINE icsum(C,x,d,y,nx,ny,nd,cmin,cmax)

cf2py intent(inplace) C
cf2py intent(in) x,d,y
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny,nd
cf2py threadsafe

      DOUBLE PRECISION C(nx,ny)
      DOUBLE PRECISION x(nd,nx),d(nd),y(nd,ny)
      INTEGER nx, ny, nd, i, j, k, cmin, cmax

      EXTERNAL DSCAL

      if (cmax.EQ.-1) then
          cmax = ny
      end if


        do j=cmin+1,cmax
            do i=1,nx
                do k=1,nd
                    C(i,j) = C(i,j) + x(k,i)*d(k)*y(k,j)
                end do
            end do
        enddo



      RETURN
      END


      SUBROUTINE iasadd(C,a,nx,ny,cmin,cmax)
cf2py intent(inplace) C
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny
cf2py threadsafe

      DOUBLE PRECISION C(nx,ny)
      DOUBLE PRECISION a
      INTEGER nx, ny, i, j, cmin, cmax
 
      EXTERNAL DSCAL
 
      if (cmax.EQ.-1) then
          cmax = ny
      end if
 
 
        do j=cmin+1,cmax
            do i=1,nx
                C(i,j) = C(i,j) + a
            end do
!          CALL DSCAL(nx,a,C(1,j),1)
        enddo
 
 
 
      RETURN
      END
 

      SUBROUTINE iaadd(C,A,nx,ny,cmin,cmax)

cf2py intent(inplace) C
cf2py intent(in) A
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny
cf2py threadsafe

      DOUBLE PRECISION C(nx,ny)
      DOUBLE PRECISION A
      INTEGER nx, ny, i, j, cmin, cmax

      EXTERNAL DSCAL

      if (cmax.EQ.-1) then
          cmax = ny
      end if


        do j=cmin+1,cmax
            do i=1,nx
                C(i,j) = C(i,j) + A(i,j)
            end do
 !          CALL DSCAL(nx,a,C(1,j),1)
        enddo



      RETURN
      END
      

      SUBROUTINE iasq(C,nx,ny,cmin,cmax)

cf2py intent(inplace) C
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny
cf2py threadsafe

      DOUBLE PRECISION C(nx,ny), cn
      INTEGER nx, ny, i, j, cmin, cmax

      EXTERNAL DSCAL

      if (cmax.EQ.-1) then
          cmax = ny
      end if


        do j=cmin+1,cmax
            do i=1,nx
                cn = C(i,j)
                C(i,j) = cn * cn
            end do
 !          CALL DSCAL(nx,a,C(1,j),1)
        enddo


      RETURN
      END

      SUBROUTINE iamul(C,A,nx,ny,cmin,cmax)

cf2py intent(inplace) C
cf2py intent(in) A
cf2py integer intent(in), optional :: cmin = 0
cf2py integer intent(in), optional :: cmax = -1
cf2py intent(hide) nx,ny
cf2py threadsafe

      DOUBLE PRECISION C(nx,ny)
      DOUBLE PRECISION A(nx,ny)
      INTEGER nx, ny, i, j, cmin, cmax

      EXTERNAL DSCAL

      if (cmax.EQ.-1) then
          cmax = ny
      end if


        do j=cmin+1,cmax
            do i=1,nx
                C(i,j) = C(i,j) * A(i,j)
            end do
 !          CALL DSCAL(nx,a,C(1,j),1)
        enddo


      RETURN
      END
