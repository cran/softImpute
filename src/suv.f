      subroutine suv (nrow,ncol,nrank,u,v,irow,jcol,nomega,r)
      integer nrow, ncol, nrank,nomega, irow(nomega),jcol(nomega)
      double precision u(nrow,nrank),v(ncol,nrank),r(nomega)
c      double precision rtemp
c      computes uv'[i,j]=<u[i,],v[j,]>
      integer ii,jj
!HPF$ INDEPENDENT
      do 10 i=1,nomega
         ii=irow(i)
         jj=jcol(i)
      r(i) = DOT_PRODUCT(u(ii,:), v(jj,:))
c          r(i)=rtemp
10    continue
      return
      end




c     do 10 i=1,nomega
c         rtemp=0d0
c         ii=irow(i)
C         jj=jcol(i)
c         do 20 k=1,nrank
c            rtemp=rtemp+u(ii,k)*v(jj,k)
c 20          continue
c          r(i)=rtemp
c 10       continue
c         return
c          end
