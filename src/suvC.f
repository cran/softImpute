      subroutine suvC (nrow,ncol,nrank,u,v,irow,pcol,nomega,r)
      integer nrow, ncol, nrank,nomega, irow(nomega),pcol(ncol+1)
      double precision u(nrow,nrank),v(ncol,nrank),r(nomega)
c      double precision rtemp
      integer ii,ni,jstart,jend,i,j
!HPF$ INDEPENDENT
      ni=0
      do 10 j=1,ncol
         jstart=pcol(j)+1
         jend=pcol(j+1)
         if(jstart > jend) goto 10
         do 20 ii =jstart,jend
            ni=ni+1
            i=irow(ii)+1
            r(ni) = DOT_PRODUCT(u(i,:), v(j,:))
 20   continue      
 10   continue
      return
      end




