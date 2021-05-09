      subroutine plusregC (nrow,ncol,nrank,u,v,irow,pcol,nomega,x,zy,zz)
      integer nrow, ncol, nrank,nomega, irow(nomega),pcol(ncol+1)
      double precision u(nrow,nrank),v(ncol,nrank),zy(nrank)
      double precision zz(nrank,nrank),x(nomega)
      integer ii,ni,jstart,jend,i,j,m,n
!HPF$ INDEPENDENT
      do 1 i=1,nrank
         zy(i)=0d0
         do 2 j=1,nrank
            zz(i,j)=0d0
 2       continue
 1       continue
      ni=0
      do 10 j=1,ncol
         jstart=pcol(j)+1
         jend=pcol(j+1)
         if(jstart > jend) goto 10
         do 20 ii =jstart,jend
            ni=ni+1
            i=irow(ii)+1
            do 30 m=1,nrank
               zy(m)=zy(m)+x(ni)*u(i,m)*v(j,m)
               zz(m,m)=zz(m,m)+(u(i,m)*v(j,m))**2
               if(m .eq. nrank) goto 30
               do 40 n=m+1,nrank
                  zz(m,n)=zz(m,n)+u(i,m)*u(i,n)*v(j,m)*v(j,n)
 40               continue
 30               continue
 20   continue      
 10   continue
      return
      end




