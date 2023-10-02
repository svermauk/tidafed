!---------------------------------------------------------------------!
!Written by Shalini Awasthi (ashalini@iitk.ac.in)
!---------------------------------------------------------------------!
      PROGRAM WSMTD_rw_2D
      IMPLICIT NONE
      REAL*8 gridmin1, gridmax1, griddif1, dummy2,dummy3,v,
     &       gridmin2, gridmax2, griddif2, data1,
     &       gridmin3, gridmax3, griddif3,
     &       gridmin4, gridmax4, griddif4,
     &       prob,den,alpha,fes,fes1,grid,Ro,
     &       cv1,cv2,cv3,cv4,dhdl,kt0,kt,ktb,
     &       dum,s1,s2,s3,s4,dummy11,av_dhdl,phi_z,num_phi_z
      ALLOCATABLE cv1(:),cv2(:),cv3(:),cv4(:),prob(:,:,:,:),
     &            fes(:,:,:,:),fes1(:),grid(:)
      ALLOCATABLE dhdl(:),av_dhdl(:,:,:,:)
      INTEGER md_steps,dummy1,i,j,index1,index2,index3,index4,k,
     &        t_min,t_max,nbin1,nbin2,nbin3,nbin4,i_md,i_s1,i_s2,i_s3,
     &        i_s4,narg,w_cv
      LOGICAL pmf,inpgrid
      CHARACTER*120 :: arg 
      REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
      REAL*8, PARAMETER :: au_to_kcal = 627.51 
      REAL*8, PARAMETER :: kj_to_kcal = 0.239006 

!      character(len=50) f1

      OPEN(11,FILE='COLVAR_1',STATUS='unknown')
      OPEN(12,FILE='cv.dat',STATUS='unknown')
      OPEN(13,FILE='ti001.out',STATUS='unknown')
!      OPEN(14,FILE='cvmdck_mtd',STATUS='unknown')
!
      CALL get_steps(11,md_steps)
!

      kt0=300.D0
      kt=300.D0
      t_min=1
      t_max=md_steps
      pmf=.FALSE.
      inpgrid=.false.
      t_max=md_steps
      narg = IARGC()
      DO i=1,narg
        CALL GETARG(i,arg)
        IF(INDEX(arg,'-T0').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)kt0
        ELSEIF(INDEX(arg,'-T').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)kt
        ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_min
        ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_max
           IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
        ELSE IF(INDEX(arg,'-grid').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)gridmin1
            CALL GETARG(i+2,arg)
            READ(arg,*)gridmax1
            CALL GETARG(i+3,arg)
            READ(arg,*)griddif1
            CALL GETARG(i+4,arg)                               
            READ(arg,*)gridmin2                       
            CALL GETARG(i+5,arg)                      
            READ(arg,*)gridmax2                       
            CALL GETARG(i+6,arg)                      
            READ(arg,*)griddif2      
            CALL GETARG(i+7,arg)                 
            READ(arg,*)gridmin3                       
            CALL GETARG(i+8,arg)                      
            READ(arg,*)gridmax3                       
            CALL GETARG(i+9,arg)                      
            READ(arg,*)griddif3                       
            CALL GETARG(i+10,arg)                 
            READ(arg,*)gridmin4                       
            CALL GETARG(i+11,arg)                      
            READ(arg,*)gridmax4                       
            CALL GETARG(i+12,arg)                      
            READ(arg,*)griddif4                       
            
            inpgrid=.true.
        ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)w_cv
        END IF
      END DO

      md_steps=t_max
      WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
      WRITE(*,'(A,I10)')'No: of max MD  steps        =',t_max
      WRITE(*,'(A,I10)')'No: of min MD  steps        =',t_min

      ALLOCATE(cv1(md_steps),cv2(md_steps),cv3(md_steps),cv4(md_steps))
      ALLOCATE(dhdl(md_steps))

      DO i_md=1,md_steps
        READ(11,*)dummy11,cv1(i_md),cv2(i_md),cv3(i_md),cv4(i_md)
        
        READ(13,*)dhdl(i_md)
!        IF( cv1(i_md) .gt.  3.14d0)  cv1(i_md) = cv1(i_md) - 6.28d0
!        IF( cv1(i_md) .lt. -3.14d0 ) cv1(i_md) = cv1(i_md) + 6.28d0
        WRITE(12,*)dummy11,cv1(i_md),cv2(i_md),cv3(i_md),cv4(i_md)
      END DO
      WRITE(*,'(A)')'CV values written in cv.dat'

      nbin1 = NINT((gridmax1-gridmin1)/griddif1)+1
      nbin2 = NINT((gridmax2-gridmin2)/griddif2)+1
      nbin3 = NINT((gridmax3-gridmin3)/griddif3)+1
      nbin4 = NINT((gridmax4-gridmin4)/griddif4)+1
      WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
      WRITE(*,'(A10,3F8.4,I10)')'CV1  COORD:',
     &          gridmin1,gridmax1,griddif1,nbin1
      WRITE(*,'(A10,3F8.4,I10)')'CV2  COORD',
     &          gridmin2,gridmax2,griddif2,nbin2
      WRITE(*,'(A10,3F8.4,I10)')'CV3  COORD',
     &          gridmin3,gridmax3,griddif3,nbin3
      WRITE(*,'(A10,3F8.4,I10)')'CV4  COORD',
     &          gridmin4,gridmax4,griddif4,nbin4

      ALLOCATE(prob(nbin1,nbin2,nbin3,nbin4))
      ALLOCATE(fes(nbin1,nbin2,nbin3,nbin4))
      ALLOCATE(av_dhdl(nbin1,nbin2,nbin3,nbin4))

!calculate prob        
      den=0.d0
      prob=0.d0
      av_dhdl=0.d0
      DO i_md=1,md_steps
      if ((i_md.ge.t_min).and.(i_md.le.t_max)) then
          index1 = nint((cv1(i_md)-gridmin1)/griddif1) +1
          index2 = nint((cv2(i_md)-gridmin2)/griddif2) +1
          index3 = nint((cv3(i_md)-gridmin3)/griddif3) +1
          index4 = nint((cv4(i_md)-gridmin4)/griddif4) +1
!          if(index1.eq.nbin1)index1=1 ! WARNING! this is to avoid -pi and +pi bins to be counted seperately
          if ((gridmin1.le.cv1(i_md)).and.(gridmin2.le.cv2(i_md)).and.
     &    (gridmin3.le.cv3(i_md)).and.(gridmin3.le.cv3(i_md)).and.
     &    (gridmax1.ge.cv1(i_md)).and.(gridmax2.ge.cv2(i_md)).and.
     &    (gridmax3.ge.cv3(i_md)).and.(gridmax4.ge.cv4(i_md))) then
            prob(index1,index2,index3,index4)=prob(index1,index2,index3,
     &                                             index4)+ 1.d0
            av_dhdl(index1,index2,index3,index4)=av_dhdl(index1,index2,
     &                                        index3,index4)+dhdl(i_md)
          end if
      end if
      END DO
      
      DO index1=1,nbin1
        DO index2=1,nbin2
          DO index3=1,nbin3
            DO index4=1,nbin4
              den=den+prob(index1,index2,index3,index4)
            end do
          end do
        end do
      end do   
      
!      dum=dfloat(t_max-t_min+1)!*griddif1
!      write(*,*) den, dum
!      dum=dfloat(den)*griddif1
      OPEN(2,FILE='Pu.dat',STATUS='unknown',FORM='unformatted')
        DO i_s1=1,nbin1
          s1=DFLOAT(i_s1-1)*griddif1+gridmin1
          DO i_s2=1,nbin2
            s2=DFLOAT(i_s2-1)*griddif2+gridmin2
            DO i_s3=1,nbin3
              s3=DFLOAT(i_s3-1)*griddif3+gridmin3
              DO i_s4=1,nbin4
                s4=DFLOAT(i_s4-1)*griddif4+gridmin4
                dum=prob(i_s1,i_s2,i_s3,i_s4) !dum has the number of visits on a given bin
                data1=0.d0
                if (dum > 0.0) data1 = av_dhdl(i_s1,i_s2,i_s3,i_s4)/dum ! Local average <du/dl> in bin
                if(dum+1.eq.dum) dum=1.0d-16  !remove infinity
                phi_z =-(kb*kt)*dlog(dum/(den*griddif1*griddif2*
     &                                        griddif3*griddif4))
                num_phi_z=dexp(-kb*kt0*phi_z)                     ! Numerator of A_\lambda(Z) in our paper
!             WRITE(2,'(2E16.8)')s1,s2,prob(i_s1,i_s2)
                WRITE(2)num_phi_z,data1
              END DO
            END DO
         END DO
      END DO
      WRITE(*,'(A)')'Unbiased distribution written in Pu.dat'

      END PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
      SUBROUTINE get_steps(iunit,nsteps)
      IMPLICIT NONE
      INTEGER iunit, nsteps
      INTEGER ios
      nsteps=0
      REWIND(iunit)
      Read_Loop: DO
         READ(iunit,*,IOSTAT=ios)
         IF(ios.ne.0)EXIT Read_Loop
         nsteps=nsteps+1
      END DO Read_Loop 
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      SUBROUTINE check_files(iunit,dt)
      IMPLICIT NONE
      INTEGER iunit, dt
      INTEGER ii, jj,i,ios
      dt=0
      i=2
      REWIND(iunit)
      READ(iunit,*)ii
      READ(iunit,*)jj
      dt=jj-ii
      ii=jj
      RLoop: DO 
        i=i+1
        READ(iunit,*,IOSTAT=ios)jj
        IF(ios.ne.0)EXIT RLoop
        IF(jj.ne.ii+dt)THEN
           print *, '!!ERROR: Steps are not at constant stride!!'
           print *, '!!       Unit No:',iunit,'!!'
           print *, '!!       Line No:',i,'!!'
           print *, '!! Expected stride =', dt,'!!'
           print *, '!! Actual stride =', jj-ii,'!!'
           STOP
        END IF
        ii=jj
      END DO RLoop
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      ! We are not using this subroutine.
      ! There is no need to modify this.
      SUBROUTINE get_gridmin_max(iunit,gridmin1,gridmax1,griddif1,
     &                           gridmin2,gridmax2,griddif2)
      IMPLICIT NONE 
      INTEGER :: iunit 
      REAL*8  :: gridmin1, gridmax1, griddif1, 
     &           gridmin2, gridmax2, griddif2
      INTEGER :: ii, ios
      REAL*8  :: cv1, cv2
      INTEGER, PARAMETER :: Def_Grid_Size=101
      REWIND(iunit)
      READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
      if(ios.ne.0)stop 'ERROR reading CV.dat'
      gridmin1=cv1
      gridmax1=cv1
      gridmin2=cv2
      gridmax2=cv2
      RLoop: DO 
        READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
        if(ios.ne.0)EXIT RLoop
        gridmin1=MIN(gridmin1,cv1)
        gridmin2=MIN(gridmin2,cv2)
        gridmax1=MAX(gridmax1,cv1)
        gridmax2=MAX(gridmax2,cv2)
      END DO RLoop
      griddif1=(gridmax1-gridmin1)/DFLOAT(Def_Grid_Size)
      griddif2=(gridmax2-gridmin2)/DFLOAT(Def_Grid_Size)
      END
!---------------------------------------------------------------------!

