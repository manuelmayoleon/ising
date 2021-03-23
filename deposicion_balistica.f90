!deposicion_balistica.f90
PROGRAM deposicion_balistica
    implicit none
REAL(kind=8),ALLOCATABLE,DIMENSION(:):: h,w,hmedia,hmedia2
INTEGER,ALLOCATABLE,DIMENSION(:)::ti
INTEGER::l,tmax, nm , i_dran,iseed,n,m,run,nruns,i,t,tini,nmlim
REAL(kind=8)::dl,alfa
l=600
dl=REAL(l)
tmax=100000
nruns=100
tini=8
nmlim=INT(dlog(dble(tmax/tini))/dlog(dble((1+tini))/dble(tini)))+1
nm=25
alfa=exp(dlog(dble(tmax)/dble(tini))/(nm-1))
iseed=1121
print*,'nm:',nm,alfa 
call dran_ini(iseed)

allocate(h(0:l+1),w(nm),hmedia(nm),hmedia2(nm),ti(nm))

w=0.d0


print*, size(h)
IF(nm>nmlim) THEN
 PRINT*, ' Has cogido mal el nmedidas'
ELSE  
    DO run=1,nruns


        h=0.d0
        i=1 !primera medida
        ti(i)=tini !tiempo de primera medida


        DO t=1,tmax

            

            DO n=1,l
            m=i_dran(l)
            !print*, 'random number', m 
            h(m)=max(h(m-1),h(m)+1,h(m+1))
            !print*,h(m)
            h(0)=h(l)
            h(l+1)=h(1)
            END DO

            IF (t == ti(i)) THEN 
            hmedia(i)=0.d0
            hmedia2(i)=0.d0
                !print*, 'suma',sum(h)
            DO m=1,l
            hmedia(i)= hmedia(i)+h(m)    
            hmedia2(i)=hmedia2(i)+h(m)**2
            !print*, hmedia(i)
            END DO
            hmedia(i)= hmedia(i)/dl    
            hmedia2(i)=hmedia2(i)/dl
            !print*, hmedia2(i)

            w(i)=w(i)+dsqrt(hmedia2(i)-hmedia(i)**2)
            i=i+1
            ti(i)=INT(alfa**(i-1)*tini)
            !print*,  ti(i)

            ENDIF


        END DO 
        OPEN(9,FILE='db.txt',STATUS='unknown')
        DO n=1,nm
            WRITE(9,*) ti(n), hmedia(n)/dble(nruns),hmedia2(n)/dble(nruns),w(n)/dble(nruns)
        END DO  
        CLOSE(9)



    END DO

END IF

END PROGRAM deposicion_balistica