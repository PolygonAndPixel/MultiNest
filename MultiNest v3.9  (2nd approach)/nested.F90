! Do nested sampling algorithm to calculate Bayesian evidence
! Aug 2013
! Farhan Feroz

module Nested
  use utils1
  use kmeans_clstr
  use xmeans_clstr
  use posterior
  use priors
  use omp_lib ! For OpenMP
  implicit none

  integer maxCls,maxeCls
  integer min_pt,totPar
  integer nCdims !total no. of parameters on which clustering should be done
  integer nlive ! Number of live points
  integer ndims ! Number of dimensions
  integer nsc_def !no. of iterations per every sub-clustering step
  integer updInt !update interval
  double precision Ztol !lowest local evidence for which samples to produce
  double precision tol ! tolerance at end
  double precision ef
  integer numlike,globff
  double precision logZero
  integer maxIter
  integer minIter ! Minimal number of multinest iterations required. Set to
                  ! negative value to have same behaviour as vanilla MultiNest
  logical fback,resumeFlag,dlive,genLive,dino
  !output files name
  character(LEN=100)physname,broot,rname,resumename,livename,evname
  !output file units
  integer u_ev,u_resume,u_phys,u_live
  double precision gZ,ginfo !total log(evidence) & info
  integer sCount
  logical, dimension(:), allocatable :: pWrap
  logical mWrap,aWrap !whether to do wraparound for mode separation
  logical debug, prior_warning, resume, outfile
  integer totalIterations, totalLikelihoodCalls, procs
 
contains
  
subroutine nestRun(nest_IS,nest_mmodal,nest_ceff,nest_nlive,nest_tol,nest_ef,nest_ndims,nest_totPar,nest_nCdims,maxClst, &
  nest_updInt,nest_Ztol,nest_root,seed,nest_pWrap,nest_fb,nest_resume,nest_outfile,&
  initMPI,nest_logZero,nest_maxIter, nest_minIter, &
  loglike,dumper,context, nest_totalLikelihoodCalls, nest_totalIterations, nest_procs)
        
    implicit none
        
    integer nest_ndims,nest_nlive,nest_updInt,context,seed,i,&
        nest_totalLikelihoodCalls,nest_totalIterations
    integer maxClst,nest_nsc,nest_totPar,nest_nCdims,nest_pWrap(*),&
        nest_maxIter, nest_minIter,nest_procs
    logical nest_IS,nest_mmodal,nest_fb,nest_resume,nest_ceff,nest_outfile,initMPI
    character(LEN=100) nest_root
    double precision nest_tol,nest_ef,nest_Ztol,nest_logZero

    INTERFACE
            !the likelihood function
        subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
            integer n_dim,nPar,context_pass
            double precision lnew,Cube(nPar)
        end subroutine loglike
    end INTERFACE

    INTERFACE
        !the user dumper function
        subroutine dumper(nSamples, nlive, nPar, physLive, posterior, &
            paramConstr, maxLogLike, logZ, logZerr, context_pass)
            integer nSamples, nlive, nPar, context_pass
            double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
            double precision maxLogLike, logZ, logZerr
        end subroutine dumper
    end INTERFACE

    procs = nest_procs
    nest_totalLikelihoodCalls = 0
    nest_totalIterations = 0
    totalIterations = 0
    totalLikelihoodCalls = 0
    nest_nsc=50
    nlive=nest_nlive
    Ztol=nest_Ztol
    updInt=nest_updInt
    logZero=nest_logZero
    maxIter=nest_maxIter
    minIter=nest_minIter
    if(maxIter<=0) maxIter=huge(1)

    ndims=nest_ndims
    totPar=nest_totPar
    nCdims=nest_nCdims
    debug=.false.
    prior_warning=.true.
    resume=nest_resume
    outfile=nest_outfile
    if(.not.outfile) resume=.false.
    
    if(nCdims>ndims) then
        write(*,*)"ERROR: nCdims can not be greater than ndims."
        write(*,*)"Aborting"
        stop
    endif
	
    sCount=0
    
    broot=nest_root
    rname = trim(broot)
    fback = nest_fb

    !output file info
    !setup the output files
    resumename = trim(rname)//'resume.dat'
    physname = trim(rname)//'phys_live.points'
    livename = trim(rname)//'live.points'
    evname = trim(rname)//'ev.dat'
    u_ev=55
    u_phys=57
    u_live=59
    u_resume=61
    
    allocate(pWrap(ndims))
    mWrap=.false.
    aWrap=.true.
    do i=1,ndims
        if(nest_pWrap(i)==0) then
            pWrap(i)=.false.
            if(i<=nCdims) aWrap=.false.
        else
            pWrap(i)=.true.
            if(i<=nCdims) mWrap=.true.
        endif
    enddo
    tol=nest_tol
    ef=nest_ef
    if(ef>=1d6) then
        dino=.false.
        ef=1d0
        min_pt=ndims+1
    else
        dino=.true.
        min_pt=2
    endif
    nsc_def=nest_nsc

    maxCls=maxClst*2
      
    maxeCls=nlive/min_pt

    if(seed<0) then
        !take the seed from system clock
        call InitRandomNS(procs)
    else
        call InitRandomNS(procs,seed)
    endif

    !set the resume flag to true if the resume file exists else to false
    if(resume) then
        inquire(file=resumename,exist=resumeFlag)
        if(.not.resumeFlag) write(*,*)"MultiNest Warning: no resume file found, starting from scratch"
    else
        resumeFlag=.false.
    endif
    
    write(*,*)"*****************************************************"
    write(*,*)"MultiNest v3.1"
    write(*,*)"Copyright Farhan Feroz & Mike Hobson"
    write(*,*)"Release Aug 2013"
    write(*,*)'Modified by Maicon Hieronymus'
    write(*,*)
    write(*,'(a,i4)')" no. of live points = ",nest_nlive
    write(*,'(a,i4)')" dimensionality = ",nest_ndims
    if(resumeFlag) write(*,'(a)')" resuming from previous job"
    write(*,*)"*****************************************************"
    write(*,*) 'WARNING: This version version does not support constant'
    write(*,*) '         efficiency mode and importance sampling.'
    write(*,*) 'This is an optimized version for multimodal regions.'
    write(*,*) 'Number of threads: ', procs
    write(*,*) '2nd parallelization approach'
    if (fback) write (*,*) 'Starting MultiNest'

    !create the output files
    if(.not.resumeFlag .and. outfile) then
        open(unit=u_ev,file=evname,status='replace')
        close(u_ev)
    endif

    call Nestsample(loglike, dumper, context)
    deallocate(pWrap)
    call killRandomNS()

    nest_totalIterations = totalIterations
    nest_totalLikelihoodCalls = totalLikelihoodCalls

  end subroutine nestRun

!----------------------------------------------------------------------

  subroutine Nestsample(loglike, dumper, context)

    implicit none

    integer context
    double precision, allocatable :: p(:,:), phyP(:,:) !live points
    double precision, allocatable :: l(:) !log-likelihood
    double precision vnow1!current vol
    double precision ltmp(totPar+2)
    character(len=100) fmt
    integer np,i,j,k,ios


    INTERFACE
        !the likelihood function
        subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
            integer n_dim,nPar,context_pass
            double precision lnew,Cube(nPar)
        end subroutine loglike
    end INTERFACE
    
    INTERFACE
        !the user dumper function
        subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
            integer nSamples, nlive, nPar, context_pass
            double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
            double precision maxLogLike, logZ, logZerr
        end subroutine dumper
    end INTERFACE


    allocate( p(ndims,nlive+1), phyP(totPar,nlive+1), l(nlive+1) )

    np=ndims
    globff=0
    numlike=0
    vnow1=1.d0

    write(fmt,'(a,i5,a)')  '(',np+1,'E28.18)'

    genLive=.true.

    if(resumeflag) then
        !check if the last job was aborted during the live points generation
        open(unit=u_resume,file=resumename,status='old')
        read(u_resume,*)genLive
        
        if( .not.genLive ) then
            read(u_resume,*)i,j,j,j
    
            if( j /= nlive ) then
                write(*,*)"ERROR: no. of live points in the resume file is not equal to the the no. passed to nestRun."
                write(*,*)"Aborting"
                stop
            endif
        endif
        close(u_resume)
    
        if( .not.genLive ) then
            j = 0
            open(unit=u_ev,file=evname,status='old') 
            write(fmt,'(a,i2.2,a)')  '(',totPar+2,'E28.18,i3)'
            do
                read(55,*,IOSTAT=ios) ltmp(1:totPar+2),k
                !end of file?
                if(ios<0) exit
                j = j + 1
            enddo
            close(u_ev)
            if( j + nlive /= i ) then
                write(*,*)"ERROR: no. of points in ev.dat file is not equal to the no. specified in resume file."
                write(*,*)"Aborting"
                stop
            endif
        endif
    
    endif

    gZ=logZero
    ginfo=0.d0

    if(genLive) then
        if(fback) write(*,*) 'generating live points'
        
        call gen_initial_live(p,phyP,l,loglike,dumper,context)

        globff=nlive
        numlike=nlive
        if(fback) write(*,*) 'live points generated, starting sampling'
    endif

    call clusteredNest(p,phyP,l,loglike,dumper,context)

    write(*,*)"ln(ev)=",gZ,"+/-",sqrt(ginfo/dble(nlive))
    write(*,'(a,i12)')' Total Likelihood Evaluations: ', numlike
    write(*,*)"Sampling finished. Exiting MultiNest"
    setBlk=.false.

    deallocate( p, phyP, l )

  end subroutine Nestsample

!----------------------------------------------------------------------

  subroutine gen_initial_live(p,phyP,l,loglike,dumper,context)
    
    implicit none
    
    integer i,j,iostatus,idum,k,m,nptPerProc,nGen,nstart,nend,context
    double precision, allocatable :: pnewP(:,:), phyPnewP(:,:), lnewP(:)
    double precision p(ndims,nlive+1), phyP(totPar,nlive+1), l(nlive+1)
    integer id
    integer nLikelihoodCalls
    character(len=100) fmt,fmt2
    
    INTERFACE
        !the likelihood function
        subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
            integer n_dim,nPar,context_pass
            double precision lnew,Cube(nPar)
        end subroutine loglike
    end INTERFACE

    INTERFACE
        !the user dumper function
        subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
            integer nSamples, nlive, nPar, context_pass
            double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
            double precision maxLogLike, logZ, logZerr
        end subroutine dumper
    end INTERFACE

    allocate( pnewP(ndims,10), phyPnewP(totPar,10), lnewP(10) )

    if(outfile) then
            open(unit=u_resume,file=resumename,form='formatted',status='replace')
            write(u_resume,'(l2)')genLive
            close(u_resume)
            write(fmt,'(a,i5,a)')  '(',ndims+1,'E28.18)'
            write(fmt2,'(a,i5,a)')  '(',totPar+1,'E28.18,i4)'
    endif
    i=0

    !resume from previous live points generation?
    if(outfile) then
        if(resumeflag) then
            !read hypercube-live file
            open(unit=u_live,file=livename,status='old')
            do
                i=i+1
                read(u_live,*,IOSTAT=iostatus) p(:,i),l(i)
                if(iostatus<0) then
                    i=i-1
                    if(i>nlive) then
                        write(*,*)"ERROR: more than ",nlive," points in the live points file."
                        write(*,*)"Aborting"
                        stop
                    endif
                    exit
                endif
            enddo
            close(u_live)

            if(i>0) then
                !read physical live file
                open(unit=u_phys,file=physname,status='old')
                do j=1,i
                    read(u_phys,*) phyP(:,j),l(j),idum
                enddo
                close(u_phys)
            endif
                open(unit=u_live,file=livename,form='formatted',status='old',position='append')
                open(unit=u_phys,file=physname,form='formatted',status='old',position='append')
            else
                open(unit=u_live,file=livename,form='formatted',status='replace')
                open(unit=u_phys,file=physname,form='formatted',status='replace')
            
        endif
    endif

    j=i
    nend=i
    
    nGen = nlive - j
    nptPerProc = ceiling( dble(nGen) / dble(1) )
      
    if(nGen==0) then
        genLive=.false.
        resumeFlag=.false.
    
        deallocate( pnewP, phyPnewP, lnewP )

        if( outfile ) then
            close(u_live)
            close(u_phys)
        endif
        return
    
    elseif(nGen<0) then
        write(*,*)"ERROR: live points files have more live points than required."
        write(*,*)"Aborting"
        stop
    endif

    k=0
    j=0
    id=0
    do
        k=k+1
        j=j+1
        nLikelihoodCalls=0
        do
            call getrandom(ndims,pnewP(:,j),id)  ! start points
            phyPnewP(1:ndims,j)=pnewP(1:ndims,j)
            lnewP(j)=logZero
            call loglike(phyPnewP(:,j),ndims,totPar,lnewP(j),context)
            nLikelihoodCalls=nLikelihoodCalls+1
            if(lnewP(j)>logZero) exit
            if(nLikelihoodCalls==100) then
                write(*,*) "During initialization a large number of likelihood calls &
                            &is being made."
                write(*,*) "Please verify your configuration file. &
                            &This message will only be printed once."
                write(*,*) "     Last calculated likelihood       :",lnewP(j)
                write(*,*) "     Limit of tolerance for likelihood:",logZero
            endif
        enddo
        if(k==nptPerProc .or. j==10) then
            if(k==nptPerProc) then
                i=mod(nptPerProc,10)
                if(i==0) i=10
            else
                i=10
                j=0
            endif
            
            !first write the points generated by the root node
            nstart=nend+1
            p(1:ndims,nstart:nstart+i-1)=pnewP(1:ndims,1:i)
            phyP(1:totPar,nstart:nstart+i-1)=phyPnewP(1:totPar,1:i)
            l(nstart:nstart+i-1)=lnewP(1:i)
            nend=nstart+i-1
            
            if( outfile ) then
                !now write this batch to the files
                do m=nstart,nend
                    write(u_live,fmt) p(1:ndims,m),l(m)
                    write(u_phys,fmt2) phyP(:,m),l(m),1
                enddo
            endif
        endif
        if(k==nptPerProc) exit
    enddo

    deallocate( pnewP, phyPnewP, lnewP )

    if( outfile ) then
        close(u_live)
        close(u_phys)
    endif
    genLive=.false.
    resumeFlag=.false.
    
  end subroutine gen_initial_live

!----------------------------------------------------------------------
  
  subroutine getrandom(n,x,id)
    
    implicit none
    
    integer, intent(in) :: n
    double precision, intent(out) :: x(:)
    integer i,id
    
    ! --- uniform prior ----
    do i = 1,n

       x(i)=ranmarNS(id)
    enddo
    return
  end subroutine getrandom
  
!----------------------------------------------------------------------

!read the resume file
subroutine readResumeFile(funit1, peswitch, phyP, l, ic_n, ic_npt, ic_done, &
        ic_fNode, ic_nBrnch, ic_brnch, ic_Z, ic_info, ic_vnow, ic_eff, &
        ic_reme, ic_rFlag, eswitch, eswitchff, totVol, p)
    implicit none
    
    double precision p(ndims,nlive+1)
    integer funit1, i, iostatus, j, k
    logical peswitch
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood
    double precision d1
    logical eswitch
    integer, allocatable :: eswitchff(:)
    double precision, dimension(:), allocatable :: totVol
    
    !isolated cluster info
    integer ic_n !no. of nodes
    integer, allocatable :: ic_npt(:)
    logical, allocatable :: ic_done(:)
    integer, dimension(:), allocatable :: ic_fNode, ic_nBrnch
    double precision, dimension(:,:,:),  allocatable :: ic_brnch   
    double precision, dimension(:),  allocatable :: ic_Z, ic_info, ic_vnow
    double precision, dimension(:,:),  allocatable :: ic_eff
    logical, dimension(:),  allocatable :: ic_reme, ic_rFlag
    
    funit1=u_resume

    open(unit=funit1,file=resumename,status='old')
                
    read(funit1,*)genLive
    read(funit1,*)globff,numlike,ic_n,nlive
    read(funit1,*)gZ,ginfo
    ic_rFlag(1:ic_n)=.true.
    read(funit1,*)eswitch
    peswitch=eswitch
    if(eswitch) totVol=1.d0
    
    !read branching info
    do i=1,ic_n
        read(funit1,*)ic_nBrnch(i)
        if(ic_nBrnch(i)>0) read(funit1,*)ic_brnch(i,1:ic_nBrnch(i),1),ic_brnch(i,1:ic_nBrnch(i),2)
    enddo
    
    !read the node info
    do i=1,ic_n
        read(funit1,*)ic_done(i),ic_reme(i),ic_fNode(i),ic_npt(i)
        read(funit1,*)ic_vnow(i),ic_Z(i),ic_info(i)
    enddo
    if(.not.eswitch .or. ic_n==1) ic_npt(1)=nlive
    close(funit1)
    eswitchff=0
    
    !read hypercube-live file
    open(unit=u_live,file=livename,status='old')
    i=0
    do
        i=i+1
        read(u_live,*,IOSTAT=iostatus) p(1:ndims,i),l(i)
        if(iostatus<0) then
            i=i-1
            if(i<nlive) then
                write(*,*)"ERROR: live points file has less than ",nlive," points."
                write(*,*)"Aborting"
                stop
            endif
            exit
        endif
        if(i>nlive) then
            write(*,*)"ERROR: live points file has greater than ",nlive," points."
            write(*,*)"Aborting"
            stop
        endif
    enddo
    close(u_live)

    !read physical-live file
    open(unit=u_phys,file=physname,status='old')
    i=0
    do
        i=i+1
        read(u_phys,*,IOSTAT=iostatus) phyP(1:totPar,i),d1,j
        if(iostatus<0) then
            i=i-1
            if(i<nlive) then
                write(*,*)"ERROR: phys live points file has less than ",nlive," points."
                write(*,*)"Aborting"
                stop
            endif
            exit
        endif
        if(i>nlive) then
            write(*,*)"ERROR: phys live points file has greater than ",nlive," points."
            write(*,*)"Aborting"
            stop
        endif
    enddo
    close(u_phys)

    return
end subroutine readResumeFile

!----------------------------------------------------------------------

!find the highest likelihood for each node
subroutine getHighestLLHperNode(ic_n, ic_npt, ic_done, ic_llimits, ic_plimits, &
        ic_climits, ic_Z, ic_vnow, ic_hilike, ic_inc, l, phyP, lowlike, p)
    implicit none
    
    integer j, i, k
    double precision l(nlive+1) !log-likelihood
    double precision phyP(totPar,nlive+1) !physical live points
    double precision lowlike !lowest log-like
    double precision p(ndims,nlive+1)
     
    !isolated cluster info
    integer ic_n !no. of nodes
    integer, allocatable :: ic_npt(:)
    logical, allocatable :: ic_done(:)
    double precision, dimension(:,:,:),  allocatable :: ic_llimits, ic_plimits
    double precision, dimension(:,:,:) :: ic_climits
    double precision, dimension(:),  allocatable :: ic_vnow, ic_hilike, ic_inc
    double precision, dimension(:),  allocatable :: ic_Z
    
    j=0
    ic_done(0)=.true.
    do i=1,ic_n
        ic_hilike(i)=maxval(l(j+1:j+ic_npt(i)))
        lowlike=minval(l(j+1:j+ic_npt(i)))
        ic_inc(i)=ic_hilike(i)+log(ic_vnow(i))-ic_Z(i)
        if(ic_npt(i)<ndims+1 .or. ((abs(lowlike-ic_hilike(i))<= 0.0001d0 .or. &
            (ic_inc(i)<log(tol) .and. globff-nlive>50)) .and. minIter<1)) then
            ic_done(i)=.true.
        else
            ic_done(i)=.false.
            ic_done(0)=.false.
        endif
        
        !set the limits
        if(.not.ic_done(i)) then
            ic_llimits(i,:,1)=p(:,j+1)
            ic_llimits(i,:,2)=p(:,j+1)
            ic_plimits(i,1:nCdims,1)=phyP(1:nCdims,j+1)
            ic_plimits(i,1:nCdims,2)=phyP(1:nCdims,j+1)
            do k=j+2,j+ic_npt(i)
            call setLimits(ndims,nCdims,ic_llimits(i,:,:), &
                ic_plimits(i,:,:),p(:,k),phyP(1:nCdims,k),ic_climits(i,:,:))
                
            enddo
        endif
        j=j+ic_npt(i)
    enddo
    return
end subroutine getHighestLLHperNode

!----------------------------------------------------------------------

subroutine writeResumeFile(fName1, ic_n, eswitch, ic_nBrnch, ic_brnch, &
        ic_done, ic_reme, ic_fNode, ic_npt, ic_vnow, ic_Z, ic_info, ic_eff)
    implicit none
    
    character(len=100) fName1
    integer funit1
    integer ic_n, i
    logical eswitch
    integer, dimension(:), allocatable :: ic_nBrnch
    double precision, dimension(:,:,:),  allocatable :: ic_brnch
    logical, allocatable :: ic_done(:)
    logical, dimension(:),  allocatable :: ic_reme
    integer, dimension(:), allocatable :: ic_fNode
    integer, allocatable :: ic_npt(:)
    double precision, dimension(:),  allocatable :: ic_vnow, ic_Z, ic_info
    double precision, dimension(:,:),  allocatable :: ic_eff
    character(len=100) fmt

    funit1=u_resume
    fName1=resumename
    write(fmt,'(a,i5,a)')  '(',totPar+1,'E28.18,i4)'
    open(unit=funit1,file=fName1,form='formatted',status='replace')
    write(funit1,'(l2)').false.
    write(funit1,'(4i12)')globff,numlike,ic_n,nlive
    write(funit1,'(2E28.18)')gZ,ginfo
    write(funit1,'(l2)')eswitch
    !write branching info
    do i=1,ic_n
        write(funit1,'(i4)')ic_nBrnch(i)
        if(ic_nBrnch(i)>0) then
            write(fmt,'(a,i5,a)')  '(',2*ic_nBrnch(i),'E28.18)'
            write(funit1,fmt)ic_brnch(i,1:ic_nBrnch(i),1),ic_brnch(i,1:ic_nBrnch(i),2)
        endif
    enddo
    !write the node info
    do i=1,ic_n
        write(funit1,'(2l2,2i6)')ic_done(i),ic_reme(i),ic_fNode(i),ic_npt(i)
        write(funit1,'(3E28.18)')ic_vnow(i),ic_Z(i),ic_info(i)
    enddo
    close(funit1)
    return
end subroutine writeResumeFile    

!----------------------------------------------------------------------

subroutine addContribution(ic_n, ic_vnow, ic_npt, logX, l, ic_Zold, ic_Z, &
        ic_info, gZold)
    implicit none
    
    integer j, i1, ic_n, i
    double precision d1, logX
    double precision l(nlive+1) !log-likelihood
    double precision, dimension(:),  allocatable :: ic_vnow, ic_Z, ic_info
    double precision, dimension(:),  allocatable :: ic_Zold
    integer, allocatable :: ic_npt(:)
    double precision gZold !global evidence & info
    
    j=0
    do i1=1,ic_n
        !live point's contribution to the evidence
        logX=log(ic_vnow(i1)/dble(ic_npt(i1)))
        do i=j+1,j+ic_npt(i1)
            d1=l(i)+logX
            !global evidence
            gZold=gZ
            gZ=LogSumExp(gZ,d1)
            !ginfo=exp(d1-gZ)*l(i)+exp(gZold-gZ)*(ginfo+gZold)-gZ
            ginfo=ginfo*exp(gZold-gZ)+exp(d1-gZ)*l(i)
            !local evidence
            !gZold=ic_Z(i1)
            ic_Zold(i1)=ic_Z(i1)
            ic_Z(i1)=LogSumExp(ic_Z(i1),d1)
            !ic_info(i1)=exp(d1-ic_Z(i1))*l(i)+ &
            !exp(gZold-ic_Z(i1))*(ic_info(i1)+gZold)-ic_Z(i1)
            ic_info(i1)=ic_info(i1)*exp(ic_Zold(i1)-ic_Z(i1))+exp(d1-ic_Z(i1))*l(i)
        enddo
        j=j+ic_npt(i1)
    enddo

    ginfo=ginfo-gZ
    return
end subroutine addContribution

!----------------------------------------------------------------------

subroutine modeSeparation(l, phyP, p, ic_n, ic_npt, ic_done, &
        ic_reme, ic_vnow, ic_rFlag, ic_climits, ic_sc, sc_npt, sc_n, ic_nsc, &
        ic_fNode, ic_nBrnch, ic_brnch, nptk, ic_chk, modeFound, nodek, q, k, &
        ic_Zold, ic_Z, ic_info, ic_hilike, ic_eff, eVolFrac, ic_plimits, &
        ic_llimits, ic_volFac, ic_inc, aux, naux)
    implicit none

    integer naux !dimensionality of aux points 
    integer i, j, m, k, q
    integer nCdim
    double precision, dimension(:,:), allocatable :: aux
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood
    
    !isolated cluster info
    integer ic_n !no. of nodes
    integer, allocatable :: ic_sc(:), ic_npt(:)
    logical, allocatable :: ic_done(:)
    integer, dimension(:), allocatable :: ic_fNode, ic_nsc, ic_nBrnch
    double precision, dimension(:,:,:),  allocatable :: ic_brnch, ic_llimits, ic_plimits
    double precision, allocatable :: ic_volFac(:)
    double precision, dimension(:,:,:) :: ic_climits
    double precision, dimension(:),  allocatable :: ic_Z, ic_Zold, ic_info, ic_vnow, ic_hilike, ic_inc
    double precision, dimension(:,:),  allocatable :: ic_eff
    logical, dimension(:),  allocatable :: ic_reme, ic_rFlag, ic_chk
    logical modeFound

    !sub-cluster properties
    integer sc_n !no. of sub-clusters
    integer, dimension(:), allocatable :: sc_npt, nptk, nodek
    
    double precision, dimension(:,:,:), allocatable :: eVolFrac
    
    nCdim=nCdims
    !aux information to be re-arranged with the live points
    naux=ndims+totPar+1-nCdim
    aux(1,1:nlive)=l(1:nlive)
    aux(2:totPar+1,1:nlive)=phyP(1:totPar,1:nlive)
    aux(totPar+2:naux,1:nlive)=p(nCdim+1:ndims,1:nlive)
    
    !save old no. of modes
    i=ic_n
    nptk(1:i)=ic_npt(1:i)
    
    !decide which nodes to check for mode separation
    !All isolated clusters, which are not done, may be used for
    !mode separation. If ic_inc(j) < log(0.5) then separation is 
    !not necessary.
    ic_chk(1:ic_n)=.not.(ic_done(1:ic_n))
    do j=1,ic_n
        if(ic_chk(j)) then
            if(ic_inc(j)<log(0.5d0)) ic_chk(j)=.false.
        endif
    enddo
    modeFound=isolateModes2(nlive,ndims,nCdim,p(1:nCdim,1:nlive),naux,aux(1:naux,:), &
        ic_n,ic_fnode,ic_npt,ic_reme,ic_chk(1:ic_n),ic_vnow(1:ic_n),ic_rFlag,ic_climits(:,1:nCdim,:))
    if(modeFound) then			
        !re-arrange info
        l(1:nlive)=aux(1,1:nlive)
        phyP(1:totPar,1:nlive)=aux(2:totPar+1,1:nlive)
        p(nCdim+1:ndims,1:nlive)=aux(totPar+2:naux,1:nlive)
        if(nCdims<ndims .or. .true.) then
            ic_sc(i+1:ic_n)=1
            sc_npt(sc_n+1:sc_n+ic_n-i)=ic_npt(i+1:ic_n)
            sc_n=sc_n+ic_n-i
        endif
        !new nodes should inherit nsc from their parents
        do j=i+1,ic_n
            ic_nsc(j)=nsc_def
            ic_nsc(ic_fNode(j))=nsc_def
            !branching info
            ic_nBrnch(ic_fNode(j))=ic_nBrnch(ic_fNode(j))+1
            ic_brnch(ic_fNode(j),ic_nBrnch(ic_fNode(j)),1)=dble(j)
            ic_brnch(ic_fNode(j),ic_nBrnch(ic_fNode(j)),2)=dble(ic_npt(j))/dble(nptk(ic_fNode(j)))
        enddo
        
        !TODO: Check what is meant with new nodes and father nodes.
        ! Are these the clusters? If yes, why do we need the father 
        ! nodes?
        !assign prior volumes
        nodek(1:i)=0
        !first the new nodes
        q=sum(ic_npt(1:i))
        do j=i+1,ic_n
            if(ic_rFlag(j)) then
                !parent node
                k=ic_fNode(j)
                !divide prior volume
                ic_vnow(j)=ic_vnow(k)*dble(ic_npt(j))/dble(nptk(k))
                !divide the evidences & info
                ic_zold(j)=ic_z(k)
                ic_Z(j)=log(dble(ic_npt(j))/dble(nptk(k)))+ic_Z(k)
                ic_info(j)=ic_info(k)*dble(ic_npt(j))/dble(nptk(k))*exp(ic_zold(j)-ic_z(j))
                nodek(k)=1 !note the parent node
                ic_hilike(j)=maxval(l(q+1:q+ic_npt(j)))
                if(ic_hilike(j)-minval(l(q+1:q+ic_npt(j)))<= 0.0001) ic_done(j)=.true.
                !ellipsoidal decomposition prediction variables
                eVolFrac(j,:,:)=0d0
                eVolFrac(k,:,:)=0d0
                
                !set the limits
                ic_climits(j,:,:)=ic_climits(k,:,:)
                ic_llimits(j,:,:)=ic_llimits(k,:,:)
                ic_plimits(j,:,:)=ic_plimits(k,:,:)
                ic_volFac(j)=ic_volFac(k)
            endif
            q=q+ic_npt(j)
        enddo
        !now the father nodes
        q=0
        m=0
        do j=1,i
            if(nodek(j)==1) then
                ic_vnow(j)=ic_vnow(j)*dble(ic_npt(j))/dble(nptk(j))
                if(ic_npt(j)==0) then
                    ic_Z(j)=logZero
                    ic_info(j)=0.d0
                else
                    ic_zold(j)=ic_z(j)
                    ic_Z(j)=log(dble(ic_npt(j))/dble(nptk(j)))+ic_Z(j)
                    ic_info(j)=ic_info(j)*dble(ic_npt(j))/dble(nptk(j))*exp(ic_zold(j)-ic_z(j))
                endif
                ic_hilike(j)=maxval(l(q+1:q+ic_npt(j)))
                if(ic_hilike(j)-minval(l(q+1:q+ic_npt(j)))<= 0.0001) ic_done(j)=.true.
                
                sc_npt(m+1)=ic_npt(j)
                sc_npt(m+2:m+ic_sc(j))=0
            endif
            
            q=q+ic_npt(j)
            m=m+ic_sc(j)
        enddo
        ic_rFlag(1:ic_n)=.true.
    endif
    return
end subroutine modeSeparation

!----------------------------------------------------------------------
! Each cluster got a subcluster by G-means. If the cluster
! has no points, we can set its volume to 0 and skip the 
! corresponding subclusters.
! If the isolated cluster is done, we skip the subclusters
! as well.
subroutine subclustering(p, phyP, l, eswitch, neVol, totVol, x1, x2, &
        y1, y2, slope, intcpt, cVolFrac, pVolFrac, eVolFrac, ic_n, ic_sc, &
        ic_npt, ic_done, ic_nsc, ic_llimits, ic_climits, ic_volFac, ic_vnow, &
        ic_hilike, ic_eff, ic_rFlag, sc_n, sc_npt, nptk, sc_node, nodek, sck, &
        meank, sc_eval, evalk, sc_invcov, invcovk, sc_evec, eveck, tMatk, &
        kfack, volk, effk, sc_mean, sc_tMat, sc_kfac, sc_eff, sc_vol, naux, &
        aux, pt, lowlike, indx, flag, flag2, eswitchff, iteration, &
        sc_newCluster, d1, d4)
    implicit none
    
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood
    logical eswitch !whether to do ellipsoidal sampling or not
    
    !diagnostics for determining when to do eigen analysis
    integer neVol
    double precision, dimension(:), allocatable :: totVol, x1, x2, y1, y2, &
        slope, intcpt, cVolFrac, pVolFrac
    double precision, dimension(:,:,:), allocatable :: eVolFrac
    
    !isolated cluster info
    integer ic_n !no. of nodes
    integer, allocatable :: ic_sc(:), ic_npt(:)
    logical, allocatable :: ic_done(:)
    integer, dimension(:), allocatable :: ic_nsc !no. of iterations per every sub-clustering step
    double precision, dimension(:,:,:),  allocatable :: ic_llimits
    double precision, allocatable :: ic_volFac(:)
    double precision, dimension(:,:,:) :: ic_climits
    double precision, dimension(:),  allocatable :: ic_vnow, ic_hilike
    double precision, dimension(:,:),  allocatable :: ic_eff
    logical, dimension(:),  allocatable :: ic_rFlag
    
    !sub-cluster properties
    integer sc_n !no. of sub-clusters
    integer sc_newCluster ! No. of sub-clusters created in the last subclustering
    integer, dimension(:), allocatable :: sc_npt, nptk, sc_node, nodek, sck
    double precision, dimension(:,:), allocatable :: meank, sc_eval, evalk
    double precision, dimension(:,:,:), allocatable :: sc_invcov, invcovk, &
        sc_evec, eveck, tMatk
    double precision, dimension(:), allocatable :: kfack, volk, effk
    double precision, allocatable :: sc_mean(:,:),  sc_tMat(:,:,:), &
        sc_kfac(:), sc_eff(:), sc_vol(:)
    
    !auxiliary points (to be re-arranged with main points during clustering)
    integer naux !dimensionality of aux points 
    double precision, dimension(:,:), allocatable :: aux, pt
    
    !rejected point info
    double precision lowlike !lowest log-like
    integer indx(1) !point no. of lowlike TODO: Check, why we have an array here.
    
    logical flag, flag2
    integer, allocatable :: eswitchff(:)
    integer iteration, loopi
    integer noTraversedSC, noTraversedPoints, i1, i2, i3, i4, i, n2
    double precision d1, d2, d3, d4, d5
    
    logical giveMoreInfo
    
    noTraversedSC=0 !no. of sub-clusters traversed q
    noTraversedPoints=0 !no. of points traversed m
    sc_newCluster=0 !no. of sub-clusters created so far n
    do i1=1,ic_n
        if(ic_npt(i1)==0) then
            totvol(i1) = 0d0
            sck(i1)=0
            noTraversedSC=noTraversedSC+ic_sc(i1) !no. of sub-clusters traversed
            cycle
        endif
        
        if(ic_done(i1)) then
            sck(i1)=1
            noTraversedSC=noTraversedSC+ic_sc(i1) !no. of sub-clusters traversed
            noTraversedPoints=noTraversedPoints+ic_npt(i1)
            nptk(sc_newCluster+1)=ic_npt(i1)
            nodek(sc_newCluster+1)=i1
            sc_newCluster=sc_newCluster+1
            cycle
        endif
        if(ic_npt(i1)>=ndims+1 .and. totVol(i1)==0.d0) ic_rFlag(i1)=.true.
        
        if(.not.eswitch) then
            if(mod(iteration-1-eswitchff(i1),ic_nsc(i1))==0) then
                flag=.true.
            else
                flag=.false.
            endif
            !ic_rFlag might be false because of something earlier.
        elseif(ic_npt(i1)<ndims+1 .and. .not.ic_rFlag(i1)) then
            flag=.false.
        elseif(ic_rFlag(i1)) then
            flag=.true.
        else
            flag=.false.
            d4=ef
            
            !current vol fraction
            cVolFrac(i1)=totVol(i1)*d4/((ic_vnow(i1)*ic_volFac(i1)))
    
            !predicted vol fraction
            
            pVolFrac(i1)=0d0
            i3=0
            do i2=1,neVol
                if(eVolFrac(i1,i2,1)<=0d0) exit
                pVolFrac(i1)=pVolFrac(i1)+eVolFrac(i1,i2,1)
                i3=i3+1
            enddo
            if(pVolFrac(i1)==0d0) then
                pVolFrac(i1)=1d0
            else
                pVolFrac(i1)=pVolFrac(i1)/dble(i3)
            endif

            if(.not.flag .and. (cVolFrac(i1)>1.1 .or. .not.dino) .and.  &
            (pVolFrac(i1)<cVolFrac(i1) .or. mod(iteration-1-eswitchff(i1),ic_nsc(i1))==0)) flag=.true.
    
            if(flag .and. dino) then
                do i=noTraversedSC+1,noTraversedSC+ic_sc(i1)
                    if(sc_eff(i)<=1.00001) exit
                    if(i==sc_n) flag=.false.
                enddo
            endif

            if(.not.flag .and. eVolFrac(i1,2*neVol,1)==0.d0 .and. mod(iteration-1-eswitchff(i1),ic_nsc(i1))==0)    flag=.true.
        endif
        ! flag has been set just before.
        if(flag) then
            !target vol
            d4=ef

            d5=1d0
            do i2=1,ndims
                d5=d5/(ic_llimits(i1,i2,2)-ic_llimits(i1,i2,1))
            enddo
            d1=ic_vnow(i1)*d5/d4
    
            !volume threshold
            if(ic_rFlag(i1)) then
                d2=1000.*d1
            else
                d2=totVol(i1)*d5/ic_volFac(i1)
            endif
        
            !volume tolerance
            d3=0.d0
        
            eswitchff(i1)=iteration-1
            
            flag=.false.
            
            !aux information to be re-arranged with the live points
            naux=totPar+2
            !rescaling
            do i3=1,ndims
                !rescale back into unit hypercube
                d4=ic_climits(i1,i3,2)-ic_climits(i1,i3,1)
                pt(i3,1:ic_npt(i1))=ic_climits(i1,i3,1)+d4*p(i3,noTraversedPoints+1:noTraversedPoints+ic_npt(i1))
                
                if(pWrap(i3)) then
                    do i4=1,ic_npt(i1)
                        call wraparound(pt(i3,i4),pt(i3,i4))
                    enddo
                endif
                
                !scale to the new limits
                d4=ic_llimits(i1,i3,2)-ic_llimits(i1,i3,1)
                pt(i3,1:ic_npt(i1))=(pt(i3,1:ic_npt(i1))-ic_llimits(i1,i3,1))/d4
            enddo
            aux(1,1:ic_npt(i1))=l(noTraversedPoints+1:noTraversedPoints+ic_npt(i1))
            aux(2:totPar+1,1:ic_npt(i1))=phyP(1:totPar,noTraversedPoints+1:noTraversedPoints+ic_npt(i1))
            lowlike=minval(l(noTraversedPoints+1:noTraversedPoints+ic_npt(i1)))
            aux(naux,1:ic_npt(i1))=(l(noTraversedPoints+1:noTraversedPoints+ic_npt(i1))-lowlike)/(ic_hilike(i1)-lowlike)
            
            !max no. of sub-clusters allowed
            n2=min(ceiling(dble(ic_npt(i1))/dble(min_pt)),maxeCls-sc_newCluster)
            if(dino) then
                d5=d1
                do i3=1,50
                    if(Dinosaur(pt(:,1:ic_npt(i1)),ic_npt(i1),ndims,sck(i1),nptk(sc_newCluster+1:sc_newCluster+n2), &
                        naux,aux(1:naux,1:ic_npt(i1)),min_pt,n2,meank(sc_newCluster+1:sc_newCluster+n2,:),invcovk(sc_newCluster+1:sc_newCluster+n2,:,:), &
                        tMatk(sc_newCluster+1:sc_newCluster+n2,:,:),evalk(sc_newCluster+1:sc_newCluster+n2,:),eveck(sc_newCluster+1:sc_newCluster+n2,:,:),kfack(sc_newCluster+1:sc_newCluster+n2), &
                        effk(sc_newCluster+1:sc_newCluster+n2),volk(sc_newCluster+1:sc_newCluster+n2),d1,d2,neVol,eVolFrac(i1,:,:),globff,d3,.false., &
                        ic_rFlag(i1),.false.,nCdims)) then
                        if(eswitch) then
                            ic_nsc(i1)=max(1,ic_nsc(i1)-10)
                        else
                            ! Only for the first time subclustering is ever done.
                            eswitch=.true.
                            ic_sc(1)=sck(i1)
                        endif
                    
                        scount=scount+1
                        totVol(i1)=sum(volk(sc_newCluster+1:sc_newCluster+sck(i1)))
                    
                        flag=.true.
                        flag2=.true.
                        
                        !set the limits
                        ic_climits(i1,:,:)=ic_llimits(i1,:,:)
                        
                        ic_volFac(i1)=1d0
                        do i2=1,ndims
                            ic_volFac(i1)=ic_volFac(i1)/(ic_climits(i1,i2,2)-ic_climits(i1,i2,1))
                        enddo
                        
                        exit
                    elseif(.not.ic_rFlag(i1)) then
                        if(eswitch) ic_nsc(i1)=ic_nsc(i1)+10
                        exit
                    endif
                    eVolFrac(i1,1,1)=eVolFrac(i1,1,1)*d1/d5
                
                    if(mod(i3,5)==0) then
                        d2=d2*2.d0
                        d1=d1*2.d0
                    endif
                enddo
            endif
            if(.not.Flag .and. ic_rFlag(i1)) then
                sck(i1)=1
                nptk(sc_newCluster+1)=ic_npt(i1)
                kfack(sc_newCluster+1)=0.d0
                volk(sc_newCluster+1)=0.d0
                effk(sc_newCluster+1)=1.d0
                totVol(i1)=0.d0
                flag=.true.
                ic_done(i1)=.true.
                ic_done(0)=.true.
                do i3=1,ic_n
                    if(.not.ic_done(i3)) then
                        ic_done(0)=.false.
                        exit
                    endif
                enddo
            else
                !predict the current vol fraction
                x1(i1)=sum(eVolFrac(i1,neVol+1:2*neVol,2))/neVol
                y1(i1)=sum(eVolFrac(i1,neVol+1:2*neVol,1))/neVol
                x2(i1)=sum(eVolFrac(i1,1:neVol,2))/neVol
                y2(i1)=sum(eVolFrac(i1,1:neVol,1))/neVol
                slope(i1)=(y2(i1)-y1(i1))/(x2(i1)-x1(i1))
                intcpt(i1)=(x2(i1)*y1(i1)-x1(i1)*y2(i1))/(x2(i1)-x1(i1))
            endif
        endif
        if(ic_rFlag(i1) .and. .not.eswitch) exit
            
        if(flag) then
            !aux information to be re-arranged with the live points
            p(:,noTraversedPoints+1:noTraversedPoints+ic_npt(i1))=pt(:,1:ic_npt(i1))
            l(noTraversedPoints+1:noTraversedPoints+ic_npt(i1))=aux(1,1:ic_npt(i1))
            phyP(1:totPar,noTraversedPoints+1:noTraversedPoints+ic_npt(i1))=aux(2:totPar+1,1:ic_npt(i1))
        else
            sck(i1)=ic_sc(i1)
            meank(sc_newCluster+1:sc_newCluster+sck(i1),:)=&
                    sc_mean(noTraversedSC+1:noTraversedSC+sck(i1),:)
            invcovk(sc_newCluster+1:sc_newCluster+sck(i1),:,:)=&
                    sc_invcov(noTraversedSC+1:noTraversedSC+sck(i1),:,:)
            tMatk(sc_newCluster+1:sc_newCluster+sck(i1),:,:)=&
                    sc_tMat(noTraversedSC+1:noTraversedSC+sck(i1),:,:)
            evalk(sc_newCluster+1:sc_newCluster+sck(i1),:)=&
                    sc_eval(noTraversedSC+1:noTraversedSC+sck(i1),:)
            eveck(sc_newCluster+1:sc_newCluster+sck(i1),:,:)=&
                    sc_evec(noTraversedSC+1:noTraversedSC+sck(i1),:,:)
            kfack(sc_newCluster+1:sc_newCluster+sck(i1))=&
                    sc_kfac(noTraversedSC+1:noTraversedSC+sck(i1))
            effk(sc_newCluster+1:sc_newCluster+sck(i1))=&
                    sc_eff(noTraversedSC+1:noTraversedSC+sck(i1))
            volk(sc_newCluster+1:sc_newCluster+sck(i1))=&
                    sc_vol(noTraversedSC+1:noTraversedSC+sck(i1))
            nptk(sc_newCluster+1:sc_newCluster+sck(i1))=&
                    sc_npt(noTraversedSC+1:noTraversedSC+sck(i1))
            slope(i1)=slope(i1)*1.01
        endif
        nodek(sc_newCluster+1:sc_newCluster+sck(i1))=i1
        
        noTraversedSC=noTraversedSC+ic_sc(i1) !no. of sub-clusters traversed
        noTraversedPoints=noTraversedPoints+ic_npt(i1) !no. of points traversed
        sc_newCluster=sc_newCluster+sck(i1) !no. of sub-clusters created so far
    enddo
    
    ! True if Dinosaur clustering has been used to create more subcluster.
    if(flag2) then
        !update mode info
        ic_sc(1:ic_n)=sck(1:ic_n)
        !update sub-cluster info
        sc_n=sum(ic_sc(1:ic_n))
        sc_mean(1:sc_n,:)=meank(1:sc_n,:)
        sc_invcov(1:sc_n,:,:)=invcovk(1:sc_n,:,:)
        sc_tMat(1:sc_n,:,:)=tMatk(1:sc_n,:,:)
        sc_eval(1:sc_n,:)=evalk(1:sc_n,:)
        sc_evec(1:sc_n,:,:)=eveck(1:sc_n,:,:)
        sc_kfac(1:sc_n)=kfack(1:sc_n)
        sc_eff(1:sc_n)=effk(1:sc_n)
        sc_vol(1:sc_n)=volk(1:sc_n)
        sc_npt(1:sc_n)=nptk(1:sc_n)
        sc_node(1:sc_n)=nodek(1:sc_n)
        !find the index of the low-like points because of re-arrangement
        indx=minloc(l(1:nlive))
    endif
    
    return
end subroutine subclustering

!----------------------------------------------------------------------

subroutine adjustPriorVol(p, phyP, l, i3, nd, nd_j, d1, d4, h, vprev, vnext, &
        shrink, ic_npt, ic_llimits, ic_plimits, ic_climits, ic_Z, ic_vnow, &
        ic_hilike, ic_inc, lowlike, lowp, lowPhyP, indx, pt)
    implicit none
    
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood
    
    !misc
    integer i3, nd, nd_j
    double precision d1, d4
    double precision h, vprev, vnext, shrink !prior volume

    !isolated cluster info
    integer, allocatable :: ic_npt(:)
    double precision, dimension(:,:,:),  allocatable :: ic_llimits, ic_plimits
    double precision, dimension(:,:,:) :: ic_climits
    double precision, dimension(:),  allocatable :: ic_Z, ic_vnow, ic_hilike, ic_inc
    
    !rejected point info
    double precision lowlike !lowest log-like
    double precision, allocatable :: lowp(:), lowPhyP(:) !point with the lowlike
    integer indx(1) !point no. of lowlike
    
    double precision, dimension(:,:), allocatable :: pt
    
    vprev=ic_vnow(nd)/shrink
    vnext=ic_vnow(nd)*shrink
    h=(vprev-vnext)/2.d0
    ic_inc(nd)=ic_hilike(nd)+log(ic_vnow(nd))-ic_Z(nd)

    !find lowlike. nd_j is the number of points in finished cluster.
    indx=minloc(l(nd_j+1:nd_j+ic_npt(nd)))
    indx(1)=indx(1)+nd_j
    lowlike=l(indx(1))
    lowp(:)=p(:,indx(1))
    lowphyP(:)=phyP(:,indx(1))
    
    !set the limits
    do i3=1,ndims
        d4=ic_climits(nd,i3,2)-ic_climits(nd,i3,1)
        d1=ic_climits(nd,i3,1)+d4*lowp(i3)
        
        if( abs( d1 - ic_llimits(nd,i3,1) ) < d4 * 1d-5 ) then
            pt(1,1:ic_npt(nd))=ic_climits(nd,i3,1)+d4*p(i3,nd_j+1:nd_j+ic_npt(nd))
            ic_llimits(nd,i3,1)=minval(pt(1,1:ic_npt(nd)),MASK=pt(1,1:ic_npt(nd))>d1)
        elseif( abs( d1 - ic_llimits(nd,i3,2) ) < d4 * 1d-5 ) then
            pt(1,1:ic_npt(nd))=ic_climits(nd,i3,1)+d4*p(i3,nd_j+1:nd_j+ic_npt(nd))
            ic_llimits(nd,i3,2)=maxval(pt(1,1:ic_npt(nd)),MASK=pt(1,1:ic_npt(nd))<d1)
        endif
    enddo
    ! multimodal.
    do i3=1,nCdims
        d4=ic_plimits(nd,i3,2)-ic_plimits(nd,i3,1)
        if( abs( lowPhyP(i3) - ic_plimits(nd,i3,1) ) < d4 * 1d-5 ) then
            ic_plimits(nd,i3,1)=minval(phyP(i3,nd_j+1:nd_j+ic_npt(nd)), &
                MASK=phyP(i3,nd_j+1:nd_j+ic_npt(nd))>lowPhyP(i3))
        elseif( abs( lowPhyP(i3) - ic_plimits(nd,i3,2) ) < d4 * 1d-5 ) then
            ic_plimits(nd,i3,2)=maxval(phyP(i3,nd_j+1:nd_j+ic_npt(nd)), &
                MASK=phyP(i3,nd_j+1:nd_j+ic_npt(nd))<lowPhyP(i3))
        endif
    enddo
    return
end subroutine adjustPriorVol

!----------------------------------------------------------------------
! Find a point inside the hard constraint, set the limits.
subroutine findPointHardConstr(context, p, phyP, l, i, n1, n2, nd, i2, i3, n, &
        eswitch, remFlag, acpt, ic_llimits, ic_plimits, ic_climits, ic_hilike, &
        sc_mean, sc_tMat, lowlike, indx, lnew, pnew, phyPnew, pnewaList, &
        phyPnewaList, lnewaList, loglike, ic_npt, nd_j)
    implicit none
    
    INTERFACE
        !the likelihood function
        subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
            integer n_dim,nPar,context_pass
            double precision lnew,Cube(nPar)
        end subroutine loglike
    end INTERFACE
    
    integer context
    integer id
    integer numberUsedThreads
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood

    !work variables

    !misc
    integer i, n1, n2, nd, i2, i3, n
    double precision d1
    logical eswitch !whether to do ellipsoidal sampling or not
    logical remFlag, acpt
    
    !isolated cluster info
    double precision, dimension(:,:,:),  allocatable ::  ic_llimits, ic_plimits
    double precision, dimension(:,:,:) :: ic_climits
    double precision, dimension(:),  allocatable :: ic_hilike

    !sub-cluster properties
    double precision, allocatable :: sc_mean(:,:), sc_tMat(:,:,:)

    !rejected point info
    double precision lowlike !lowest log-like
    integer indx(1) !point no. of lowlike

    !new point
    double precision lnew
    double precision, allocatable :: pnew(:), phyPnew(:) ! new point

    ! 2nd approach for parallelizing
    double precision, dimension(:,:,:), allocatable :: pnewaList, phyPnewaList
    double precision, dimension(:,:), allocatable :: lnewaList
    integer :: loopi, loopj, nd_j, abort, maxNumLike
    integer, allocatable :: ic_npt(:)
    logical :: logic(ic_npt(nd))
    
    integer :: counti, countf, count_rate, calls
    real :: dt
    calls = 0
    acpt=.false.
    !$OMP PARALLEL &
    !$OMP& NUM_THREADS(procs), &
    !$OMP& PRIVATE(id, logic, n, abort, loopi, loopj, counti, countf, count_rate, dt), &
    !$OMP& FIRSTPRIVATE(lowlike, indx), REDUCTION(MAX:maxNumLike), REDUCTION(+:calls)
    calls = 0
    maxNumLike = 0
    id = OMP_get_thread_num() +1
    ! Each thread shall get a different lowest likelihood, 
    ! e.g. thread 0 (with id = 1) gets the lowest likelihood, thread 1 the
    ! second lowest likelihood and so on. Unfortunately Fortran has a bug when
    ! using minloc() with a Mask created on the fly. We have to create an array
    ! of logicals beforehand or else the mask is used randomised.
    ! See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=19016 for
    ! further information. (26.11.2015 at 15:40 UMT+1)
    if(id > 1) then
        do loopi=1, id-1
            do loopj = nd_j+1, nd_j+ic_npt(nd)
                if(l(loopj) .gt. lowlike) then
                    logic(loopj-nd_j) = .true.
                else
                    logic(loopj-nd_j) = .false.
                endif
            enddo
            indx=minloc(l(nd_j+1:nd_j+ic_npt(nd)), MASK = logic)
            indx(1)=indx(1)+nd_j
            if(lowlike == l(indx(1))) then
                write(*,*) '***** - Failure - ***** '
                write(*,*) 'Thread', id, 'failed to get another point'
                write(*,*) 'Round', loopi
                write(*,*) 'nd_j: ', nd_j
                write(*,*) 'ic_npt(nd): ', ic_npt(nd)
                write(*,*) 'nd: ', nd
                write(*,*) 'Active points: ', l(nd_j+1:nd_j+ic_npt(nd))
                write(*,*) '*********************** '
                stop
            endif
            lowlike=l(indx(1))
        enddo
    endif
    abort = 0
    ! Now each thread samples a new point.
    do
        abort = abort + 1
        call samp(pnewaList(id,nd, :), phyPnewaList(id,nd, :), &
            lnewaList(id, nd),sc_mean(1,:),d1, sc_tMat(1,:,:), &
            ic_climits(nd,:,:), loglike,eswitch,lowlike,n,context, id-1)
        calls = calls + n
        maxNumLike = maxNumLike + 1
        if(lnewaList(id, nd)>lowlike) then
            exit
        endif
        if(mod(abort, 10000000) == 0) then
            write(*,*) '----- - Failure? - -----'
            write(*,*) 'Tried a lot of samplings: ', abort
            write(*,*) 'lowlike', lowlike
            write(*,*) 'lnew: ', lnewaList(id, nd)
            write(*,*) 'pnew: ', pnewaList(id, nd, :)
            write(*,*) 'Thread: ', id
            write(*,*) 'ic_npt(nd): ', ic_npt(nd)
            do loopi=1,ic_npt(nd)
                if(l(loopi+nd_j) < lowlike) write(*,*) 'at', loopi, 'lowlike', l(loopi+nd_j)
            enddo
            write(*,*) '------------------------'
            stop
        endif
    enddo
    !$OMP CRITICAL
    p(:,indx(1))=pnewaList(id, nd, :)
    phyP(:,indx(1))=phyPnewaList(id, nd, :)
    l(indx(1))=lnewaList(id, nd)
    if(lnewaList(id, nd)>ic_hilike(nd)) ic_hilike(nd)=lnewaList(id, nd)
    !$OMP END CRITICAL
    !$OMP END PARALLEL
    do i = 1, procs
        pnew(:)=pnewaList(i, nd, :)
        phyPnew(:)=phyPnewaList(i, nd, :)
        lnew=lnewaList(i, nd)
        call setLimits(ndims,nCdims,ic_llimits(nd,:,:), &
            ic_plimits(nd,:,:),pnew,phyPnew(1:nCdims),ic_climits(nd,:,:))
    enddo
    ! Since the +-reduction in OpenMP initializes variables to 0
    ! we have to add the number of samplings like this.
    numlike =  maxNumLike + numlike
    totalLikelihoodCalls = totalLikelihoodCalls + calls
    return
end subroutine findPointHardConstr

!----------------------------------------------------------------------
! Find a point inside the hard constraint, accept it with probability 
! 1/#clusters the point lies in and set the limits
subroutine findPointHardConstrEllipsoids(context, p, phyP, l, i, j, m, n, j1, nd, nd_i, &
        num_old, d1, urv, eswitch, acpt, ic_sc, ic_llimits, ic_plimits, &
        ic_climits, ic_hilike, sc_npt, sc_invcov, sc_mean, sc_tMat, sc_kfac, &
        sc_eff, sc_vol, lowlike, lnew, pnew, phyPnew, pnewaList, phyPnewaList, &
        lnewaList, loglike, ic_done, iList, lnewList, &
        pnewList, phyPnewList, indxList, lowpList, lowp, nd_j, indx, ic_npt, acptList)
    implicit none
    
    INTERFACE
        !the likelihood function
        subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
            integer n_dim,nPar,context_pass
            double precision lnew,Cube(nPar)
        end subroutine loglike
    end INTERFACE
    
    integer context
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood

    !misc
    ! nd_i is the number of traversed ellipsoids. If no cluster is finished
    ! nd_i is zero. Else it is the number of subclusters in a finished
    ! isolated cluster.
    integer i, j, m, n, j1, nd, nd_i, id, loopi, nd_j
    integer num_old
    double precision d1, urv
    logical eswitch !whether to do ellipsoidal sampling or not
    logical acpt
    
    !isolated cluster info
    integer, allocatable :: ic_sc(:)
    double precision, dimension(:,:,:),  allocatable :: ic_llimits, ic_plimits
    double precision, allocatable :: ic_climits(:,:,:)
    double precision, dimension(:),  allocatable :: ic_hilike

    !sub-cluster properties
    integer, dimension(:), allocatable :: sc_npt
    double precision, dimension(:,:,:), allocatable :: sc_invcov
    double precision, allocatable :: sc_mean(:,:), sc_tMat(:,:,:), sc_kfac(:), sc_eff(:), sc_vol(:)

    !rejected point info
    double precision lowlike !lowest log-like

    !new point
    double precision lnew
    double precision, allocatable :: pnew(:), phyPnew(:) ! new point
    
    ! 1st approach for parallelizing
    double precision, dimension(:,:,:), allocatable :: pnewaList, phyPnewaList
    double precision, dimension(:,:), allocatable :: lnewaList
    logical, allocatable, dimension(:) :: acptList
    double precision, allocatable :: pnewList(:,:), phyPnewList(:,:) ! new point for every thread
    integer, allocatable, dimension(:) :: iList
    double precision, allocatable, dimension(:) :: lnewList
    double precision, allocatable, dimension(:,:) :: lowpList
    integer, dimension(:,:), allocatable :: indxList
    double precision, allocatable :: lowp(:)
    integer abort, abortTwo, auxCounter, auxI
    logical, allocatable :: ic_done(:)
    integer :: indx(1), maxNumLike
    integer, allocatable :: ic_npt(:)
    integer :: test(1), loopj, loopk
    double precision :: tmpLowest
    logical :: logic(ic_npt(nd))
    
    integer :: counti, countf, count_rate, calls
    real :: dt
    
    num_old=numlike
    abortTwo = 0
    !first pick an ellipsoid according to the vol
    do

        abortTwo = abortTwo + 1
        !pick a sub-cluster according to the vol
        call selectEll(ic_sc(nd),sc_vol(nd_i+1:nd_i+ic_sc(nd)),i)
        !nd_i is the number of finished clusters.
        i=i+nd_i
        if(sc_kfac(i)==0.d0 .or. sc_vol(i)==0.d0) then
            cycle
        else
            exit
        endif
    enddo
    calls=0
    !$OMP PARALLEL &
    !$OMP& NUM_THREADS(procs), &
    !$OMP& PRIVATE (id, d1, abort, n,  loopi, logic, loopj, counti, countf, dt, count_rate), &
    !$OMP& FIRSTPRIVATE(lowlike, indx), REDUCTION(MAX:maxNumLike), REDUCTION(+:calls)
    calls=0
    maxNumLike = 0
    id = OMP_get_thread_num() + 1
    abort = 0
    acptList(id) = .false.
    !find lowlike. nd_j is the number of points in finished cluster.
    if(id < ic_npt(nd) + 1) then
        if(id > 1) then
        
            do loopi=1, id-1
                ! Each thread shall get a different lowest likelihood, 
                ! ! e.g. thread 0 (with id = 1) gets the lowest likelihood, thread 1 the
                ! second lowest likelihood and so on. Unfortunately Fortran has a bug when
                ! using minloc() with a Mask created on the fly. We have to create an array
                ! of logicals beforehand or else the mask is used randomised.
                ! See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=19016 for
                ! further information. (26.11.2015 at 15:40 UMT+1)
                do loopj = nd_j+1, nd_j+ic_npt(nd)
                    if(l(loopj) .gt. lowlike) then
                        logic(loopj-nd_j) = .true.
                    else
                        logic(loopj-nd_j) = .false.
                    endif
                enddo
                indx=minloc(l(nd_j+1:nd_j+ic_npt(nd)), MASK = logic)
                indx(1)=indx(1)+nd_j
                if(lowlike == l(indx(1))) then
                    write(*,*) 'Thread', id, 'failed to get another point'
                    write(*,*) 'Round', loopi
                    stop
                endif
                lowlike=l(indx(1))
                lowpList(id, :) = p(:, indx(1))
            enddo
        endif
    
        do
            abort = abort + 1
            
            d1=sc_kfac(i)*sc_eff(i)
            call samp(pnewaList(id,nd, :), phyPnewaList(id,nd, :), &
                lnewaList(id, nd),sc_mean(i,:),d1, &
                sc_tMat(i,:,:),ic_climits(nd,:,:),loglike,eswitch,lowlike,n,context, id-1)
            calls = calls + n
            !check if any of them is inside the hard edge
            maxNumLike=maxNumLike+1
            if(lnewaList(id, nd)>lowlike .and. &
                (.not. acptList(id) .or. lnewaList(id, nd) > lnewList(id))) then
                if(sc_npt(i)>0) then
                    !find the no. of ellipsoids j, the new points lies in
                    j=1
                    acptList(id)=.true.
                    ! Check if the point is in more than one ellipsoid.
                    do m=nd_i+1,nd_i+ic_sc(nd)
                        if(m==i .or. sc_npt(m)==0) cycle
                        if(ptIn1Ell(ndims,pnewaList(id, nd, :),sc_mean(m,:), &
                            sc_invcov(m,:,:),sc_kfac(m)*sc_eff(m))) j=j+1
                    enddo
                endif
                                
                if(acptList(id)) then
                    !accept it with probability 1/j where j is the number of 
                    ! ellipsoids the number lies in.
                    if(j>1) then
                        urv=ranmarNS(id-1)
                        if(urv<=(1.d0/dble(j))) then
                            acptList(id)=.true.
                            lnewList(id)=lnewaList(id, nd)
                            pnewList(id, :)=pnewaList(id, nd, :)
                            phyPnewList(id, :)=phyPnewaList(id, nd, :)
                            iList(id) = i
                        else
                            acptList(id)=.false.
                        endif
                    else
                        acptList(id)=.true.
                        lnewList(id)=lnewaList(id, nd)
                        pnewList(id, :)=pnewaList(id, nd, :)
                        phyPnewList(id, :)=phyPnewaList(id, nd, :)
                        iList(id)=i
                    endif
                endif
            endif

            ! If one point is accepted (due to probability, there might
            ! be no point accepted), we set the
            ! limits according to the new point.
            if(acptList(id)) then   
                !$OMP CRITICAL
                if(lnewList(id)>ic_hilike(nd)) ic_hilike(nd)=lnewList(id)
                !$OMP END CRITICAL
                indxList(id,:) = indx(:)
                exit
            endif
        enddo
    endif
    !$OMP END PARALLEL
    totalLikelihoodCalls = totalLikelihoodCalls + calls
    ! Since the +-reduction in OpenMP initializes variables to 0
    ! We have to add the number of samplings like this.
    numlike = maxNumLike + numlike
    do loopi=1, procs
        if(acptList(loopi)) then
            !set the limits
            call setLimits(ndims,nCdims,ic_llimits(nd,:,:), &
                ic_plimits(nd,:,:),pnewList(loopi,:),phyPnewList(loopi,1:nCdims),ic_climits(nd,:,:))
        endif
    enddo
    return
end subroutine findPointHardConstrEllipsoids

!----------------------------------------------------------------------
! Remove the rejected point and evolve the ellipsoid with the rejected and the 
! new point.
subroutine evolveEllipsoid(p, phyP, l, i, j, m, i1, n1, n2, q, nd, nd_i, nd_j, &
        d1, d4, shrink, ic_sc, ic_npt, ic_fnode, ic_volFac, ic_vnow, ic_inc, &
        ic_eff, sc_npt, sc_eval, sc_invcov, sc_mean, sc_kfac, sc_eff, sc_vol, &
        lowp, indx, pnew)
    implicit none
    
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood

    !misc
    integer i, j, m, i1, n1, n2, q, nd, nd_i, nd_j
    double precision d1, d4
    double precision shrink !prior volume
    
    !isolated cluster info
    integer, allocatable :: ic_sc(:), ic_npt(:)
    integer, dimension(:), allocatable :: ic_fnode
    double precision, allocatable :: ic_volFac(:)
    double precision, dimension(:),  allocatable ::  ic_vnow, ic_inc
    double precision, dimension(:,:),  allocatable :: ic_eff

    !sub-cluster properties
    integer, dimension(:), allocatable :: sc_npt 
    double precision, dimension(:,:), allocatable :: sc_eval
    double precision, dimension(:,:,:), allocatable :: sc_invcov
    double precision, allocatable :: sc_mean(:,:), sc_kfac(:), sc_eff(:), sc_vol(:)

    !rejected point info
    double precision, allocatable :: lowp(:) !point with the lowlike
    integer indx(1) !point no. of lowlike

    !new point
    double precision, allocatable :: pnew(:) ! new point
    
    !find the sub-cluster in which the rejected point lies
    ! nd_j is the current amount of points traversed. It is 0 if no cluster
    ! is finished, else it is the number of points of all finished cluster.
    n1=nd_j	
    do q=nd_i+1,nd_i+ic_sc(nd)
        ! If the current subcluster has no points in it, then we cycle.
        if(sc_npt(q)==0) cycle
        n1=n1+sc_npt(q)
        ! indx(1) is the pointer (actually just the index) to the current point
        ! in the list with all points. If the number of points saved in n1 is
        ! higher than the index of the current point, we just got the right 
        ! subcluster and we exit. The right subcluster has index q.
        if(indx(1)<=n1) then
            n1=n1-sc_npt(q)
            exit
        endif
    enddo

    !remove the rejected point & its likelihood & update the no. of points in the cluster
    if(indx(1)<nlive) then
        p(:,indx(1):nd_j+ic_npt(nd)-1)=p(:,indx(1)+1:nd_j+ic_npt(nd))
        phyP(:,indx(1):nd_j+ic_npt(nd)-1)=phyP(:,indx(1)+1:nd_j+ic_npt(nd))
        l(indx(1):nd_j+ic_npt(nd)-1)=l(indx(1)+1:nd_j+ic_npt(nd))
    endif
    sc_npt(q)=sc_npt(q)-1
    
    ! If this subcluster has no more points in it, we set its volume to zero 
    ! and ignore it in later iterations.
    if(sc_npt(q)==0) then
        sc_vol(q)=0.d0
        sc_kfac(q)=0.d0
        sc_eff(q)=1.d0
        if(ic_sc(nd)==1) ic_inc(nd)=log(tol)
    ! Else check if the current subcluster has a volume
    ! and kfac bigger than zero. 
    ! Also see if this is the ellipsoid i from sEll from sampling. If yes, it
    ! should have at least one more point.
    elseif(sc_vol(q)>0.d0 .and. sc_kfac(i)>0.d0 &
        .and. (i/=q .or. (i==q .and. sc_npt(q)>0))) then

        !min vol this ellipsoid should occupy
        if(dino) then
            d4=ef
            d1=(ic_vnow(nd)*ic_volFac(nd)*shrink/d4)*(sc_npt(q))/dble(ic_npt(nd))
        else
            d1=tiny(1.d0)
        endif
        !now evolve the ellipsoid with the rejected point
        call evolveEll(0,sc_npt(q),ndims,lowp,p(:,n1+1:n1+sc_npt(q)),sc_mean(q,:), &
        sc_eval(q,:),sc_invcov(q,:,:),sc_kfac(q),sc_eff(q),sc_vol(q),d1)
    endif
    
    n1=sum(sc_npt(1:i-1))
    
    if(i/=q .or. (i==q .and. sc_npt(q)>0)) then
        !min vol this ellipsoid should occupy
        if(dino) then
            d4=ef
            d1=(ic_vnow(nd)*ic_volFac(nd)*shrink/d4)*(dble(sc_npt(i))+1.d0)/dble(ic_npt(nd))
        else
            d1=tiny(1.d0)
        endif
        !now evolve the ellipsoid with the inserted point
        call evolveEll(1,sc_npt(i),ndims,pnew,p(:,n1+1:n1+sc_npt(i)),sc_mean(i,:), &
        sc_eval(i,:),sc_invcov(i,:,:),sc_kfac(i),sc_eff(i),sc_vol(i),d1)
    endif
    return
end subroutine evolveEllipsoid

!----------------------------------------------------------------------

subroutine addNewPoint(p, phyP, l, n1, nd_j, i, nd, sc_npt, ic_npt, lnew, &
        pnew, phyPnew)
    implicit none
    
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood
    
    integer n1, nd_j, i, nd
    
    integer, dimension(:), allocatable :: sc_npt
    integer, allocatable :: ic_npt(:)
    
    !new point
    double precision lnew
    double precision, allocatable :: pnew(:), phyPnew(:) ! new point
    
    !add the new point & its likelihood & update the no. of points 
    !in the cluster new point's index
    n1=n1+sc_npt(i)+1
    !make room
    p(:,n1+1:nd_j+ic_npt(nd))=p(:,n1:nd_j+ic_npt(nd)-1)
    phyP(:,n1+1:nd_j+ic_npt(nd))=phyP(:,n1:nd_j+ic_npt(nd)-1)
    l(n1+1:nd_j+ic_npt(nd))=l(n1:nd_j+ic_npt(nd)-1)
    !insert the new point
    p(:,n1)=pnew(:)
    phyP(:,n1)=phyPnew(:)
    l(n1)=lnew
    !increment no. of points in sub-cluster i
    sc_npt(i)=sc_npt(i)+1
    return
end subroutine addNewPoint

!----------------------------------------------------------------------

  !TODO: Check, if all these variables are necessary. Some might be for
  ! subroutines only.
  subroutine clusteredNest(p,phyP,l,loglike,dumper,context)
  	
    implicit none

    !input variables

    integer context
    double precision p(ndims,nlive+1) !live points
    double precision phyP(totPar,nlive+1) !physical live points
    double precision l(nlive+1) !log-likelihood

    !work variables

    !misc
    integer i, j, k, m, n, j1, i1, i2, i3, i4, ff, sff, n1, n2, q, nd, nd_i
    integer nd_j, iostatus, loopi
    integer num_old
    integer, allocatable :: eswitchff(:), escount(:), dmin(:)
    double precision d1, d2, d3, d4, d5, urv
    double precision h, logX, vprev, vnext, shrink !prior volume
    double precision mar_r !marginal acceptance rate
    double precision gZold !global evidence & info
    logical eswitch,peswitch !whether to do ellipsoidal sampling or not
    logical remFlag, acpt, flag, flag2
    integer funit1, funit2 !file units
    character(len=100) fName1, fName2 !file names
    character(len=100) fmt,fmt1
    
    !diagnostics for determining when to do eigen analysis
    integer neVol
    parameter(neVol=4)
    double precision, dimension(:), allocatable :: totVol, x1, x2, y1, y2, slope, intcpt, cVolFrac, pVolFrac
    double precision, dimension(:,:,:), allocatable :: eVolFrac

    !info for output file update
    double precision, allocatable :: evData(:,:), evDataAll(:), evDataTemp(:)

    !isolated cluster info
    integer ic_n !no. of nodes
    integer, allocatable :: ic_sc(:), ic_npt(:)
    logical, allocatable :: ic_done(:)
    integer, dimension(:), allocatable :: ic_fNode, ic_nsc, ic_nBrnch
    double precision, dimension(:,:,:),  allocatable :: ic_brnch, ic_llimits, ic_plimits
    double precision, allocatable :: ic_climits(:,:,:), ic_volFac(:)
    double precision, dimension(:),  allocatable :: ic_Z, ic_Zold, ic_info, ic_vnow, ic_hilike, ic_inc
    double precision, dimension(:,:),  allocatable :: ic_eff
    logical, dimension(:),  allocatable :: ic_reme, ic_rFlag, ic_chk
    logical modeFound

    !means & standard deviations of the live points (for prior edge detection)
    double precision, dimension(:,:), allocatable :: ic_mean, ic_sigma
    double precision lPts(ndims)

    !sub-cluster properties
    integer sc_n !no. of sub-clusters
    integer, dimension(:), allocatable :: sc_npt, nptk, nptx ,sc_node, nodek, sck
    double precision, dimension(:,:), allocatable :: meank, sc_eval, evalk
    double precision, dimension(:,:,:), allocatable :: sc_invcov, invcovk, sc_evec, eveck, tMatk
    double precision, dimension(:), allocatable :: kfack, volk, effk
    double precision, allocatable :: sc_mean(:,:), sc_tmat(:,:,:), sc_kfac(:), sc_eff(:), sc_vol(:)

    !auxiliary points (to be re-arranged with main points during clustering)
    integer naux !dimensionality of aux points 
    double precision, dimension(:,:), allocatable :: aux, pt

    !rejected point info
    double precision lowlike !lowest log-like
    double precision, allocatable, dimension(:) :: lowp, lowphyP !point with the lowlike
    integer indx(1) !point no. of lowlike

    !new point
    double precision lnew
    double precision, allocatable :: pnew(:), phyPnew(:) ! new point
    double precision, dimension(:,:), allocatable :: pnewa, phyPnewa
    double precision, dimension(:), allocatable :: lnewa
    integer, dimension(:), allocatable :: sEll

    !mode separation
    integer nCdim
    
    ! 1st approach for parallelizing
    double precision, allocatable, dimension(:) :: llhList
    logical llhFound
    double precision, dimension(:,:,:), allocatable :: pnewaList, phyPnewaList
    double precision, dimension(:,:), allocatable :: lnewaList
    double precision, dimension(:), allocatable :: lnewList
    logical stopLoop
    double precision, allocatable :: ic_climitsList(:,:,:,:)
    integer, dimension(:,:), allocatable :: sEllList
    
    logical, allocatable, dimension(:) :: acptList
    double precision, allocatable :: pNewList(:,:), phyPnewList(:,:) ! new point for every thread
    integer, allocatable, dimension(:) :: iList
    double precision, allocatable, dimension(:,:) :: lowpList
    integer, dimension(:,:), allocatable :: indxList
    
    integer :: counti, countf, count_rate
    real :: dt
                
    INTERFACE
        !the likelihood function
        subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
            integer n_dim,nPar,context_pass
            double precision lnew,Cube(nPar)
        end subroutine loglike
    end INTERFACE

    INTERFACE
            !the likelihood function
        subroutine loglikeDEBUG(Cube,n_dim,nPar,lnew,context_pass, id)
            integer n_dim,nPar,context_pass, id
            double precision lnew,Cube(nPar)
        end subroutine loglikeDEBUG
    end INTERFACE
    
    INTERFACE
        !the user dumper function
        subroutine dumper(nSamples, nlive, nPar, physLive, posterior, paramConstr, maxLogLike, logZ, logZerr, context_pass)
            integer nSamples, nlive, nPar, context_pass
            double precision, pointer :: physLive(:,:), posterior(:,:), paramConstr(:)
            double precision maxLogLike, logZ, logZerr
        end subroutine dumper
    end INTERFACE
    allocate( eswitchff(maxCls), escount(maxCls), dmin(maxCls) )
    allocate( evData(updInt,totPar+3) )
    allocate( ic_sc(maxCls), ic_npt(maxCls) )
    allocate( ic_done(0:maxCls) )
    allocate( ic_climits(maxCls,ndims,2), ic_volFac(maxCls) )
    allocate( sc_mean(maxeCls,ndims), sc_tmat(maxeCls,ndims,ndims), sc_kfac(maxeCls), sc_eff(maxeCls), sc_vol(maxeCls) )
    allocate( lowp(ndims), lowphyP(totPar) )
    allocate( pnew(ndims), phyPnew(totPar) )
    allocate( pNewList(procs, ndims), phyPnewList(procs, totPar) )
    allocate( lnewList(procs), llhList(procs) )
    allocate( ic_climitsList(procs, maxCls,ndims,2) )
    allocate( acptList(procs) )
    allocate( lowpList(procs, ndims) )
    allocate( iList(procs), indxList(procs, 1) )

    
    !initializations
    ic_done=.false.
    ic_npt=nlive
    ic_climits(:,:,2)=1d0
    ic_climits(:,:,1)=0d0
    ic_volFac(:)=1d0
	
	
    ! Allocate memory and initialize variables. If there is a resume, read
    ! the resumefile and find the highest Likelihood for each node.

    !memory allocation
    allocate(evDataAll(1))
    allocate(sc_npt(maxeCls), nptk(maxeCls), nptx(nlive), meank(maxeCls,ndims), &
    sc_eval(maxeCls,ndims), evalk(maxeCls,ndims), sc_invcov(maxeCls,ndims,ndims), &
    invcovk(maxeCls,ndims,ndims), tMatk(maxeCls,ndims,ndims), &
    sc_evec(maxeCls,ndims,ndims), eveck(maxeCls,ndims,ndims), kfack(maxeCls), &
    effk(maxeCls), volk(maxeCls), &
    sc_node(maxeCls),nodek(maxeCls),sck(maxCls))
    allocate(pt(ndims,nlive), aux(ndims+totPar+4-nCdims,nlive))
    allocate(ic_fNode(maxCls),ic_nsc(maxCls),ic_nBrnch(maxCls), &
    ic_brnch(maxCls,maxCls,2),ic_reme(maxCls),ic_rFlag(maxCls),ic_z(maxCls),ic_zold(maxCls),ic_info(maxCls), &
    ic_vnow(maxCls),ic_hilike(maxCls),ic_inc(maxCls),ic_chk(maxCls),ic_llimits(maxCls,ndims,2))
    allocate(ic_plimits(maxCls,nCdims,2))
    allocate(pnewa(maxCls,ndims), &
        phyPnewa(maxCls,totPar), &
        lnewa(maxCls), &
        sEll(maxCls))
    allocate(pnewaList(procs, maxCls, ndims),&
        phyPnewaList(procs, maxCls, totPar), &
        lnewaList(procs, maxCls) )
    allocate(totVol(maxCls),eVolFrac(maxCls,neVol*2,neVol),x1(maxCls),x2(maxCls),y1(maxCls), &
        y2(maxCls),slope(maxCls),intcpt(maxCls),cVolFrac(maxCls),pVolFrac(maxCls))
    if( prior_warning ) allocate( ic_mean(maxCls, ndims), ic_sigma(maxCls, ndims) )

    !global logZ = log(0)
    gZ=logZero
    gZOld=logZero
    ic_Z=logZero
    ic_zold=logZero
    ginfo=0.d0
    ic_info=0.d0
    
    ic_llimits(:,:,2)=1d0
    ic_llimits(:,:,1)=0d0
    ic_plimits(:,:,2)=huge(1d0)
    ic_plimits(:,:,1)=-huge(1d0)
    
    !just one node to start with
    ic_n=1
    ic_fNode(1)=0
    ic_sc(1)=0
    ic_nsc=nsc_def
    ic_reme=.false.
    ic_nBrnch=0

    !set the prior volume
    ic_vnow=0.d0
    ic_vnow(1)=1.d0

    !no ellipsoidal sampling to start with
    eswitch=.false.
    peswitch=.false.
    escount=0
    
    !no sub-clusters at the beginning
    sc_n=0
    totVol=1.d0
    sc_node=1

    !eigen analysis diagnostics
    eVolFrac=0.d0

    dmin=ndims+1
    ic_rFlag=.false.
    sff=0

    if(resumeFlag) then
        call readResumeFile(funit1, peswitch, phyP, l, ic_n, ic_npt, &
                ic_done, ic_fNode, ic_nBrnch, ic_brnch, ic_Z, ic_info, &
                ic_vnow, ic_eff, ic_reme, ic_rFlag, eswitch, eswitchff, totVol, p)
    endif

    call getHighestLLHperNode(ic_n, ic_npt, ic_done, ic_llimits, &
            ic_plimits, ic_climits, ic_Z, ic_vnow, ic_hilike, ic_inc, l, &
            phyP, lowlike, p)
    do ff=1,maxIter
        !stopping condition reached
        if(ic_done(0)) then
            if(outfile) then
                call writeResumeFile(fName1, ic_n, eswitch, ic_nBrnch, &
                    ic_brnch, ic_done, ic_reme, ic_fNode, ic_npt, ic_vnow, &
                    ic_Z, ic_info, ic_eff)
            endif
            ! Add the contribution to the global evidence and stuff.
            !fback
            if(fback) call gfeedback(gZ,numlike,globff,.false.)
            call pos_samp(Ztol,globff,broot,nlive,ndims,nCdims,totPar,&
                outfile,gZ,ginfo,ic_n,ic_Z(1:ic_n), &
                ic_info(1:ic_n),ic_reme(1:ic_n),ic_vnow(1:ic_n),&
                ic_npt(1:ic_n),ic_nBrnch(1:ic_n),ic_brnch(1:ic_n,:,1),&
                phyP(:,1:nlive),l(1:nlive),evDataAll,dumper,context)
            !if done then add in the contribution to the global evidence 
            !from live points
            call addContribution(ic_n, ic_vnow, ic_npt, logX, l, ic_Zold, &
                ic_Z, ic_info, gZold)
            !memory deallocation
            deallocate(sc_npt, sc_invcov, sc_node, sc_eval, sc_evec)
            deallocate(nptk, nptx, meank, evalk, invcovk, tMatk, eveck, kfack, effk, volk, &
                nodek ,sck)
            deallocate(pt, aux)
            deallocate(ic_fNode,ic_nsc,ic_nBrnch,ic_brnch,ic_reme,ic_rFlag,ic_Z,ic_zold, &
                ic_info,ic_vnow,ic_hilike,ic_inc,ic_chk,ic_llimits)
            deallocate(ic_plimits)
            deallocate(pnewa,phyPnewa,lnewa,sEll)
            deallocate(totVol,eVolFrac,x1,x2,y1,y2,slope,intcpt,cVolFrac,pVolFrac)
            if( prior_warning ) deallocate( ic_mean, ic_sigma )
            deallocate( evDataAll )
            deallocate(phyPnewaList, phyPnewList, pNewaList, pNewList, llhlist, &
                lnewaList, lnewList)
            deallocate( eswitchff, escount, dmin, evData, ic_sc, ic_npt, ic_done, ic_climits, &
            ic_volFac, sc_mean, sc_tmat, sc_kfac, sc_eff, sc_vol, lowp, lowphyP, pnew, phyPnew )
            deallocate( ic_climitsList )
            totalIterations = ff
            return
        endif
    
        modeFound=.false.
        
        !mode separation
        if(ff/=1 .and. eswitch .and. mod(ff,15)==0 .and. ic_n<maxCls) then
            call modeSeparation(l, phyP, p, ic_n, ic_npt, ic_done, &
                    ic_reme, ic_vnow, ic_rFlag, ic_climits, ic_sc, sc_npt, &
                    sc_n, ic_nsc, ic_fNode, ic_nBrnch, ic_brnch, nptk, &
                    ic_chk, modeFound, nodek, q, k, ic_Zold, ic_Z, &
                    ic_info, ic_hilike, ic_eff, eVolFrac, ic_plimits, &
                    ic_llimits, ic_volFac, ic_inc, aux, naux)
        endif
        
        num_old=numlike
    
        !ellipsoidal decomposition
        flag2=.false.
        ! Check, if one (isolated) cluster can be divided into more 
        ! (sub) clusters. If yes, then the number of subclusters sc_n > 0
        ! and we let eswitch be true.
        if(peswitch) then
            call subclustering(p, phyP, l, eswitch, neVol, &
                totVol, x1, x2, y1, y2, slope, intcpt, cVolFrac, pVolFrac, &
                eVolFrac, ic_n, ic_sc, ic_npt, ic_done, ic_nsc, &
                ic_llimits, ic_climits, ic_volFac, ic_vnow, ic_hilike, &
                ic_eff, ic_rFlag, sc_n, sc_npt, nptk, sc_node, nodek, sck, &
                meank, sc_eval, evalk, sc_invcov, invcovk, sc_evec, eveck, &
                tMatk, kfack, volk, effk, sc_mean, sc_tMat, sc_kfac, &
                sc_eff, sc_vol, naux, aux, pt, lowlike, indx, &
                flag, flag2, eswitchff, ff, &
                n, d1, d4)
        endif
        ic_rFlag(1:ic_n)=.false.
        if(sc_n==0) eswitch=.false.

        nd_i=0 !no. of ellipsoids traversed
        nd_j=0 !no. of points traversed
        do nd=1,ic_n
        
            !Check if any of the clusters is done. If one of the clusters is
            !not done, then the algorithm isn't finished <=> ic_done(0) = false.
            if(ic_npt(nd)<ndims+1 .or. (ic_sc(nd)==1 .and. sc_vol(nd_i+1)==0.d0)) then
                ic_done(nd)=.true.
                ic_done(0)=.true.
                do i3=1,ic_n
                    if(.not.ic_done(i3)) then
                        ic_done(0)=.false.
                        exit
                    endif
                enddo
            endif
            ! Add the number of points and subclusters to nd_j and nd_i if
            ! the isolated cluster is done and not all clusters are done.
            if(ic_done(nd) .and. .not.ic_done(0)) then
                nd_i=nd_i+ic_sc(nd)
                nd_j=nd_j+ic_npt(nd)
                cycle
            endif
            ! For every cluster, which is not finished, adjust the prio volumes.
            ! Then find the point with the lowest likelihood and set the 
            ! corresponding limits.
            if(.not.ic_done(nd)) then
                !adjust the prior volumes & inc
                shrink=exp(-1.d0/dble(ic_npt(nd)))
                call adjustPriorVol(p, phyP, l, i3, nd, nd_j, d1, d4, h, &
                    vprev, vnext, shrink, ic_npt, ic_llimits, ic_plimits, &
                    ic_climits, ic_Z, ic_vnow, &
                    ic_hilike, ic_inc, lowlike, lowp, lowPhyP, indx, pt)
                !now find a new point inside the hard constraint
                if(.not.eswitch) then
                    call findPointHardConstr(context, p, phyP, l, i, n1, n2, nd, i2, i3, n, eswitch, &
                        remFlag, acpt, ic_llimits, ic_plimits, ic_climits, ic_hilike, &
                        sc_mean, sc_tMat, lowlike, indx, lnew, pnew, phyPnew, pnewaList, phyPnewaList, &
                        lnewaList, loglike, ic_npt, nd_j)
                else
                    call findPointHardConstrEllipsoids(context, p, phyP, l, i, &
                        j, m, n, j1, nd, nd_i, num_old, d1, urv, eswitch, &
                        acpt, ic_sc, ic_llimits, ic_plimits, ic_climits, &
                        ic_hilike, sc_npt, sc_invcov, sc_mean, sc_tMat, &
                        sc_kfac, sc_eff, sc_vol, lowlike, lnew, pnew, phyPnew, &
                        pnewaList, phyPnewaList, lnewaList, loglike, &
                        ic_done, iList, lnewList, pnewList, &
                        phyPnewList, indxList, lowpList, lowp, nd_j, indx, ic_npt, acptList)
                    do loopi=1, procs
                        lowp(:) = lowpList(loopi, :)
                        pnew(:) = pnewList(loopi, :)
                        phyPnew(:) = phyPnewList(loopi, :)
!                         write(*,*) 'Evolve ellipsoid'
                        call evolveEllipsoid(p, phyP, l, iList(loopi), j, m, i1, n1, n2, q, &
                            nd, nd_i, nd_j, d1, d4, shrink, ic_sc, ic_npt, &
                            ic_fNode, ic_volFac, ic_vnow, ic_inc, ic_eff, sc_npt, &
                            sc_eval, sc_invcov, sc_mean, sc_kfac, sc_eff, sc_vol, &
                            lowp, indxList(loopi, :), pnew)
                    
                        !now evolve the rest of the sub-clusters
                        !only if sub-clustering wasn't done in the current iteration
                        do j=nd_i+1,nd_i+ic_sc(nd)
                            !sub-clusters with rejected/inserted point already evolved
                            if(j==iList(loopi) .or. j==q) cycle
                            if(sc_eff(j)>1.d0 .and. sc_vol(j)>0.d0 .and. sc_kfac(iList(1))>0.d0) then
                                d1=sc_eff(j)
                                sc_eff(j)=max(1.d0,sc_eff(j)*(shrink**(2.d0/dble(ndims))))
                                d1=sc_eff(j)/d1
                                sc_vol(j)=sc_vol(j)*(d1**(dble(ndims)/2.d0))
                            endif						
                        enddo
                        call addNewPoint(p, phyP, l, n1, nd_j, iList(loopi), nd, &
                            sc_npt, ic_npt, lnewList(loopi), pnew, &
                            phyPnew)
                    enddo
                endif
                !update evidence, info, prior vol, sampling statsitics
                globff=globff+1
                sff=sff+1
                gZold=gZ
                ic_zold(nd)=ic_z(nd)
                d1=lowlike+log(h)
                gZ=LogSumExp(gZ,d1)
                ic_Z(nd)=LogSumExp(ic_Z(nd),d1)
                ginfo=ginfo*exp(gzold-gz)+exp(d1-gz)*lowlike
                ic_info(nd)=ic_info(nd)*exp(ic_zold(nd)-ic_z(nd))+exp(d1-ic_z(nd))*lowlike
    
                !data for ev.dat file
                j1=mod(sff-1,updInt)+1
                evData(j1,1:totPar)=lowPhyP(1:totPar)
                evData(j1,totPar+1)=lowlike
                evData(j1,totPar+2)=log(h)
                evData(j1,totPar+3)=dble(nd)
                
                lowlike=minval(l(nd_j+1:nd_j+ic_npt(nd)))
            
                !Check if one of the terminating conditions is reached.
                if(((abs(lowlike-ic_hilike(nd))<= 0.0001 .or. (ic_inc(nd)<log(tol) .and. &
                    globff-nlive>50)) .and. minIter<ff) .or. ff==maxIter) then
                    ic_done(nd)=.true.
                        
                    !Why do we check this here? We could do this out of 
                    ! the do-loop I believe.
                    !check if all done
                    ic_done(0)=.true.
                    do i=1,ic_n
                        if(.not.ic_done(i)) then
                            ic_done(0)=.false.
                            exit
                        endif
                    enddo
                endif
            endif
    
            if(sff>0 .and. (mod(sff,updInt)==0 .or. ic_done(0))) then
            
                if( .not.outfile ) then
                    k=0
                    if( sff > updInt ) then
                        k=size(evDataAll)
                        allocate( evDataTemp(k) )
                        evDataTemp=evDataAll
                        deallocate( evDataAll )
                        allocate( evDataAll(k+j1*(totPar+3)) )
                        evDataAll(1:k)=evDataTemp(1:k)
                        deallocate( evDataTemp )
                    else
                        deallocate( evDataAll )
                        allocate( evDataAll(j1*(totPar+3)) )
                    endif
                    do i=1,j1
                        evDataAll(k+1:k+totPar+3) = evData(i,1:totPar+3)
                        k=k+totPar+3
                    enddo
                else
                    !write the evidence file
                    funit1=u_ev
                    fName1=evname
                    open(unit=funit1,file=fName1,form='formatted',status='old', position='append')
                    write(fmt,'(a,i5,a)')  '(',totPar+2,'E28.18,i5)'	
                    do i=1,j1
                        write(funit1,fmt) evData(i,1:totPar+2),int(evData(i,totPar+3))
                    enddo
                    !close the files
                    close(funit1)
                
                    !write the live file
                    funit1=u_phys
                    fName1=physname
                    funit2=u_live
                    fName2=livename
                    open(unit=funit1,file=fName1,form='formatted',status='replace')
                    open(unit=funit2,file=fName2,form='formatted',status='replace')
                        
                    write(fmt,'(a,i5,a)')  '(',totPar+1,'E28.18,i4)'
                    write(fmt1,'(a,i5,a)')  '(',ndims+1,'E28.18)'
                    k=0
                    do i=1,ic_n
                        do j=1,ic_npt(i)
                            k=k+1
                            write(funit1,fmt) phyP(1:totPar,k),l(k),i
                            lPts(1:ndims) =ic_climits(i,1:ndims,1)+ &
                                (ic_climits(i,1:ndims,2)-ic_climits(i,1:ndims,&
                                1))*p(1:ndims,k)
                            write(funit2,fmt1) lPts(1:ndims),l(k)
                        enddo
                    enddo
                    !close the files
                    close(funit1)
                    close(funit2)
                    
                    !write the resume file
                    funit1=u_resume
                    fName1=resumename
                    open(unit=funit1,file=fName1,form='formatted',status='replace')
                        
                    write(funit1,'(l2)')genLive
                    write(funit1,'(4i12)')globff,numlike,ic_n,nlive
                    write(funit1,'(2E28.18)')gZ,ginfo
                    write(funit1,'(l2)')eswitch
                
                    !write branching info
                    do i=1,ic_n
                        write(funit1,'(i4)')ic_nBrnch(i)
                        if(ic_nBrnch(i)>0) then
                            write(fmt,'(a,i5,a)')  '(',2*ic_nBrnch(i),'E28.18)'
                            write(funit1,fmt)ic_brnch(i,1:ic_nBrnch(i),1), ic_brnch(i,1:ic_nBrnch(i),2)
                        endif
                    enddo
            
                    !write the node info
                    do i=1,ic_n
                        write(funit1,'(2l2,2i6)')ic_done(i),ic_reme(i),ic_fNode(i), ic_npt(i)
                        write(funit1,'(3E28.18)')ic_vnow(i),ic_Z(i),ic_info(i)
                    enddo
                    close(funit1)
                    
                    
                endif
                
                !check if the parameters are close the prior edges
                if( .not.ic_done(0) .and. prior_warning .and. mod(ff,50)== 0 ) then
                    flag = .false.
                    k=0
                    do i=1,ic_n
                        if( ic_npt(i) == 0 .or. ic_done(i) ) cycle
                            
                        ic_mean(i,1:ndims) = 0d0
                        ic_sigma(i,1:ndims) = 0d0
                        
                        do j=1,ic_npt(i)
                            k=k+1
                            lPts(1:ndims) = ic_climits(i,1:ndims,1)+&
                                (ic_climits(i,1:ndims,2)-ic_climits(i,1:ndims,&
                                1))*p(1:ndims,k)
                            if( prior_warning .and. mod(ff,50)== 0 ) then
                                ic_mean(i,1:ndims) = ic_mean(i,1:ndims) + lPts(1:ndims)
                                ic_sigma(i,1:ndims) = ic_sigma(i,1:ndims) + lPts(1:ndims) * lPts(1:ndims)
                            endif
                        enddo
                        
                        ic_mean(i,1:ndims) = ic_mean(i,1:ndims) / dble(ic_npt(i))
                        ic_sigma(i,1:ndims) = sqrt(max(0d0, ic_sigma(i,1:ndims) &
                            / dble(ic_npt(i)) + ic_mean(i,1:ndims) &
                            * ic_mean(i,1:ndims)))
                        
                        do j = 1, ndims
                            if( ic_sigma(i,j) <= 0.05 .and. ( ic_sigma(i,j) <= 0.05 .or. ic_sigma(i,j) >= 0.95 ) ) then
                                if( .not. flag ) then
                                    write(*,*)
                                    write(*,*)"MultiNest Warning!"
                                    flag = .true.
                                endif
                                write(*,*)"Parameter ", j, " of mode ", i, " is converging towards the edge of the prior."
                            endif
                        enddo
                    enddo
                endif
                
                if(mod(sff,updInt*10)==0 .or. ic_done(0)) then
                    call pos_samp(Ztol,globff,broot,nlive,ndims,nCdims,totPar, &
                        outfile,gZ,ginfo,ic_n,ic_Z(1:ic_n),ic_info(1:ic_n),ic_reme(1:ic_n),ic_vnow(1:ic_n),ic_npt(1:ic_n), &
                        ic_nBrnch(1:ic_n),ic_brnch(1:ic_n,:,1),phyP(:,1:nlive),l(1:nlive),evDataAll,dumper,context)
                endif
            endif
            nd_i=nd_i+ic_sc(nd)
            
            !prior volume shrinks
            if(.not.ic_done(nd)) ic_vnow(nd)=ic_vnow(nd)*shrink
            nd_j=nd_j+ic_npt(nd)
            if( ic_done(0) ) exit
        enddo
        
        !update the total volume
        if(eswitch) then
            k = 0
            do i = 1, ic_n
                if( ic_done(i) ) then
                    totVol(i) = 0d0
                else
                    d1 = totVol(i)
                    totVol(i) = sum(sc_vol(k+1:k+ic_sc(i)))
                endif
                k = k + ic_sc(i)
            enddo
        endif
        
        !calculate the global evidence & info
        if(mod(ff,50)==0) then
            
            if(fback) then
                call gfeedback(gZ,numlike,globff,.false.)
            
                if(debug) then
                    d1=0.d0
                    d2=0.d0
                    j=0
                    do i=1,ic_n
                        if(ic_done(i)) then
                            j=j+ic_sc(i)
                            cycle
                        endif
                        d1=d1+ic_vnow(i)*ic_volFac(i)
                        d2=d2+sum(sc_vol(j+1:j+ic_sc(i)))
                        j=j+ic_sc(i)
                    enddo
                    write(*,*)ic_n, sc_n, d2/d1, scount
                endif
            endif
        endif
        
        if(.not.peswitch) then
            !marginal acceptance rate
            if(.not.peswitch) mar_r=1.d0/dble(numlike-num_old)
!             write(*,*) 'ef', ef, ' mar_r', mar_r
            if(mar_r<ef/2) then
                escount(1)=escount(1)+1
                if(escount(1)==5) then
                    peswitch=.true.
                    eswitchff(1)=ff
                    totVol=1.d0
                    escount(1)=0
                endif
            else
                escount(1)=0
            endif
        endif
    enddo

    deallocate( eswitchff, escount, dmin, evData, ic_sc, ic_npt, ic_done, &
    ic_climits, ic_volFac, sc_mean, sc_tmat, sc_kfac, sc_eff, sc_vol, lowp, lowphyP, pnew, phyPnew )
    
    

  end subroutine clusteredNest
  
!----------------------------------------------------------------------
 
 logical function ellIntersect(ndim,mean1,mean2,eval1,evec1,ef1,ef2,inv_cov1,inv_cov2)
 	
      implicit none
      
      integer ndim !dimensionality
      double precision eval1(:), evec1(:,:), mean1(:), mean2(:), inv_cov1(:,:), inv_cov2(:,:)
      ! matrices for ellipsoid interaction detection
      double precision, allocatable :: matA(:,:), matB(:,:), matR(:,:), matAinvB(:,:)
      ! variables for calculating eigenvalues of N*N real non-sym matrix
      double precision, allocatable :: evalR(:), evalI(:), VL(:,:), VR(:,:), WORK(:), delMean(:)
      ! effective enlargement factor
      double precision ef1,ef2
      integer inf,k,i
      integer i1,i2,i3,i4,i5
      
      
      
	if( ef1 == 0.d0 .and. ef2 == 0.d0 ) then
      		ellIntersect=.false.
            	return
      	else if( ef1 == 0.d0 ) then
      		ellIntersect = ptIn1Ell(ndim, mean1, mean2, inv_cov2, ef2)
            	return
	else if( ef2==0.d0 ) then
      		ellIntersect = ptIn1Ell(ndim, mean2, mean1, inv_cov1, ef1)
            	return
	endif
	
	allocate( matA(ndim+1,ndim+1), matB(ndim+1,ndim+1), matR(ndim+1,ndim+1), matAinvB(ndim+1,ndim+1) )
	allocate( evalR(ndim+1),evalI(ndim+1),VL(ndim+1,ndim+1),VR(ndim+1,ndim+1), WORK(4*ndim+4), delMean(ndim) )
      
      	delMean(1:ndim) = mean1(1:ndim) - mean2(1:ndim)
      	matA = 0.d0

	do i = 1, ndim
		matA(i,i)=eval1(i)
		matB(i,1:ndim)=inv_cov2(i,1:ndim)
		matB(ndim+1,i)=sum(-delMean(:)*inv_cov2(:,i))
		matB(i,ndim+1)=matB(ndim+1,i)
      	enddo
      	
	matA(ndim+1,ndim+1)=-1.d0/(ef1)
      	matB(ndim+1,ndim+1)=sum(matB(ndim+1,1:ndim)*(-delMean(:)))-ef2
      	matR=0.d0
      	
	do i=1,ndim
		matR(i,1:ndim)=evec1(:,i)
	enddo
      
      	matR(ndim+1,ndim+1)=1.d0
      	matB=MATMUL(MATMUL(matR,matB),TRANSPOSE(matR))
      	matAinvB=MATMUL(matA,matB)
      	i1=ndim+1
      	i2=ndim+1
      	i3=ndim+1
      	i4=ndim+1
      	i5=4*ndim+4
      	call DGEEV('N','N',i1,matAinvB,i2,evalR,evalI,VL,i3,VR,i4,WORK,i5,INF)
      	k=0
      
      	do i=1,ndim+1
		if(evalI(i)/=0.d0) then
			ellIntersect=.true.
			deallocate( matA, matB, matR, matAinvB, evalR, evalI, VL,VR, WORK, delMean )
			return
		else if(evalR(i)<0.d0) then
            		k=k+1
		endif
      	enddo
      
      	if(k<2) then
		ellIntersect=.true.
		deallocate( matA, matB, matR, matAinvB, evalR, evalI, VL,VR, WORK, delMean )
		return
      	endif
            
      	ellIntersect=.false.
	deallocate( matA, matB, matR, matAinvB, evalR, evalI, VL,VR, WORK, delMean )
 end function ellIntersect

!---------------------------------------------------------------------- 

!sample a point inside the given ellipsoid with log-likelihood>lboundary
subroutine samp(pnew,phyPnew,lnew,mean,ekfac,TMat,limits,loglike,eswitch,lowlike,n,context, id)

    implicit none
    double precision lnew
    double precision pnew(ndims),spnew(ndims),phyPnew(totPar),ekfac
    double precision mean(ndims),TMat(ndims,ndims)
    double precision limits(ndims,2)
    double precision lowlike	!likelihood threshold
    integer n			!no. of points drawn
    logical eswitch
    integer id,i,context
    
    INTERFACE
        !the likelihood function
        subroutine loglike(Cube,n_dim,nPar,lnew,context_pass)
        integer n_dim,nPar,context_pass
        double precision lnew,Cube(nPar)
    end subroutine loglike
    end INTERFACE
    

    n = 0
    do
        
        if(.not.eswitch) then
            !generate a random point inside unit hypercube
            call getrandom(ndims,pnew(1:ndims),id)
            phyPnew(1:ndims)=pnew(1:ndims)
            lnew=lowlike
            call loglike(phyPnew,ndims,totPar,lnew,context)
            n = n+1
        else
            !generate a point uniformly inside the given ellipsoid
            call genPtInEll(ndims,mean,ekfac,TMat,id,pnew(1:ndims))
            spnew(:)=limits(:,1)+(limits(:,2)-limits(:,1))*pnew(:)
            do i=1,ndims
                if(pWrap(i)) then
                    call wraparound(spnew(i),spnew(i))
                endif
            enddo
            if(.not.inprior(ndims,spnew(1:ndims))) cycle
            phyPnew(1:ndims)=spnew(1:ndims)
            lnew=lowlike
            call loglike(phyPnew,ndims,totPar,lnew,context)
            n = n+1
        endif
        if(lnew>logZero) exit
    enddo
        
  end subroutine samp
  
  !----------------------------------------------------------------------
  
  subroutine ApplyLimits(flag,limits,pt,transpt)

    implicit none

    !input parameters
    integer flag !0 => given point in unit hypercube, apply the limits
                 !1 => point to be transformed to unit hypercube, reverse the limits
    double precision limits(ndims,2)	!current limits
    double precision pt(ndims)		!point

    !output variables
    double precision transpt(ndims)		!final result


    if( flag == 0 ) then
        !apply the limits to point in unit hypercube
        transpt(:)=(pt(:)-limits(:,1))/(limits(:,2)-limits(:,1))
    elseif( flag == 1 ) then
        !reverse the limits such that the point is in unit hypercube
        transpt(:)=limits(:,1)+(limits(:,2)-limits(:,1))*pt(:)
    else
        write(*,*)'Incorrect value of flag passed to ApplyLimits'
        stop
    endif

  end subroutine ApplyLimits
  
  !----------------------------------------------------------------------
 
  subroutine selectEll(n,volx,sEll)

    implicit none
    
    integer n!total no. of ellipsoids
    double precision volx(n)!no. points in each ellipsoid
    integer sEll!answer, the ellipsoid to sample from
    double precision volTree(n),totvol,urv
    integer i
 
    totvol=0.d0
    volTree=0.d0
    do i=1,n
        volTree(i)=totvol+volx(i)
        totvol=volTree(i)
    enddo
    volTree=volTree/totvol
    urv=ranmarNS(0)
    do i=1,n
        if(urv<=volTree(i)) exit
    enddo
    sEll=i
  end subroutine selectEll
  
!----------------------------------------------------------------------
   
   !provide fback to the user
  subroutine gfeedback(logZ,nlike,nacc,dswitch)
    
    implicit none
    !input variables
    double precision logZ !log-evidence
    integer nlike !no. of likelihood evaluations
    integer nacc !no. of accepted samples
    logical dswitch !dynamic live points

    write(*,'(a,F14.6)')	     'Acceptance Rate:                  ',dble(nacc)/dble(nlike)
    write(*,'(a,i14)')   	     'Replacements:                     ',nacc
    write(*,'(a,i14)')   	     'Total Samples:                    ',nlike
    write(*,'(a,F14.6)')	     'Nested Sampling ln(Z):            ',logZ
    if(dswitch) write(*,'(a,i5)')'Total No. of Live Points:         ',nlive
    
  end subroutine gfeedback
  
!----------------------------------------------------------------------
 
 !isolate the modes
 logical function isolateModes2(npt,ndim,nCdim,pt,naux,aux,ic_n,ic_fnode,ic_npt,ic_reme,ic_chk,ic_vnow,reCluster,limits)
    implicit none

    !input parameters
    integer npt !total no. of points
    integer ndim !dimensionality
    integer nCdim !clustering dimensionality
    double precision pt(nCdim,npt)
    integer naux !no. of dimensions of the aux array
    double precision aux(naux,npt)
    logical ic_chk(ic_n) !whether to check a node for mode separation
    double precision ic_vnow(ic_n) !prior volumes
    double precision limits(maxCls,nCdim,2) !physical parameter ranges

    !input/output parameters
    !mode info
    integer ic_n !input: no. of existing modes. output: updated no. of existing modes
    integer ic_fnode(maxCls),ic_npt(maxCls)
    logical ic_reme(maxCls)

    !output variables
    logical reCluster(maxCls)

    !work variables
    integer i,j,k,i1,j2,j3,j4,i2,i3,i4,i5,n,n1,n2,n_mode,npt_mode,m,nLost
    integer, allocatable :: order(:), nptx(:), nodex(:)
    double precision d1,d2,d4,ef0, ef1, ef2
    integer nN,sc_n
    logical, allocatable :: gList(:), lList(:), toBeChkd(:), overlapk(:,:)
    logical flag,intFlag
    double precision, allocatable :: ptk(:,:), auxk(:,:), ptx(:,:), auxx(:,:)
    double precision, allocatable :: mMean(:), lMean(:), mStdErr(:), lStdErr(:)
    double precision, allocatable :: mean1(:), mean2(:), mean1w(:), mean2w(:)
    double precision, allocatable :: eval1(:), evec1(:,:), invcov1(:,:), invcov2(:,:)
    logical, allocatable :: wrapEll(:), wrapDim(:,:,:)
    integer, allocatable :: wrapN(:)
    double precision, allocatable :: meanw(:,:),meank(:,:),evalk(:,:),eveck(:,:,:),invcovk(:,:,:),tmatk(:,:,:),kfack(:)


    allocate( order(nCdim), nptx(npt/(ndim+1)+1), nodex(npt/(ndim+1)+1) )
    allocate( gList(npt/(ndim+1)+1), lList(npt/(ndim+1)+1), &
        toBeChkd(npt/(ndim+1)+1), overlapk(npt/(ndim+1)+1,npt/(ndim+1)+1) )
    allocate( ptk(nCdim,npt), auxk(naux,npt), ptx(nCdim,npt), auxx(naux,npt), mMean(nCdim), lMean(nCdim), mStdErr(nCdim), &
        lStdErr(nCdim), mean1(nCdim), mean2(nCdim), mean1w(nCdim), mean2w(nCdim), eval1(nCdim), evec1(nCdim,nCdim), &
        invcov1(nCdim,nCdim), invcov2(nCdim,nCdim) )

    nN=ic_n
    isolateModes2=.false.
    reCluster=.false.

    i1=0 !no. of points traversed
    sc_n=0 !no. of clusters created so far

    do i=1,nN
        !enlargement factor
        ef0 = max( ( ic_vnow(i) * 100d0 / 111d0 ) + ( 11d0 / 111d0 ) , 0.8d0 )

        i3=ic_n
        
        if(ic_npt(i)==0) cycle
        
        if(.not.ic_chk(i)) then
            do j=1,nCdim
                ptk(j,i1+1:i1+ic_npt(i))=pt(j,i1+1:i1+ic_npt(i))
            enddo
            auxk(1:naux,i1+1:i1+ic_npt(i))=aux(1:naux,i1+1:i1+ic_npt(i))
            nptx(sc_n+1)=ic_npt(i)
            nodex(sc_n+1)=i
            i1=i1+ic_npt(i)
            sc_n=sc_n+1
            cycle
        endif
        
        nLost=0
        overlapk=.true.
        gList=.false.
        
        !cluster using G-means
        do j=1,nCdim
            d4=limits(i,j,2)-limits(i,j,1)
            ptk(j,i1+1:i1+ic_npt(i))=pt(j,i1+1:i1+ic_npt(i))
            
            if(pWrap(j)) then
                do i2=i1+1,i1+ic_npt(i)
                    !scale to unit hypercube
                    d2=limits(i,j,1)+d4*ptk(j,i2)
                    
                    call wraparound(d2,d1)
                    
                    !scale back to the limits
                    if(d1 /= d2) ptk(j,i2)=(d1-limits(i,j,1))/d4
                enddo
            endif
        enddo
        
        auxk(1:naux,i1+1:i1+ic_npt(i))=aux(1:naux,i1+1:i1+ic_npt(i))
        n1=max(nCdim+1,3) !min no. of points allowed in a cluster
        n2=ic_npt(i)/n1+1 !max no. of clusters possible
        
        call doGmeans(ptk(1:nCdim,i1+1:i1+ic_npt(i)),ic_npt(i),nCdim,k,nptx(sc_n+1:sc_n+n2), &
            naux,auxk(:,i1+1:i1+ic_npt(i)),n1,n2)

        allocate(wrapN(k),wrapEll(k),wrapDim(k,nCdim,2),meank(k,nCdim),evalk(k,nCdim), &
            eveck(k,nCdim,nCdim),invcovk(k,nCdim,nCdim),tmatk(k,nCdim,nCdim),kfack(k),meanw(k,nCdim))
        wrapEll=.false.
        wrapDim=.false.
        
        nodex(sc_n+1:sc_n+k)=i

        !enclose sub-clusters in bounding ellipsoids
        n=i1
        wrapN=1
        do j=1,k
            call CalcBEllInfo(nptx(sc_n+j),nCdim,ptk(1:nCdim,n+1:n+nptx(sc_n+j)),meank(j,:), &
            evalk(j,:),eveck(j,:,:),invcovk(j,:,:),tmatk(j,:,:),kfack(j),n1)
            
            if(mWrap) then
                ef1=kfack(j)*max(2d0,((1.d0+ef0*sqrt(50.d0/dble(nptx(sc_n+j))))**(1d0/nCdim)))
                call wrapEllCheck(nCdim,meank(j,:),tmatk(j,:,:),ef1,limits(i,:,:), &
                wrapEll(j),wrapDim(j,:,:))
                
                !transform the mean to hypercube
                call Scaled2Cube(nCdim,limits(j,:,:),meank(j,:),meanw(j,:))
                
                if(wrapEll(j)) then
                    do i2=1,nCdim
                        if(wrapDim(j,i2,1)) then
                            wrapN(j)=wrapN(j)+1
                            meanw(j,i2)=1d0+meanw(j,i2)
                        elseif(wrapDim(j,i2,2)) then
                            wrapN(j)=wrapN(j)+1
                            meanw(j,i2)=meanw(j,i2)-1d0
                        endif
                    enddo
                endif
                
                !scale the mean
                call Cube2Scaled(nCdim,limits(j,:,:),meanw(j,:),meanw(j,:))
            endif

            n=n+nptx(sc_n+j)
        enddo
        
        !calculate the standard error, required for localization
        if(.not.aWrap) then
            mMean=0.d0
            mStdErr=0.d0
            do j=i1+1,i1+ic_npt(i)
                mMean(:)=mMean(:)+ptk(:,j)
                mStdErr(:)=mStdErr(:)+ptk(:,j)**2
            enddo
            mMean=mMean/dble(ic_npt(i))
            mStdErr=mStdErr/dble(ic_npt(i))
            mStdErr=sqrt(mStdErr-mMean**2)
        endif
        
        do
            npt_mode=0
            lList=.false.
            toBeChkd=.false.
            n_mode=0 !no. of clusters in the mode
        
            !find a starting sub-cluster
            do j=1,k
                if(gList(j)) then
                    cycle
                else
                    n_mode=1
                    npt_mode=nptx(sc_n+j)
                    lList(j)=.true.
                    gList(j)=.true.
                    toBeChkd(j)=.true.
                    m=j
                    exit
                endif
            enddo
        
            !didn't find a starting position?
            if(n_mode==0) exit
            
            do
                flag=.false.
                do j=1,k
                    if(.not.toBeChkd(j)) cycle
                    flag=.true.
                    toBeChkd(j)=.false.
                    mean1(:)=meank(j,:)
                    eval1(:)=evalk(j,:)
                    ef1=kfack(j)*((1.d0+ef0*sqrt(60.d0/dble(nptx(sc_n+j))))**(1d0/nCdim))
                    ef1=kfack(j)*1.5d0
                    invcov1(:,:)=invcovk(j,:,:)
                    evec1(:,:)=eveck(j,:,:)
                    exit
                enddo
                
                if(.not.flag) exit
                        
                do n=1,k
                    if(lList(n) .or. n==j .or. .not.overlapk(n,j)) cycle
                    mean2(:)=meank(n,:)
                    ef2=kfack(n)*((1.d0+ef0*sqrt(60.d0/dble(nptx(sc_n+n))))**(1d0/nCdim))
                    ef2=kfack(n)*1.5d0
                    invcov2(:,:)=invcovk(n,:,:)
                    
                    intFlag=.false.
                    
                    if(mWrap .and. (wrapEll(j) .or. wrapEll(n))) then
                        do i2=0,2**(wrapN(j)-1)-1
                            call returnOrder(wrapN(j)-1,i2,order(1:wrapN(j)-1))
                            
                            i4=0
                            do i5=1,nCdim
                                if(wrapDim(j,i5,1) .or. wrapDim(j,i5,2)) then
                                    i4=i4+1
                                    if(order(i4)==0) then
                                        mean1w(i5)=meanw(j,i5)
                                    else
                                        mean1w(i5)=mean1(i5)
                                    endif
                                else
                                    mean1w(i5)=mean1(i5)
                                endif
                            enddo
                            
                            
                            do j2=0,2**(wrapN(n)-1)-1
                                call returnOrder(wrapN(n)-1,j2,order(1:wrapN(n)-1))
                                
                                j4=0
                                do j3=1,nCdim
                                    if(wrapDim(n,j3,1) .or. wrapDim(n,j3,2)) then
                                        j4=j4+1
                                        if(order(j4)==0) then
                                            mean2w(j3)=meanw(n,j3)
                                        else
                                            mean2w(j3)=mean2(j3)
                                        endif
                                    else
                                        mean2w(j3)=mean2(j3)
                                    endif
                                enddo
                            
                                if(ellIntersect(nCdim,mean1w,mean2w,eval1,evec1,ef1,ef2,invcov1,invcov2)) then
                                    intFlag=.true.
                                    exit
                                endif
                            enddo
                            
                            if(intFlag) exit
                        enddo
                    elseif(ellIntersect(nCdim,mean1,mean2,eval1,evec1,ef1,ef2,invcov1,invcov2)) then
                        intFlag=.true.
                    endif
                    
                    if(intFlag) then
                        lList(n)=.true.
                        gList(n)=.true.
                        toBeChkd(n)=.true.
                        n_mode=n_mode+1
                        npt_mode=npt_mode+nptx(sc_n+n)
                    else
                        overlapk(n,j)=.false.
                        overlapk(j,n)=.false.
                    endif
                enddo
            enddo
            
            !found a candidate?
            if((n_mode<k .or. ic_reme(i)) .and. npt_mode>=2*(ndim+1) .and. ((ic_npt(i)-npt_mode)>=2*(ndim+1) &
            .or. (ic_reme(i) .and. ic_npt(i)-npt_mode==0))) then
                flag=.true.
            else
                flag=.false.
            endif

            !now check for localization in ndim parameters
            !calculate the standard error
            if(flag .and. n_mode<k .and. npt_mode<ic_npt(i) .and. .not.aWrap) then
                n=i1
                lMean=0.d0
                lStdErr=0.d0
                do j=sc_n+1,sc_n+k
                    if(lList(j)) then
                        do i2=n+1,n+nptx(j)
                            lMean(:)=lMean(:)+ptk(:,i2)
                            lStdErr(:)=lStdErr(:)+ptk(:,i2)**2
                        enddo
                    endif
                    n=n+nptx(j)
                enddo
                lMean=lMean/dble(npt_mode)
                lStdErr=lStdErr/dble(npt_mode)
                lStdErr=sqrt(lStdErr-lMean**2)
                
                do j=1,nCdim
                    if(lStdErr(j)<mStdErr(j)) exit
                    if(j==nCdim) flag=.false.
                enddo
            endif
            
            if(flag) then
                nLost=nLost+npt_mode
                ic_n=ic_n+1
                if(ic_n>maxCls) then
                    write(*,*)"ERROR: More modes found than allowed memory."
                    write(*,*)"Increase maxmodes in the call to nestrun and run MultiNest again."
                    stop
                endif
                ic_fnode(ic_n)=i
                ic_reme(ic_n)=.false.
                ic_reme(i)=.true.
                isolateModes2=.true.
                reCluster(i)=.true.
                reCluster(ic_n)=.true.
                do j=1,k
                    if(lList(j)) nodex(sc_n+j)=ic_n
                enddo
            endif
        enddo
        
        !make the leftover points into a new node
        if(ic_reme(i) .and. ic_npt(i)-nLost>=ndim+1 .and. .false.) then
            ic_n=ic_n+1
            ic_fnode(ic_n)=i
            ic_reme(ic_n)=.false.
            reCluster(ic_n)=.true.
            do j=1,k
                if(nodex(sc_n+j)==i) nodex(sc_n+j)=ic_n
            enddo
        endif
        
        if(ic_n==i3) then
            do j=1,nCdim
                ptk(j,i1+1:i1+ic_npt(i))=pt(j,i1+1:i1+ic_npt(i))
            enddo
            auxk(1:naux,i1+1:i1+ic_npt(i))=aux(1:naux,i1+1:i1+ic_npt(i))
        endif
        
        i1=i1+ic_npt(i)
        sc_n=sc_n+k
        deallocate(wrapN,wrapEll,wrapDim,meank,evalk,eveck,invcovk,tmatk,kfack,meanw)
    enddo

    !if modes found then re-arrange everything
    if(isolateModes2) then
        i1=0 !no. of points re-arranged
        do i=1,ic_n
            ic_npt(i)=0
            k=0 !no. of points traversed
            do j=1,sc_n
                if(nodex(j)==i) then
                    !arrange the points
                    do i2=1,nCdim
                        ptx(i2,i1+1:i1+nptx(j))=ptk(i2,k+1:k+nptx(j))
                    enddo
                    auxx(1:naux,i1+1:i1+nptx(j))=auxk(1:naux,k+1:k+nptx(j))
                    i1=i1+nptx(j)
                    ic_npt(i)=ic_npt(i)+nptx(j)
                endif
                k=k+nptx(j)
            enddo
        enddo
        pt=ptx
        aux=auxx
    endif

    deallocate( order, nptx, nodex )
    deallocate( gList, lList, toBeChkd, overlapk )
    deallocate( ptk, auxk, ptx, auxx, mMean, lMean, mStdErr, lStdErr, mean1, mean2, mean1w, &
    mean2w, eval1, evec1, invcov1, invcov2 )
 
 end function isolateModes2
  
!----------------------------------------------------------------------

 subroutine setLimits(ndim,nCdim,llimits,plimits,pnew,phyPnew,climits)
 
    implicit none

    !input variables
    integer ndim !dimensionality
    integer nCdim !clustering dimension
    double precision pnew(ndim) !new point
    double precision phyPnew(nCdim) !new physical point
    double precision climits(ndim,2) !current scaling limits

    !input/output variables
    double precision llimits(ndim,2) !current limits
    double precision plimits(nCdim,2) !current clustering limits

    !work variables
    integer i
    double precision pt(ndim)


    !first scale the point
    do i=1,ndim
        pt(i)=climits(i,1)+(climits(i,2)-climits(i,1))*pnew(i)
        
        !wraparound
        if(pWrap(i)) call wraparound(pt(i),pt(i))
    enddo
                    
    do i=1,ndim
        if(pt(i)<llimits(i,1)) then
            llimits(i,1)=pt(i)
        elseif(pt(i)>llimits(i,2)) then
            llimits(i,2)=pt(i)
        endif
    enddo

    do i=1,nCdim
        if(phyPnew(i)<plimits(i,1)) then
            plimits(i,1)=phyPnew(i)
        elseif(phyPnew(i)>plimits(i,2)) then
            plimits(i,2)=phyPnew(i)
        endif
    enddo
 
 end subroutine setLimits
  
!----------------------------------------------------------------------

 subroutine wraparound(oPt,wPt)

    implicit none
    double precision oPt !actual point
    double precision wPt !wrapped-around point

    wPt=oPt
    do
        if(wPt<0.d0 .or. wPt>1.d0) then
            wPt=wPt-floor(wPt)
        else
            exit
        endif
    enddo

 end subroutine wraparound
  
!----------------------------------------------------------------------

 subroutine wrapEllCheck(ndim,mean,TMat,ef,limits,wrapEll,wrapDim)
                        
    implicit none

    !input variables
    integer ndim !dimensionality
    double precision mean(ndim)
    double precision TMat (ndim,ndim) !transformation matrix
    double precision ef !enlargement factor
    double precision limits(ndim,2) !current scale limits

    !output variable
    logical wrapEll
    logical wrapDim(ndim,2) !which dimensions in which directions to wrap

    !work variable
    integer i,j,k
    double precision cubeEdge(ndim),sCubeEdge
    double precision pnewM(1,ndim),u(1,ndim)


    wrapEll=.false.
    wrapDim=.false.

    !loop over the principle directions
    do i=1,ndim
        !loop over the two edges
        do k=1,2
            u(1,:)=0d0
            if(k==1) then
                u(1,i)=1d0
            else
                u(1,i)=-1d0
            endif
                
            pnewM=MatMul(u,TMat)
            cubeEdge(:)=sqrt(ef)*pnewM(1,:)+mean(:)
        
            !loop over dimensions
            do j=1,ndim
                if(pWrap(j) .and. (.not.wrapDim(j,1) .or. .not.wrapDim(j,2))) then
                    !transform back to unit hypercube
                    sCubeEdge=limits(j,1)+(limits(j,2)-limits(j,1))*cubeEdge(j)
                
                    if(sCubeEdge<0d0) then
                        wrapDim(j,1)=.true.
                        wrapEll=.true.
                    endif
                    if(sCubeEdge>1d0) then
                        wrapDim(j,2)=.true.
                        wrapEll=.true.
                    endif
                endif
            enddo
        enddo
    enddo

    !sanity check
    !no wraparound if both edges are outside the hyper cube boundary
    if(wrapEll) then
        wrapEll=.false.
        do i=1,ndim
            if(wrapDim(i,1) .and. wrapDim(i,2)) then
                wrapDim(i,1)=.false.
                wrapDim(i,2)=.false.
            else
                wrapEll=.true.
            endif
        enddo
    endif  


 end subroutine wrapEllCheck
  
!----------------------------------------------------------------------

 subroutine scaled2Cube(ndim,limits,sP,cP)
 
    implicit none

    !input variables
    integer ndim
    double precision limits(ndim,2) !current limits
    double precision sP(ndim) !scaled point

    !output variable
    double precision cP(ndim) !point in unit hypercube

    !work variables
    integer i


    do i=1,ndim
        cP(i)=limits(i,1)+(limits(i,2)-limits(i,1))*sP(i)
    enddo
 
 end subroutine scaled2Cube
  
!----------------------------------------------------------------------

 subroutine Cube2Scaled(ndim,limits,sP,cP)
 
    implicit none

    !input variables
    integer ndim
    double precision limits(ndim,2) !current limits
    double precision sP(ndim) !scaled point

    !output variable
    double precision cP(ndim) !point in unit hypercube

    !work variables
    integer i


    do i=1,ndim
        sP(i)=(cP(i)-limits(i,1))/(limits(i,2)-limits(i,1))
    enddo
 
 end subroutine Cube2Scaled
  
!----------------------------------------------------------------------

 subroutine returnOrder(nBits,n,order)
 
    implicit none

    !input variables
    integer nBits !no. of bits
    integer n !number under consideration

    !output variable
    integer order(nBits) !list with the order

    !work variables
    integer i,j


    order=0

    do i=1,nBits
        j=2**(i-1)
        if(mod(n/j,2)==0) then
            order(i)=0
        else
            order(i)=1
        endif
    enddo
 
 end subroutine returnOrder
  
!----------------------------------------------------------------------
 
 !check if node1 is ancestor of node2
 logical function isAncestor(node1, node2, fnode)
 
    implicit none

    !input variables
    integer node1			!ancestor node to be checked
    integer node2			!child node to be checked
    integer fnode(node2)		!array with parent nodes

    !work variables
    integer i, n2

    isAncestor = .false.
    if( node1 > node2 ) return

    if( node1 == 1 ) then
        isAncestor = .true.
        return
    endif

    n2 = node2
    do
        if( n2 == node1 ) then
            isAncestor = .true.
            return
        else
            n2 = fnode(n2)
            if( node1 > n2 ) return
        endif
    enddo
 
 end function isAncestor
  
!----------------------------------------------------------------------

end module Nested
