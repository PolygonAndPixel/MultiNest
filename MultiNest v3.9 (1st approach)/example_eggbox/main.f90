program main

    use params
    use nestwrapper
    use omp_lib
    implicit none
    
      
    integer i
    double precision pi
    integer :: counti, countf, count_rate, funit, threads
    integer :: rounds = 100
    integer, allocatable :: calls(:), iterations(:)
    double precision, allocatable :: dt(:)
    double precision meanCalls, meanIterations, meanTime
    allocate(iterations(rounds))
    allocate(calls(rounds))
    allocate(dt(rounds))
    
    do threads = 2, OMP_get_num_procs() 
        do i=1,rounds
            write(*,*) 'Round', i
    
            pi = 4d0*atan(1d0)
                
            !setting priors
            spriorran(1:sdim,1)=0d0
            spriorran(1:sdim,2)=10d0*pi
            calls(i) = 0
            iterations(i) = 0
            !no parameters to wrap around
            nest_pWrap=0
            call system_clock(counti, count_rate)
            call nest_Sample(threads, calls(i), iterations(i))
            call system_clock(countf)
            dt(i) = DBLE(countf-counti)/DBLE(count_rate)
            
        enddo
        write(*,*) '######################################################'
        meanTime = 0
        meanCalls = 0
        meanIterations = 0
        do i=1,rounds
            meanTime = meanTime + dt(i)
            meanCalls = meanCalls + calls(i)
            meanIterations = meanIterations + iterations(i)
        enddo
        meanCalls = meanCalls/rounds
        meanIterations = meanIterations/rounds
        meanTime = meanTime/rounds

        open(unit=funit,file='benchmarks/benchmark_eggbox_1st', access='sequential',&
        status='unknown',position='append')
        write(funit,*) 'Number of threads: ', threads
        write(funit,*) 'total likelihood calls: ', meanCalls
        write(funit,*) 'total iterations: ', meanIterations
        write(funit,*) 'Elapsed time: ', meanTime
        write(funit,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
        close(funit)
    enddo
    stop
end
