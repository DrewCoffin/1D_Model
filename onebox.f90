!#########################################
!           Matthew Copper
! One box model of io plasma torus
! Based on IDL program by Andrew Steffl 
! Started on 2/25/2013
!#########################################

PROGRAM Onebox

  USE DEBUG
  USE FUNCTIONS
  USE TIMESTEP 
  USE ReadEmis
  USE INPUTS
  USE FTMIX
  USE PARALLELVARIABLES
  USE OUTPUTDATA
  USE MPI

  IMPLICIT NONE

  character(len=8)    ::x1

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierr)
    mygrid=mype+1 
  open(unit=21,file="/home/dcoffin/1D_Model/values.dat")
  write (x1, '(I3.3)') mygrid !write the integer 'mygrid' to non-existent file
 
  num_char=trim(x1)  !trim non-existent file and store as num_char
  !num_char can be used to name output files by grid space using "output"//num_char//".dat" 

  if ( npes .ne. LNG_GRID) then
    print *, "The current version only supports ", LNG_GRID, " processors."   
  else 
   call model()
  endif
 
call MPI_FINALIZE(ierr)

CONTAINS 

subroutine model()
  integer             ::nit
  real                ::lontemp, day
  real                ::tm, tm0
  type(density)       ::n, ni, np
  real                ::const
  type(temp)          ::T, Ti, Tp
  double precision    ::Te0, Ti0, Teh0
  type(height)        ::h, hi
  type(r_ind)         ::ind
  type(nT)            ::nrg, nTi, nTp
  integer             ::i, j, k
  real                ::var, var2, n_height
  type(nu)            ::v, vi
  type(r_dep)         ::dep, depi
  type(lat_dist)      ::lat, lati
  type(ft_int)        ::ft
  type(energy)        ::nrgy
  type(ft_mix)        ::plot
  real                ::output_it
  character(len=8)    ::x1
  character(len=4)    ::day_char
  integer             ::file_num
  real                ::longitude, elecHot_multiplier, intensity, n_ave, T_ave, test_multiplier, volume, dr

!  call initNu(v)
 
  write(21,*) "~~~~FROM ONEBOX~~~~"
  write(21,*) "mype: ",mype !!!!!
  longitude = (mype * 360.0 / LNG_GRID)

  write(21,*) "longitude: ",longitude !!!!!

  do i=1, LAT_SIZE
    lat%z(i)= (i-1) * Rj / 5.0   !Initializing lat%z
  end do

  call readInputs()  !call to input.f90 module to read initial variables from 'input.dat'
!print *, source
  call read_rec_tables()

!set trans_ variables (user prompt or formatted file migh be used in the future)
  trans_exp=1.0
  trans_type=.false.

!set dt (2000)
  dt=750.0 
  write(21,*) "dt: ",dt !!!!! 
! source = source *2000.0/dt

!set run time
  write(21,*) "run_days: ", run_days !!!!!
  runt=run_days*8.64e4 !one day = 86400 seconds
  write(21,*) "runt: ",runt !!!!! 
  nit=(runt/dt)+1 ! number of iterations to reach run_days
  !nit=5
 write(21,*) "nit: ",nit !!!!! 

!set radial distance
  rdist= 6.0   !in Rj
  dr=Rj
  volume=PI*((((rdist*Rj+1.5*Rj)*1.0e5))**2 - ((rdist*Rj)*1.0e5)**2)*0.5*ROOTPI*Rj*1.0e5
  write(21,*) "volume: ",volume !!!!! 
  torus_circumference = Rj * rdist * 2.0 * PI
  dx = torus_circumference / LNG_GRID
  write(21,*) "dx: ",dx !!!!! 
  numerical_c_neutral = v_neutral*dt/dx
  write(21,*) "numerical_c_neutral: ",numerical_c_neutral !!!!! 
  numerical_c_ion = v_ion*dt/dx
  write(21,*) "numerical_c_ion: ",numerical_c_ion !!!!! 

!set sys3 longitude of box
  lon3=200

!set zoff
  write(21,*) "longitude: ", longitude  !!!!! 
  !theta_offset = (6.4*cos((lon3-longitude)*dTOr)*dTOr)
  !write(21,*) "theta_offset: ", theta_offset  !!!!! 
  write(21,*) "Rj: ", Rj
!!!!!  zoff= abs((6.4*cos((lon3-longitude)*dTOr)*dTOr) * rdist * Rj) !in km (This calculation depends on mype, before switch to mathc idl model

zoff = abs((6.4*cos((110-200)*dTOr)*dTOr)*rdist*Rj) !!!!! Matches IDL model
  write(21,*) "zoff: ",zoff !!!!! 
  n_height = Rj/2.0

  tm0=0.01 
  write(21,*) "tm0: ",tm0 !!!!! 

!set density values
  const=1800.0

  test_multiplier=1.0
  if( test_pattern ) then
!!    if( rdist .gt. 7.0 .and. rdist .lt. 9.0) test_multiplier=4.0
    test_multiplier=1.0+0.2*cos(2.0*longitude*dTOr)
    write(21,*) "test_multiplier: ",test_multiplier !!!!! 
  endif

  n%sp = 0.060 * const * test_multiplier!* (rdist/6.0)**(-8.0)
  write(21,*) "n%sp: ",n%sp !!!!! 
  n%s2p= 0.212 * const * test_multiplier!* (rdist/6.0)**(-8.0)
  write(21,*) "n%s2p: ",n%s2p !!!!!
  n%s3p= 0.034 * const * test_multiplier!* (rdist/6.0)**(-3.0)
  write(21,*) "n%s3p: ",n%s3p !!!!!
  n%op = 0.242 * const * test_multiplier!* (rdist/6.0)**(-8.0)
  write(21,*) "n%op: ",n%op !!!!!
  n%o2p= 0.123 * n%op * test_multiplier!* (rdist/6.0)**(-3.0)
  write(21,*) "n%o2p: ",n%o2p !!!!!

!  n%sp = 150.0 * test_multiplier
!  n%s2p= 600.0 * test_multiplier
!  n%s3p=  100.0 * test_multiplier
!  n%op = 400.0 * test_multiplier
!  n%o2p=  40.0 * test_multiplier

  n%s=25.0 * test_multiplier
  write(21,*) "n%s: ",n%s !!!!!
  n%o=50.0 * test_multiplier
  write(21,*) "n%o: ",n%o !!!!!


  Te0 = 5.0
  Ti0 = 70.0
  Teh0= tehot
  write(21,*) "Teh0: ",Teh0 !!!!!
  if(Teh0 .gt. 400.0) Teh0=400.0
  n%fh=fehot_const
  write(21,*) "n%fh: ",n%fh !!!!!
  trans = 1.646851e-7 !4.62963e-7
  write(21,*) "trans: ",trans !!!!!
  !net_source=source/volume
  net_source = source !!!!!
  !net_source = 6.3e6 !!!!! set to match idl calculation independant of volume
  write(21,*) "net_source: ",net_source !!!!!
  print *, net_source
  n%elec = (n%sp + n%op + (2.0 * (n%s2p + n%o2p)) + (3.0 * n%s3p)) * (1.0 - protons)
  write(21,*) "n%elec: ",n%elec !!!!!
  n%elecHot = n%fh * n%elec / (1.0-n%fh)
  write(21,*) "n%elecHot: ",n%elecHot !!!!!

  n%fc = 1.0 - n%fh
  write(21,*) "n%fc: ",n%fc !!!!!

!set temp values
  T%sp      = Ti0
  T%s2p     = Ti0
  T%s3p     = Ti0
  T%op      = Ti0
  T%o2p     = Ti0
  T%elec    = Te0
  T%elecHot = Teh0
  
  write(21,*) "~~~SETTING INITIAL TEMPS~~~" !!!!!
  write(21,*) "T%sp: ",T%sp !!!!!  
  write(21,*) "T%s2p: ",T%s2p !!!!!
  write(21,*) "T%s3p: ",T%s3p !!!!!
  write(21,*) "T%op: ",T%op !!!!!
  write(21,*) "T%o2p: ",T%o2p !!!!!
  write(21,*) "T%elec: ",T%elec !!!!!
  write(21,*) "T%elecHot: ",T%elecHot !!!!!
  

!get scale heights 
  call get_scale_heights(h, T, n)
  write(21,*) "h: ",h !!!!!
  write(21,*) "T: ",T !!!!!
  write(21,*) "n: ",n !!!!! 
  if (protons > 0.0) then
    n%protons = protons
  endif
  write(21,*) "n%protons: ",n%protons !!!!!
  ind%o_to_s= o_to_s
 ! write(21,*) "ind%o: ",ind%o !!!!!
  ind%o2s_spike=2.0
 ! write(21,*) "ind%o2s: ",ind%o2s !!!!!

  tau0=transport !1.0/(trans*8.64e4)
  write(21,*) "tau0: ",tau0 !!!!!
  print *, tau0
  net_source0=net_source 
  !fh0 = fehot_const

  h%s=n_height
  write(21,*) "h%s: ",h%s !!!!!
  h%o=n_height
  write(21,*) "h%o: ",h%o !!!!!

  call InitIndependentRates(ind)

  T%pu_s = Tpu(32.0, rdist*1.0)
  write(21,*) "T%pu_s: ",T%pu_s !!!!!
  T%pu_o = Tpu(16.0, rdist*1.0)
  write(21,*) "T%pu_o: ",T%pu_o !!!!!

  T%elecHot=Teh0
  write(21,*) "T%elecHot: ",T%elecHot !!!!!

  call independent_rates(ind, T, h)

  !n%fh = protons
  !n%fc= 1.0 - n%fh 
  write(21,*) "n%fc: ",n%fc !!!!!

  write(21,*) "n%sp: ",n%sp !!!!!   
  write(21,*) "n%s2p: ",n%s2p !!!!!   
  write(21,*) "n%s3p: ",n%s3p !!!!!   
  write(21,*) "n%op: ",n%op !!!!!   
  write(21,*) "n%o2p: ",n%o2p !!!!!   
  write(21,*) "protons: ", protons !!!!!
  n%elec = ((n%sp + n%op) + 2.0*(n%s2p + n%o2p) + 3.0 * n%s3p)*(1.0-n%protons)
  write(21,*) "n%elec: ",n%elec !!!!!
  n%elecHot = n%elec * n%fh/n%fc
  write(21,*) "n%elecHot: ",n%elecHot !!!!!
  nrg%elec = n%elec * T%elec
  write(21,*) "nrg%elec: ",nrg%elec !!!!!
  nrg%elecHot = n%elecHot * T%elecHot
  write(21,*) "nrg%elecHot: ",nrg%elecHot !!!!!
  nrg%sp = n%sp * T%sp
  write(21,*) "nrg%sp: ",nrg%sp !!!!!
  nrg%s2p = n%s2p * T%s2p
  write(21,*) "nrg%s2p: ",nrg%s2p !!!!!
  nrg%s3p = n%s3p * T%s3p
  write(21,*) "nrg%s3p: ",nrg%s3p !!!!!
  nrg%op = n%op * T%op
  write(21,*) "nrg%op: ",nrg%op !!!!!
  nrg%o2p = n%o2p * T%o2p
  write(21,*) "nrg%o2p: ",nrg%o2p !!!!!

  ni=n
  np=n

  Ti=T
  Tp=T

  hi=h
 
  nTi=nrg
  nTp=nrg

  vi=v

  lati=lat

  call get_scale_heights(h, T, n)

  output_it = 0 !This variable determine when data is output. 
  Io_loc=0      !Io's location in the torus
  sys4_loc=0    !The location of the sys4 hot electron population
  file_num=0    !Output files are numbered so they can be assembled as a animated visualization (refer to scripts)

  mass_loading=1.0 !1e-33 !set to arbitarily small value
  ave_loading=1.0 !1e-33

!-----------------------Interation loop-------------------------------------------------------------------------------------------------
  do i=1, nit
    tm = tm0 + (i-1) * dt / 86400.0
     write(21,*) "iteration: ",i !!!!!   
!----------------------time dependent neutral source rate-------------------------------------------------------------------------------
    var =exp(-((tm-neutral_t0)/neutral_width)**2)
    var2 =exp(-((tm-hote_t0)/neutral_width)**2)

  !  net_source = net_source0*(1.0 + neutral_amp*var) !Ubiquitous source
    if( moving_Io ) then
      if( mype .eq. int(Io_loc*LNG_GRID/torus_circumference) )then
        net_source = LNG_GRID*net_source0*(1.0+neutral_amp*var)
      else
        if( i .eq. 1 ) then
         write(21,*) "~~~net_source is being called at if i=1 ~~~"
          net_source = net_source0*(1.0+neutral_amp*var)
        else
          write(21,*) "~~~net_source is being called at net_source=0 ~~~"
          net_source=0
        endif
      endif
    endif

    if( .not. moving_Io ) then
       write(21,*) "~~~net_source is being called at .not. moving_io~~~"
       write(21,*) "net_source0: ",net_source0 !!!!! 
       write(21,*) "neutral_amp: ",neutral_amp !!!!!
       write(21,*) "var: ",var !!!!! 
      net_source = (net_source0*(1.0 + neutral_amp*var))!/LNG_GRID !ubiquitous
    endif
!----------------------time dependent neutral source rate-------------------------------------------------------------------------------

    ind%o_to_s = o_to_s
!    ind%o_to_s = (otos + o2s_spike * neutral_amp * var) & !o2s_spike
!               / (1.0 + neutral_amp * var)
!    n%fh  = fehot_const * (1.0 + hote_amp * var)

    elecHot_multiplier=1.0

!----------------------hot electrons-------------------------------------------------------------------------------
    if( sys3hot ) then
      elecHot_multiplier=elecHot_multiplier*(1.0) !+sys3_amp*(cos((290.0-longitude)*dTOr)))
    endif
    

    if( sys4hot ) then
      elecHot_multiplier=elecHot_multiplier&
             *(1.0+(sys4_amp+(hote_amp*var2))*cos(((mype/(LNG_GRID-1.0))-(sys4_loc/torus_circumference))*2.0*PI))
    endif

 !   elecHot_multiplier=elecHot_multiplier*(1.0+0.4*(mass_loading/ave_loading))
   
    write(21,*) "elecHot_multiplier: ",elecHot_multiplier !!!!!  

    n%fh  = fehot_const * (1.0 + hote_amp * var)*elecHot_multiplier

    write(21,*) "next calculation of n%fh: ",n%fh !!!!!  
    ni%fh = n%fh
    write(21,*) "ni%fh: ",ni%fh !!!!!  
    np%fh = n%fh
    write(21,*) "np%fh: ",np%fh !!!!!  
    n%fc  = 1.0 - n%fh
    write(21,*) "n%fc: ",n%fc !!!!!  
    ni%fc = n%fc
    write(21,*) "ni%fc: ",ni%fc !!!!!  
    np%fc = n%fc
    write(21,*) "np%fc: ",np%fc !!!!!  

   !!!!! n%elecHot = n%elec * n%fh/n%fc!!!!!
    n%elecHot = (n%fh * n%elec)/(n%fc) !!!!! Added this instead of above calculation
    nrg%elecHot = n%elecHot * T%elecHot
    write(21,*) "nrg%elecHot: ",nrg%elecHot !!!!!  

    do j=1, LAT_SIZE
      lat%elec(j) = n%elec!*exp(-(lat%z(j)/h%elec)**2)
      lati%elecHot(j) = n%elecHot!*exp(-(lat%z(j)/h%elec)**2)
    end do

    if ( DEBUG_GEN ) then !this variable set in debug.f90
      call DebugOutput(i, n, h, T, v, nrg)
    endif

    call cm3_latavg_model(n, T, nrg, h, v, ni, Ti, nTi, hi &
                         ,vi, np, Tp, nTp, ind, dep, depi, lat, lati, ft, zoff) 

    call MPI_ALLREDUCE(mass_loading, ave_loading, 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
    ave_loading=ave_loading/(npes*1.0)
!    call MPI_ALLREDUCE(mass_loading, min_loading, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, ierr)
!    call MPI_ALLREDUCE(mass_loading, max_loading, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
!    v_ion=2.0
!    print *, longitude, mass_loading/ave_loading
!    v_ion=0.0+((mass_loading/ave_loading)**1.50)*2.05 !ion lag in km/s where v_ion is scaled on mass_loading
!    n%fh=n%fh+n%fh*2.0*(((mass_loading/ave_loading)**1.0)-1.0) !Should replace sys4 population.
!    v_ion=(((mass_loading*volume*mp)*57.0*(rdist)**5)/(0.3*sqrt(1.0-(1.0/(rdist)))*4.0*PI*1.5*((Rj*1.0e3)**2)*(4.2e-4)**2))
!    v_ion=v_ion+0.5
    if( mype .eq. 0) then
      open(unit=199, file='FFT.dat', status='unknown', position='append')
      write(199, *) tm, n%s3p
      close(199)
!      print *, v_ion, mass_loading*volume, mass_loading2*volume, ave_loading*volume*mp
    endif
    numerical_c_ion = v_ion*dt/dx

    call update_temp(n, nrg, T)

    call get_scale_heights(h, T, n)

    call energyBudget(n, h, T, dep, ind, ft, lat, v, nrgy)

    if (nint(output_it)+1 .eq. i .and. (OUTPUT_MIXR .or. OUTPUT_DENS .or. OUTPUT_TEMP .or. OUTPUT_INTS)) then !Output at set intervals when OUTPUT_MIX is true (from debug.f90)
        day = (i-1.0)*dt/86400
        write (x1, '(I4.4)') file_num
        day_char=trim(x1)  !trim non-existent file and store as day_char
        !if( mype .eq. 0 ) then
        !endif
        call dens_ave(n_ave, n) 
        call temp_ave(T_ave, T) 
        do j=0, LNG_GRID - 1
          if( mype .eq. j) then
            if(OUTPUT_DENS) call IonElecOutput(n%sp, n%s2p, n%s3p, n%op, n%o2p, n%elec, longitude, day_char, 'DENS')
            if(OUTPUT_MIXR) then  
              plot = ftint_mix(n, h) !calculate values to be plotted
              call IonOutput(plot%sp, plot%s2p, plot%s3p, plot%op, plot%o2p, longitude, day_char, 'MIXR')
            endif
            if(OUTPUT_TEMP) call IonElecOutput(T%sp, T%s2p, T%s3p, T%op, T%o2p, T%elec, longitude, day_char, 'TEMP')
            if(OUTPUT_INTS) then !Intensity
              call IonOutput(n%sp*T%sp, n%s2p*T%s2p, n%s3p*T%s3p, n%op*T%op, n%o2p*T%o2p, longitude, day_char, 'INTS')
              intensity= n%sp*T_ave/(n_ave*T%sp)
              open(unit=120, file='intensity'//day_char//'.dat', status='unknown', position='append')
                write(120,*) longitude, intensity 
              close(120)
            end if
            call OtherOutput(100.0*mass_loading/ave_loading, longitude, day_char, 'LOAD')
          endif
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        end do
        if( mype .eq. 0 ) then
          if(OUTPUT_DENS) call IonElecOutput(n%sp, n%s2p, n%s3p, n%op, n%o2p, n%elec, 360.0, day_char, 'DENS')
          if(OUTPUT_MIXR) then  
            plot = ftint_mix(n, h) !calculate values to be plotted
            !print *, "ftmix = ", plot 
            call IonOutput(plot%sp, plot%s2p, plot%s3p, plot%op, plot%o2p, 360.0, day_char, 'MIXR')
          endif
          if(OUTPUT_TEMP) call IonElecOutput(T%sp, T%s2p, T%s3p, T%op, T%o2p, T%elec, 360.0, day_char, 'TEMP')
          if(OUTPUT_INTS) then !Intensity
            call IonOutput(n%sp*T%sp, n%s2p*T%s2p, n%s3p*T%s3p, n%op*T%op, n%o2p*T%o2p, 360.0, day_char, 'INTS')
            intensity= n%sp*T_ave/(n_ave*T%sp)
            write(21,*) "intensity: ", intensity !!!!!  
            open(unit=120, file='intensity'//day_char//'.dat', status='unknown', position='append')
              write(120,*) 360.0, intensity 
            close(120)
          end if
          call OtherOutput(100.0*mass_loading/ave_loading, 360.00, day_char, 'LOAD')
        endif

        output_it=output_it + (86400.0/(dt*per_day)) !Determines when data is output. Set for once each run day (86400/dt).
        write(21,*) "output_it: ",output_it !!!!!  
        file_num = file_num + 1
    endif        
 
    !call Grid_transport(n, nrg)
!    isNaN=NaNcatch(n%sp, 100, mype) !fix
    Io_loc = mod(Io_loc+(dt*v_Io), torus_circumference)
    write(21,*) "Io_loc: ",Io_loc !!!!!  
    sys4_loc = mod(sys4_loc+(dt*v_sys4), torus_circumference)
    write(21,*) "sys4_loc: ",sys4_loc !!!!!  

  end do
!----------------------------------------------------------------------------------------------------------------------------

call FinalOutput(nrgy)

end subroutine model

subroutine dens_ave(n_ave, n)!, i)
  real                ::n_tot, n_ave, space_ave
  type(density)       ::n
!  integer             ::i

  n_tot=n%sp !+n%s2p+n%s3p+n%op+n%o2p
 ! write(21,*) "n_tot: ",n_tot !!!!!  
  call MPI_REDUCE(n_tot, n_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  n_ave=n_ave/LNG_GRID
 ! write(21,*) "n_ave: ",n_ave !!!!! 
  call MPI_BCAST(n_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

!  n_tot=n%sp+n%s2p+n%s3p+n%op+n%o2p
!  call MPI_REDUCE(n_tot, space_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!  space_ave=space_ave/LNG_GRID
!  n_ave=(n_ave*(i-1.0)+space_ave)/(i*1.0)
!  call MPI_BCAST(n_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

end subroutine dens_ave

subroutine temp_ave(T_ave, T)!, i)
  real                ::T_tot, T_ave, space_ave
  type(temp)          ::T
!  integer             ::i

  T_tot=T%sp !+T%s2p+T%s3p+T%op+T%o2p
  call MPI_REDUCE(T_tot, T_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  T_ave=T_ave/LNG_GRID
  call MPI_BCAST(T_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

!  T_tot=T%sp+T%s2p+T%s3p+T%op+T%o2p
!  call MPI_REDUCE(T_tot, space_ave, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
!  space_ave=space_ave/LNG_GRID
!  T_ave=(T_ave*(i-1.0)+space_ave)/(i*1.0)
!  call MPI_BCAST(T_ave, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

end subroutine temp_ave

subroutine FinalOutput(nrgy)
  type(energy)        ::nrgy, avg
  integer             ::j

  call MPI_REDUCE(nrgy%s_ion, avg%s_ion, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%s_ion=avg%s_ion/LNG_GRID

  call MPI_REDUCE(nrgy%s_cx, avg%s_cx, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%s_cx=avg%s_cx/LNG_GRID

  call MPI_REDUCE(nrgy%o_ion, avg%o_ion, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%o_ion=avg%o_ion/LNG_GRID

  call MPI_REDUCE(nrgy%o_cx, avg%o_cx, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%o_cx=avg%o_cx/LNG_GRID

  call MPI_REDUCE(nrgy%elecHot_eq, avg%elecHot_eq, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%elecHot_eq=avg%elecHot_eq/LNG_GRID

  call MPI_REDUCE(nrgy%tot_eq, avg%tot_eq, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%tot_eq=avg%tot_eq/LNG_GRID

  call MPI_REDUCE(nrgy%P_in, avg%P_in, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%P_in=avg%P_in/LNG_GRID

  call MPI_REDUCE(nrgy%Puv, avg%Puv, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Puv=avg%Puv/LNG_GRID

  call MPI_REDUCE(nrgy%Pfast, avg%Pfast, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Pfast=avg%Pfast/LNG_GRID

  call MPI_REDUCE(nrgy%Ptrans, avg%Ptrans, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Ptrans=avg%Ptrans/LNG_GRID

  call MPI_REDUCE(nrgy%Ptrans_elecHot, avg%Ptrans_elecHot, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%Ptrans_elecHot=avg%Ptrans_elecHot/LNG_GRID

  call MPI_REDUCE(nrgy%P_out, avg%P_out, LNG_GRID, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  avg%P_out=avg%P_out/LNG_GRID

  if( mype .eq. 0 ) then
    print *, "AVERAGE VALUES"
    call FinalTable(avg)
    print *, ""
  endif

!  do j=1, LNG_GRID
!    if( mygrid .eq. j ) then
!      print *, "mygrid = ", mygrid
!      call FinalTable(nrgy)
!      print *, ""
!    endif
!    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!  enddo

end subroutine FinalOutput

subroutine FinalTable(nrgy)
  type(energy)        ::nrgy

  print *, '$$--------------------------------'
  print *, '$$ INPUT PARAMETERS'
  print *, '$$--------------------------------'
  print *, '$$ Source Rate..........', source 
  print *, '$$ Hot Elec Fraction....', fehot_const 
  print *, '$$ Transport............', transport 
  print *, '$$ Hot Elec Temp........', tehot  
  print *, '$$ Subcorotation........', v_ion
  print *, '$$ Sys 3 Hot Elec Frac..', sys3_amp 
  print *, '$$ Sys 4 Hot Elec Amp...', sys4_amp
  print *, '$$ Sys 4 Speed..........', v_sys4

  print *, ''
  print *, '$$ GAUSSIAN SOURCE CHANGE VARIABLES'  
  print *, '$$ Neutral Variation....', neutral_amp 
  print *, '$$ Neutral Temp.........', neutral_t0 
  print *, '$$ Variation Duration...', neutral_width
  print *, '$$ Hot Elec Variation...', hote_amp  
  print *, '$$ Hot Elec Temp........', hote_t0
  print *, '$$ Variation Duration...', hote_width
 
  print *, ''
  print *, '$$ Run Length(days)....',  run_days
  print *, '$$ Outputs per day.....',  per_day

  print *, ''
  print *, '$$--------------------------------'
  print *, '$$ IN-CODE ENERGY BUDGET'
  print *, '$$--------------------------------'
  print *, '$$ ionized S............', nrgy%s_ion
  print *, '$$ ionized O............', nrgy%o_ion
  print *, '$$ charge exchange S....', nrgy%s_cx
  print *, '$$ charge exchange O....', nrgy%o_cx
  print *, '$$ equil with ehot......', nrgy%elecHot_eq + nrgy%tot_eq
  print *, '$$ total in.............', nrgy%P_in + nrgy%tot_eq
  print *, '$$ puv..................', nrgy%Puv
  print *, '$$ fast/ena.............', nrgy%pfast - nrgy%tot_eq
  print *, '$$ transport............', nrgy%ptrans + nrgy%ptrans_elecHot
  print *, '$$ total out............', nrgy%P_out - nrgy%tot_eq
  print *, '$$ in/out...............', (nrgy%P_in + nrgy%tot_eq )/(nrgy%P_out - nrgy%tot_eq )
!  print *, ""
!  print *, '++++++++++++++++++++++++++++++++++++'
!  print *, 'Final Variable Values'
!  print *, '++++++++++++++++++++++++++++++++++++'
!  print *, 'O/S.........................', o_to_s
!  print *, 'Fraction of Hot Electrons...', fehot_const
!  print *, 'Transport...................', transport
!  print *, 'Hot Electron Temp...........', tehot
!  print *, 'Lag Constant................', lag_const
!  print *, 'Neutral Amplitude...........', neutral_amp
!  print *, 'Inital Neutral Temperature..', neutral_t0
!  print *, 'Neutral Width...............', neutral_width
!  print *, 'Hot Electron Amplitude......', hote_amp
!  print *, 'Hot Electron Initial Temp...', hote_t0
!  print *, 'Hot Electron Width..........', hote_width
close(21)
end subroutine FinalTable

subroutine DebugOutput(i, n, h, T, v, nrg)
  integer             ::i
  type(density)       ::n
  type(height)        ::h
  type(temp)          ::T
  type(nu)            ::v
  type(nT)            ::nrg

  print *,  "||||||||||||||||||||||||||||||||||||||||||||||"
  print *,  "grid = ", mygrid
  print *,  "i = ", i-1
  print *,  "||||||||||||||||||||||||||||||||||||||||||||||"
  print *, "~~~~~~~~~~~~~DENSITY~~~~~~~~~~~~~"
  call output(n)
  print *, "~~~~~~~~~~~~~HEIGHT~~~~~~~~~~~~~~"
  call output(h)
  print *, "~~~~~~~~~~~TEMPERATURE~~~~~~~~~~~"
  call output(T)
  print *, "~~~~~~~~~~~~~~~NU~~~~~~~~~~~~~~~~"
  call output(v)
  print *, "~~~~~~~~~~~~~ENERGY~~~~~~~~~~~~~~"
  call output(nrg)
 

end subroutine DebugOutput

subroutine Grid_transport(n, nrg)
  type(density)       ::n, dens_source
  type(nT)            ::nrg, nrg_source

  call az_transport(n, nrg)

end subroutine Grid_transport

subroutine Communicate(dens_source, nrg_source)
  type(density)       ::dens_source
  type(nT)            ::nrg_source

  call MPI_SEND(dens_loss%s, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(dens_source%s, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(dens_loss%sp, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(dens_source%sp, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(dens_loss%s2p, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(dens_source%s2p, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(dens_loss%s3p, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(dens_source%s3p, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(dens_loss%o, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(dens_source%o, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(dens_loss%op, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(dens_source%op, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(dens_loss%o2p, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(dens_source%o2p, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg_loss%sp, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nrg_source%sp, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg_loss%s2p, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nrg_source%s2p, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg_loss%s3p, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nrg_source%s3p, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg_loss%op, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nrg_source%op, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg_loss%o2p, 1, MPI_DOUBLE_PRECISION, mod(mype+1, npes), 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nrg_source%o2p, 1, MPI_DOUBLE_PRECISION, mod(mype-1, npes), 22, MPI_COMM_WORLD, stat, ierr)

end subroutine Communicate

subroutine az_transport(n, nrg)
  type(density)     ::n
  type(nT)          ::nrg
  integer           ::left, right
  real              ::c_left, c_right
  real              ::cleft, cright, cn, c

  left = mype - 1 
  right= mype + 1

  if( left  < 0 )           left = left+LNG_GRID
  if( right .ge. LNG_GRID ) right=right-LNG_GRID

  cn=numerical_c_neutral
  c=numerical_c_ion
 
  call MPI_SEND(numerical_c_ion, 1, MPI_REAL, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(c_left, 1, MPI_REAL, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(numerical_c_ion, 1, MPI_REAL, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(c_right, 1, MPI_REAL, right, 22, MPI_COMM_WORLD, stat, ierr)

  cleft=c_left
  cright=c_right

  if(UseLaxWendroff) then
    call GetNeighbors(n, nrg)
    n%s=LaxWendroff(nleft%s, n%s, nright%s, cn, cn, cn)
    n%sp=LaxWendroff(nleft%sp, n%sp, nright%sp, cleft, c, cright)
    n%s2p=LaxWendroff(nleft%s2p, n%s2p, nright%s2p, cleft, c, cright)
    n%s3p=LaxWendroff(nleft%s3p, n%s3p, nright%s3p, cleft, c, cright)
    n%o=LaxWendroff(nleft%o, n%o, nright%o, cn, cn, cn)
    n%op=LaxWendroff(nleft%op, n%op, nright%op, cleft, c, cright)
    n%o2p=LaxWendroff(nleft%o2p, n%o2p, nright%o2p, cleft, c, cright)
!double precision function LaxWendroff(left, center, right, uleft, u, uright)

    nrg%sp=LaxWendroff(nTleft%sp, nrg%sp, nTright%sp, cleft, c, cright)
    nrg%s2p=LaxWendroff(nTleft%s2p, nrg%s2p, nTright%s2p, cleft, c, cright)
    nrg%s3p=LaxWendroff(nTleft%s3p, nrg%s3p, nTright%s3p, cleft, c, cright)
    nrg%op=LaxWendroff(nTleft%op, nrg%op, nTright%op, cleft, c, cright)
    nrg%o2p=LaxWendroff(nTleft%o2p, nrg%o2p, nTright%o2p, cleft, c, cright)
  endif  
 
  if( Upwind ) then
    call GetNeighbors(n, nrg)

    n%s   =UpwindTransport(nleft%s  , nright%s  , n%s  ,numerical_c_neutral, numerical_c_neutral)
    n%sp  =UpwindTransport(nleft%sp , nright%sp , n%sp ,numerical_c_ion, c_left)
    n%s2p =UpwindTransport(nleft%s2p, nright%s2p, n%s2p,numerical_c_ion, c_left)
    n%s3p =UpwindTransport(nleft%s3p, nright%s3p, n%s3p,numerical_c_ion, c_left)
    n%o   =UpwindTransport(nleft%o  , nright%o, n%o  ,numerical_c_neutral ,numerical_c_neutral)
    n%op  =UpwindTransport(nleft%op , nright%op , n%op ,numerical_c_ion, c_left)
    n%o2p =UpwindTransport(nleft%o2p, nright%o2p, n%o2p,numerical_c_ion, c_left)
  
    nrg%sp  =UpwindTransport(nTleft%sp , nTright%sp , nrg%sp ,numerical_c_ion, c_left)
    nrg%s2p =UpwindTransport(nTleft%s2p, nTright%s2p, nrg%s2p,numerical_c_ion, c_left)
    nrg%s3p =UpwindTransport(nTleft%s3p, nTright%s3p, nrg%s3p,numerical_c_ion, c_left)
    nrg%op  =UpwindTransport(nTleft%op , nTright%op , nrg%op ,numerical_c_ion, c_left)
    nrg%o2p =UpwindTransport(nTleft%o2p, nTright%o2p, nrg%o2p,numerical_c_ion, c_left)
  endif

  if( Euler ) then
    n%s   = EulerTransport(n%s  , v_neutral)
    n%sp  = EulerTransport(n%sp , v_ion)
    n%s2p = EulerTransport(n%s2p, v_ion)
    n%s3p = EulerTransport(n%s3p, v_ion)
    n%o   = EulerTransport(n%o  , v_neutral)
    n%op  = EulerTransport(n%op , v_ion)
    n%o2p = EulerTransport(n%o2p, v_ion)

    nrg%sp  = EulerTransport(nrg%sp , v_ion)
    nrg%s2p = EulerTransport(nrg%s2p, v_ion)
    nrg%s3p = EulerTransport(nrg%s3p, v_ion)
    nrg%op  = EulerTransport(nrg%op , v_ion)
    nrg%o2p = EulerTransport(nrg%o2p, v_ion)
  endif

end subroutine az_transport

double precision function UpwindTransport(left, right, center, c, c_left)
  double precision    ::left, right, center
  real                ::c, c_left

!  UpwindTransport = (numerical_s + c)*left + (1 - 2*numerical_s - c)*center + numerical_s*right
  UpwindTransport =  c_left*left + (1.0 - c)*center

end function UpwindTransport

double precision function EulerTransport(old, v) !improved euler method applied to azimuthal transport
  double precision    ::old, intermediate, loss
  real                ::v

  loss = getLoss(v, old)

  intermediate = old - loss

  EulerTransport = old - .5 * (loss + getLoss(v, intermediate))

end function EulerTransport

double precision function LaxWendroff(left, center, right, cleft, c, cright)
  double precision    :: left, center, right
  real                :: cleft, c, cright
  
  LaxWendroff=center+0.5*(cleft*left - cright*right)+0.5*((cleft**2.0)*left-2.0*c*c*center+(cright**2.0)*right)

end function LaxWendroff

double precision function getLoss(v, val)
  real                ::v
  double precision    ::val, source
  integer             ::left, right

  getLoss = val * dt * v * LNG_GRID / torus_circumference
  
  call MPI_SEND(getLoss, 1, MPI_DOUBLE_PRECISION,  right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(source, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  getLoss = getLoss - source

end function getLoss
subroutine GetNeighbors(n, nrg)
  type(density)       ::n
  type(nT)            ::nrg
  integer             ::left, right

  left = mype - 1 
  right= mype + 1

  if( left  < 0 )           left = left+LNG_GRID
  if( right .ge. LNG_GRID ) right=right-LNG_GRID
! ALGORITHM
! Send to right
! Receive from left
! Send to left
! Receive from right
!!!!!!!!!!!!!!!!!!!!DENSITY!!!!!!!!!!!!!!!!!!!!!
  call MPI_SEND(n%s, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nleft%s, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(n%s, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nright%s, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(n%sp, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nleft%sp, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(n%sp, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nright%sp, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(n%s2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nleft%s2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(n%s2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nright%s2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(n%s3p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nleft%s3p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(n%s3p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nright%s3p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(n%o, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nleft%o, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(n%o, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nright%o, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(n%op, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nleft%op, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(n%op, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nright%op, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(n%o2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nleft%o2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(n%o2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nright%o2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)

!!!!!!!!!!!!!!!!!!!!ENERGY!!!!!!!!!!!!!!!!!!!!!!
  call MPI_SEND(nrg%sp, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTleft%sp, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg%sp, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTright%sp, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(nrg%s2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTleft%s2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg%s2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTright%s2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(nrg%s3p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTleft%s3p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg%s3p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTright%s3p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(nrg%op, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTleft%op, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg%op, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTright%op, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)
!!
  call MPI_SEND(nrg%o2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTleft%o2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, stat, ierr)

  call MPI_SEND(nrg%o2p, 1, MPI_DOUBLE_PRECISION, left, 22, MPI_COMM_WORLD, ierr)
  call MPI_RECV(nTright%o2p, 1, MPI_DOUBLE_PRECISION, right, 22, MPI_COMM_WORLD, stat, ierr)

end subroutine GetNeighbors



END PROGRAM Onebox

