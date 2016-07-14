     module brom_transport
    
    use fabm
    use fabm_config
    use fabm_types, only: attribute_length
    use io_netcdf
    use io_ascii
    
    implicit none
    private
    public init_brom_transport, do_brom_transport, clear_brom_transport
    
!FABM model with all data and procedures related to biogeochemistry
    type (type_model)                                               :: model
    integer                                                         :: i, i_min, i_water, i_max ! x-axis related
    integer                                                         :: k, k_min, k_wat_bbl, k_bbl_sed,  k_max, k_points_below_water  ! z-axis related
    integer                                                         :: ip, par_max              ! parameter`s index, max number of parameters
    integer                                                         :: i_day,  year             ! time related

!grids for depth(z), gradients of depth(dz), bioturbations coefficients, hice -  "time-averaged average ice thickness in cell"
    real(8), pointer, dimension(:)                                  :: z, dz, kz_bio, pco2atm, windspeed, hice
    real(8), pointer, dimension(:, :, :)                            :: t, s, kz   !grids for temperature, salinity, coefficient of turbulence(kz)
    real(8), pointer, dimension(:, :, :)                            :: vv, dvv  
    real(8), pointer, dimension(:)                                  :: z_w, dz_w      ! in water column only
    real(8), pointer, dimension(:, :, :)                            :: t_w, s_w, kz_w ! in water column only
    real(8), allocatable, dimension(:, :, :)                        :: cc, dcc, fick, wbio   !cc - array of concentrations
    real(8), allocatable, dimension(:, :)                           :: surf_flux, bound_up, bound_low
    logical, allocatable, dimension(:, :)                           :: use_bound_up, use_bound_low
    character(len = attribute_length), allocatable, dimension(:)    :: par_name
    real(8), allocatable, dimension(:, :)                           :: Izt, density, pressure  !Izt - irradiance at z
    real(8), allocatable, dimension(:, :)                           :: kz_tot  ! to delete 
    real(8)                                                         :: Io, Iot, wind_speed, pco2_atm, &
                                                                       water_layer_thickness, bbl_thickness , dz_bbl_min, &
                                                                       dz_sed_min, dz_sed_max                  

!indexes of state variables that are needed in brom-transport and subroutines (i.e. boundary conditions description)
    integer         :: id_O2, id_Mn2, id_Mn3, id_Mn4, id_H2S, id_Fe2, id_Fe3, id_FeS,   &
                       id_MnS, id_DON, id_PON, id_NH4, id_NO2, id_NO3, id_S0, id_S2O3,  &
                       id_SO4, id_DIC, id_Alk, id_PO4, id_Si, id_Sipart, id_Phy, id_Het,&
                       id_Baae, id_Bhae, id_Baan, id_Bhan, id_Hplus, id_CaCO3, id_FeS2, id_MnCO3
    
    contains
!=================================================================================    
    subroutine init_brom_transport()
    
    use smooth

    implicit none
    
!reading brom.yaml
    call init_common()
!getting variable values from brom.yaml
    water_layer_thickness = get_brom_par("water_layer_thickness")    
    bbl_thickness = get_brom_par("bbl_thickness")   
    dz_bbl_min = get_brom_par("dz_bbl_min") 
    dz_sed_min = get_brom_par("dz_sed_min") 
    dz_sed_max = get_brom_par("dz_sed_max") 
    k_min = get_brom_par("k_min") 
    k_wat_bbl = get_brom_par("k_wat_bbl") 
    k_points_below_water = get_brom_par("k_points_below_water")   
    i_min = get_brom_par("i_min") 
    i_water = get_brom_par("i_water") 
    i_max = get_brom_par("i_max")  
    year = get_brom_par("year")
! Uncomment  1 of 3 inputs of physical scenarios:  
! (1) - primitive physics  sinusoidal seasonal changes  
!     call input_primitive_physics (z_w, dz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_water, i_max)   
!         allocate(hice(365))
! (2)- input  physics  from ascii  
     call input_ascii_physics (z_w, dz_w, k_wat_bbl, water_layer_thickness, t_w, s_w, kz_w, i_water, i_max)   
         allocate(hice(365))
! (3) - input  physics  from netcdf i.e. for Kara Sea, kz - AKs, hice -  "time-averaged average ice thickness in cell" generates grid according to input data
!    call input_netcdf('KaraSea.nc', z_w, dz_w, k_wat_bbl, t_w, s_w, kz_w, hice, year, i_water, i_max)
!______________________________________________________________________________________________________________________
         k_max = k_wat_bbl + k_points_below_water !total number of vertical gridpoints
! generates vertical grid according to input data: z(), dz(), grid numbers for interfaces   
    call make_vert_grid (z, dz, z_w, dz_w, k_min,  k_wat_bbl, k_bbl_sed, k_max, bbl_thickness, dz_bbl_min, dz_sed_min, dz_sed_max)
! provide hydrophysical data (t, s, kz) and kz_bio for BBL and sediments
    call make_physics_bbl_sed (t, s, kz, t_w, s_w, kz_w, kz_bio, z, dz,  i_max, i_water, k_min, k_wat_bbl, k_bbl_sed, k_max)      
!  for Berre
     do i_day = 1, 365  
        kz(:,k_wat_bbl-1,i_day) = 0.5*(kz(:,k_wat_bbl-1,i_day)+kz(:,k_wat_bbl,i_day))
     enddo
! /for Berre    
! initialize FABM model from fabm.yaml
    call fabm_create_model_from_yaml_file(model)
    par_max = size(model%state_variables)
! set model domain    
    call fabm_set_domain(model,i_max,k_max)
! Specify vertical index of surface and bottom
    call model%set_surface_index(k_min+1)
    call model%set_bottom_index (k_max-1)        

    allocate(vv(i_max, k_max, 1))
    allocate(dvv(i_max, k_max, 1))
    allocate(cc(i_max,k_max, par_max))
    allocate(dcc(i_max,k_max, par_max))
    allocate(fick(i_max,k_max, par_max))
    allocate(wbio(i_max,k_max, par_max))    !vertical velocity (m/s, negative for sinking)
    allocate(surf_flux(i_max, par_max))     !surface flux (tracer unit * m/s, positive for tracer entering column)
    allocate(bound_up(i_max, par_max))
    allocate(bound_low(i_max,par_max))
    allocate(use_bound_up(i_max, par_max))
    allocate(use_bound_low(i_max, par_max))
    allocate(par_name(par_max))  
    
    use_bound_up  = .false.
    use_bound_low = .false.
    
    do ip = 1, par_max
        par_name(ip) = model%state_variables(ip)%name
    end do
! retrieve indexes of parameters needed in brom-transport and subroutines (i.e. boundary conditions description)
    id_O2  = find_index(par_name, 'niva_brom_bio_O2')              
    id_NO3 = find_index(par_name, 'niva_brom_bio_NO3')             
    id_NO2 = find_index(par_name, 'niva_brom_redox_NO2')            
    id_NH4 = find_index(par_name, 'niva_brom_bio_NH4')         
    id_DON = find_index(par_name, 'niva_brom_bio_DON')         
    id_PON = find_index(par_name, 'niva_brom_bio_PON')            
    id_PO4 = find_index(par_name, 'niva_brom_bio_PO4')      
    id_Si = find_index(par_name, 'niva_brom_redox_Si')      
    id_Sipart = find_index(par_name, 'niva_brom_redox_Sipart')          
    id_Mn4 = find_index(par_name, 'niva_brom_redox_Mn4')           
    id_Mn2 = find_index(par_name, 'niva_brom_redox_Mn2')          
    id_Mn3 = find_index(par_name, 'niva_brom_redox_Mn3')          
    id_MnS = find_index(par_name, 'niva_brom_redox_MnS')
    id_MnCO3 = find_index(par_name, 'niva_brom_redox_MnCO3')       
    id_Fe3 = find_index(par_name, 'niva_brom_redox_Fe3')          
    id_Fe2 = find_index(par_name, 'niva_brom_redox_Fe2')          
    id_FeS = find_index(par_name, 'niva_brom_redox_FeS')          
    id_SO4 = find_index(par_name, 'niva_brom_redox_SO4')          
    id_S2O3 = find_index(par_name, 'niva_brom_redox_S2O3')          
    id_S0 = find_index(par_name, 'niva_brom_redox_S0')           
    id_H2S = find_index(par_name, 'niva_brom_redox_H2S')         
    id_Phy = find_index(par_name, 'niva_brom_bio_Phy')          
    id_Het = find_index(par_name, 'niva_brom_bio_Het')          
    id_Baae = find_index(par_name, 'niva_brom_redox_Baae')           
    id_Bhae = find_index(par_name, 'niva_brom_redox_Bhae')          
    id_Baan = find_index(par_name, 'niva_brom_redox_Baan')          
    id_Bhan = find_index(par_name, 'niva_brom_redox_Bhan')          
    id_DIC = find_index(par_name, 'niva_brom_carb_DIC')          
    id_Alk = find_index(par_name, 'niva_brom_carb_Alk')          
    id_Hplus = find_index(par_name, 'niva_brom_redox_Hplus')          
    id_CaCO3 = find_index(par_name, 'niva_brom_redox_CaCO3')          
    id_FeS2 = find_index(par_name, 'niva_brom_redox_FeS2')
     
!Initial volumes of layers:
    do k = 1, k_max 
         vv(i_water, k, 1)=1. ! 10000.-50.*k
    end do       
    
    !Boudary conditions:
    !Where constant (Dirichlet) bcs are assumed, the constants are imported from brom.yaml
    !and the logical use_bound_up/low is set to .true.
    !Otherwise no-gradient (Neumann) bcs are assumed.
    !Note: logicals use_bound_up/low could be included in brom.yaml for more flexibility
    !      present code assumes that these choices will be infrequently changed
    if (id_O2 /= -1) then
        use_bound_low(i_water,id_O2) = .true.
        bound_low(i_water,id_O2) = get_brom_par("bound_lowO2") 
    end if
    if (id_SO4 /= -1) then
        use_bound_up(i_water,id_SO4) = .true.
        bound_up(i_water,id_SO4) = get_brom_par("bound_upSO4")     
        use_bound_low(i_water,id_SO4) = .true.
        bound_low(i_water,id_SO4) = get_brom_par("bound_lowSO4")
    end if
    if (id_Mn4/=-1) then
        use_bound_up(i_water,id_Mn4) = .true.
        bound_up(i_water,id_Mn4) = get_brom_par("bound_upMn4")
    end if
    if (id_Fe3/=-1) then
        use_bound_up(i_water,id_Fe3) = .true.
        bound_up(i_water,id_Fe3) = get_brom_par("bound_upFe3")
    end if
    if (id_ALK/=-1) then
        use_bound_up(i_water,id_ALK) = .true.
        bound_up(i_water,id_ALK) = get_brom_par("bound_upALK")
        use_bound_low(i_water,id_ALK) = .true.
        bound_low(i_water,id_ALK) = get_brom_par("bound_lowALK")
    end if
    if (id_DIC/=-1) then
        use_bound_low(i_water,id_DIC) = .true.
        bound_low(i_water,id_DIC) = get_brom_par("bound_lowDIC")
    end if
    if (id_NH4/=-1) then
        use_bound_low(i_water,id_NH4) = .true.
        bound_low(i_water,id_NH4) = get_brom_par("bound_lowNH4")    
    end if
    if (id_NO3/=-1) then
        use_bound_low(i_water,id_NO3) = .true.
        bound_low(i_water,id_NO3) = get_brom_par("bound_lowNO3")
    end if
    if (id_PO4/=-1) then
        use_bound_low(i_water,id_PO4) = .true.
        bound_low(i_water,id_PO4) = get_brom_par("bound_lowPO4")
    end if
    if (id_Si/=-1) then
        use_bound_low(i_water,id_Si) = .true.
        bound_low(i_water,id_Si) = get_brom_par("bound_lowSi")    
    end if
    if (id_H2S/=-1) then
        use_bound_low(i_water,id_H2S) = .true.
        bound_low(i_water,id_H2S) = get_brom_par("bound_lowH2S")    
    end if    
    if (id_NO3 /= -1) use_bound_up(i_water,id_NO3) = .true.
    if (id_PO4 /= -1) use_bound_up(i_water,id_PO4) = .true.
    if (id_Si /= -1)  use_bound_up(i_water,id_Si)  = .true.   
    
!point FABM to array slices with biogeochemical state.
    do ip = 1, par_max
        call fabm_link_bulk_state_data(model, ip, cc(:, :, ip))
    end do
    
    allocate(Izt(i_max,k_max))
    allocate(density(i_max,k_max))
    allocate(pressure(i_max,k_max))
    allocate(kz_tot(i_max,k_max))
    allocate(pco2atm(i_max))
    allocate(windspeed(i_max))
!provide initial array slices with temperature and salinity.
!these will be resent every time julianday is updated, below.
    call fabm_link_bulk_data(model, standard_variables%temperature, t(:,:,1))
    call fabm_link_bulk_data(model, standard_variables%practical_salinity, s(:,:,1))
!provide initial array slices with irraduance.
    call fabm_link_bulk_data(model, standard_variables%downwelling_photosynthetic_radiative_flux, Izt)           !W m-2
    !!!!!call fabm_link_horizontal_data(model, standard_variables%surface_downwelling_shortwave_flux, io)            !W m-2
!additional environmental variables needed by e.g., PML carbonate chemistry.
!!!!!    call fabm_link_bulk_data(model, standard_variables%density, density)                                        !kg m-3
    call fabm_link_bulk_data(model, standard_variables%pressure, pressure)                                      !dbar
    call fabm_link_horizontal_data(model, standard_variables%wind_speed, windspeed)                            !m s-1
    call fabm_link_horizontal_data(model, standard_variables%mole_fraction_of_carbon_dioxide_in_air, pco2atm)  !ppm
    
    call fabm_link_horizontal_data(model, standard_variables%ice_thickness, hice)                            !m
    call fabm_link_bulk_data(model, standard_variables%volume_of_cell, vv(:, :, 1)) 
    
    call fabm_check_ready(model)
    
!default environment in absence of actual data (could be retrievd from GOTM output)
    density    = get_brom_par("density")       ! kg m-3
    wind_speed = get_brom_par("wind_speed")    ! wind speed, m s-1
    pco2_atm   = get_brom_par("pco2_atm")      ! CO2 partical pressure,  ppm
    Io         = get_brom_par("Io")            ! W m-2 maximum surface irradiance at latitudes <= 23.5N,S
                                               !Note: should account for local cloudiness
    
!allow FABM models to use their default initialization (this sets cc)
    do k=1,k_max
        call fabm_initialize_state(model, 1, i_max, k)  !!    
    enddo
    do i=1,i_max
        pco2atm(i)   = pco2_atm       ! put pCO2 values to all the surface gridpoints
        windspeed(i) = wind_speed     ! put wind speed values to all the surface gridpoints
    enddo        
!read initials from file and save it inside cc massive (reset cc from start.dat)
 !   call porting_initial_state_variables_data(k_max, par_name, cc)
    call porting_initial_state_variables_n_volumes ('start_01.dat',year, i_day, i_max, k_max, par_max, par_name, z, dz, Cc, Vv, t, s, kz)
!        call saving_state_variables_n_volumes(model_year, julianday, i_max, k_max, par_max, par_name, z, dz, cc, vv, t, s, kz)i_day,  year    
    

!smoothing
    !call tem_sal_kz_smoothing(t, s, kz)
 !   call cc_smoothing(cc, k_max, par_max)
    
    !initializing output
    call init_netcdf("BROM_out.nc", i_max, k_max, model)
    
    end subroutine init_brom_transport
!=================================================================================    
    subroutine do_brom_transport()
    
    use calculate, only: calculate_phys, calculate_sed
    
    implicit none
    
    integer     :: id, idt, idf                  !time related parameters
    integer     :: julianday,last_day, model_year = 0
    real(8)     :: lat_light
    real(8)     :: kc                       !attenuation constant for the self shading effect
    real(8)     :: k_erlov                  !extinction coefficient
    real(8)     :: turbid                   ! parameter characterizing turbidity due to particled presence (for irradiance(z) calculations)
    real(8)     :: dt                       !time step in [1/day]
    integer     :: freq_turb                  !vert.turb. / bhc frequency
    integer     :: freq_sed                 !sinking / bhc frequency
    real(8), parameter :: pi=3.141592653589793_8
    
    
    lat_light = get_brom_par("lat_light")
    kc = get_brom_par("kc")
    k_erlov = get_brom_par("k_erlov")
    dt = get_brom_par("dt")
    freq_turb = get_brom_par("freq_turb")
    freq_sed  = get_brom_par("freq_sed ")
    last_day = get_brom_par("last_day")
    !open(5, file = 'read_from_nc.dat')
    !    do i=1, 365
    !        write (5,('(i5, 4(1x,f8.3))')) i,hice(i),t(i_water,i), s(i_water,i), kz(i_water,i)
    !    enddo
    !close (5)
    i=i_water

    idt = int(1. / dt)                                      !number of cycles per day
    
    do  i_day = 0, last_day - 1                          !BIG Cycle ("i_day"-days)
        julianday = max(1, i_day - int(i_day / 365) * 365 + 1)
        if (julianday == 365) then
            model_year = model_year + 1
  !          call saving_state_variables_data(model_year, julianday, k_max, par_max, par_name, z, cc)
 !ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ           
            !here save together with vv
 !ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ            
        end if
        !compute surface irradiance
 !!       io = max (0., 10. * cos((lat_light - (23.5 * sin(2. * 3.14 * (julianday - 81.) /365.))) * 3.14 / 180.) ) - from eya
        !compute irradiance at depth  80* cos((40. - (23.5 * sin(2. * 3.14 * (X - 81.) /365.))) * 3.14 / 180.)
        Iot = max(0., 10.*cos((lat_light-(23.5*sin(2.*pi*(julianday-81.)/365.)))*pi/180.))
 !!       Iot = max(0., Io*cos((lat_light-(23.5*sin(2.*pi*(julianday-81.)/365.)))*pi/180.)) - from phil
 !!       !This is the approximation used in Yakushev and Sorensen (2013) for OXYDEP
        
        !compute irradiance at depth
        turbid = 0.
        do k = 1, (k_max - 1)
            !we check that we are not in the sediments
            !turbidity caused by particles in the water above
            !irradiance changing with depth
            if (kz(i_water,k, julianday) > 5E-9) then
                turbid = turbid + (cc(i,k,id_Phy) + cc(i,k,id_PON) + cc(i,k,id_Baae) + cc(i,k,id_Bhae)   &
                    + cc(i,k,id_Baan) + cc(i,k,id_Bhan) + cc(i,k,id_Mn4) + cc(i,k,id_FeS) + cc(i,k,id_S0) &
                    + cc(i,k,id_Sipart) + cc(i,k,id_CaCO3) + cc(i,k,id_MnCO3)) * dz(k)
                Izt(i_water,k) = Iot*exp(-k_erlov*z(k))  &
                    * exp(-kc * turbid)   ! phytoplankton autoshading  +  particles shading
!                   * max(0.,(0.6-hice(julianday))) !  restricted by ice with thickness 0.5 m
            else
                Izt(i_water,k) = 0.
            end if
        end do
        
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ            
!            Calaulate changes of volumes vv as a function of changing salts, i.e.
            !i=3  !solid
            ! dvv(i)= dcc(i,k, ip_NaCl)/dens(ip_NaCl)
            !I=2  !water
            ! dvv(i)= dcc(i,k, ip_NaCl)
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ 
        Vv=0.         
!        Vv(1,k_min+1, 1)= max(0.,hice ( julianday))   ! i.e. Vv() is an amount of solid matter in the cell
 !       Vv(i_water,k_min+1, 1)= 0. !-Vv(1,k_min+1, 1)
        
        if (id_PO4 /= -1) bound_up(i_water,id_PO4) = 0. + (1. + sin(2 * pi * (julianday - 60.) / 365.)) * 1.  !1.25 !; % max.conc=1 uM 
        if (id_NO3 /= -1) bound_up(i_water,id_NO3) = 0. + (1. + sin(2 * pi * (julianday - 60.) / 365.)) * 3.   !; % max.conc=1 uM       
        if (id_Si  /= -1) bound_up(i_water,id_Si)  = 0. + (1. + sin(2 * pi * (julianday - 60.) / 365.)) * 3.   
!resend data that depend on julianday to FABM
        call fabm_link_bulk_data(model, standard_variables%temperature, t(:,:,julianday))
        call fabm_link_bulk_data(model, standard_variables%practical_salinity, s(:,:,julianday))
!call fabm_link_horizontal_data(model, standard_variables%ice_thickness, hice(julianday))

!-------NETCDF-----------------------------------------------------------------------------------------    
        write (*,'(a, i4, a, i4)') " model year:", model_year, "; julianday:", julianday   
        call saving_state_variables_n_volumes('finish_01.dat',model_year, julianday, i_max, k_max, par_max, par_name, z, dz, cc, vv, t, s, kz)    
        call save_netcdf(i_max, k_max, julianday, cc, t, s, kz, model, -z, hice)
        !if (i == last_day - 1) then
        !    call close_netcdf()
        !end if
!---END-of-NETCDF--------------------------------------------------------------------------------------

        do id = 1, idt  ! timesteps in the course of a day
 !_________________________________           
 !_______biogeochemistry___________!    
            
                 dcc = 0.
            do  k = 2, (k_max - 1)
                 call fabm_do(model,i_water, i_water, k, dcc(i_water:i_water,k, :))                            
            end do	
   
            do ip= 1, par_max 
                do  k = 2, (k_max - 1)  
                    cc(i_water, k, ip) = max(0.00000000001, (cc(i_water, k, ip) + dt * dcc(i_water, k, ip) * 86400.))
                end do	
            end do		
! in the solid column           
!                 dcc = 0.
!                 k=2
!!            do  k = 2, (k_max - 1)
!                 call fabm_do(model,2, 2, k, dcc(1:1,k, :))                            
!!            end do	
!   
!            do ip= 1, par_max 
!!                do  k = 2, (k_max - 1)  
!                    cc(1, k, ip) = max(0.00000000001, (cc(1, k, ip) + dt * dcc(1, k, ip) * 86400.))
! !               end do	
!            end do		            
 !_________________________________           
 !_______vertical diffusion________!

            do  idf = 1, freq_turb
                 surf_flux=0.
                call fabm_do_surface(model,  i_water, i_water,  surf_flux(i_water:i_water,:))                    
                 surf_flux = surf_flux*86400. !           NB FABM rates are in s-1
                call calculate_phys(i_water, k_max, par_max,  cc, fick, dcc, use_bound_up, use_bound_low, bound_up, bound_low, &
                                surf_flux, k_bbl_sed, dz, kz, kz_bio, kz_tot, julianday,  id_O2,  dt, freq_turb)     
            end do  
            
 !_________________________________           
 !_______Horizontal excange_________!      
                  fick(1, k_min+1, :) = 0.           
!            do ip= 1, par_max 
!!                do  k = 2, (k_max - 1)  
! !                      k=2
!                    if (ip.ne.id_Het.and.ip.ne.id_Phy) then
!                        fick(1, k_min+1, ip)     =  1.e-1 *Vv(1,k_min+1, 1)*(-cc(1, k_min+1, ip) + cc(i_water, k_min+1, ip))
!                        cc(1, k_min+1, ip)       =  cc(1, k_min+1, ip) + dt*fick(1, k_min+1, ip) 
!                        
!                   !     dcc(i_water, k_min+1, ip) = 1.e-9 *Vv(1,k_min+1, 1)*( cc(1, k_min+1, ip) - cc(i_water, k_min+1, ip))                        
!                        cc(i_water, k_min+1, ip) =  cc(i_water, k_min+1, ip) -dt*fick(1, k_min+1, ip)                       
!                    endif
! !!               end do	
!            end do	          
 !_________________________________           
 !_______Particles sinking_________!             
            !compute residual vertical velocity (sinking/floating) in FABM.
            do k=k_min+1,k_max-1
                call fabm_get_vertical_movement(model,  i_water, i_water,  k, wbio(i_water:i_water,k,:))  !  use or store the returned vertical velocities    
            enddo
            wbio= wbio*86400.             !NB FABM rates are in s-1
            !sedimentation of particulate matter
            do  ip = 1, freq_sed
                call calculate_sed (i_water, k_max, par_max, cc, dcc, dz, kz, wbio, julianday, dt, freq_sed) 
            end do
!22 continue            
            if (any(isnan(cc(:,2:k_max-1,:)))) then
                stop
            end if
        end do
    end do
        
    end subroutine do_brom_transport
!=================================================================================    
    subroutine clear_brom_transport()
    
    deallocate(cc)
    deallocate(dcc)
    deallocate(vv)
    deallocate(dvv)        
    deallocate(fick)
    deallocate(wbio) 
    deallocate(surf_flux)
    deallocate(bound_up)
    deallocate(bound_low)
    deallocate(use_bound_up)
    deallocate(use_bound_low)
    deallocate(par_name)
    deallocate(Izt)
    deallocate(density)
    deallocate(pressure)
    deallocate(kz_tot)
    
    end subroutine clear_brom_transport
!=================================================================================    
    end module brom_transport