    module io_ascii
    
    use fabm_driver
    use fabm_config_types
    use fabm_yaml,yaml_parse=>parse,yaml_error_length=>error_length

    implicit none
    private
    
    type brom_par ! Type for sparse matrix
        real(8) :: realvalue
        character(len=64) :: initname
        type(brom_par), pointer :: next
    end type brom_par
    
    type (brom_par), pointer :: first
    
    public find_index, porting_initial_state_variables_data, init_common, saving_state_variables_data, &
           get_brom_par, svan, saving_state_variables_n_volumes, porting_initial_state_variables_n_volumes, &
           make_vert_grid, input_primitive_physics,input_ascii_physics, make_physics_bbl_sed

    contains
!=================================================================================    
    subroutine init_common()
        
    class (type_node), pointer       :: node
    character(len=yaml_error_length) :: yaml_error
    integer                          :: unit_eff
    character(len=256)               :: path_eff
         
    unit_eff = 1
    path_eff = 'brom.yaml'
    ! Parse YAML file.
    node => yaml_parse(trim(path_eff),unit_eff,yaml_error) 
    if (yaml_error/='') call fatal_error('fabm_create_model_from_yaml_file',trim(yaml_error))
    if (.not.associated(node)) call fatal_error('fabm_create_model_from_yaml_file', &
         'No configuration information found in '//trim(path_eff)//'.')
    select type (node)
         class is (type_dictionary)
            call tree(node)
         class is (type_node)
            call fatal_error('brom', trim(path_eff)//' must contain a dictionary &
               &at the root (non-indented) level, not a single value. Are you missing a trailing colon?')
      end select
    
    end subroutine init_common
!=================================================================================    
    subroutine tree(mapping)
    
    class (type_dictionary), intent(in) :: mapping
    class (type_node), pointer          :: node
    type (type_key_value_pair), pointer :: pair
    character(len=64)                   :: instancename
    
    node => mapping%get('instances')
    if (.not.associated(node)) &
        call fatal_error('create_model_tree_from_dictionary', 'No "instances" dictionary found at root level.')
    select type (node)
        class is (type_dictionary)
            pair => node%first
        class is (type_node)
            nullify(pair)
            call fatal_error('create_model_tree_from_dictionary',trim(node%path)// &
               ' must be a dictionary with (model name : information) pairs, not a single value.')
    end select
        
    ! only brom acceptance, iterate throw other models with no action
    do while (associated(pair))
        instancename = trim(pair%key)
        if (instancename == "brom") then
            select type (dict=>pair%value)
                class is (type_dictionary)
                    call from_tree(instancename, dict)
                class is (type_null)
                    call fatal_error('create_model_tree_from_dictionary','Configuration information for model "'// &
                        trim(instancename)//'" must be a dictionary, not a single value or nothing')
                class is (type_node)
                    call fatal_error('create_model_tree_from_dictionary','Configuration information for model "'// &
                        trim(instancename)//'" must be a dictionary, not a single value.')
                end select
        end if
        pair => pair%next
    end do
    
    end subroutine tree
!=================================================================================    
    subroutine from_tree(instancename, node)
    
    character(len=*),       intent(in)          :: instancename
    class (type_dictionary),intent(in)          :: node
    class (type_dictionary),pointer             :: childmap
    type (type_error),pointer                   :: config_error
    
    childmap => node%get_dictionary('initialization',required=.false.,error=config_error)
    if (associated(config_error)) call fatal_error('create_model_from_dictionary',config_error%message)
    if (associated(childmap)) call parse_initialization(childmap)
      
    end subroutine from_tree
!=================================================================================    
    subroutine parse_initialization(node)
    
    class (type_dictionary),intent(in)  :: node

    type (type_key_value_pair), pointer :: pair
    type (brom_par), pointer :: current
    
    logical                             :: success
    real(8)                             :: realvalue
    character(len=64)                   :: initname
        
    ! Transfer user-specified initial state to the model.
    nullify(first)
    pair => node%first
    do while (associated(pair))
        select type (value=>pair%value)
            class is (type_scalar)
                initname = trim(pair%key)
                realvalue = value%to_real(default=real(0,real_kind),success=success)
                if (.not.success) call fatal_error('parse_initialization', &
                    trim(value%path)//': "'//trim(value%string)//'" is not a real number.')
            class is (type_null)
                call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to null.')
            class is (type_dictionary)
                call fatal_error('parse_initialization',trim(value%path)//' must be set to a real number, not to a dictionary.')
            end select
        allocate (current)
        current = brom_par(realvalue, initname, first)
        first => current
        pair => pair%next
    end do

    end subroutine parse_initialization
!=================================================================================    
    function get_brom_par(initname) result (realvalue)
    
    real(8)                             :: realvalue
    character(len=*), intent(in)        :: initname
    type (brom_par), pointer            :: current
    
    current => first
    do
        if (.not.associated(current)) exit
        if (trim(initname) == trim(current%initname)) then
            realvalue = current%realvalue
            return
        end if
        current => current%next
    end do
    
    print *, "Check brom.yaml or name of the input variable please"
    stop
    
    end function get_brom_par
!=================================================================================    
    function find_index(names, name) result(index)
    
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: name
    integer :: index
    
    do index = 1, size(names)
        if (names(index) == name) return
    end do
    index = -1
    
    end function
!=================================================================================    
    subroutine porting_initial_state_variables_n_volumes (file_name, dYear, julianday, i_max, k_max, par_max, par_name, z, dz, Cc, Vv, t, s, kz)

    character (len = *), intent(in)             :: file_name    
    integer, intent(in)                         :: dYear, julianday, i_max, k_max, par_max
    real(8), dimension(:, :, :), intent(out)    :: Cc
    real(8), dimension(:, :, :), intent(out)    :: Vv
    real(8), dimension (:, :, :), pointer, intent(out)      :: t, s, Kz    
    character(len = *), intent(in)              :: par_name(:)
    real(8), intent(in)                         :: z(:), dz(:)
    integer                                 :: d_year_start, julianday_start, k_max_start, par_max_start
    integer, allocatable                    :: column2state(:)
    character(20000)                        :: labels
    character(200)                          :: comments
    integer                                 :: i, j, k, ip, istart, istop, foo, m
    real(8)                                 :: value    

    open(9, file = file_name) !'start_u_01.dat')
    read(9, '( 1x,i4,1x,i5,1x,i5,1x,i5 )') d_year_start, julianday_start, k_max_start, par_max_start
    if (k_max_start /= k_max) stop 'number of levels in start-up file does not match number of levels in model'
!    allocate(column2state(par_max_start))
    allocate(column2state(par_max))    
    read(9, '( 1x, a )') labels
    istart = 1
    do ip = 1, par_max_start
        istop = istart - 1 + index(labels(istart:), ' ')
        column2state(ip) = find_index(par_name, labels(istart:istop-1))
        istart = istop + 1
    end do
    read (9, *) comments
    do i=1, i_max
      do k = 1, k_max
        read(9, '( i5 )', advance = 'no') foo
        read(9, '( i5 )', advance = 'no') foo
        read(9, '( 1x,f8.4 )', advance = 'no') value !we don't read again z()
       do m=1,4  
        read(9, '( 1x,f15.9 )', advance = 'no') value !we don't read again t,s,Kz, dz
       enddo    
            read(9, '( 1x,f15.9 )', advance = 'no') Vv(i,k,1)
        do ip = 1, par_max
            read(9, '( 1x, f15.9 )', advance = 'no') value
            if (ip /= -1) then
                if (value /= 0.) then
                    Cc(i,k,ip) = value                    
                else
                    cc(i,k,ip) = 0.0001
                end if
            end if
        end do
        read(9, *)
      end do
    enddo
    close (9)
    end subroutine
!=================================================================================    
    subroutine porting_initial_state_variables_data(k_max, par_name, cc)
    
    integer, intent(in)                     :: k_max
    character(len = *), intent(in)          :: par_name(:)
    real(8), dimension(:, :), intent(out)   :: cc
    integer                                 :: d_year_start, julianday_start, k_max_start, par_max_start
    integer, allocatable                    :: column2state(:)
    character(20000)                        :: labels
    character(200)                          :: comments
    integer                                 :: i, j, istart, istop, foo
    real(8)                                 :: value
    
    open(9, file = 'start10.dat')
    read(9, '( 1x,i4,1x,i5,1x,i5,1x,i5 )') d_year_start, julianday_start, k_max_start, par_max_start
    if (k_max_start /= k_max) stop 'number of levels in start-up file does not match number of levels in model'
    allocate(column2state(par_max_start))

    !"labels" is a space-separated string with state variable names.
    !find indiviudal names and look up their state variable index.
    read(9, '( 1x, a )') labels
    istart = 1
        do i = 1, par_max_start
        istop = istart - 1 + index(labels(istart:), ' ')
        column2state(i) = find_index(par_name, labels(istart:istop-1))
        istart = istop + 1
    end do
    read (9, *) comments
    do i = 1, k_max_start
        read(9, '( i5 )', advance = 'no') foo
        read(9, '( 1x,f8.4 )', advance = 'no') value
        do j = 1, par_max_start
            read(9, '( 1x, f15.9 )', advance = 'no') value
            if (column2state(j) /= -1) then
                if (value /= 0.) then
                    cc(i, column2state(j)) = value
                else
                    cc(i, column2state(j)) = 0.0001
                end if
            end if
        end do
        read(9, *)
    end do
    
    close (9)
    
    end subroutine porting_initial_state_variables_data
!=================================================================================    
    subroutine saving_state_variables_data(dYear, julianday, k_max, par_max, ParName, z, Cc)
    
    integer, intent(in)                         :: dYear, julianday, k_max, par_max
    real(8), dimension(:, :), intent(in)        :: Cc
    character(len = *), intent(in)              :: ParName(:)
    real(8), intent(in)                         :: z(:)
    
    integer :: ip, k
    
    open(10,FILE = 'finish10.dat') 

    ! First line with year, julian day, number of levels, number of state variables
    write(10,'( i4,1x,i5,1x,i5,1x,i5 )') dYear, julianday, k_max, par_max
    do ip=1, par_max
        write(10,'(1x,a)',advance='NO') trim(ParName(ip))
    end do
    write (10,*)
    
    write(10,'(3h n ,6hDepth )',advance='NO')
    do ip=1, par_max
    !       write(10,'(1x,a)',advance='NO') trim(ParName(ip)(len_trim(ParName(ip))-3:))
        write(10,'(1x,a)',advance='NO') trim(ParName(ip)(15:))
    end do
    write (10,*)

    ! Subsequent lines: one for each depth level
    do k=1, k_max
        write(10,'( i5 )',advance='NO') k
        if (k.eq.1) then
            write(10,'(1x,f8.4 )',advance='NO') 0.
        else
            write(10,'(1x,f8.4 )',advance='NO') z(k-1)
        end if
        do ip=1,par_max
            write(10,'( 1x,f15.9 )',advance='NO') Cc(k,ip)
        end do
        write (10,*)
    end do
    close (10)
    
    end subroutine saving_state_variables_data
 !=================================================================================    
    subroutine saving_state_variables_n_volumes (file_name, dYear, julianday, i_max, k_max, par_max, ParName, z, dz, Cc, Vv, t, s, kz)

    character (len = *), intent(in)             :: file_name      
    integer, intent(in)                         :: dYear, julianday, i_max, k_max, par_max
    real(8), dimension(:, :, :), intent(in)     :: Cc
    real(8), dimension(:, :, :), intent(in)        :: Vv
    real(8), dimension (:, :, :), pointer, intent(out)      :: t, s, Kz    
    character(len = *), intent(in)              :: ParName(:)
    real(8), intent(in)                         :: z(:), dz(:)
    
    integer :: ip, k, i
    
    open(10,FILE = file_name) 

    ! First line with year, julian day, number of levels, number of state variables
    write(10,'( i4,1x,i5,1x,i5,1x,i5 )') dYear, julianday, k_max, par_max
    do ip=1, par_max
        write(10,'(1x,a)',advance='NO') trim(ParName(ip))
    end do
    write (10,*)
    
    write(10,'(3h i ,3h k ,6hDepth ,4h t  ,4h s   ,4h Kz ,4h dz ,4hVol )',advance='NO')
    
    do ip=1, par_max
    !       write(10,'(1x,a)',advance='NO') trim(ParName(ip)(len_trim(ParName(ip))-3:))
        write(10,'(1x,a)',advance='NO') trim(ParName(ip)(15:))
    end do
    write (10,*)

    ! Subsequent lines: one for each depth level
    do i=1, i_max
        do k=1, k_max
            write(10,'( i5 )',advance='NO') i
        
            write(10,'( i5 )',advance='NO') k
            if (k.eq.1) then
                write(10,'(1x,f8.4 )',advance='NO') 0.
            else
                write(10,'(1x,f8.4 )',advance='NO') z(k)
            end if
            
                write(10,'( 1x,f15.9 )',advance='NO') t(i,k,julianday)            
                write(10,'( 1x,f15.9 )',advance='NO') s(i,k,julianday)            
                write(10,'( 1x,f15.9 )',advance='NO') kz(i,k,julianday)    
                write(10,'( 1x,f15.9 )',advance='NO') dz(k)         
                write(10,'( 1x,f15.9 )',advance='NO') Vv(i,k,1)
            do ip=1,par_max
                write(10,'( 1x,f15.9 )',advance='NO') Cc(i,k,ip)
            end do
            write (10,*)
        end do
    enddo
    
    close (10)
    
    end subroutine saving_state_variables_n_volumes   
!=================================================================================    
    subroutine make_vert_grid (z, dz, z_w, dz_w, k_min,  k_wat_bbl, k_bbl_sed, k_max, BBL_thickness, & 
                               dz_bbl_min, dz_sed_min, dz_sed_max ) 
    
    real(8), dimension (:), pointer, intent(out) :: z, dz 
    integer, intent(out)                         :: k_bbl_sed
    real(8), dimension (:), pointer, intent(in)  :: z_w, dz_w 
    integer, intent(in)                          :: k_min, k_max, k_wat_bbl 
    integer                                      :: k
    real(8)                                      :: water_layer_thickness != 12.   ! water column heigh, m
    real(8)                                      :: BBL_thickness != 12.   ! sediment layer thickness, m
    real(8)                                      :: dz_bbl_min != 0.03   ! minimum dz in the BBL, m
    real(8)                                      :: dz_sed_min != 0.0005 ! minimum dz in the sediments   
    real(8)                                      :: dz_sed_max != 0.01   ! maximum dz in the ediment        
    real(8)                                      :: z_correction

    allocate(dz(k_max))    
    allocate(z(k_max))    
!____________________________________________________________________       
!  define grid points depths in [m] in the water column...    
    do k=1,k_wat_bbl  
      dz(k) = dz_w(k) 
      z(k)  = z_w(k)
    enddo
!  define grid points depths in [m] in the BBL ... 
    do  k=k_wat_bbl+1,k_max 
      dz(k-1) = BBL_thickness/(2.**(k-k_wat_bbl))
      z (k) = z(k-1) + dz(k-1)    
      if (dz(k-1)<dz_bbl_min) then ! we decrease twice dz() since it will  be 3 cm ....  
         k_bbl_sed=k
         goto 1 
      endif
    enddo
1   continue    
!  define grid points depths in [m] in the sediment...     
    do  k=k_bbl_sed+1,k_max 
       dz(k-1)= min(dz_sed_min*(2**(k-k_bbl_sed-1)),dz_sed_max)
       z(k)  =  z(k-1) + dz(k-1)          
    enddo    
!---------------------------      
!correct depths for a beautiful number at the SWI
        z_correction = z(k_bbl_sed)-int(z(k_bbl_sed))
     do k=1,k_max  
        z(k)  =  z(k) - z_correction
     enddo
!____________________________________________________________________    
     open (10,FILE = 'Vertical_grid.dat') 
           do k=1,k_max
                    write (10,'(1x,i4, 2(1x,f8.4))') k,dz(k),z(k)
           enddo
     close (10)    
    
    end subroutine make_vert_grid 
!=================================================================================    
    subroutine input_primitive_physics (z, dz, k_max, water_layer_thickness, t, s, Kz, i_water, i_max)
    
!!- THIS IS THE FOR THE SIMPLIEST "HYDROPHYSICAL SCENARIO" WITH SINUSOIDAL SEASONAL CHANGES
    real(8), dimension (:, :, :), pointer, intent(out)      :: t, s, Kz
    real(8), dimension (:), pointer, intent(out)            :: z, dz
    integer, intent(in)                                     :: k_max
    integer, intent(in)                                     :: i_water, i_max
    real(8), intent(in)                                    :: water_layer_thickness != 12.   ! 12 meters for maximum lake depth
   
    real(8), allocatable, dimension(:)                      :: sal_sum, sal_win, tem_sum, tem_win, dens
    integer                                                 :: i, k, iday, k_min, k_wat_bbl !, k_bbl_sed
    integer, parameter                                      :: days_in_yr=365   ! k_max= 12
         
    allocate(sal_sum(k_max))
    allocate(sal_win(k_max))
    allocate(tem_sum(k_max))
    allocate(tem_win(k_max))    
    allocate(dens(k_max))      
    allocate( z(k_max))   
    allocate( dz(k_max))       
    allocate( t(i_max,k_max,days_in_yr))      
    allocate( s(i_max,k_max,days_in_yr))      
    allocate( Kz(i_max,k_max,days_in_yr)) 
    
    dz(1)=water_layer_thickness/k_max
    z(1)=0.
     do k = 2, k_max    
       z(k)=z(k-1)+dz(k-1)
       dz(k)=dz(k-1)
     enddo
!____________________________________________________________________       
! summer and winter vert. distributions     
    do k = 1, k_max !k_wat_bbl
            sal_sum(k)= 35.+0.01*k
            sal_win(k)= 35.+0.001*k       
            tem_sum(k)= 25.-0.1*k
            tem_win(k)= 10.-0.001*k               
    enddo

! seasonal variability of T and S       
    do iday = 1, days_in_yr
        do k = 1, k_max !k_wat_bbl
            s(i_water,k,iday)=sal_win(k)+0.5*(1.+sin(2*3.14*(iday+240.)/365.))*(sal_sum(k)-sal_win(k))
            t(i_water,k,iday)=tem_win(k)+0.5*(1.+sin(2*3.14*(iday+240.)/365.))*(tem_sum(k)-tem_win(k))
        enddo
    enddo
! Calculate  vertical turbulent coefficient Kz from denity (Gargett)
     ! K_z = 1.62*10-3 * 1./( (-g/ro) * (d ro/d z) )*0.5 ! sudya po vsemu cm2/s  (Samodurov, Lyubitsky, Panteleev, 1994). 
    do iday = 1, days_in_yr
        do k = 1, k_max !k_wat_bbl   
             CALL svan( s(i_water,k,iday),  t(i_water,k,iday),   z(k),   dens(k)) !calculate density as f(p,t,s)
        end do   
        do k = 1, k_max-1 !k_wat_bbl-1    
             Kz(i_water,k,iday)=0.5E-6 /&		!0.25 !0.30         
		        ( &
		        (9.81/(1000.+(dens(k)+dens(k+1))/2.)& 
		        *(abs(dens(k+1)-dens(k))/dz(k)) &
		        )**0.5)  ! m2/s due to /10000
        end do   
    enddo
    
     open (12,FILE = 'Hydrophysics4.dat') 
       do iday=1,days_in_yr
           do k=1,k_max
                    write (12,'(2(1x,i4))',advance='NO') iday, k
               do i=1,i_max
                    write (12,'(1x,i4, 2(1x,f8.4),1x,f15.11,2f7.3)',advance='NO') i, t(i_water,k,iday),s(i_water,k, iday),Kz(i_water,k,iday),dz(k),z(k)
               enddo
                    write (12,*)
           enddo
       enddo
     close (12)    
    
    end subroutine input_primitive_physics 
!=================================================================================    
    subroutine input_ascii_physics (z, dz, k_max, water_layer_thickness, t, s, Kz, i_water, i_max)
   
!!- THIS IS THE FOR INPUT FROM ASCII FILES for T,S,(and Kz)
    real(8), dimension (:, :, :), pointer, intent(out)      :: t, s, Kz
    real(8), dimension (:), pointer, intent(out)            :: z, dz
    integer, intent(in)                                     :: k_max
    integer, intent(in)                                     :: i_water, i_max
    real(8), intent(in)                                    :: water_layer_thickness != 12.   ! 12 meters for maximum lake depth
   
    real(8), allocatable, dimension(:)                      :: sal_sum, sal_win, tem_sum, tem_win, dens
    integer                                                 :: i, j, k, iday, k_min, k_wat_bbl , k_max_TEL
    integer, parameter                                      :: days_in_yr=365   
    real(8)     :: tem_TEL(55,365), sal_TEL(55,365), kz_TEL(55,365) ! this is for TELEMARK model outputs for Berre Lagoon
    
    allocate(sal_sum(k_max))
    allocate(sal_win(k_max))
    allocate(tem_sum(k_max))
    allocate(tem_win(k_max))    
    allocate(dens(k_max))      
    allocate( z(k_max))   
    allocate( dz(k_max))       
    allocate( t(i_max,k_max,days_in_yr))      
    allocate( s(i_max,k_max,days_in_yr))      
    allocate( Kz(i_max,k_max,days_in_yr)) 
    
    dz(1)=water_layer_thickness/k_max
    z(1)=0.
     do k = 2, k_max    
       z(k)=z(k-1)+dz(k-1)
       dz(k)=dz(k-1)
     enddo
!____________________________________________________________________       
! read data for Berre 
     k_max_TEL= 35
    open (19, file='tem2_new.dat')   ! open (19, file='tem2_in.dat')
    do i = 1, days_in_yr
        do j = 1, k_max_TEL
            read (19, *) tem_TEL(j, i) ! measured and interpolated data from TELEMARC model
        end do
    end do
    close (19)

    open (20, file='sal2_new.dat')     !  open (20, file='sal2_in.dat')
    do i = 1, days_in_yr
        do j = 1, k_max_TEL
            read (20, *) sal_TEL(j, i) ! measured and interpolated data from TELEMARC model
        end do
    end do
    close (20)
    
    !open (21, file='kz2_in.dat')
    !do i = 1, days_in_yr
    !    do j = 1, k_max_TEL
    !        read (21, *) kz_TEL(j, i) ! measured and interpolated data from TELEMARC model
    !    end do
    !end do
! seasonal variability of T and S       
    do iday = 1, days_in_yr
        do k = 1, k_max !k_wat_bbl
            s(i_water,k,iday)= sal_TEL(k,iday) 
            t(i_water,k,iday)= tem_TEL(k,iday) 
        enddo
    enddo    
! Calculate  vertical turbulent coefficient Kz from denity (Gargett)
     ! K_z = 1.62*10-3 * 1./( (-g/ro) * (d ro/d z) )*0.5 ! sudya po vsemu cm2/s  (Samodurov, Lyubitsky, Panteleev, 1994). 
    do iday = 1, days_in_yr
        do k = 1, k_max    
             CALL svan( s(i_water,k,iday),  t(i_water,k,iday),   z(k),   dens(k)) !calculate density as f(p,t,s)
        end do   
        do k = 1, k_max-1    
             Kz(i_water,k,iday)=0.5E-6 /&		!0.25 !0.30         
		        ( &
		        (9.81/(1000.+(dens(k)+dens(k+1))/2.)& 
		        *(abs(dens(k+1)-dens(k))/dz(k)) &
		        )**0.5)  ! m2/s due to /10000
        end do   
    enddo
    
     open (12,FILE = 'Hydrophysics_Berre_TEL.dat') 
       do iday=1,days_in_yr
           do k=1,k_max
                    write (12,'(2(1x,i4))',advance='NO') iday, k
               do i=1,i_max
                    write (12,'(1x,i4, 2(1x,f8.4),1x,f15.11,2f7.3)',advance='NO') i, t(i_water,k,iday),s(i_water,k, iday),Kz(i_water,k,iday),dz(k),z(k)
               enddo
                    write (12,*)
           enddo
       enddo
     close (12)    

    end subroutine input_ascii_physics 
!=================================================================================    
    subroutine make_physics_bbl_sed (t, s, kz, t_w, s_w, kz_w, kz_bio, z, dz,  i_max, i_water, k_min, k_wat_bbl, k_bbl_sed, k_max)
!!- THIS IS THE FOR THE SIMPLIEST "HYDROPHYSICAL SCENARIO" WITH SINUSOIDAL SEASONAL CHANGES
    real(8), dimension (:, :, :), pointer, intent(out)      :: t, s, kz, t_w, s_w, kz_w
    real(8), dimension (:), pointer, intent(out)            :: z, dz, kz_bio
    real(8), allocatable, dimension(:)                      :: sal_sum, sal_win, tem_sum, tem_win, dens
    integer                                                 :: i, k, iday, k_min, k_wat_bbl, k_bbl_sed, k_max, i_max, i_water
    real(8), parameter                                      :: depth_range= 12.   ! 12 meters for maximum lake depth
    real(8), parameter                                      :: dz_bbl_min= 0.03   ! minimum dz in the BBL
    real(8), parameter                                      :: dz_sed_min= 0.0005 ! minimum dz in the sediments   
    real(8), parameter                                      :: dz_sed_max= 0.01   ! maximum dz in the ediment        
    integer, parameter                                      :: days_in_yr=365   ! k_max= 12
         
    allocate(sal_sum(k_max))
    allocate(sal_win(k_max))
    allocate(tem_sum(k_max))
    allocate(tem_win(k_max))    
    allocate(dens(k_max))      
    allocate(kz_bio(k_max))      
    allocate( t(i_max,k_max,days_in_yr))      
    allocate( s(i_max,k_max,days_in_yr))      
    allocate(kz(i_max,k_max,days_in_yr)) 
!____________________________________________________________________   
    
    i=i_water
    
    do iday = 1, days_in_yr
        !do k = 1, k_wat_bbl   
        !     CALL svan( s(i_water,k,iday),  t(i_water,k,iday),   z(k),   dens(k)) !calculate density as f(p,t,s)
        !end do   
        do k = 1, k_wat_bbl  
            s(i_water,k,iday)  = s_w(i_water,k,iday)
            t(i_water,k,iday)  = t_w(i_water,k,iday)  
            kz(i_water,k,iday) = kz_w(i_water,k,iday)       
            
            !s(i_water,k,iday)  = s_w(1,k,iday)
            !t(i_water,k,iday)  = t_w(1,k,iday)  
            !kz(i_water,k,iday) = kz_w(1,k,iday)                   
             kz_bio=0.
        enddo   
!            kz(i_water,k_wat_bbl,iday) = (kz_w(i_water,k_wat_bbl-1,iday)+1.e-6)/2.        
            
        do k = k_wat_bbl,k_bbl_sed-1   
            s(i_water,k,iday)=s(i_water,k_wat_bbl,iday)
            t(i_water,k,iday)=t(i_water,k_wat_bbl,iday)  
            kz(i_water,k,iday)=1.e-6  
             kz_bio=0.
        enddo   
        do k = k_bbl_sed,k_max   
            s(i_water,k,iday)=s(i_water,k_wat_bbl,iday)
            t(i_water,k,iday)=t(i_water,k_wat_bbl,iday)  
            kz(i_water,k,iday)=1.e-11  
            kz_bio=1.e-10  ! bioturbation exist in all the layer
        enddo           
    enddo
    do iday = 1, days_in_yr    
        do k = 1,k_max      
            s(1,k,iday)  = s(i_water,k,iday)
            t(1,k,iday)  = t(i_water,k,iday)  
            kz(1,k,iday) = kz(i_water,k,iday)
        enddo       
    enddo
    
     open (12,FILE = 'Hydrophysics3.dat') 
       do iday=1,days_in_yr
           do k=1,k_max
                    write (12,'(2(1x,i4))',advance='NO') iday, k
               do i=1,i_max
                    write (12,'(1x,i4, 2(1x,f8.4),1x,f15.11,2f7.3)',advance='NO') i, t(i_water,k,iday),s(i_water,k, iday),kz(i_water,k,iday),dz(k),z(k)
               enddo
                    write (12,*)
           enddo
       enddo
     close (12)    
    
    end subroutine make_physics_bbl_sed 
!=================================================================================    
	subroutine svan( s,  t,   po,   sigma)
!/*
!c------specific volume anomaly based on 1980 equation of state for
!c      seawater and 1978 practical salinity scale
!c      pressure          PO     decibars
!c      temperature        T     degree celsius (IPTS-68)
!c      salinitty          S     (PSS-78)
!c      spec.vol.anom.  SVAN     1.0E-8 m**3/Kg
!c      density anom.   SIGMA    Kg/m**3
!*/

		real(8)  p,sig,sr,r1,r2,r3,s,t,po,sigma
		real(8)  a,b,c,d,e,a1,b1,aw,bw,k,ko,kw,k35,v350p,sva,dk,gam,pk,dr35p,dvan

	 real(8) r3500, r4 ,dr350
	 data  r3500 /1028.1063/, r4/4.8314E-4/,dr350/28.106331/

	p=po/10.
	sr=sqrt(abs(s))

        r1= ((((6.536332E-9*t-1.120083E-6)*t+1.001685E-4)*t &
	       -9.095290E-3)*t+6.793952E-2)*t-28.263737
        r2= (((5.3875E-9*t-8.2467E-7)*t+7.6438E-5)*t-4.0899E-3)*t &
	      +8.24493E-1
		r3= (-1.6546E-6*t+1.0227E-4)*t-5.72466E-3

	sig=(r4*s + r3*sr + r2)*s +r1

	v350p=1.0/r3500
	sva=-sig*v350p/(r3500+sig)
	sigma= sig + dr350

	if (p.eq.0.0) return

	e = (9.1697E-10*t+2.0816E-8)*t-9.9348E-7
       bw = (5.2787E-8*t-6.12293E-6)*t+3.47718E-5
	b = bw + e*s

	d= 1.91075E-4
	c = (-1.6078E-6*t-1.0981E-5)*t+2.2838E-3
	aw = ((-5.77905E-7*t+1.16092E-4)*t+1.43713E-3)*t-0.1194975
	a = (d*sr + c)*s + aw

	b1 = (-5.3009E-4*t+1.6483E-2)*t+7.944E-2
	a1 = ((-6.1670E-5*t+1.09987E-2)*t-0.603459)*t+54.6746
	kw = (((-5.155288E-5*t+1.360477E-2)*t-2.327105)*t &
	       +148.4206)*t-1930.06
	ko = (b1*sr + a1)*s + kw

	dk = (b*p+a)*p+ko
	k35 = (5.03217E-5*p+3.359406)*p+21582.27
	gam=p/k35
	pk=1.0-gam
	sva = sva * pk + (v350p+sva)*p*dk/(k35*(k35+dk))

	v350p= v350p*pk
	dr35p=gam/v350p
	dvan= sva/(v350p*(v350p+sva))
	sigma = dr350 + dr35p -dvan
    return
    end
     !ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    end module io_ascii