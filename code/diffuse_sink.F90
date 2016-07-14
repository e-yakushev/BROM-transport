    module diffuse_sink
    
    implicit none
    private
    public calculate_phys, calculate_sed
    
   
    contains
    
    
    subroutine calculate_phys(i_water, k_max, par_max,  cc, fick, dcc, use_bound_up, use_bound_low, bound_up, bound_low, &
                            surf_flux, boundary_bbl_sediments, dz, kz, kz_bio, kz_tot, julianday,  i_O2,  dt, freq_az)
    
    !subroutine calculate_phys(cc, fick, dcc, k_max, par_max, use_bound_up, use_bound_low, bound_up, bound_low, &
    !                        surf_flux, boundary_bbl_sediments, kz_tot, kz, julianday, kz_bio, i_O2, dz, dt, freq_az)
    
    implicit none
    real(8), dimension(:, :, :), intent(inout)          :: cc, fick, dcc
    integer, intent(in)                                 :: i_water
    integer, intent(in)                                 :: k_max
    integer, intent(in)                                 :: par_max
    logical, dimension(:,:), intent(in)                   :: use_bound_up, use_bound_low
    real(8), dimension(:,:), intent(in)                   :: surf_flux, bound_up, bound_low
    integer, intent(in)                                 :: boundary_bbl_sediments
    real(8), dimension(:), intent(inout)                :: kz_tot
    real(8), intent(in)                                 :: kz(:, :, :)
    integer, intent(in)                                 :: julianday
    real(8), intent(in)                                 :: kz_bio(:)
    integer, intent(in)                                 :: i_O2
    real(8), intent(in)                                 :: dz(:)
    real(8), intent(in)                                 :: dt
    integer, intent(in)                                 :: freq_az
    integer                                             :: i, k, ip 
    
    i=2
    
    do k = 2, (k_max - 1)
        !boundary conditions
        !upper boundary
        if (k .eq. 2) then
            do ip = 1, par_max
                if (use_bound_up(i,ip)) then
                    cc(i,k - 1, ip) = bound_up(i,ip)
                else
                    cc(i,k - 1, ip) = cc(i,k, ip)
                end if
            enddo
        end if
        !low Boundary
        if (k .eq. (k_max - 1)) then
            do ip = 1, par_max
                if (use_bound_low(i,ip)) then
                    cc(i,k + 1, ip) = bound_low(i,ip)
                else
                    cc(i,k + 1, ip) = cc(i,k, ip)
                end if
            end do
        end if
        !turbulence and advection
        do ip = 1, par_max
            !here we calculate turbulent exchange (difference of turb.fluxes (Fick) for all the layers)
            !if (k > boundary_bbl_sediments - 1) then
            !    kz_tot(k) = (kz(k,  julianday) + kz_bio(k)                                             &
            !        * cc(boundary_bbl_sediments - 1, i_O2) / (cc(boundary_bbl_sediments - 1, i_O2) + 5.))!for plots
            !else
            !    kz_tot(k) = 0.
            !end if 
            !fluxes of parameters, bioturbation (kz_bio) depend on Zoo or O2 availability at the surface  
            fick(i,k, ip) = (kz(i,k, julianday) + kz_bio(k)                                                &
                * cc(i,boundary_bbl_sediments - 1, i_O2) / (cc(i,boundary_bbl_sediments - 1, i_O2) + 5.))   &
                * (cc(i,k + 1, ip) - cc(i,k, ip)) / dz(k)
            fick(i,k - 1, ip) = (kz(i,k - 1, julianday) + kz_bio(k - 1)                                    &
                * cc(i,boundary_bbl_sediments - 1, i_O2) / (cc(i,boundary_bbl_sediments - 1, i_O2) + 5.))   &
                * (cc(i,k, ip) - cc(i,k - 1, ip)) / dz(k - 1) 
          
            dcc(i,k, ip) = 86400. * (fick(i,k, ip) - fick(i,k - 1, ip)) / ((dz(k) + dz(k - 1)) / 2.) 
        end do
        !fluxes of O2 and nutrients from the air/rivers:
        if (k .eq. 2) then 
            do ip = 1,par_max
                dcc(i,k, ip) = dcc(i,k, ip) + surf_flux(i,ip) / dz(k)
            end do
        end if
        !fluxes from the deep sediments (i.e. CO2 leaks):
    end do
    !time integration
    do k = 2, (k_max - 1)
        do ip = 1, par_max
            cc(i,k, ip) = cc(i,k, ip) + (dt / freq_az) * dcc(i,k, ip)
        end do
    end do
    
    end subroutine calculate_phys
                            
!    subroutine calculate_sed(par_max, k_max, kz, julianday, wbio, cc, dcc, dz, dt, freq_sed)
    subroutine calculate_sed(i_water, k_max, par_max, cc, dcc, dz, kz, wbio, julianday, dt, freq_sed)      
    
    implicit none
    integer, intent(in)                                 :: par_max
    integer, intent(in)                                 :: i_water
    integer, intent(in)                                 :: k_max
    real(8), intent(in)                                 :: kz(:, :, :)
    integer, intent(in)                                 :: julianday
    real(8), dimension(:, :, :), intent(in)                :: wbio
    real(8), dimension(:, :, :), intent(inout)             :: cc, dcc
    real(8), intent(in)                                 :: dz(:)
    real(8), intent(in)                                 :: dt
    integer, intent(in)                                 :: freq_sed
    
    real(8)                                             :: w_u(par_max)
    real(8)                                             :: w_d(par_max)
    real(8), parameter                                  :: w_buruing = 15.e-6 ! rate of burying 5*0.1 cm/y
    integer                                             :: i, k, ip
    
    i=2
    
    do  k = 2, k_max - 1
        !boundary conditions
        !from "up": we check either we are in water or sediment
        if (kz(i, k, julianday) > 5E-9) then
            !in FABM: negative w means sinking!
            w_u(:) = -wbio(i,k, :)
        else
            w_u(:) = w_buruing
        end if
        !to "down"           
        if (kz(i,k + 1, julianday)  > 1.E-9 ) then
            w_d(:) = -wbio(i,k, :)
        else
            w_d(:) = w_buruing
        end if
        !upper boundary 
        if (k < 3) then
            w_u = 0.0
        end if   
        !f surface ----------
        if  (k .eq. 2) then
            do ip = 1, par_max
                cc(i,k - 1, ip)=0.
            end do
        end if        
        dcc(i,k, :) = 0.5 * (w_u * cc(i,k - 1, :) / ((dz(k - 1) + dz(k)) /2.) &
        -w_d * cc(i,k, :)  / ((dz(k - 1) + dz(k)) /2.))
    end do
    !time integration
    do  k = 2, (k_max - 1)
        cc(i,k, :) = cc(i,k, :) + (dt / freq_sed) * dcc(i,k, :)
    end do

    end subroutine calculate_sed

    end module diffuse_sink