    module smooth
    
    implicit none
    private
    public tem_sal_kz_smoothing, cc_smoothing
    
    contains
    
    subroutine tem_sal_kz_smoothing(tem2, sal2, Kz2)
    
    implicit none
    real(8), dimension(:, :), intent(inout) :: tem2, sal2, Kz2
    integer                                 :: ijk, iday, k, inext, ilast
    
    do ijk=1,5
        !with time    
        do iday=1,365
            do k=1,74
                if (iday.ge.365) then
                    inext=1
                else
                    inext=iday+1
                end if     

                if (iday.le.1) then
                    ilast=365
                else
                    ilast=iday-1
                end if       
                tem2(k,iday) = (tem2(k,iday)+tem2(k,ilast)+tem2(k,inext))/3.
                sal2(k,iday) = (sal2(k,iday)+sal2(k,ilast)+sal2(k,inext))/3.
                Kz2(k,iday)  = (Kz2(k,iday)+Kz2(k,ilast)+Kz2(k,inext))/3.
            end do    
        end do  
        !!with depth
        !do iday=1,365
        !    do k=2,18
        !        tem2(k,iday) = (tem2(k,iday)+tem2(k,ilast)+tem2(k,inext))/3.
        !        sal2(k,iday) = (sal2(k,iday)+sal2(k,ilast)+sal2(k,inext))/3.
        !        Kz2(k,iday)  = (Kz2(k,iday)+Kz2(k-1,iday)+Kz2(k+1,iday))/3.
        !    enddo    
        !enddo   
 
    end do
    
    end subroutine tem_sal_kz_smoothing
    
    subroutine cc_smoothing(Cc, k_max, par_max)
    
    implicit none
    real(8), dimension(:, :), intent(inout) :: Cc
    integer, intent(in)                     :: k_max, par_max
    integer                                 :: k, ip
    
    do k=2, k_max-1
        do ip=1,par_max
            Cc(k,ip)=(2.*Cc(k,ip)+Cc(k-1,ip)+Cc(k+1,ip))/4.
        end do
    end do
    
    end subroutine cc_smoothing
    
    end module smooth