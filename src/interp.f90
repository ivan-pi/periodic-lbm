module interp

   use precision, only: wp
   implicit none

contains 


   !   u1 ------- u2
   !   |          |
   !   |          |
   !   |          |
   !   |          |
   !   u3 ---o--- u4
   !   |    /| ry |
   !   |   @--    |
   !   |    rx    |
   !   |          |
   !   u5 ------- u6

   function interp_vbox6(u1,u2,u3,u4,u5,u6,rx,ry) result(uinterp)
      real(wp), intent(in) :: u1, u2, u3, u4, u5, u6
      real(wp), intent(in) :: rx, ry

      real(wp) :: uinterp

      real(wp) :: tx
      real(wp) :: v1, v2, v3
      real(wp) :: h1, h2, h3

      tx = 0.5_wp - rx

      ! interpolate linearly in first direction
      v1 = u1 + (u2 - u1)*tx
      v2 = u3 + (u4 - u3)*tx
      v3 = u5 + (u6 - u5)*tx

      ! interpolate quadratically along second direction
      
      h1 = 0.5_wp*ry*(ry + 1)
      h2 = (ry + 1)*(ry - 1)
      h3 = 0.5_wp*ry*(ry - 1)

      uinterp = h1*v1 + h2*v2 + h3*v3

   end function


   !   u5 ------- u3 ------- u1
   !   |          |          |
   !   |          |          |
   !   |          o          |
   !   |        / |          |
   !   |       @--|          |
   !   u6 ------- u4 ------- u2

   function interp_hbox6(u1,u2,u3,u4,u5,u6,rx,ry) result(uinterp)
      real(wp), intent(in) :: u1, u2, u3, u4, u5, u6
      real(wp), intent(in) :: rx, ry

      real(wp) :: uinterp

      real(wp) :: ty
      real(wp) :: v1, v2, v3
      real(wp) :: h1, h2, h3

      ty = 0.5_wp - ry

      ! interpolate linearly in first direction
      v1 = u2 + (u1 - u2)*ty
      v2 = u4 + (u3 - u4)*ty
      v3 = u6 + (u5 - u6)*ty

      ! interpolate quadratically along second direction
      h1 = 0.5_wp*rx*(rx + 1)
      h2 = (rx + 1)*(rx - 1)
      h3 = 0.5_wp*rx*(rx - 1)

      uinterp = h1*v1 + h2*v2 + h3*v3

   end function

end module
