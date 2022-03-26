module precision

   implicit none
   private

   public :: wp

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)

   integer, parameter :: wp = dp

end module