module precision

   implicit none
   private

   public :: wp

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)

#if PRECISION_SINGLE
   integer, parameter :: wp = sp
#else
   integer, parameter :: wp = dp
#endif
   
end module