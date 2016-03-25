      use iso_fortran_env, only: stdout=>output_unit, stderr=>error_unit
      use, intrinsic :: iso_c_binding, only: sp=>C_FLOAT, dp=>C_DOUBLE,
     &              cp=>C_FLOAT_COMPLEX, zp=>C_DOUBLE_COMPLEX,
     &                  i64=>C_LONG_LONG, sizeof=>c_sizeof
      use, intrinsic :: ieee_arithmetic, only: qnan=>ieee_quiet_nan,
     &  ieee_value
