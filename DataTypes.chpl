module DataTypes {
  private use SysCTypes;

  private use LongDouble;

  config param USE_DOUBLE  = true;
  config param USE_FLOAT   = false;
  config param USE_COMPLEX = true;

  type Index_type = int;
  type Int_type = c_int;
  type Real_type = if USE_DOUBLE then real(64)
                   else if USE_FLOAT then real(32)
                   else compilerError("Real_type is undefined!");
  type Complex_type = if USE_COMPLEX then complex(2*numBits(Real_type))
                      else nothing;
  type Elapsed_type = real;
  type Checksum_type = longdouble;
}
