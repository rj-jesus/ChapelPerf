module DataTypes {
  private use LongDouble;

  type Checksum_type = longdouble;
  type Complex_type = complex;
  type Elapsed_type = real;
  type Index_type = int;
  type Real_type = real;

  //
  // This is a vector wrapper that uses 0-based indexing (only primitive types
  // supported).
  //
  class vector {
    type eltType;
    var A: list(eltType);

    //
    // This vector supports 0-based indexing.
    //
    proc this(i: int) ref {
      return A[i];
    }
    proc push_back(e: eltType) {
      A.append(e);
    }
    proc size {
      return A.size;
    }
    iter these() {
      for a in A do yield a;
    }
  }
}
