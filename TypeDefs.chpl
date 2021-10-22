module TypeDefs {
  use LongDouble;

  type Index_type = int;
  type Real_type = real;
  type Checksum_type = real;
  //type Checksum_type = longdouble;
  type Elapsed_type = real;

  //
  // This is a vector wrapper that uses 0-based indexing.
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
