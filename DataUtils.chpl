module DataUtils {
  var data_init_count: uint = 0;

  /*
   * Reset counter for data initialization.
   */
  proc resetDataInitCount() { data_init_count = 0; }

  /*
   * Increment counter for data initialization.
   */
  proc incDataInitCount() { data_init_count += 1; }


  /*
   * Allocate and initialize aligned integer, real, or complex data arrays.
   */
  proc allocAndInitData(type t, len: int /*, VariantID vid */)
  {
    if !(t == int || t == real || t == complex) then
      compilerError('Invalid type ' + t:string);

    var a: [0..<len] t;

    return a;

    //// Should we do this differently for alignment?? If so, change dealloc()
    //ptr = new Int_type[len];
    //initData(ptr, len, vid);
  }

//  /*
//   * Allocate and initialize aligned data arrays.
//   */
//  void allocAndInitData(Real_ptr& ptr, int len, VariantID vid )
//  {
//    ptr =
//      RAJA::allocate_aligned_type<Real_type>(RAJA::DATA_ALIGN,
//          len*sizeof(Real_type));
//    initData(ptr, len, vid);
//  }
//
//  void allocAndInitDataConst(Real_ptr& ptr, int len, Real_type val,
//      VariantID vid)
//  {
//    (void) vid;
//
//    ptr =
//      RAJA::allocate_aligned_type<Real_type>(RAJA::DATA_ALIGN,
//          len*sizeof(Real_type));
//    initDataConst(ptr, len, val, vid);
//  }
//
//  void allocAndInitDataRandSign(Real_ptr& ptr, int len, VariantID vid)
//  {
//    ptr =
//      RAJA::allocate_aligned_type<Real_type>(RAJA::DATA_ALIGN,
//          len*sizeof(Real_type));
//    initDataRandSign(ptr, len, vid);
//  }
//
//  void allocAndInitDataRandValue(Real_ptr& ptr, int len, VariantID vid)
//  {
//    ptr =
//      RAJA::allocate_aligned_type<Real_type>(RAJA::DATA_ALIGN,
//          len*sizeof(Real_type));
//    initDataRandValue(ptr, len, vid);
//  }
//
//  void allocAndInitData(Complex_ptr& ptr, int len, VariantID vid)
//  {
//    // Should we do this differently for alignment?? If so, change dealloc()
//    ptr = new Complex_type[len];
//    initData(ptr, len, vid);
//  }
}
