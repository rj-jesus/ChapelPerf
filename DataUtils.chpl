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
   * Allocate and initialize integer, real, or complex data arrays.
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

  proc allocAndInitDataConst(len: int, val: real /*, VariantID vid */)
  {
    var a: [0..<len] real;

    return a;
    initDataConst(ptr, len, val, vid);
  }

  proc allocAndInitDataRandSign(len: int /*, VariantID vid */)
  {
    var a: [0..<len] real;
    initDataRandSign(ptr, len, vid);
  }

  proc allocAndInitDataRandValue(len: int /*, VariantID vid */)
  {
    var a: [0..<len] real;
    initDataRandValue(ptr, len, vid);
  }
}
