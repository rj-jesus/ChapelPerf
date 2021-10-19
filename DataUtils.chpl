module DataUtils {
  use Random;

  var data_init_count: uint = 0;

  type VariantID = nothing;  // temporary...

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
  proc allocAndInitData(type t, len: int, vid: VariantID = none)
  {
    if !(t == int || t == real || t == complex) then
      compilerError('Invalid type ' + t:string);

    var a: [0..<len] t;
    initData(a, len, vid);
    return a;
  }

  proc allocAndInitDataConst(len: int, val: real, vid: VariantID)
  {
    var a: [0..<len] real;
    initDataConst(a, len, val, vid);
    return a;
  }

  proc allocAndInitDataRandSign(len: int, vid: VariantID)
  {
    var a: [0..<len] real;
    initDataRandSign(a, len, vid);
    return a;
  }

  proc allocAndInitDataRandValue(len: int, vid: VariantID)
  {
    var a: [0..<len] real;
    initDataRandValue(a, len, vid);
    return a;
  }

  /*
   * \brief Initialize int array to
   * randomly signed positive and negative values.
   */
  proc initData(a: [] int, len: int, vid: VariantID)
  {
//    // First touch...
//#if defined(RAJA_ENABLE_OPENMP) && defined(RUN_OPENMP)
//    if (vid == Base_OpenMP ||
//        vid == Lambda_OpenMP ||
//        vid == RAJA_OpenMP ) {
//#pragma omp parallel for
//      for (int i = 0; i < len; ++i) {
//        ptr[i] = 0;
//      };
//    }
//#endif

    var randStream = new RandomStream(real, 4793);

    for (i,r) in zip(a.domain, randStream) do
      a[i] = (r < 0.5 ? -1 : 1);

    var ilo = (len * randStream.getNext()):int;
    a[ilo] = -58;

    var ihi = (len * randStream.getNext()):int;
    ptr[ihi] = 19;

    incDataInitCount();
  }

  /*
   * \brief Initialize real array to non-random
   * positive values (0.0, 1.0) based on their array position
   * (index) and the order in which this method is called.
   */
  proc initData(a: [] real, len: int, vid: VariantID)
  {
    var factor = if data_init_count % 2 then 0.1 else 0.2;

//    // first touch...
//#if defined(RAJA_ENABLE_OPENMP) && defined(RUN_OPENMP)
//    if ( vid == Base_OpenMP ||
//        vid == Lambda_OpenMP ||
//        vid == RAJA_OpenMP ) {
//#pragma omp parallel for
//      for (int i = 0; i < len; ++i) {
//        ptr[i] = factor*(i + 1.1)/(i + 1.12345);
//      };
//    }
//#endif

    for i in a.domain do
      a[i] = factor*(i + 1.1)/(i + 1.12345);

    incDataInitCount();
  }

//  /*
//   * Initialize Real_type data array to constant values.
//   */
//  void initDataConst(Real_ptr& ptr, int len, Real_type val,
//      VariantID vid)
//  {
//
//    // first touch...
//#if defined(RAJA_ENABLE_OPENMP) && defined(RUN_OPENMP)
//    if ( vid == Base_OpenMP ||
//        vid == Lambda_OpenMP ||
//        vid == RAJA_OpenMP ) {
//#pragma omp parallel for
//      for (int i = 0; i < len; ++i) {
//        ptr[i] = 0;
//      };
//    }
//#else
//    (void) vid;
//#endif
//
//    for (int i = 0; i < len; ++i) {
//      ptr[i] = val;
//    };
//
//    incDataInitCount();
//  }
//
//  /*
//   * Initialize Real_type data array with random sign.
//   */
//  void initDataRandSign(Real_ptr& ptr, int len, VariantID vid)
//  {
//    (void) vid;
//
//    // First touch...
//#if defined(RAJA_ENABLE_OPENMP) && defined(RUN_OPENMP)
//    if ( vid == Base_OpenMP ||
//        vid == Lambda_OpenMP ||
//        vid == RAJA_OpenMP ) {
//#pragma omp parallel for
//      for (int i = 0; i < len; ++i) {
//        ptr[i] = 0.0;
//      };
//    }
//#endif
//
//    Real_type factor = ( data_init_count % 2 ? 0.1 : 0.2 );
//
//    srand(4793);
//
//    for (int i = 0; i < len; ++i) {
//      Real_type signfact = Real_type(rand())/RAND_MAX;
//      signfact = ( signfact < 0.5 ? -1.0 : 1.0 );
//      ptr[i] = signfact*factor*(i + 1.1)/(i + 1.12345);
//    };
//
//    incDataInitCount();
//  }
//
//  /*
//   * Initialize Real_type data array with random values.
//   */
//  void initDataRandValue(Real_ptr& ptr, int len, VariantID vid)
//  {
//    (void) vid;
//
//    // First touch...
//#if defined(RAJA_ENABLE_OPENMP) && defined(RUN_OPENMP)
//    if ( vid == Base_OpenMP ||
//        vid == Lambda_OpenMP ||
//        vid == RAJA_OpenMP ) {
//#pragma omp parallel for
//      for (int i = 0; i < len; ++i) {
//        ptr[i] = 0.0;
//      };
//    }
//#endif
//
//    srand(4793);
//
//    for (int i = 0; i < len; ++i) {
//      ptr[i] = Real_type(rand())/RAND_MAX;
//    };
//
//    incDataInitCount();
//  }
//
//  /*
//   * Initialize Complex_type data array.
//   */
//  void initData(Complex_ptr& ptr, int len, VariantID vid)
//  {
//    (void) vid;
//
//    Complex_type factor = ( data_init_count % 2 ?  Complex_type(0.1,0.2) :
//        Complex_type(0.2,0.3) );
//
//#if defined(RAJA_ENABLE_OPENMP) && defined(RUN_OPENMP)
//    if ( vid == Base_OpenMP ||
//        vid == Lambda_OpenMP ||
//        vid == RAJA_OpenMP ) {
//#pragma omp parallel for
//      for (int i = 0; i < len; ++i) {
//        ptr[i] = factor*(i + 1.1)/(i + 1.12345);
//      };
//    }
//#endif
//
//    for (int i = 0; i < len; ++i) {
//      ptr[i] = factor*(i + 1.1)/(i + 1.12345);
//    }
//
//    incDataInitCount();
//  }
//
//  /*
//   * Initialize scalar data.
//   */
//  void initData(Real_type& d, VariantID vid)
//  {
//    (void) vid;
//
//    Real_type factor = ( data_init_count % 2 ? 0.1 : 0.2 );
//    d = factor*1.1/1.12345;
//
//    incDataInitCount();
//  }
}
