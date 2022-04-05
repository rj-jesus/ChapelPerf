module DataUtils {
  private use DataTypes;
  private use Enums;
  private use Utils;

  // Note: Rectangular domain indices in Chapel are ordered according to the
  // lexicographic order of their values, i.e. the index with the highest rank
  // is listed first and changes most slowly (as in row-major ordering)
  // https://chapel-lang.org/docs/language/spec/domains.html#rectangular-domain-values

  /*
   * Chapel helper functions.
   */
  pragma "fn returns aliasing array"
  inline proc reindex(ref arr, newRange) return arr.reindex(newRange#arr.size);

  inline proc setv(ref args...?n, rhs) { for arg in args do arg = rhs; }

  inline proc vcopy(ref dst, const ref src) {
    assert(src.size == dst.size);
    for (_, s, d) in zip(0.., src, dst) do d = s;
  }

  proc makeDomain(shape: int ...?rank) {
    var dims: rank*range;
    for param i in 0..<rank do dims[i] = 0..<shape[i];
    return {(...dims)};
  }

  proc makeArrayFromArray(const ref A, shape) {
    var dom = makeDomain((...shape));
    var B: [dom] A.eltType;

    vcopy(B, A);

    return B;
  }

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
  proc allocAndInitData(type t: numeric, lens: int ...?n, vid:VariantID)
    return allocAndInitData(t, makeDomain((...lens)), vid);

  proc allocAndInitData(type t: numeric, d: domain, vid:VariantID) {
    var a: [d] t;
    initData(a, vid);
    return a;
  }

  proc allocAndInitDataConst(type t: numeric, lens: int ...?n, val, vid:VariantID)
    return allocAndInitDataConst(t, makeDomain((...lens)), val, vid);

  proc allocAndInitDataConst(type t: numeric, d: domain, val, vid:VariantID) {
    var a: [d] t;
    initDataConst(a, val:t, vid);
    return a;
  }

  proc allocAndInitDataRandSign(type t: numeric, len: int, vid: VariantID) {
    var a: [0..<len] t;
    initDataRandSign(a, len, vid);
    return a;
  }

  proc allocAndInitDataRandValue(type t: numeric, len: int, vid: VariantID) {
    var a: [0..<len] Real_type;
    initDataRandValue(a, len, vid);
    return a;
  }

  /*
   * Initialize scalar data.
   */
  proc initData(type t: Real_type, vid: VariantID): t {
    const factor: t = if data_init_count % 2 then 0.1 else 0.2;
    incDataInitCount();
    return (factor*1.1/1.12345):t;
  }

  /*
   * Initialize int array to randomly signed positive and negative values.
   */
  proc initData(A: [] Int_type, vid: VariantID) where A.isRectangular() {
    srand(4793);

    inline proc signfact return rand():Real_type/RAND_MAX;

    for a in A do
      a = (if signfact < 0.5 then -1 else 1):Int_type;

    A[(A.size*signfact):int] = -58;
    A[(A.size*signfact):int] =  19;

    incDataInitCount();
  }

  /*
   * Initialize real array to non-random positive values (0.0, 1.0) based on
     their array position (index) and the order in which this method is called.
   */
  proc initData(A: [] Real_type, vid: VariantID) where A.isRectangular() {
    const factor = (if data_init_count % 2 then 0.1 else 0.2):Real_type;

    for (a,i) in zip(A, 0..<A.size) do
      a = (factor*(i + 1.1)/(i + 1.12345)):Real_type;

    incDataInitCount();
  }

  /*
   * Initialize complex array.
   */
  proc initData(A: [] Complex_type, vid: VariantID) where A.isRectangular() {
    const factor = (if data_init_count % 2
                    then 0.1+0.2i
                    else 0.2+0.3i):Complex_type;

    for (a,i) in zip(A, 0..<A.size) do
      a = (factor*(i + 1.1)/(i + 1.12345)):Complex_type;

    incDataInitCount();
  }

  /*
   * Initialize array to constant values.
   */
  proc initDataConst(A: [] ?t, val: t, vid: VariantID) {
    A = val;
    incDataInitCount();
  }

  /*
   * Initialize real array with random sign.
   */
  proc initDataRandSign(A: [] Real_type, len: int=A.size, vid: VariantID) where A.isRectangular() {
    srand(4793);

    const factor = if data_init_count % 2 then 0.1 else 0.2;
    inline proc signfact return if rand():Real_type/RAND_MAX < 0.5 then -1.0 else 1.0;

    for (a, i) in zip(A, 0..<len) do
      a = (signfact*factor*(i + 1.1)/(i + 1.12345)):Real_type;

    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random values.
   */
  proc initDataRandValue(A: [] Real_type, len: int=A.size, vid: VariantID) where A.isRectangular() {
    srand(4793);

    for (a, i) in zip(A, 0..<len) do
      a = rand():Real_type/RAND_MAX;

    incDataInitCount();
  }

  proc calcChecksum(const A: [] Real_type, len: int=A.size, scale_factor: Real_type=1.0): Checksum_type where A.isRectangular() {
    var s: Checksum_type = 0;
    for (a, i) in zip(A, 1..len) do s += a*i*scale_factor;
    return s;
  }

  proc calcChecksum(const A: [] Complex_type, len: int=A.size, scale_factor: Real_type=1.0): Checksum_type where A.isRectangular() {
    var s: Checksum_type = 0;
    for (a, i) in zip(A, 1..len) do s += (a.re+a.im)*i*scale_factor;
    return s;
  }
}
