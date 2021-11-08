module DataUtils {
  private use DataTypes;
  private use Enums;
  private use Utils;

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
  proc allocAndInitData(type t: numeric, len: int, vid:VariantID) {
    return allocAndInitData(t, {0..<len}, vid);
  }

  proc allocAndInitData(type t: numeric, d: domain, vid:VariantID) {
    var a: [d] t;
    initData(a, vid);
    return a;
  }

  proc allocAndInitDataConst(type t: numeric, len: int, val, vid:VariantID) {
    return allocAndInitDataConst(t, {0..<len}, val, vid);
  }

  proc allocAndInitDataConst(type t: numeric, d: domain, val, vid:VariantID) {
    var a: [d] t;
    initDataConst(a, val:t, vid);
    return a;
  }

  proc allocAndInitDataRandSign(len: int, vid: VariantID) {
    var a: [0..<len] real;
    initDataRandSign(a, len, vid);
    return a;
  }

  proc allocAndInitDataRandValue(type t: numeric, len: int, vid: VariantID) {
    var a: [0..<len] Real_type;
    initDataRandValue(a, len, vid);
    return a;
  }

  // Note: Rectangular domain indices are ordered according to the
  // lexicographic order of their values, i.e. the index with the highest rank
  // is listed first and changes most slowly (as in row-major ordering)
  // https://chapel-lang.org/docs/language/spec/domains.html#rectangular-domain-values

  proc isCLike(d: domain) where isRectangularDom(d) && isIntegralType(d.idxType) && d.rank == 1 {
    return d.low == 0 && d.stride == 1;
  }

  proc isCLike(d: domain) where isRectangularDom(d) && isIntegralType(d.idxType) && d.rank  > 1 {
    return && reduce(for i in 0..<d.rank do d.low(i) == 0 && d.stride(i) == 1);
  }

  /*
   * \brief Initialize scalar data.
   */
  proc initData(type t: numeric = Real_type, vid: VariantID) {
    const factor = if data_init_count % 2 then 0.1 else 0.2;
    incDataInitCount();
    return (factor*1.1/1.12345):t;
  }

  /*
   * \brief Initialize int array to randomly signed positive and negative
   * values.
   */
  proc initData(A: [?d] Int_type, vid: VariantID) where isRectangularDom(d) && A.rank == 1 {
    srand(4793);

    proc signfact return rand():Real_type/RAND_MAX;

    for a in A do
      a = (if signfact < 0.5 then -1 else 1):Int_type;

    A[(A.size * signfact):int] = -58;
    A[(A.size * signfact):int] =  19;

    incDataInitCount();
  }

  /*
   * \brief Initialize real array to non-random positive values (0.0, 1.0)
   * based on their array position (index) and the order in which this method
   * is called.
   */
  proc initData(A: [?d] Real_type, vid: VariantID) where isRectangularDom(d) {
    const factor = (if data_init_count % 2 then 0.1 else 0.2):Real_type;

    for (a,i) in zip(A, 0..<A.size) do
      a = (factor*(i + 1.1)/(i + 1.12345)):Real_type;

    incDataInitCount();
  }

  /*
   * \brief Initialize complex array.
   */
  proc initData(A: [?d] Complex_type, vid: VariantID) where isRectangularDom(d) {
    const factor = (if data_init_count % 2 then 0.1+0.2i
                                           else 0.2+0.3i):Complex_type;

    for (a,i) in zip(A, 0..<A.size) do
      a = (factor*(i + 1.1)/(i + 1.12345)):Complex_type;

    incDataInitCount();
  }

  /*
   * \brief Initialize array to constant values.
   */
  proc initDataConst(A: [] ?t, val: t, vid: VariantID) {
    A = val;
    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random sign.
   */
  proc initDataRandSign(A: [] Real_type, vid: VariantID) where isRectangularArr(A) {
    srand(4793);

    const factor = if data_init_count % 2 then 0.1 else 0.2;
    proc signfact return if rand():Real_type/RAND_MAX < 0.5 then -1.0 else 1.0;

    for (a,i) in zip(A, 0..<A.size) do
      a = (signfact*factor*(i + 1.1)/(i + 1.12345)):Real_type;

    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random values.
   */
  proc initDataRandValue(A: [?d] Real_type, len: int=A.size, vid: VariantID) where isRectangularArr(A) {
    srand(4793);

    for (a, i) in zip(A, 0..<len) do
      a = rand():Real_type/RAND_MAX;

    incDataInitCount();
  }

  proc calcChecksum(const A: [] Real_type, len:int=A.size, scale_factor:Real_type=1.0): Checksum_type where isRectangularArr(A) {
    var s:Checksum_type = 0;
    for (a,i) in zip(A, 1..len) do s += a*i*scale_factor;
    return s;
  }

  proc calcChecksum(const A: [] Complex_type, len:int=A.size, scale_factor:Real_type=1.0): Checksum_type where isRectangularArr(A) {
    var s:Checksum_type = 0;
    for (a,i) in zip(A, 1..len) do s += (a.re+a.im)*i*scale_factor;
    return s;
  }
}
