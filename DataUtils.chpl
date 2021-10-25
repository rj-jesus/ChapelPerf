module DataUtils {
  private use Random;

  private use DataTypes;
  private use KernelBase;

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
  proc allocAndInitData(type t: numeric, len: int, vid: VariantID = VariantID.NONE)
  {
    return allocAndInitData(t, {0..<len});
  }

  proc allocAndInitData(type t: numeric, d: domain, vid: VariantID = VariantID.NONE)
  {
    var a: [d] t;
    initData(a, vid);
    return a;
  }

  proc allocAndInitDataConst(type t: numeric, len: int, val, vid: VariantID = VariantID.NONE) where isCoercible(val.type, t)
  {
    return allocAndInitDataConst(t, {0..<len}, val);
  }

  proc allocAndInitDataConst(type t: numeric, d: domain, val, vid: VariantID = VariantID.NONE) where isCoercible(val.type, t)
  {
    var a: [d] t;
    initDataConst(a, val:t, vid);
    return a;
  }

  proc allocAndInitDataRandSign(len: int, vid: VariantID = VariantID.NONE)
  {
    var a: [0..<len] real;
    initDataRandSign(a, len, vid);
    return a;
  }

  proc allocAndInitDataRandValue(len: int, vid: VariantID = VariantID.NONE)
  {
    var a: [0..<len] real;
    initDataRandValue(a, len, vid);
    return a;
  }

  /*
   * \brief Initialize int array to randomly signed positive and negative
   * values.
   */
  proc initData(A: [] int, len: int, vid: VariantID = VariantID.NONE)
  {
    var randStream = new RandomStream(real, 4793);

    for (i,r) in zip(0..<len, randStream) do
      A[i] = (r < 0.5 ? -1 : 1);

    A[(len * randStream.getNext()):int] = -58;
    A[(len * randStream.getNext()):int] =  19;

    incDataInitCount();
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
  proc initData(vid: VariantID = VariantID.NONE): real
  {
    const factor = if data_init_count % 2 then 0.1 else 0.2;
    incDataInitCount();
    return factor*1.1/1.12345;
  }

  /*
   * \brief Initialize real array to non-random positive values (0.0, 1.0)
   * based on their array position (index) and the order in which this method
   * is called.
   */
  proc initData(A: [?d] real, vid: VariantID = VariantID.NONE) where isRectangularDom(d)
  {
    const factor = if data_init_count % 2 then 0.1 else 0.2;
    for (a,i) in zip(A, 0..<A.size) do a = factor*(i + 1.1)/(i + 1.12345);
    incDataInitCount();
  }

  /*
   * \brief Initialize complex array.
   */
  proc initData(A: [?d] complex, vid: VariantID = VariantID.NONE) where isRectangularDom(d)
  {
    const factor = if data_init_count % 2 then 0.1+0.2i else 0.2+0.3i;
    for (a,i) in zip(A, 0..<A.size) do a = factor*(i + 1.1)/(i + 1.12345);
    incDataInitCount();
  }

  /*
   * \brief Initialize array to constant values.
   */
  proc initDataConst(A: [] ?t, val, vid: VariantID = VariantID.NONE) where isCoercible(val.type, t)
  {
    A = val;
    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random sign.
   */
  proc initDataRandSign(A: [] real, vid: VariantID = VariantID.NONE) where isRectangularArr(A)
  {
    const factor = if data_init_count % 2 then 0.1 else 0.2;
    var randStream = new RandomStream(real, 4793);

    for (a,i,r) in zip(A, 0..<A.size, randStream) {
      const signfact = if r < 0.5 then -1.0 else 1.0;
      a = signfact*factor*(i + 1.1)/(i + 1.12345);
    }

    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random values.
   */
  proc initDataRandValue(A: [], vid: VariantID = VariantID.NONE) where isRectangularArr(A)
  {
    fillRandom(A, seed=4793);
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
