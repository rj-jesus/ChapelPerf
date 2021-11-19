module DataUtils {
  private use DataTypes;
  private use Enums;
  private use Utils;

  // Note: Rectangular domain indices are ordered according to the
  // lexicographic order of their values, i.e. the index with the highest rank
  // is listed first and changes most slowly (as in row-major ordering)
  // https://chapel-lang.org/docs/language/spec/domains.html#rectangular-domain-values

  var data_init_count: uint = 0;

  pragma "fn returns aliasing array"
  inline proc reindex(ref arr, newRange) return arr.reindex(newRange#arr.size);

  inline proc setv(ref args...?n, rhs) { for arg in args do arg = rhs; }

  /*
   * Reset counter for data initialization.
   */
  proc resetDataInitCount() { data_init_count = 0; }

  /*
   * Increment counter for data initialization.
   */
  proc incDataInitCount() { data_init_count += 1; }

  proc makeDomain(shape: int ...?rank) {
    var dims: rank*range;
    for param i in 0..<rank do dims[i] = 0..<shape[i];
    return {(...dims)};
  }

  //
  // A simple array view, see
  // https://matrix.to/#/!PBYDSerrfYujeStENM:gitter.im/$gjozTZJcBgZbThV9d7UPIGTOcz7JUMyxItjxXo1Id34
  // and previous messages
  //
  class SimpleArrayView {
    param rank;
    const shape;
    const arr;
    const dom;

    proc init(arr, shape) where isRectangularArr(arr) && arr.rank == 1 {
      const size = * reduce shape;

      if arr.size < size then
        halt("number of elements in the view must not " +
             "exceed the size of the underlying array");

      this.rank = shape.size;
      this.shape = shape;
      this.arr = arr._value;
      this.dom = makeDomain((...shape));
    }

    pragma "reference to const when const this"
    inline proc ref this(args: int ...rank) ref {
      var idx = 0;
      for param i in 0..#rank-1 do idx = (idx+args[i])*shape[i+1];
      idx += args[rank-1];

      return arr.dsiAccess(idx);
    }

    override proc writeThis(f) throws {
      //var tmpArr: [dom] arr.eltType;
      //for (i, x) in zip(0..<dom.size, tmpArr) do x = arr.dsiAccess(i);
      //f.write(tmpArr);

      f <~> shape <~> ": ";
      arr.dsiSerialWrite(f);
    }
  }

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
  proc initData(type t: numeric = Real_type, vid: VariantID) {
    const factor = if data_init_count % 2 then 0.1 else 0.2;
    incDataInitCount();
    return (factor*1.1/1.12345):t;
  }

  /*
   * Initialize int array to randomly signed positive and negative values.
   */
  proc initData(A: [] Int_type, vid: VariantID) where isRectangularArr(A) {
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
  proc initData(A: [] Real_type, vid: VariantID) where isRectangularArr(A) {
    const factor = (if data_init_count % 2 then 0.1 else 0.2):Real_type;

    for (a,i) in zip(A, 0..<A.size) do
      a = (factor*(i + 1.1)/(i + 1.12345)):Real_type;

    incDataInitCount();
  }

  /*
   * Initialize complex array.
   */
  proc initData(A: [] Complex_type, vid: VariantID) where isRectangularArr(A) {
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
  proc initDataRandSign(A: [] Real_type, len: int=A.size, vid: VariantID) where isRectangularArr(A) {
    srand(4793);

    const factor = if data_init_count % 2 then 0.1 else 0.2;
    inline proc signfact return if rand():Real_type/RAND_MAX < 0.5 then -1.0 else 1.0;

    for (a,i) in zip(A, 0..<len) do
      a = (signfact*factor*(i + 1.1)/(i + 1.12345)):Real_type;

    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random values.
   */
  proc initDataRandValue(A: [] Real_type, len: int=A.size, vid: VariantID) where isRectangularArr(A) {
    srand(4793);

    for (a, i) in zip(A, 0..<len) do
      a = rand():Real_type/RAND_MAX;

    incDataInitCount();
  }

  proc calcChecksum(const A: [] Real_type, len:int=A.size, scale_factor: Real_type=1.0): Checksum_type where isRectangularArr(A) {
    var s:Checksum_type = 0;
    for (a,i) in zip(A, 1..len) do s += a*i*scale_factor;
    return s;
  }

  proc calcChecksum(const A: [] Complex_type, len:int=A.size, scale_factor: Real_type=1.0): Checksum_type where isRectangularArr(A) {
    var s:Checksum_type = 0;
    for (a,i) in zip(A, 1..len) do s += (a.re+a.im)*i*scale_factor;
    return s;
  }
}
