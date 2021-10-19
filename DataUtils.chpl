module DataUtils {
  use Random;

  var data_init_count: uint = 0;

  type VariantID = int;  // temporary...

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
  proc allocAndInitData(type t: numeric, len: int, vid: VariantID = 0)
  {
    return allocAndInitData(t, {0..<len});
  }

  proc allocAndInitData(type t: numeric, d: domain, vid: VariantID = 0)
  {
    var a: [d] t;
    initData(a, vid);
    return a;
  }

  proc allocAndInitDataConst(type t: numeric, len: int, val, vid: VariantID = 0) where isCoercible(val.type, t)
  {
    return allocAndInitDataConst(t, {0..<len}, val);
  }

  proc allocAndInitDataConst(type t: numeric, d: domain, val, vid: VariantID = 0) where isCoercible(val.type, t)
  {
    var a: [d] t;
    initDataConst(a, val:t, vid);
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
   * \brief Initialize int array to randomly signed positive and negative
   * values.
   */
  proc initData(a: [] int, len: int, vid: VariantID)
  {
    var randStream = new RandomStream(real, 4793);

    for (i,r) in zip(0..<len, randStream) do
      a[i] = (r < 0.5 ? -1 : 1);

    a[(len * randStream.getNext()):int] = -58;
    a[(len * randStream.getNext()):int] =  19;

    incDataInitCount();
  }

  ///*
  // * \brief Initialize real array to non-random positive values (0.0, 1.0)
  // * based on their array position (index) and the order in which this method
  // * is called.
  // */
  //proc initData(a: [] real, len: int, vid: VariantID)
  //{
  //  const factor = if data_init_count % 2 then 0.1 else 0.2;
  //  for i in 0..<len do a[i] = factor*(i + 1.1)/(i + 1.12345);
  //  incDataInitCount();
  //}

  /*
   * \brief Initialize real array to non-random positive values (0.0, 1.0)
   * based on their array position (index) and the order in which this method
   * is called.
   */
  proc initData(a: [?d] real, vid: VariantID)
  {
    const factor = if data_init_count % 2 then 0.1 else 0.2;

    writeln(d);

    for c in a {
      writeln(c);
      return;
    }
  }

  /*
   * \brief Initialize complex array.
   */
  proc initData(a: [] complex, len: int, vid: VariantID)
  {
    const factor = if data_init_count % 2 then 0.1+0.2i else 0.2+0.3i;
    for i in 0..<len do a[i] = factor*(i + 1.1)/(i + 1.12345);
    incDataInitCount();
  }

  /*
   * \brief Initialize array to constant values.
   */
  proc initDataConst(a: [?d] ?t, val, vid: VariantID) where isCoercible(val.type, t)
  {
    a = val;
    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random sign.
   */
  proc initDataRandSign(a: [] real, len: int, vid: VariantID)
  {
    const factor = if data_init_count % 2 then 0.1 else 0.2;

    var randStream = new RandomStream(real, 4793);

    for (i,r) in zip(0..<len, randStream) {
      const signfact = if r < 0.5 then -1.0 else 1.0;
      a[i] = signfact*factor*(i + 1.1)/(i + 1.12345);
    }

    incDataInitCount();
  }

  /*
   * \brief Initialize real array with random values.
   */
  proc initDataRandValue(a: [], len: int, vid: VariantID)
  {
    fillRandom(a, seed=4793);
    incDataInitCount();
  }

  ///*
  // * \brief Initialize scalar data.
  // */
  //proc initData(out d: real, vid: VariantID)
  //{
  //  const factor = if data_init_count % 2 then 0.1 else 0.2;
  //  d = factor*1.1/1.12345;
  //  incDataInitCount();
  //}
}
