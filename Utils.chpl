module Utils {
  require "utils.h";

  extern proc elapsed_time(): real;

  inline proc sizeof(type t) param { return numBytes(t); };

  iter count(n, low:int=0) where isIntegralValue(n) {
    for i in low..#n do
      yield i;
  }

  iter count(shape, low:int=0) {
    for i in count(*reduce shape, low) do
      yield i;
  }

  iter count(param tag: iterKind, low:int=0, followThis) where tag == iterKind.follower && followThis.size == 1 {
    for i in followThis(0).translate(low) do
      yield i;
  }

  iter count(param tag: iterKind, const shape, low:int=0, followThis) where tag == iterKind.follower && followThis.size > 1 {
    assert(shape.size == followThis.size);

    var off = shape;

    for param i in 0..<followThis.size-1 by -1 do
      off(i) = shape(i)*off(i+1);

    for i in _follower_counter(off, low, followThis) do
      yield i;
  }

  iter _follower_counter(const offsets, in low:int=0, followThis, param _idx:int=0) {
    if _idx == followThis.size-1 then
      for i in followThis(_idx).translate(low) do
        yield i;
    else
      for r in followThis(_idx) do
        for c in _follower_counter(offsets, r*offsets(_idx+1)+low, followThis, _idx+1) do
          yield c;
  }

  iter _array.flatIdx(low:int=0) where isRectangularArr(this) {
    for i in low..#this.size do
      yield i;
  }

  iter _array.flatIdx(param tag: iterKind, low:int=0, followThis) where isRectangularArr(this) && tag == iterKind.follower {
    for v in Utils.count(tag, this.shape, low, followThis) do
      yield v;
  }

  extern proc printf(fmt: c_string, vals...?numvals): int;
}
