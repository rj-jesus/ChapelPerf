module Utils {
  require "utils.h";

  extern proc elapsed_time(): real;

  //iter count(n: int, low: int=1) {
  //  for i in low..#n do
  //    yield i;
  //}

  iter count(shape, low: int=0) {
    for i in low..#*reduce(shape) do
      yield i;
  }

  iter count(param tag: iterKind, const ref shape, low: int=0, param _idx:int=0, followThis) where tag == iterKind.follower && _idx == followThis.size-1 {
    for i in followThis(_idx).translate(low) do
      yield i;
  }

  iter count(param tag: iterKind, const ref shape, low: int=0, param _idx:int=0, followThis) where tag == iterKind.follower && _idx  < followThis.size-1 {
    assert(shape.size == followThis.size);

    var off = shape;

    if _idx == 0 then
      for param i in 0..<followThis.size-1 by -1 do
        off(i) = shape(i)*off(i+1);

    for r in followThis(_idx) do
      for c in count(iterKind.follower, off, r*off(_idx+1)+low, _idx+1, followThis) do
        yield c;
  }
}
