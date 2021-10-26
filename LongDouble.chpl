/* A module to provide a "long double" type to Chapel that matches the
   underlying C compiler's "long double" type.

   https://raw.githubusercontent.com/chapel-lang/chapel/main/test/release/examples/benchmarks/lcals/LongDouble.chpl
 */

private use CPtr;
private use SysCTypes;

require "longdouble.h";

extern type longdouble;

operator  :(ld: longdouble, type t) where isRealType(t) || isIntegralType(t) return __primitive("cast", t, ld);

operator  :(d:?tt, type t:longdouble) where isRealType(tt) || isIntegralType(tt) return __primitive("cast", t, d);

operator  +(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("+", ld, d);
operator  -(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("-", ld, d);
operator  *(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("*", ld, d);
operator  /(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("/", ld, d);

operator +=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld + d; }
operator -=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld - d; }
operator *=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld * d; }
operator /=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld / d; }

operator  =(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { __primitive("=", ld, d); }

private extern proc snprintf(str:c_ptr(c_char), size:size_t, fmt: c_string, args...): int;

proc longdouble.writeThis(writer) throws {
  var buf = new c_array(c_char, 255);

  writer.lock();
  const style = writer._style();
  writer.unlock();

  // e.g., printf("%#-*.*Lf\n", width, prec, 123.0);
  var fmt = "%";
  fmt += if style.showpoint != 0 then "#" else "";
  fmt += if style.leftjustify != 0 then "-" else "";
  fmt += "*.*Lg";

  snprintf(buf:c_ptr(c_char), buf.size:size_t, fmt.localize().c_str(),
           style.min_width_columns:c_int, style.precision:c_int, this);

  writer <~> (buf:c_string:string);
}

proc c_snprintf(fmt, args...) {
  var buf = new c_array(c_char, 255);

  snprintf(buf:c_ptr(c_char), buf.size:size_t, fmt.localize().c_str(), (...args));

  return buf:c_string:string;
}

proc longdouble.cvt(fmt, args...): string {
  var buf = new c_array(c_char, 255);

  snprintf(buf:c_ptr(c_char), buf.size:size_t, fmt.localize().c_str(),
           (...args), this);

  return buf:c_string:string;
}

proc cprintf(writer, fmt, args...) {
  var buf = new c_array(c_char, 255);

  var ret = snprintf(buf:c_ptr(c_char), buf.size:size_t,
                     fmt.localize().c_str(), (...args));

  writer <~> buf:c_string:string;

  return ret;
}
