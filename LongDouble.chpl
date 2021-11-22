/* A module to provide a "long double" type to Chapel that matches the
   underlying C compiler's "long double" type.

   https://raw.githubusercontent.com/chapel-lang/chapel/main/test/release/examples/benchmarks/lcals/LongDouble.chpl
 */

private use CPtr;

private use Utils;

require "longdouble.h";

extern type longdouble;

operator :(ld: longdouble, type t)   where isRealType(t)  || isIntegralType(t)  return __primitive("cast", t, ld);
operator :(d:?tt, type t:longdouble) where isRealType(tt) || isIntegralType(tt) return __primitive("cast", t,  d);

operator :(s: string, type t: longdouble) {
  var ld: longdouble;
  assert(sscanf(s.localize().c_str(), "%Lf", c_ptrTo(ld)) == 1);
  return ld;
}

operator +(d: ?t, ld: longdouble): longdouble where isRealType(t) || isIntegralType(t) return __primitive("+", ld, d);
operator -(d: ?t, ld: longdouble): longdouble where isRealType(t) || isIntegralType(t) return __primitive("-", ld, d);
operator *(d: ?t, ld: longdouble): longdouble where isRealType(t) || isIntegralType(t) return __primitive("*", ld, d);
operator /(d: ?t, ld: longdouble): longdouble where isRealType(t) || isIntegralType(t) return __primitive("/", ld, d);

operator +(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("+", ld, d);
operator -(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("-", ld, d);
operator *(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("*", ld, d);
operator /(ld: longdouble, d: ?t): longdouble where t == longdouble || isRealType(t) || isIntegralType(t) return __primitive("/", ld, d);

operator +=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld + d; }
operator -=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld - d; }
operator *=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld * d; }
operator /=(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { ld = ld / d; }

operator  =(ref ld: longdouble, d: ?t) where t == longdouble || isRealType(t) || isIntegralType(t) { __primitive("=", ld, d); }

proc longdouble.writeThis(f) throws return cprintf(f, "%Lf", this);
