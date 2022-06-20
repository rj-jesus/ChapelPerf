//----------------------------------------------------------------------------
// Copyright (C) 2022 Ricardo Jesus, EPCC, United Kingdom.
// All rights reserved.
//
// Redistribution and use of this software, with or without modification, is
// permitted provided that the following conditions are met:
//
// 1. Redistributions of this software must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//----------------------------------------------------------------------------

/* A module to provide a "long double" type to Chapel that matches the
   underlying C compiler's "long double" type. Based on:
   https://raw.githubusercontent.com/chapel-lang/chapel/1.25.0/test/release/examples/benchmarks/lcals/LongDouble.chpl
 */

private use CTypes;

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
