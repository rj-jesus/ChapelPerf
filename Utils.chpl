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

module Utils {
  inline proc sizeof(type t) param { return numBytes(t); };

  private use CPtr;
  private use SysCTypes;

  extern const RAND_MAX: c_int;
  extern proc rand(): c_int;
  extern proc srand(seed: c_uint);

  extern proc sscanf(str: c_string, fmt: c_string, args...): c_int;

  private extern proc printf(fmt: c_string, vals...?numvals): c_int;

  private extern proc snprintf(str: c_ptr(c_char), size: size_t, fmt: c_string, args...): c_int;

  proc cprintf(writer, fmt, args...) throws {
    var buf = new c_array(c_char, 255);

    var ret = snprintf(buf:c_ptr(c_char), buf.size:size_t,
                       fmt.localize().c_str(), (...args));

    writer <~> buf:c_string:string;

    return ret;
  }
}
