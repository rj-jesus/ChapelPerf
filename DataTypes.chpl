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

module DataTypes {
  private use CTypes;

  private use LongDouble;

  config param USE_DOUBLE  = true;
  config param USE_FLOAT   = false;
  config param USE_COMPLEX = true;

  type Index_type = int;
  type Int_type = c_int;
  type Real_type = if USE_DOUBLE then real(64)
                   else if USE_FLOAT then real(32)
                   else compilerError("Real_type is undefined!");
  type Complex_type = if USE_COMPLEX then complex(2*numBits(Real_type))
                      else nothing;
  type Elapsed_type = real;
  type Checksum_type = longdouble;
}
