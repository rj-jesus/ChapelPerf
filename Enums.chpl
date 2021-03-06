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

module Enums {

  enum CSVRepMode {
    Timing = 0,
    Speedup,
  };

  /*!
   ***************************************************************************
   *
   * \brief Enumeration defining unique id for each group of kernels in suite.
   *
   ***************************************************************************
   */
  enum GroupID {
    Basic = 0,
    Lcals,
    Polybench,
    Stream,
    Apps,
    Algorithm,
  };

  /* Return group name associated with GroupID enum value. */
  proc getGroupName(gid:GroupID) { return gid:string; }

  /*!
   ***************************************************************************
   *
   * \brief Enumeration defining unique id for each KERNEL in suite.
   *
   ***************************************************************************
   */
  enum KernelID {
    //
    // Basic kernels...
    //
    Basic_DAXPY = 0,
    Basic_IF_QUAD,
    Basic_INIT3,
    Basic_INIT_VIEW1D,
    Basic_INIT_VIEW1D_OFFSET,
    Basic_MAT_MAT_SHARED,
    Basic_MULADDSUB,
    Basic_NESTED_INIT,
    Basic_PI_ATOMIC,
    Basic_PI_REDUCE,
    Basic_REDUCE3_INT,
    Basic_TRAP_INT,

    //
    // Lcals kernels...
    //
    Lcals_DIFF_PREDICT,
    Lcals_EOS,
    Lcals_FIRST_DIFF,
    Lcals_FIRST_MIN,
    Lcals_FIRST_SUM,
    Lcals_GEN_LIN_RECUR,
    Lcals_HYDRO_1D,
    Lcals_HYDRO_2D,
    Lcals_INT_PREDICT,
    Lcals_PLANCKIAN,
    Lcals_TRIDIAG_ELIM,

    //
    // Polybench kernels...
    //
    Polybench_2MM,
    Polybench_3MM,
    Polybench_ADI,
    Polybench_ATAX,
    Polybench_FDTD_2D,
    Polybench_FLOYD_WARSHALL,
    Polybench_GEMM,
    Polybench_GEMVER,
    Polybench_GESUMMV,
    Polybench_HEAT_3D,
    Polybench_JACOBI_1D,
    Polybench_JACOBI_2D,
    Polybench_MVT,

    //
    // Stream kernels...
    //
    Stream_ADD,
    Stream_COPY,
    Stream_DOT,
    Stream_MUL,
    Stream_TRIAD,

    //
    // Apps kernels...
    //
    //Apps_COUPLE,
    Apps_DEL_DOT_VEC_2D,
    Apps_DIFFUSION3DPA,
    Apps_ENERGY,
    Apps_FIR,
    Apps_HALOEXCHANGE,
    Apps_HALOEXCHANGE_FUSED,
    Apps_LTIMES,
    Apps_LTIMES_NOVIEW,
    Apps_MASS3DPA,
    Apps_PRESSURE,
    Apps_VOL3D,

    //
    // Algorithm kernels...
    //
    Algorithm_SORT,
    Algorithm_SORTPAIRS,
  };

  /* Return kernel name associated with KernelID enum value. */
  proc getKernelName(kid:KernelID) {
    const pos = (kid:string).find("_");
    return (kid:string)[pos+1..];
  }

  /* Return full kernel name associated with KernelID enum value. */
  proc getFullKernelName(kid:KernelID) { return kid:string; }

  /*!
   ***************************************************************************
   *
   * \brief Enumeration defining unique id for each FEATURE used in suite.
   *
   ***************************************************************************
   */
  enum FeatureID {
    Forall = 0,
    Kernel,
    Teams,

    Sort,
    Scan,
    Workgroup,

    Reduction,
    Atomic,

    View,
  };

  /* Return feature name associated with FeatureID enum value. */
  proc getFeatureName(fid:FeatureID) { return fid:string; }

  /* Enumeration defining unique id for each VARIANT in suite. */
  enum VariantID {
    Base_Chpl = 0,   // using a for-loop
    Forall_Chpl,     //  ''   a forall-loop
    Promotion_Chpl,  //  ''   promotions
    Reduction_Chpl,  //  ''   reductions

    //Seq_2D,          // like Base but using a 2D structure (not yet supported)
  };

  /* Return variant name associated with VariantID enum value. */
  proc getVariantName(vid:VariantID) { return vid:string; }

  /* Return true if variant associated with VariantID enum value is available
   * to run; else false. */
  proc isVariantAvailable(vid:VariantID) { return true; }

  /*!
   ***************************************************************************
   *
   * \brief Enumeration indicating state of input options requested
   *
   ***************************************************************************
   */
  enum InputOpt {
    InfoRequest,  /* option requesting information */
    DryRun,       /* report summary of how suite will run w/o running */
    CheckRun,     /* run suite with small rep count to make sure everything
                     works properly */
    PerfRun,      /* input defines a valid performance run, suite will run as
                     specified */
    BadInput,     /* erroneous input given */
    Undefined,    /* input not defined (yet) */
  };

  /* make InputOpt's visible without the enum type prefix */
  /*public*/ use InputOpt;  // This is giving me a compiler bug

  /*!
   ***************************************************************************
   *
   * \brief Enumeration indicating how to interpret size input
   *
   ***************************************************************************
   */
  enum SizeMeaning {
    Unset,    /* indicates value is unset */
    Factor,   /* multiplier on default kernel iteration space */
    Direct,   /* directly use as kernel iteration space */
  };

}
