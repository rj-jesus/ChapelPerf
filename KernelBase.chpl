module Kernel {
  use Set;
  import RunParams;

  type Index_type = uint;

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
    Apps_COUPLE,
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

  /*
   * \brief Enumeration defining unique id for each (RAJA) FEATURE used in
     suite.
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

  class KernelBase {
    //
    // Static properties of kernel, independent of run
    //
    var kernel_id: KernelID;
    var name: string;

    var default_prob_size: Index_type;
    var default_reps: Index_type;

    var actual_prob_size: Index_type;

    var uses_feature: set(FeatureID);

    //
    // Properties of kernel dependent on how kernel is run
    //
    var its_per_rep: Index_type;
    var kernels_per_rep: Index_type;
    var bytes_per_rep: Index_type;
    var FLOPs_per_rep: Index_type;

    proc getTargetProblemSize(): Index_type
    {
      return
        if RunParams.size_meaning == RunParams.SizeMeaning.Factor then
          default_prob_size*RunParams.size_factor
        else if RunParams.size_meaning == RunParams.SizeMeaning.Direct then
          RunParams.size
        else
          0;
    }

    proc getRunReps(): Index_type
    {
      return
        if RunParams.input_state == RunParams.InputOpt.CheckRun then
          RunParams.checkrun_reps
        else
          default_reps*RunParams.rep_fact;
    }

    proc setUsesFeature(fid: FeatureID) { uses_feature += fid; }
    proc usesFeature(fid: FeatureID) { return uses_feature.contains(fid); };
  };

  proc calcChecksum(const arr: [] real, scale_factor: real): real {
    return + reduce ((1..arr.size)*arr*scale_factor);
  }

  proc calcChecksum(const arr: [] complex, scale_factor: real): real {
    return + reduce ((1..arr.size)*(arr.re+arr.im)*scale_factor);
  }

  proc main() {
    writeln(calcChecksum([1.1,2.1,3.1,4.1], 1.0));
  }
}
