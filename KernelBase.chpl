module KernelBase {
  private use Set;
  private use Time;

  private use DataUtils;
  private use TypeDefs;
  private import RunParams;

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
   * \brief Enumeration defining unique id for each FEATURE used in
   * suite.
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

  /*
   * \brief Enumeration defining unique id for each VARIANT in suite.
   */
  enum VariantID {
    Seq = 0,  // using a for-loop
    Seq_2D,   // like Seq but using a 2D structure

    Forall,     // using a forall-loop
    Promotion,  //  ''   promotions
    Reduction,  //  ''   reductions 

    NONE,  // temp
  };

  class Kernel {
    //
    // Static properties of kernel, independent of run
    //
    var kernel_id: KernelID;
    proc name return kernel_id:string;

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
    var flops_per_rep: Index_type;

    // Checksums
    var checksum = 0;
    //var checksum: [0..<VariantID.size] Checksum_type = 0;
    //var checksum: [0..1] Checksum_type = 0;
    var checksum_scale_factor: Checksum_type;

    // Elapsed time in seconds
    var timer: Timer;
    var min_time: Elapsed_type = 0;
    var max_time: Elapsed_type = 0;
    var tot_time: Elapsed_type = 0;

    proc init(kernel_id: KernelID)
    {
      this.kernel_id = kernel_id;
    }

    proc target_problem_size: Index_type
    {
      return
        if RunParams.size_meaning == RunParams.SizeMeaning.Factor then
          (default_prob_size*RunParams.size_factor):Index_type
        else if RunParams.size_meaning == RunParams.SizeMeaning.Direct then
          (RunParams.size):Index_type
        else
          0:Index_type;
    }

    proc run_reps: Index_type
    {
      return
        if RunParams.input_state == RunParams.InputOpt.CheckRun then
          (RunParams.checkrun_reps):Index_type
        else
          (default_reps*RunParams.rep_fact):Index_type;
    }

    proc setUsesFeature(fid: FeatureID) { uses_feature.add(fid); }

    proc usesFeature(fid: FeatureID) { return uses_feature.contains(fid); };

    proc startTimer() {
      timer.clear();
      timer.start();
    }

    proc stopTimer() {
      timer.stop();

      const exec_time = timer.elapsed():Elapsed_type;
      min_time = min(min_time, exec_time);
      max_time = max(max_time, exec_time);
      tot_time += exec_time;
    }

    proc readTimer(unit:TimeUnits=TimeUnits.seconds) { return timer.elapsed(unit); }

    proc log(checksum) {
      //printf("%s: done in %f seconds (%f)\n", name.c_str(), readTimer(), checksum);
      writef("%s: done in %dr seconds (%dr)\n", kernel_id, readTimer(), checksum:real);
    }
  };

  class FIRST_MIN: Kernel {
    var m_x: [1..0] real;
    var m_xmin_init: Real_type;
    var m_initloc: Index_type;
    var m_minloc: Index_type;

    var m_N: Index_type;

    proc init() { super.init(KernelID.Lcals_FIRST_MIN); }

    proc setUp(vid:VariantID)
    {
      m_x = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      m_x[m_N/2] = -1.0e+10;
      m_xmin_init = m_x[0];
      m_initloc = 0;
      m_minloc = -1;
    }

    proc updateChecksum(vid:VariantID) {
      checksum[vid] += m_minloc:Checksum_type;
    }
  }
}
