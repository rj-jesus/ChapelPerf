module lcals {
  use Reflection;

  private use DataTypes;
  private use DataUtils;
  private use KernelBase;
  private use Utils;
  private use Enums;

  class DIFF_PREDICT : KernelBase {
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Lcals_DIFF_PREDICT);

      setDefaultProblemSize(1000000);
      setDefaultReps(200);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());

      setKernelsPerRep(1);
      setBytesPerRep((10*sizeof(Real_type) + 10*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(9 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      m_array_length = getActualProblemSize() * 14;
      const offset = getActualProblemSize();

      var px = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var cx = allocAndInitData(Real_type, m_array_length, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {
        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              var ar, br, cr: Real_type;

              ar                  =      cx[i + offset *  4];
              br                  = ar - px[i + offset *  4];
              px[i + offset *  4] = ar;
              cr                  = br - px[i + offset *  5];
              px[i + offset *  5] = br;
              ar                  = cr - px[i + offset *  6];
              px[i + offset *  6] = cr;
              br                  = ar - px[i + offset *  7];
              px[i + offset *  7] = ar;
              cr                  = br - px[i + offset *  8];
              px[i + offset *  8] = br;
              ar                  = cr - px[i + offset *  9];
              px[i + offset *  9] = cr;
              br                  = ar - px[i + offset * 10];
              px[i + offset * 10] = ar;
              cr                  = br - px[i + offset * 11];
              px[i + offset * 11] = br;
              px[i + offset * 13] = cr - px[i + offset * 12];
              px[i + offset * 12] = cr;
            }
          }

          stopTimer();
        }
      }

      // update checksum
      checksum[vid] += calcChecksum(px, m_array_length);
    }
  }

  proc diff_predict_2() {
    var kernel = new KernelBase(KernelID.Lcals_DIFF_PREDICT);

    kernel.default_prob_size = 1000000;
    kernel.default_reps = 200;
    kernel.actual_prob_size = kernel.target_problem_size;
    kernel.its_per_rep = kernel.actual_prob_size;

    kernel.kernels_per_rep = 1;
    kernel.bytes_per_rep = 20 * sizeof(real) * kernel.actual_prob_size;
    kernel.flops_per_rep = 9 * kernel.actual_prob_size;

    kernel.setUsesFeature(FeatureID.Forall);

    // Setup
    const array_length = kernel.actual_prob_size * 14;

    var px = allocAndInitDataConst(real, {0..<14, 0..<kernel.actual_prob_size}, 0);
    var cx = allocAndInitData(real, {0..<14, 0..<kernel.actual_prob_size});

    const run_reps = kernel.run_reps;

    // Run
    kernel.startTimer();

    for 0..#run_reps {
      forall j in cx.domain.dim(1) {
        var ar, br, cr: real;

        ar        =      cx[ 4, j];
        br        = ar - px[ 4, j];
        px[ 4, j] = ar;
        cr        = br - px[ 5, j];
        px[ 5, j] = br;
        ar        = cr - px[ 6, j];
        px[ 6, j] = cr;
        br        = ar - px[ 7, j];
        px[ 7, j] = ar;
        cr        = br - px[ 8, j];
        px[ 8, j] = br;
        ar        = cr - px[ 9, j];
        px[ 9, j] = cr;
        br        = ar - px[10, j];
        px[10, j] = ar;
        cr        = br - px[11, j];
        px[11, j] = br;
        px[13, j] = cr - px[12, j];
        px[12, j] = cr;
      }
    }

    kernel.stopTimer();

    const checksum = calcChecksum(px);
    kernel.log(checksum);
  }

  class EOS : KernelBase {
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Lcals_EOS);

      setDefaultProblemSize(1000000);
      setDefaultReps(500);

      setActualProblemSize(getTargetProblemSize());

      m_array_length = getActualProblemSize() + 7;

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize() +
                      (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_array_length );
      setFLOPsPerRep(16 * getActualProblemSize());

      checksum_scale_factor = 0.0001:Checksum_type * (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var y = allocAndInitData(Real_type, m_array_length, vid);
      var z = allocAndInitData(Real_type, m_array_length, vid);
      var u = allocAndInitData(Real_type, m_array_length, vid);

      var q = initData(vid);
      var r = initData(vid);
      var t = initData(vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {
        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              x[i] = u[i] + r*( z[i] + r*y[i] ) +
                            t*( u[i+3] + r*( u[i+2] + r*u[i+1] ) +
                                         t*( u[i+6] + q*( u[i+5] + q*u[i+4] ) ) );
            }
          }

          stopTimer();
        }
      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize(), checksum_scale_factor:Real_type);
    }
  }

  class FIRST_DIFF: KernelBase {
    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Lcals_FIRST_DIFF);

      setDefaultProblemSize(1000000);
      setDefaultReps(2000);

      setActualProblemSize(getTargetProblemSize());

      m_N = getActualProblemSize()+1;

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type) + 0*sizeof(Real_type)) * getActualProblemSize() +
                      (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_N );
      setFLOPsPerRep(1 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      var y = allocAndInitData(Real_type, m_N, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {
        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              x[i] = y[i+1] - y[i];
          }

          stopTimer();
        }
      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize());
    }
  }

  class FIRST_MIN: KernelBase {
    var m_N: Index_type;
    var m_minloc = -1;

    proc init() {
      super.init(KernelID.Lcals_FIRST_MIN);

      setDefaultProblemSize(1000000);
      setDefaultReps(100);

      setActualProblemSize(getTargetProblemSize());

      m_N = getActualProblemSize();

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type ) + 1*sizeof(Real_type )) +
                      (1*sizeof(Index_type) + 1*sizeof(Index_type)) +
                      (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_N );
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.Reduction);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Reduction_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      x[m_N/2] = -1.0e+10;
      const xmin_init = x[0];
      const initloc = 0;
      m_minloc = -1;

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {
        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            var mymin = (xmin_init, initloc);

            for i in ibegin..<iend do
              if x[i] < mymin(0) then
                mymin = (x[i], i);

            m_minloc = max(m_minloc, mymin(1));
          }

          stopTimer();
        }

        when VariantID.Reduction_Chpl {
          startTimer();

          for 0..#run_reps {
            var (myminval, myminloc) = minloc reduce zip(x, x.domain);
            m_minloc = max(m_minloc, myminloc);
          }

          stopTimer();
        }
      }

      // update checksum
      checksum[vid] += m_minloc:Checksum_type;
    }
  }
}
