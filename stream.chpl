module stream {
  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private use KernelBase;
  private use Utils;

  class ADD: KernelBase {

    proc init() {
      super.init(KernelID.Stream_ADD);

      setDefaultProblemSize(1000000);
      setDefaultReps(1000);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 2*sizeof(Real_type)) *
                     getActualProblemSize());
      setFLOPsPerRep(1 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var a = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var b = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var c = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              c[i] = a[i] + b[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              c[i] = a[i] + b[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps do
            c = a + b;

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(c, getActualProblemSize());
    }
  }

  class COPY: KernelBase {

    proc init() {
      super.init(KernelID.Stream_COPY);

      setDefaultProblemSize(1000000);
      setDefaultReps(1800);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) *
                     getActualProblemSize());
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var a = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var c = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              c[i] = a[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              c[i] = a[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            c[I] = a[I];
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(c, getActualProblemSize());
    }
  }

  class DOT: KernelBase {

    proc init() {
      super.init(KernelID.Stream_DOT);

      setDefaultProblemSize(1000000);
      setDefaultReps(2000);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) +
                     (0*sizeof(Real_type) + 2*sizeof(Real_type)) *
                     getActualProblemSize());
      setFLOPsPerRep(2 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.Reduction);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Reduction_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var a = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var b = allocAndInitData(Real_type, getActualProblemSize(), vid);

      var m_dot: Real_type = 0.0;
      var m_dot_init: Real_type = 0.0;

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            var dot = m_dot_init;

            for i in ibegin..<iend do
              dot += a[i] * b[i];

            m_dot += dot;
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            var dot = m_dot_init;

            forall i in ibegin..<iend with (+ reduce dot) do
              dot += a[i] * b[i];

            m_dot += dot;
          }

          stopTimer();
        }

        when VariantID.Reduction_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            var dot = m_dot_init + +reduce(a[I]*b[I]);
            m_dot += dot;
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += m_dot;
    }
  }

  class MUL: KernelBase {

    proc init() {
      super.init(KernelID.Stream_MUL);

      setDefaultProblemSize(1000000);
      setDefaultReps(1800);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) *
                     getActualProblemSize());
      setFLOPsPerRep(1 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var b = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var c = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var alpha = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              b[i] = alpha * c[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              b[i] = alpha * c[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            b[I] = alpha * c[I];
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(b, getActualProblemSize());
    }
  }

  class TRIAD: KernelBase {

    proc init() {
      super.init(KernelID.Stream_TRIAD);

      setDefaultProblemSize(1000000);
      setDefaultReps(1000);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 2*sizeof(Real_type)) *
                     getActualProblemSize());
      setFLOPsPerRep(2 * getActualProblemSize());

      checksum_scale_factor = 0.001:Checksum_type *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var a = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var b = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var c = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var alpha = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              a[i] = b[i] + alpha * c[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              a[i] = b[i] + alpha * c[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            a[I] = b[I] + alpha * c[I];
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(a, getActualProblemSize(), checksum_scale_factor:Real_type);
    }
  }

}
