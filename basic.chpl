module basic {
  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private use KernelBase;
  private use Utils;

  class DAXPY: KernelBase {

    proc init() {
      super.init(KernelID.Basic_DAXPY);

      setDefaultProblemSize(1000000);
      setDefaultReps(500);

      setActualProblemSize( getTargetProblemSize() );

      setItsPerRep( getActualProblemSize() );
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize() );
      setFLOPsPerRep(2 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var y = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var x = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var a = initData(Real_type);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              y[i] += a * x[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              y[i] += a * x[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps do
            y += a * x;

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(y, getActualProblemSize());
    }
  }

  class REDUCE3_INT: KernelBase {

    var m_vsum: Int_type;
    var m_vsum_init: Int_type;
    var m_vmin: Int_type;
    var m_vmin_init: Int_type;
    var m_vmax: Int_type;
    var m_vmax_init: Int_type;

    proc init() {
      super.init(KernelID.Basic_REDUCE3_INT);

      setDefaultProblemSize(1000000);
      setDefaultReps(50);

      setActualProblemSize( getTargetProblemSize() );

      setItsPerRep( getActualProblemSize() );
      setKernelsPerRep(1);
      setBytesPerRep( (3*sizeof(Int_type) + 3*sizeof(Int_type)) +
                      (0*sizeof(Int_type) + 1*sizeof(Int_type)) * getActualProblemSize() );
      setFLOPsPerRep(1 * getActualProblemSize() + 1);

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.Reduction);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Reduction_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var vec = allocAndInitData(Int_type, getActualProblemSize(), vid);

      m_vsum = 0;
      m_vsum_init = 0;
      m_vmin = max(Int_type);
      m_vmin_init = max(Int_type);
      m_vmax = min(Int_type);
      m_vmax_init = min(Int_type);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            var vsum = m_vsum_init;
            var vmin = m_vmin_init;
            var vmax = m_vmax_init;

            for i in ibegin..<iend {
              vsum += vec[i];
              vmin = min(vmin, vec[i]);
              vmax = max(vmax, vec[i]);
            }

            m_vsum += vsum;
            m_vmin  = min(m_vmin, vmin);
            m_vmax  = max(m_vmax, vmax);

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            var vsum = m_vsum_init;
            var vmin = m_vmin_init;
            var vmax = m_vmax_init;

            forall i in ibegin..<iend with (+   reduce vsum,
                                            min reduce vmin,
                                            max reduce vmax) {
              vsum += vec[i];
              vmin reduce= vec[i];
              vmax reduce= vec[i];
            }

            m_vsum += vsum;
            m_vmin  = min(m_vmin, vmin);
            m_vmax  = max(m_vmax, vmax);

          }

          stopTimer();
        }

        when VariantID.Reduction_Chpl {
          startTimer();

          for 0..#run_reps {

            var vsum = +   reduce vec;
            var vmin = min reduce vec;
            var vmax = max reduce vec;

            m_vsum += vsum;
            m_vmin  = min(m_vmin, vmin);
            m_vmax  = max(m_vmax, vmax);

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += m_vsum;
      checksum[vid] += m_vmin;
      checksum[vid] += m_vmax;
    }
  }
}
