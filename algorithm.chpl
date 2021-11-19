module algorithm {
  private use Sort;

  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private use KernelBase;
  private use Utils;

  class SORT: KernelBase {

    proc init() {
      super.init(KernelID.Algorithm_SORT);

      setDefaultProblemSize(1000000);
      setDefaultReps(20);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) * getActualProblemSize());  // touched data size, not actual number of stores and loads
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Sort);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataRandValue(Real_type, getActualProblemSize()*getRunReps(), vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for irep in 0..<run_reps do
            sort(x[ibegin+iend*irep..<iend*irep+iend]);

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize()*getRunReps());
    }
  }

  class SORTPAIRS: KernelBase {

    record PairComparator { proc key(a) return abs(a[0]); }
    const pairComparator: PairComparator;

    proc init() {
      super.init(KernelID.Algorithm_SORTPAIRS);

      setDefaultProblemSize(1000000);
      setDefaultReps(20);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((2*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize());  // touched data size, not actual number of stores and loads
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Sort);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataRandValue(Real_type, getActualProblemSize()*getRunReps(), vid);
      var i = allocAndInitDataRandValue(Real_type, getActualProblemSize()*getRunReps(), vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            var vector_of_pairs = forall iemp in ibegin..<iend do
                                    (x[iend*irep + iemp], i[iend*irep + iemp]);

            sort(vector_of_pairs, comparator=pairComparator);

            forall (iemp, pair) in zip(ibegin..<iend, vector_of_pairs) do
              (x[iend*irep + iemp], i[iend*irep + iemp]) = pair;

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize()*getRunReps());
      checksum[vid] += calcChecksum(i, getActualProblemSize()*getRunReps());
    }
  }
}
