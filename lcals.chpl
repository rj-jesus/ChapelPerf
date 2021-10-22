module lcals {
  use Reflection;
  use Time;
  use Utils;

  private    use DataUtils;
  private    use KernelBase;
  private    use LongDouble;
  private    use TypeDefs;
  private import RunParams;

  use List;
  use LinkedLists;

  proc main() {
    //diff_predict();
    //diff_predict_2();
    //eos();
    //first_diff();
    //first_min();

    // void Executor::runSuite()

    var kernels = makeList(new FIRST_MIN());
    //var kernels = [new KernelBase(KernelID.NONE), new FIRST_DIFF(), new FIRST_MIN()];
    var k = new KernelBase(KernelID.NONE);
    ref k2 = k;
    var kernels: LinkedList(KernelBase);

    writeln(kernels);

    //kernels.push_back(new FIRST_MIN());

    //for ip in 0..#RunParams.getNumPasses() {
    //  if RunParams.showProgress() then
    //    writeln("\nPass through suite # " + ip:string);

    //  for kernel in kernels {
    //    if RunParams.showProgress() then
    //      writeln("\nRun kernel -- " + kernel.getName());

    //    for vid in kernel.getVariants() {
    //      if RunParams.showProgress() then
    //        writeln("   Running " + vid:string + " variant\n");

    //      kernel.execute(vid);
    //    } // loop over variants
    //  } // loop over kernels
    //} // loop over passes through suite
  }

  proc diff_predict() {
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
    const offset = kernel.actual_prob_size;

    var px = allocAndInitDataConst(real, array_length, 0.0);
    var cx = allocAndInitData(real, array_length);

    const run_reps = kernel.run_reps;
    const ibegin = 0;
    const iend = kernel.actual_prob_size;

    // Run
    kernel.startTimer();

    for 0..#run_reps {
      forall i in ibegin..<iend {
        var ar, br, cr: real;

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

    kernel.stopTimer();

    const checksum = calcChecksum(px, array_length);
    kernel.log(checksum);
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

  proc eos() {
    var kernel = new KernelBase(KernelID.Lcals_EOS);

    kernel.default_prob_size = 1000000;
    kernel.default_reps = 500;
    kernel.actual_prob_size = kernel.target_problem_size;

    const array_length = kernel.actual_prob_size + 7;

    kernel.its_per_rep = kernel.actual_prob_size;
    kernel.kernels_per_rep = 1;
    kernel.bytes_per_rep = (1*sizeof(real) + 2*sizeof(real)) * kernel.actual_prob_size +
                           (0*sizeof(real) + 1*sizeof(real)) * array_length;
    kernel.flops_per_rep = 16 * kernel.actual_prob_size;

    kernel.checksum_scale_factor = 0.0001:longdouble * kernel.default_prob_size:Checksum_type / kernel.actual_prob_size;

    kernel.setUsesFeature(FeatureID.Forall);

    // Setup
    var x = allocAndInitDataConst(real, array_length, 0.0);
    var y = allocAndInitData(real, array_length);
    var z = allocAndInitData(real, array_length);
    var u = allocAndInitData(real, array_length);

    var q = initData();
    var r = initData();
    var t = initData();

    const run_reps = kernel.run_reps;
    const ibegin = 0;
    const iend = kernel.actual_prob_size;

    // Run
    kernel.startTimer();

    for 0..#run_reps {
      for i in ibegin..<iend {
          x[i] = u[i] + r*( z[i] + r*y[i] ) +
                        t*( u[i+3] + r*( u[i+2] + r*u[i+1] ) +
                                     t*( u[i+6] + q*( u[i+5] + q*u[i+4] ) ) );
      }
    }

    kernel.stopTimer();

    const checksum = calcChecksum(x, kernel.actual_prob_size, kernel.checksum_scale_factor);
    kernel.log(checksum);
  }

  proc first_diff() {
    var kernel = new KernelBase(KernelID.Lcals_FIRST_DIFF);

    kernel.default_prob_size = 1000000;
    kernel.default_reps = 2000;

    kernel.actual_prob_size = kernel.target_problem_size;

    const N = kernel.actual_prob_size+1;

    kernel.its_per_rep = kernel.actual_prob_size;
    kernel.kernels_per_rep = 1;
    kernel.bytes_per_rep = (1*sizeof(Real_type) + 0*sizeof(Real_type)) * kernel.actual_prob_size +
                           (0*sizeof(Real_type) + 1*sizeof(Real_type)) * N;
    kernel.flops_per_rep = 1 * kernel.actual_prob_size;

    kernel.setUsesFeature(FeatureID.Forall);

    // Setup
    var x = allocAndInitDataConst(real, N, 0.0);
    var y = allocAndInitData(real, N);

    const run_reps = kernel.run_reps;
    const ibegin = 0;
    const iend = kernel.actual_prob_size;

    // Run
    kernel.startTimer();

    for 0..#run_reps {
      for i in ibegin..<iend do
        x[i] = y[i+1] - y[i];
    }

    kernel.stopTimer();

    const checksum = calcChecksum(x, kernel.actual_prob_size);
    kernel.log(checksum);
  }

  proc first_min() {
    var kernel = new KernelBase(KernelID.Lcals_FIRST_MIN);

    kernel.default_prob_size = 1000000;
    kernel.default_reps = 100;

    kernel.actual_prob_size = kernel.target_problem_size;

    const N = kernel.actual_prob_size;

    kernel.its_per_rep = kernel.actual_prob_size;
    kernel.kernels_per_rep = 1;
    kernel.bytes_per_rep = (1*sizeof(Real_type ) + 1*sizeof(Real_type )) +
                           (1*sizeof(Index_type) + 1*sizeof(Index_type)) +
                           (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * N;
    kernel.flops_per_rep = 0;

    kernel.setUsesFeature(FeatureID.Forall);
    kernel.setUsesFeature(FeatureID.Reduction);

    // Setup
    var x = allocAndInitDataConst(real, N, 0.0);
    x[N/2] = -1.0e+10;
    var xmin_init = x[0];
    var initloc = 0;
    var minloc = -1;

    const run_reps = kernel.run_reps;
    const ibegin = 0;
    const iend = kernel.actual_prob_size;

    // Run
    kernel.startTimer();

    for 0..#run_reps {
      var mymin = (xmin_init, initloc);

      for i in ibegin..<iend do
        if x[i] < mymin(0) then
          mymin = (x[i], i);

      minloc = max(minloc, mymin(1));
    }

    kernel.stopTimer();

    const checksum = minloc:Checksum_type;
    kernel.log(checksum);
  }

  class FIRST_DIFF: KernelBase {
    proc init() {
      super.init(KernelID.Lcals_FIRST_DIFF);
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

      setVariantDefined(VariantID.Seq      );
      setVariantDefined(VariantID.Reduction);
    }

    proc execute(vid:VariantID)
    {
      resetTimer();
      resetDataInitCount();

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
        when VariantID.Seq {
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

        when VariantID.Reduction {
          startTimer();

          for 0..#run_reps {
            var (myminval, myminloc) = minloc reduce zip(x, x.domain);
            m_minloc = max(m_minloc, myminloc);
          }

          stopTimer();
        }
      }

      // update checksum
      checksum[vid:int] += m_minloc:Checksum_type;
    }
  }
}
