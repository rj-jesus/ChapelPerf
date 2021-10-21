module lcals {
  use Reflection;
  use Time;
  use Utils;

  use KernelBase;
  use DataUtils;

  proc main() {
    diff_predict();
    diff_predict_2();
    eos();
  }

  proc diff_predict() {
    var kernel = new Kernel(KernelID.Lcals_DIFF_PREDICT);

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
    var kernel = new Kernel(KernelID.Lcals_DIFF_PREDICT);

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
    var kernel = new Kernel(KernelID.Lcals_EOS);

    kernel.default_prob_size = 1000000;
    kernel.default_reps = 500;
    kernel.actual_prob_size = kernel.target_problem_size;

    const array_length = kernel.actual_prob_size + 7;

    kernel.its_per_rep = kernel.actual_prob_size;
    kernel.kernels_per_rep = 1;
    kernel.bytes_per_rep = (1*sizeof(real) + 2*sizeof(real)) * kernel.actual_prob_size +
                           (0*sizeof(real) + 1*sizeof(real)) * array_length;
    kernel.flops_per_rep = 16 * kernel.actual_prob_size;

    kernel.checksum_scale_factor = 0.0001 * kernel.default_prob_size:Checksum_type / kernel.actual_prob_size;

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
}
