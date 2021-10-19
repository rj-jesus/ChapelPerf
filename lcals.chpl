module lcals {
  use Reflection;
  use Time;

  use KernelBase;
  use DataUtils;

  //config const run_reps = 200;
  //config const prob_size = 1000000;
  //config const factor = 0.1;

  //const bytes_per_rep = 20*64/8 * prob_size;
  //const flops_per_rep = 9 * prob_size;

  //const ibegin = 0;
  //const iend = prob_size;

  //var timer: Timer;

  proc main() {
    diff_predict();
    //diff_predict_2();
  }

  proc diff_predict() {
    // Define the kernel
    var kernel = new Kernel();

    kernel.default_prob_size = 1000000;
    kernel.default_reps = 200;
    kernel.actual_prob_size = kernel.getTargetProblemSize();
    kernel.its_per_rep = kernel.actual_prob_size;

    kernel.kernels_per_rep = 1;
    kernel.bytes_per_rep = 20 * numBytes(real) * kernel.actual_prob_size;
    kernel.flops_per_rep = 9 * kernel.actual_prob_size;

    //kernel.setUsesFeature(FeatureID.Forall);

    // setup
    const array_length = kernel.actual_prob_size * 14;
    const offset = kernel.actual_prob_size;

    const run_reps = kernel.getRunReps();
    const prob_size = kernel.actual_prob_size;

    //var px: [0..<array_length] real;
    var px = allocAndInitData(real, array_length);
    var cx: [0..<array_length] real;

    const factor = 0.1;

    px = 0.0;
    [i in cx.domain] cx[i] = factor*(i + 1.1)/(i + 1.12345);

    var timer: Timer; timer.start();

    for 0..#run_reps {
      forall i in 0..<prob_size {
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

    writef("%s: done in %dr seconds.\n", getRoutineName(), timer.elapsed());
  }

  //proc diff_predict_2() {
  //  // setup
  //  var px: [0..<14, 0..<prob_size] real;
  //  var cx: [0..<14, 0..<prob_size] real;

  //  px = 0.0;
  //  forall (i,j) in cx.domain {
  //    var idx = i*cx.shape[1]+j;
  //    cx[i,j] = factor*(idx + 1.1)/(idx + 1.12345);
  //  }

  //  timer.clear(); timer.start();

  //  for 0..#run_reps {
  //    forall j in cx.domain.dim(1) {
  //      var ar, br, cr: real;

  //      ar        =      cx[ 4, j];
  //      br        = ar - px[ 4, j];
  //      px[ 4, j] = ar;
  //      cr        = br - px[ 5, j];
  //      px[ 5, j] = br;
  //      ar        = cr - px[ 6, j];
  //      px[ 6, j] = cr;
  //      br        = ar - px[ 7, j];
  //      px[ 7, j] = ar;
  //      cr        = br - px[ 8, j];
  //      px[ 8, j] = br;
  //      ar        = cr - px[ 9, j];
  //      px[ 9, j] = cr;
  //      br        = ar - px[10, j];
  //      px[10, j] = ar;
  //      cr        = br - px[11, j];
  //      px[11, j] = br;
  //      px[13, j] = cr - px[12, j];
  //      px[12, j] = cr;
  //    }
  //  }

  //  //timer.stop();

  //  writef("%s: done in %dr seconds.\n", getRoutineName(), timer.elapsed());
  //}
}
