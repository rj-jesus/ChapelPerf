module lcals {
  use Time;

  config const run_reps = 200;
  config const prob_size = 1000000;
  config const factor = 0.1;

  const bytes_per_rep = 20*64/8 * prob_size;
  const flops_per_rep = 9 * prob_size;

  const array_length = prob_size * 14;
  const offset = prob_size;

  const ibegin = 0;
  const iend = prob_size;

  var timer: Timer;

  proc main() {
    diff_predict();
  }

  proc diff_predict() {
    // setup
    var px: [0..#array_length] real;
    var cx: [0..#array_length] real;

    px = 0.0;
    [i in cx.domain] cx[i] = factor*(i + 1.1)/(i + 1.12345);

    timer.start();

    for 0..#run_reps {
      forall i in 0..#prob_size {
        var ar, br, cr: real;

        ar                  = cx[i + offset * 4];       
        br                  = ar - px[i + offset * 4];  
        px[i + offset * 4]  = ar;                       
        cr                  = br - px[i + offset * 5];  
        px[i + offset * 5]  = br;                       
        ar                  = cr - px[i + offset * 6];  
        px[i + offset * 6]  = cr;                       
        br                  = ar - px[i + offset * 7];  
        px[i + offset * 7]  = ar;                       
        cr                  = br - px[i + offset * 8];  
        px[i + offset * 8]  = br;                       
        ar                  = cr - px[i + offset * 9];  
        px[i + offset * 9]  = cr;                       
        br                  = ar - px[i + offset * 10]; 
        px[i + offset * 10] = ar;                       
        cr                  = br - px[i + offset * 11]; 
        px[i + offset * 11] = br;                       
        px[i + offset * 13] = cr - px[i + offset * 12]; 
        px[i + offset * 12] = cr;
      }
    }

    //timer.stop();

    writef("Done in %dr seconds.\n", timer.elapsed());

    var sum = +reduce px;
    writef("CHECKSUM %dr.\n", sum);
  }
}
