module Executor {
  private use lcals;
  private import RunParams;

  proc main() {
    var kernels = (
        new DIFF_PREDICT(),
        new EOS(),
        new FIRST_DIFF(),
        new FIRST_MIN(),
    );

    for ip in 0..#RunParams.getNumPasses() {
      if RunParams.showProgress() then
        writeln("\nPass through suite # " + ip:string);

      for kernel in kernels {
        if RunParams.showProgress() then
          writeln("\nRun kernel -- " + kernel.getName());

        for vid in kernel.getVariants() {
          if RunParams.showProgress() then
            writeln("   Running " + vid:string + " variant\n");

          kernel.execute(vid);
        } // loop over variants
      } // loop over kernels
    } // loop over passes through suite
  }
}
