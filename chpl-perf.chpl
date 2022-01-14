module ChplPerf {
  private use IO;
  private use Executor only getCout;

  private import Executor;
  private import RunParams;

  proc main(args:[] string) throws {
    // STEP 1: Create suite executor object
    RunParams.parseCommandLineOptions(args);

    // STEP 2: Assemble kernels and variants to run
    Executor.setupSuite();

    // STEP 3: Report suite run summary
    //         (enable users to catch errors before entire suite is run)
    Executor.reportRunSummary(getCout());

    // STEP 4: Execute suite
    Executor.runSuite();

    // STEP 5: Generate suite execution reports
    Executor.outputRunData();

    getCout().writeln("\n\nDONE!!!....");

    return 0;
  }
}
