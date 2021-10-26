module RunParams {
  private use List;

  /*
   * \brief Enumeration indicating state of input options requested
   */
  enum InputOpt {
    InfoRequest,  /* option requesting information */
    DryRun,       /* report summary of how suite will run w/o running */
    CheckRun,     /* run suite with small rep count to make sure everything
                     works properly */
    PerfRun,      /* input defines a valid performance run, suite will run as
                     specified */
    BadInput,     /* erroneous input given */
    Undefined     /* input not defined (yet) */
  };

  /*
   * \brief Enumeration indicating how to interpret size input
   */
  enum SizeMeaning {
    Unset,    /* indicates value is unset */
    Factor,   /* multiplier on default kernel iteration space */
    Direct,   /* directly use as kernel iteration space */
  };

  config var input_state: InputOpt = InputOpt.Undefined;  /* state of command
                                                             line input */

  config var show_progress: bool = false;  /* true -> show run progress; false
                                              -> do not */

  config var npasses: int = 1;  /* Number of passes through suite  */

  config var rep_fact: real = 1.0;  /* pct of default kernel reps to run */

  config var size_meaning: SizeMeaning = SizeMeaning.Factor;  /* meaning of
                                                                 size value */
  config var size: real = 0;  /* kernel size to run (input option) */
  config var size_factor: real = 1.0;  /* default kernel size multipier (input
                                          option) */

  config var pf_tol: real = 0.1;  /* pct RAJA variant run time can exceed base
                                     for each PM case to pass/fail acceptance
                                     */

  config var checkrun_reps: int = 1;  /* Num reps each kernel is run in check
                                         run */

  config var reference_variant: string;  /* Name of reference variant for
                                            speedup calculations */

  //
  // Arrays to hold input strings for valid/invalid input. Helpful for
  // debugging command line args.
  //
  var kernel_input:list(string);
  var invalid_kernel_input:list(string);
  var exclude_kernel_input:list(string);
  var invalid_exclude_kernel_input:list(string);
  var variant_input:list(string);
  var invalid_variant_input:list(string);
  var exclude_variant_input:list(string);
  var invalid_exclude_variant_input:list(string);
  var feature_input:list(string);
  var invalid_feature_input:list(string);
  var exclude_feature_input:list(string);
  var invalid_exclude_feature_input:list(string);

  config const outdir = "";  /* Output directory name. */
  config const outfile_prefix = "RAJAPerf";  /* Prefix for output data file
                                                names. */

  // Methods to get/set input state

  proc getInputState(): InputOpt { return input_state; }

  proc setInputState(is:InputOpt) { input_state = is; }

  // Getters/setters for processing input and run parameters

  proc showProgress(): bool { return show_progress; }

  proc getNumPasses(): int { return npasses; }

  proc getRepFactor(): real { return rep_fact; }

  proc getSizeMeaning(): SizeMeaning { return size_meaning; }

  proc getSize(): real { return size; }

  proc getSizeFactor(): real { return size_factor; }

  proc getPFTolerance(): real { return pf_tol; }

  proc getCheckRunReps(): int { return checkrun_reps; }

  proc getOutputDirName(): string { return outdir; }
  proc getOutputFilePrefix(): string { return outfile_prefix; }

  proc print(writer) {
    compilerWarning("TODO");
  }
}
