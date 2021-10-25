module RunParams {
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

  config const input_state: InputOpt = InputOpt.Undefined;  /* state of command
                                                               line input */

  config const show_progress: bool = false;  /* true -> show run progress;
                                                false -> do not */

  config const npasses: int = 1;  /* Number of passes through suite  */

  config const rep_fact: real = 1.0;  /* pct of default kernel reps to run */

  config const size_meaning: SizeMeaning = SizeMeaning.Factor;  /* meaning of
                                                                   size value
                                                                   */
  config const size: real = 0;  /* kernel size to run (input option) */
  config const size_factor: real = 1.0;  /* default kernel size multipier
                                            (input option) */

  config const pf_tol: real = 0.1;  /* pct RAJA variant run time can exceed
                                       base for each PM case to pass/fail
                                       acceptance */

  config const checkrun_reps: int = 1;  /* Num reps each kernel is run in check
                                           run */

  config const reference_variant: string;  /* Name of reference variant for
                                              speedup calculations */

  //
  // Arrays to hold input strings for valid/invalid input. Helpful for
  // debugging command line args.
  //
  config const kernel_input: string;
  config const invalid_kernel_input: string;
  config const exclude_kernel_input: string;
  config const invalid_exclude_kernel_input: string;
  config const variant_input: string;
  config const invalid_variant_input: string;
  config const exclude_variant_input: string;
  config const invalid_exclude_variant_input: string;
  config const feature_input: string;
  config const invalid_feature_input: string;
  config const exclude_feature_input: string;
  config const invalid_exclude_feature_input: string;

  config const outdir: string;  /* Output directory name. */
  config const outfile_prefix: string = "RAJAPerf";  /* Prefix for output data
                                                        file names. */

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
}
