module RunParams {
  private use IO;
  private use List;
  private use Help;

  private use Enums;
  private use Executor;

  var input_state: InputOpt = InputOpt.Undefined;  /* state of command line
                                                      input */

  var show_progress: bool = false;  /* true -> show run progress; false -> do
                                       not */

  var npasses: int = 1;  /* Number of passes through suite  */

  var rep_fact: real = 1.0;  /* pct of default kernel reps to run */

  var size_meaning: SizeMeaning = SizeMeaning.Unset;  /* meaning of size value */
  var size: real = 0;  /* kernel size to run (input option) */
  var size_factor: real = 0.0;  /* default kernel size multipier (input option) */

  var pf_tol: real = 0.1;  /* pct RAJA variant run time can exceed base for
                              each PM case to pass/fail acceptance
                              */

  var checkrun_reps: int = 1;  /* Num reps each kernel is run in check run */

  var reference_variant: string;  /* Name of reference variant for speedup
                                     calculations */

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

  var outdir = "";  /* Output directory name. */
  var outfile_prefix = "RAJAPerf";  /* Prefix for output data file names. */

  // Methods to get/set input state

  proc getInputState()                               { return input_state; }
  proc setInputState(is:InputOpt)                    { input_state = is; }

  // Getters/setters for processing input and run parameters

  proc showProgress()                                { return show_progress; }

  proc getNumPasses()                                { return npasses; }

  proc getRepFactor()                                { return rep_fact; }

  proc getSizeMeaning(): SizeMeaning                 { return size_meaning; }

  proc getSize()                                     { return size; }

  proc getSizeFactor()                               { return size_factor; }

  proc getPFTolerance()                              { return pf_tol; }

  proc getCheckRunReps()                             { return checkrun_reps; }

  proc getReferenceVariant() const ref               { return reference_variant; }

  proc getKernelInput() const ref                    { return kernel_input; }
  proc setInvalidKernelInput(const ref svec)         { invalid_kernel_input = svec; }
  proc getInvalidKernelInput() const ref             { return invalid_kernel_input; }

  proc getExcludeKernelInput() const ref             { return exclude_kernel_input; }
  proc setInvalidExcludeKernelInput(const ref svec)  { invalid_exclude_kernel_input = svec; }
  proc getInvalidExcludeKernelInput() const ref      { return invalid_exclude_kernel_input; }

  proc getVariantInput() const ref                   { return variant_input; }
  proc setInvalidVariantInput(const ref svec)        { invalid_variant_input = svec; }
  proc getInvalidVariantInput() const ref            { return invalid_variant_input; }

  proc getExcludeVariantInput() const ref            { return exclude_variant_input; }
  proc setInvalidExcludeVariantInput(const ref svec) { invalid_exclude_variant_input = svec; }
  proc getInvalidExcludeVariantInput() const ref     { return invalid_exclude_variant_input; }

  proc getFeatureInput() const ref                   { return feature_input; }
  proc setInvalidFeatureInput(const ref svec)        { invalid_feature_input = svec; }
  proc getInvalidFeatureInput() const ref            { return invalid_feature_input; }

  proc getExcludeFeatureInput() const ref            { return exclude_feature_input; }
  proc setInvalidExcludeFeatureInput(const ref svec) { invalid_exclude_feature_input = svec; }
  proc getInvalidExcludeFeatureInput() const ref     { return invalid_exclude_feature_input; }

  proc getOutputDirName() const ref                  { return outdir; }
  proc getOutputFilePrefix() const ref               { return outfile_prefix; }

  /* Print all run params data to given output stream. */
  proc print(writer) throws {
    writer <~> "\n show_progress = " <~> show_progress;
    writer <~> "\n npasses = " <~> npasses;
    writer <~> "\n rep_fact = " <~> rep_fact;
    writer <~> "\n size_meaning = " <~> getSizeMeaning();
    writer <~> "\n size = " <~> size;
    writer <~> "\n size_factor = " <~> size_factor;
    writer <~> "\n pf_tol = " <~> pf_tol;
    writer <~> "\n checkrun_reps = " <~> checkrun_reps;
    writer <~> "\n reference_variant = " <~> reference_variant;
    writer <~> "\n outdir = " <~> outdir;
    writer <~> "\n outfile_prefix = " <~> outfile_prefix;

    writer <~> "\n kernel_input = ";
    for j in kernel_input.indices do
      writer <~> "\n\t" <~> kernel_input[j];
    writer <~> "\n invalid_kernel_input = ";
    for j in invalid_kernel_input.indices do
      writer <~> "\n\t" <~> invalid_kernel_input[j];

    writer <~> "\n exclude_kernel_input = ";
    for j in exclude_kernel_input.indices do
      writer <~> "\n\t" <~> exclude_kernel_input[j];
    writer <~> "\n invalid_exclude_kernel_input = ";
    for j in invalid_exclude_kernel_input.indices do
      writer <~> "\n\t" <~> invalid_exclude_kernel_input[j];

    writer <~> "\n variant_input = ";
    for j in variant_input.indices do
      writer <~> "\n\t" <~> variant_input[j];
    writer <~> "\n invalid_variant_input = ";
    for j in invalid_variant_input.indices do
      writer <~> "\n\t" <~> invalid_variant_input[j];

    writer <~> "\n exclude_variant_input = ";
    for j in exclude_variant_input.indices do
      writer <~> "\n\t" <~> exclude_variant_input[j];
    writer <~> "\n invalid_exclude_variant_input = ";
    for j in invalid_exclude_variant_input.indices do
      writer <~> "\n\t" <~> invalid_exclude_variant_input[j];

    writer <~> "\n feature_input = ";
    for j in feature_input.indices do
      writer <~> "\n\t" <~> feature_input[j];
    writer <~> "\n invalid_feature_input = ";
    for j in invalid_feature_input.indices do
      writer <~> "\n\t" <~> invalid_feature_input[j];

    writer <~> "\n exclude_feature_input = ";
    for j in exclude_feature_input.indices do
      writer <~> "\n\t" <~> exclude_feature_input[j];
    writer <~> "\n invalid_exclude_feature_input = ";
    for j in invalid_exclude_feature_input.indices do
      writer <~> "\n\t" <~> invalid_exclude_feature_input[j];

    writer.writeln();
    writer.flush();
  }

  /*
   ***************************************************************************
   *
   * Parse command line args to set how suite will run.
   *
   ***************************************************************************
   */
  proc parseCommandLineOptions(argv: [] string) throws {
    writeln("\n\nReading command line input...");

    var i = 1; while i < argv.size {
      var opt = argv[i];

      if opt == "--help" || opt == "-h" {

        printHelpMessage(stdout);
        input_state = InputOpt.InfoRequest;

      } else if opt == "--show-progress" || opt == "-sp" {

        show_progress = true;

      } else if opt == "--print-kernels" || opt == "-pk" {

        printFullKernelNames(stdout);
        input_state = InputOpt.InfoRequest;

      } else if opt == "--print-variants" || opt == "-pv" {

        printVariantNames(stdout);
        input_state = InputOpt.InfoRequest;

      } else if opt == "--print-features" || opt == "-pf" {

        printFeatureNames(stdout);
        input_state = InputOpt.InfoRequest;

      } else if opt == "--print-feature-kernels" || opt == "-pfk" {

        printFeatureKernels(stdout);
        input_state = InputOpt.InfoRequest;

      } else if opt == "--print-kernel-features" || opt == "-pkf" {

        printKernelFeatures(stdout);
        input_state = InputOpt.InfoRequest;

      } else if opt == "--npasses" {

        i += 1;
        if i < argv.size then
          npasses = argv[i]:int;
        else {
          writeln("\nBad input: must give --npasses a value for number of passes (int)");
          input_state = InputOpt.BadInput;
        }

      } else if opt == "--repfact" {

        i += 1;
        if i < argv.size then
          rep_fact = argv[i]:real;
        else {
          writeln("\nBad input: must give --rep_fact a value (double)");
          input_state = InputOpt.BadInput;
        }

      } else if opt == "--sizefact" {

        i += 1;
        if i < argv.size {
          if size_meaning == SizeMeaning.Direct {
            writeln("\nBad input: may only set one of --size and --sizefact");
            input_state = InputOpt.BadInput;
          } else {
            size_factor = argv[i]:real;
            if size_factor >= 0.0 then
              size_meaning = SizeMeaning.Factor;
            else {
              writeln("\nBad input: must give --sizefact a POSITIVE value (double)");
              input_state = InputOpt.BadInput;
            }
          }
        } else {
          writeln("\nBad input: must give --sizefact a value (double)");
          input_state = InputOpt.BadInput;
        }

      } else if opt == "--size" {

        i += 1;
        if i < argv.size {
          if size_meaning == SizeMeaning.Factor {
            writeln("\nBad input: may only set one of --size and --sizefact");
            input_state = InputOpt.BadInput;
          } else {
            size = argv[i]:real;
            if size >= 0.0 then
              size_meaning = SizeMeaning.Direct;
            else {
              writeln("\nBad input: must give --size a POSITIVE value (double)");
              input_state = InputOpt.BadInput;
            }
          }
        } else {
          writeln("\nBad input: must give --size a value (int)");
          input_state = InputOpt.BadInput;
        }

      } else if opt == "--pass-fail-tol" || opt == "-pftol" {

        i += 1;
        if i < argv.size then
          pf_tol = argv[i]:real;
        else {
          writeln("\nBad input: must give --pass-fail-tol (or -pftol) a value (double)");
          input_state = InputOpt.BadInput;
        }

      } else if opt == "--kernels" || opt == "-k" {

        var done = false;
        i += 1;
        while i < argv.size && !done {
          opt = argv[i];
          if opt.startsWith('-') {
            i -= 1;
            done = true;
          } else {
            kernel_input.append(opt);
            i += 1;
          }
        }

      } else if opt == "--exclude-kernels" || opt == "-ek" {

        var done = false;
        i += 1;
        while i < argv.size && !done {
          opt = argv[i];
          if opt.startsWith('-') {
            i -= 1;
            done = true;
          } else {
            exclude_kernel_input.append(opt);
            i += 1;
          }
        }

      } else if opt == "--variants" || opt == "-v" {

        var done = false;
        i += 1;
        while i < argv.size && !done {
          opt = argv[i];
          if opt.startsWith('-') {
            i -= 1;
            done = true;
          } else {
            variant_input.append(opt);
            i += 1;
          }
        }

      } else if opt == "--exclude-variants" || opt == "-ev" {

        var done = false;
        i += 1;
        while i < argv.size && !done {
          opt = argv[i];
          if opt.startsWith('-') {
            i -= 1;
            done = true;
          } else {
            exclude_variant_input.append(opt);
            i += 1;
          }
        }

      } else if opt == "--features" || opt == "-f" {

        var done = false;
        i += 1;
        while i < argv.size && !done {
          opt = argv[i];
          if opt.startsWith('-') {
            i -= 1;
            done = true;
          } else {
            feature_input.append(opt);
            i += 1;
          }
        }

      } else if opt == "--exclude-features" || opt == "-ef" {

        var done = false;
        i += 1;
        while i < argv.size && !done {
          opt = argv[i];
          if opt.startsWith('-') {
            i -= 1;
            done = true;
          } else {
            exclude_feature_input.append(opt);
            i += 1;
          }
        }

      } else if opt == "--outdir" || opt == "-od" {

        i += 1;
        if i < argv.size {
          opt = argv[i];
          if opt.startsWith('-') then
            i -= 1;
          else
            outdir = argv[i];
        }

      } else if opt == "--outfile" || opt == "-of" {

        i += 1;
        if i < argv.size {
          opt = argv[i];
          if opt.startsWith('-') then
            i -= 1;
          else
            outfile_prefix = argv[i];
        }

      } else if opt == "--refvar" || opt == "-rv" {

        i += 1;
        if i < argv.size {
          opt = argv[i];
          if opt.startsWith('-') then
            i -= 1;
          else
            reference_variant = argv[i];
        }

      } else if opt == "--dryrun" {

        if input_state != InputOpt.BadInput then
          input_state = InputOpt.DryRun;

      } else if opt == "--checkrun" {

        input_state = InputOpt.CheckRun;

        i += 1;
        if i < argv.size {
          opt = argv[i];
          if opt.startsWith('-') then
            i -= 1;
          else
            checkrun_reps = argv[i]:int;
        }

      } else {

        input_state = InputOpt.BadInput;

        var huh = argv[i];
        writeln("\nUnknown option: " + huh);
        stdout.flush();
      }

      i += 1;
    }

    // Default size and size_meaning if unset
    if size_meaning == SizeMeaning.Unset {
      size_meaning = SizeMeaning.Factor;
      size_factor = 1.0;
    }
  }

  proc printHelpMessage(writer) throws {
    printUsage();
    writer.writeln();

    writer <~> "OTHER OPTIONS:\n";
    writer <~> "==============\n";

    // Done by `printUsage' above
    //writer <~> "  --help, -h (print options with descriptions)\n\n";

    writer <~> "  --show-progress, -sp (print execution progress during run)\n\n";

    writer <~> "  --print-kernels, -pk (print names of available kernels to run)\n\n";

    writer <~> "  --print-variants, -pv (print names of available variants to run)\n\n";

    writer <~> "  --print-features, -pf (print names of RAJA features exercised in Suite)\n\n";

    writer <~> "  --print-feature-kernels, -pfk \n"
           <~> "       (print names of kernels that use each feature)\n\n";

    writer <~> "  --print-kernel-features, -pkf \n"
           <~> "       (print names of features used by each kernel)\n\n";

    writer <~> "  --npasses <int> [default is 1]\n"
           <~> "       (num passes through Suite)\n";
    writer <~> "    Example...\n"
           <~> "      --npasses 2 (runs complete Suite twice\n\n";

    writer <~> "  --repfact <double> [default is 1.0]\n"
           <~> "       (multiplier on default # reps to run each kernel)\n";
    writer <~> "    Example...\n"
           <~> "      --repfact 0.5 (runs kernels 1/2 as many times as default)\n\n";

    writer <~> "  --sizefact <double> [default is 1.0]\n"
           <~> "       (fraction of default kernel sizes to run)\n"
           <~> "       (may not be set if --size is set)\n";
    writer <~> "    Example...\n"
           <~> "      --sizefact 2.0 (kernels will run with size twice the default)\n\n";

    writer <~> "  --size <int> [no default]\n"
           <~> "       (kernel size to run for all kernels)\n"
           <~> "       (may not be set if --sizefact is set)\n";
    writer <~> "    Example...\n"
           <~> "      --size 1000000 (runs kernels with size ~1,000,000)\n\n";

    writer <~> "  --pass-fail-tol, -pftol <double> [default is 0.1; i.e., 10%]\n"
           <~> "       (slowdown tolerance for RAJA vs. Base variants in FOM report)\n";
    writer <~> "    Example...\n"
           <~> "      -pftol 0.2 (RAJA kernel variants that run 20% or more slower than Base variants will be reported as OVER_TOL in FOM report)\n\n";

    writer <~> "  --kernels, -k <space-separated strings> [Default is run all]\n"
           <~> "       (names of individual kernels and/or groups of kernels to run)\n";
    writer <~> "    Examples...\n"
           <~> "      --kernels Polybench (run all kernels in Polybench group)\n"
           <~> "      -k INIT3 MULADDSUB (run INIT3 and MULADDSUB kernels)\n"
           <~> "      -k INIT3 Apps (run INIT3 kernel and all kernels in Apps group)\n\n";

    writer <~> "  --exclude-kernels, -ek <space-separated strings> [Default is exclude none]\n"
           <~> "       (names of individual kernels and/or groups of kernels to exclude)\n";
    writer <~> "    Examples...\n"
           <~> "      --exclude-kernels Polybench (exclude all kernels in Polybench group)\n"
           <~> "      -ek INIT3 MULADDSUB (exclude INIT3 and MULADDSUB kernels)\n"
           <~> "      -ek INIT3 Apps (exclude INIT3 kernel and all kernels in Apps group)\n\n";

    writer <~> "  --variants, -v <space-separated strings> [Default is run all]\n"
           <~> "       (names of variants to run)\n";
    writer <~> "    Examples...\n"
           <~> "      --variants RAJA_CUDA (run all RAJA_CUDA kernel variants)\n"
           <~> "      -v Base_Seq RAJA_CUDA (run Base_Seq and  RAJA_CUDA variants)\n\n";

    writer <~> "  --exclude-variants, -ev <space-separated strings> [Default is exclude none]\n"
           <~> "       (names of variants to exclude)\n";
    writer <~> "    Examples...\n"
           <~> "      --exclude-variants RAJA_CUDA (exclude all RAJA_CUDA kernel variants)\n"
           <~> "      -ev Base_Seq RAJA_CUDA (exclude Base_Seq and  RAJA_CUDA variants)\n\n";

    writer <~> "  --features, -f <space-separated strings> [Default is run all]\n"
           <~> "       (names of features to run)\n";
    writer <~> "    Examples...\n"
           <~> "      --features Forall (run all kernels that use RAJA forall)\n"
           <~> "      -f Forall Reduction (run all kernels that use RAJA forall or RAJA reductions)\n\n";

    writer <~> "  --exclude-features, -ef <space-separated strings> [Default is exclude none]\n"
           <~> "       (names of features to exclude)\n";
    writer <~> "    Examples...\n"
           <~> "      --exclude-features Forall (exclude all kernels that use RAJA forall)\n"
           <~> "      -ef Forall Reduction (exclude all kernels that use RAJA forall or RAJA reductions)\n\n";

    writer <~> "  --outdir, -od <string> [Default is current directory]\n"
           <~> "       (directory path for output data files)\n";
    writer <~> "    Examples...\n"
           <~> "      --outdir foo (output files to ./foo directory\n"
           <~> "      -od /nfs/tmp/me (output files to /nfs/tmp/me directory)\n\n";

    writer <~> "  --outfile, -of <string> [Default is RAJAPerf]\n"
           <~> "       (file name prefix for output files)\n";
    writer <~> "    Examples...\n"
           <~> "      --outfile mydata (output data will be in files 'mydata*')\n"
           <~> "      -of dat (output data will be in files 'dat*')\n\n";

    writer <~> "  --refvar, -rv <string> [Default is none]\n"
           <~> "       (reference variant for speedup calculation)\n\n";
    writer <~> "    Example...\n"
           <~> "      --refvar Base_Seq (speedups reported relative to Base_Seq variants)\n\n";

    writer <~> "  --dryrun (print summary of how Suite will run without running it)\n\n";

    writer <~> "  --checkrun <int> [default is 1]\n"
           <~> "       (run each kernel a given number of times; usually to check things are working properly or to reduce aggregate execution time)\n";
    writer <~> "    Example...\n"
           <~> "      --checkrun 2 (run each kernel twice)\n\n";

    writer.writeln();
    writer.flush();
  }

  proc printFullKernelNames(writer) throws {
    writer <~> "\nAvailable kernels (<group name>_<kernel name>):";
    writer <~> "\n-----------------------------------------\n";
    for kid in KernelID do
      writer.writeln(getFullKernelName(kid));
    writer.flush();
  }

  proc printVariantNames(writer) throws {
    writer <~> "\nAvailable variants:";
    writer <~> "\n-------------------\n";
    for vid in VariantID do
      writer.writeln(getVariantName(vid));
    writer.flush();
  }

  proc printFeatureNames(writer) throws {
    writer <~> "\nAvailable features:";
    writer <~> "\n-------------------\n";
    for fid in FeatureID do
      writer.writeln(getFeatureName(fid));
    writer.flush();
  }

  proc printFeatureKernels(writer) throws {
    writer <~> "\nAvailable features and kernels that use each:";
    writer <~> "\n---------------------------------------------\n";
    for tfid in FeatureID {
      writer.writeln(getFeatureName(tfid));
      for tkid in KernelID {
        var kern = getKernelObject(tkid);
        if kern.usesFeature(tfid) then
          writer.writeln("\t" + getFullKernelName(tkid));
        delete kern;
      }  // loop over kernels
      writer.writeln();
    }  // loop over features
    writer.flush();
  }

  proc printKernelFeatures(writer) throws {
    writer <~> "\nAvailable kernels and features each uses:";
    writer <~> "\n-----------------------------------------\n";
    for tkid in KernelID {
      writer.writeln(getFullKernelName(tkid));
      var kern = getKernelObject(tkid);
      for tfid in FeatureID {
        if kern.usesFeature(tfid) then
          writer.writeln("\t" + getFeatureName(tfid));
      }  // loop over features
      delete kern;
    }  // loop over kernels
    writer.flush();
  }
}
