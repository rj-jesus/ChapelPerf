module Executor {
  private use IO;
  private use List;
  private import FileSystem;

  private use KernelBase;
  private use LongDouble;
  private import lcals;
  private import RunParams;

  enum CSVRepMode {
    Timing = 0,
    Speedup,
  };

  var variant_ids:list(VariantID);

  var reference_vid:VariantID;

  proc haveReferenceVariant() { return reference_vid:int < VariantID.size; }

  var kernels = (
      new lcals.DIFF_PREDICT(),
      new lcals.EOS(),
      new lcals.FIRST_DIFF(),
      new lcals.FIRST_MIN(),
  );

  proc init() {
    variant_ids.append(VariantID.Base_Seq);
    variant_ids.append(VariantID.Seq_2D);

    variant_ids.append(VariantID.Forall);
    variant_ids.append(VariantID.Promotion);
    variant_ids.append(VariantID.Reduction);
  }

  proc main() {
    init();
    //// STEP 1: Create suite executor object
    //rajaperf::Executor executor(argc, argv);

    //// STEP 2: Assemble kernels and variants to run
    //executor.setupSuite();

    //// STEP 3: Report suite run summary
    ////         (enable users to catch errors before entire suite is run)
    //executor.reportRunSummary(std::cout);

    // STEP 4: Execute suite
    //runSuite();

    // STEP 5: Generate suite execution reports
    outputRunData();

    writeln("\n\nDONE!!!....");

    return 0;
  }

  proc runSuite() {
    const in_state = RunParams.getInputState();

    if in_state != RunParams.InputOpt.PerfRun && in_state != RunParams.InputOpt.CheckRun then
      return;

    writeln("\n\nRun warmup kernels...");

    var warmup_kernels = (
        //new basic::DAXPY(),
        //new basic::REDUCE3_INT(),
        //new algorithm::SORT(),
        new lcals.FIRST_MIN(),
    );

    for warmup_kernel in warmup_kernels {
      writeln("Kernel : " + warmup_kernel.getName());
      for vid in warmup_kernel.getVariants() {
        if RunParams.showProgress() then
          writeln("   Running " + vid:string + " variant");

        warmup_kernel.execute(vid);
      }
    }

    writeln("\n\nRunning specified kernels and variants...");

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

  proc outputRunData() {
    const in_state = RunParams.getInputState();

    if in_state != RunParams.InputOpt.PerfRun && in_state != RunParams.InputOpt.CheckRun then
      return;

    writeln("\n\nGenerate run report files...");

    //
    // Generate output file prefix (including directory path).
    //
    const outdir:string = RunParams.getOutputDirName();
    const out_fprefix = "./" + RunParams.getOutputFilePrefix();

    if !outdir.isEmpty() then try! {
      FileSystem.mkdir(outdir, mode=0o755, parents=true);
      here.chdir(outdir);
    }

    var filename = out_fprefix + "-timing.csv";
    writeCSVReport(filename, CSVRepMode.Timing, 6 /* prec */);

    if haveReferenceVariant() {
      filename = out_fprefix + "-speedup.csv";
      writeCSVReport(filename, CSVRepMode.Speedup, 3 /* prec */);
    }

    filename = out_fprefix + "-checksum.txt";
    //writeChecksumReport(filename);

    //filename = out_fprefix + "-fom.csv";
    //writeFOMReport(filename);

    //filename = out_fprefix + "-kernels.csv";
    //ofstream file(filename.c_str(), ios::out | ios::trunc);
    //if ( !file ) {
    //  cout << " ERROR: Can't open output file " << filename << endl;
    //}

    //if ( file ) {
    //  bool to_file = true;
    //  writeKernelInfoSummary(file, to_file);
    //}
  }

  proc writeCSVReport(filename:string, mode:CSVRepMode, prec:int(32)) {
    var f = try! open(filename, iomode.cw);
    var channel = try! f.writer();
    //writeln(" ERROR: Can't open output file " + filename);

    //
    // Set basic table formatting parameters.
    //
    const kernel_col_name = "Kernel  ";
    const sepchr = " , ";

    var kercol_width = kernel_col_name.size:uint(32);
    for kernel in kernels do
      kercol_width = max(kercol_width, kernel.getName().size):kercol_width.type;
    kercol_width += 1;

    var varcol_width = for vid in variant_ids do max(prec+2, (vid:string).size):uint(32);

    //
    // Print title line.
    //
    try! channel.write(getReportTitle(mode));

    //
    // Wrtie CSV file contents for report.
    //
    try! {
      channel.write(for vid in variant_ids do sepchr);
      channel.writeln();
    }

    //
    // Print column title line.
    //
    var style:iostyle = defaultIOStyle();
    try! {
      style.leftjustify = 1;
      style.min_width_columns = kercol_width;

      channel.write(kernel_col_name, style);

      for iv in variant_ids.indices {
        style.min_width_columns = varcol_width[iv];

        channel.write(sepchr);
        channel.write(variant_ids[iv]:string, style);
      }
      channel.writeln();
    }

    //
    // Print row of data for variants of each kernel.
    //
    for kern in kernels do try! {
      style.leftjustify = 1;
      style.min_width_columns = kercol_width;
      channel.write(kern.getName(), style);
      for iv in variant_ids.indices {
        const vid = variant_ids[iv];
        channel.write(sepchr);
        style.leftjustify = 0;
        style.min_width_columns = varcol_width[iv];
        if (mode == CSVRepMode.Speedup) && (!kern.hasVariantDefined(reference_vid) || !kern.hasVariantDefined(vid)) then
          channel.write("Not run", style);
        else if (mode == CSVRepMode.Timing) && !kern.hasVariantDefined(vid) then
          channel.write("Not run", style);
        else {
          style.precision = prec;
          style.realfmt = 1;
          channel.write(getReportDataEntry(mode, kern, vid):real, style);
        }
      }
      channel.writeln();
    }

    try! channel.flush();
  } // note file will be closed when file stream goes out of scope

  proc getReportTitle(mode:CSVRepMode): string {
    var title:string;
    select mode {
      when CSVRepMode.Timing do title = "Mean Runtime Report (sec.) ";
      when CSVRepMode.Speedup do
        if haveReferenceVariant() then
          title = "Speedup Report (T_ref/T_var): ref var = "
            + reference_vid:string + " ";
      otherwise do writeln("\n Unknown CSV report mode = " + mode:string);
    };
    return title;
  }


  proc getReportDataEntry(mode:CSVRepMode, kern, vid:VariantID): longdouble {
    var retval:longdouble = 0.0;
    select mode {
      when CSVRepMode.Timing do
        retval = kern.getTotTime(vid)/RunParams.getNumPasses();
      when CSVRepMode.Speedup {
        if haveReferenceVariant() {
          if kern.hasVariantDefined(reference_vid) && kern.hasVariantDefined(vid) then
            retval = kern.getTotTime(reference_vid)/kern.getTotTime(vid);
          else
            retval = 0.0;
          if(false) {
            writeln("Kernel(iv): " + kern.getName() + "(" + vid + ")");
            writeln("\tref_time, tot_time, retval = "
              + kern.getTotTime(reference_vid) + " , "
              + kern.getTotTime(vid) + " , "
              + retval);
          }
        }
      }
      otherwise writeln("\n Unknown CSV report mode = " + mode:string);
    };
    return retval;
  }
}
