module Executor {
  private use IO;
  private use List;
  private import FileSystem;

  private use DataTypes;
  private use KernelBase;
  private use LongDouble;
  private import lcals;
  private import RunParams;

  enum CSVRepMode {
    Timing = 0,
    Speedup,
  };

  record FOMGroup {
    var base:VariantID;
    var variants:list(VariantID);
    //std::vector<VariantID> variants;
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

  proc main() throws {
    init();
    //// STEP 1: Create suite executor object
    //rajaperf::Executor executor(argc, argv);

    //// STEP 2: Assemble kernels and variants to run
    //executor.setupSuite();

    // STEP 3: Report suite run summary
    //         (enable users to catch errors before entire suite is run)
    reportRunSummary(stdout);

    // STEP 4: Execute suite
    runSuite();

    // STEP 5: Generate suite execution reports
    outputRunData();

    writeln("\n\nDONE!!!....");

    return 0;
  }

  proc runSuite() {
    const in_state = RunParams.getInputState();

    if in_state != RunParams.InputOpt.PerfRun &&
       in_state != RunParams.InputOpt.CheckRun then
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
        }  // loop over variants
      }  // loop over kernels
    }  // loop over passes through suite
  }

  proc outputRunData() throws {
    const in_state = RunParams.getInputState();

    if in_state != RunParams.InputOpt.PerfRun &&
       in_state != RunParams.InputOpt.CheckRun then
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
    writeChecksumReport(filename);

    filename = out_fprefix + "-fom.csv";
    writeFOMReport(filename);

    filename = out_fprefix + "-kernels.csv";
    writeKernelInfoSummary(filename, true);
  }

  proc writeCSVReport(filename:string, mode:CSVRepMode, prec) throws {
    var channel;

    try {
      channel = open(filename, iomode.cw).writer();
    } catch {
      writeln(" ERROR: Can't open output file " + filename);
      return;
    }

    //
    // Set basic table formatting parameters.
    //
    const kernel_col_name = "Kernel  ";
    const sepchr = " , ";

    var kercol_width = kernel_col_name.size;
    for kern in kernels do
      kercol_width = max(kercol_width, kern.getName().size);
    kercol_width += 1;

    var varcol_width = for vid in variant_ids do max(prec+2, (vid:string).size);

    //
    // Print title line.
    //
    channel.write(getReportTitle(mode));

    //
    // Wrtie CSV file contents for report.
    //
    channel.write(for vid in variant_ids do sepchr);
    channel.writeln();

    //
    // Print column title line.
    //
    channel.writef("%-*s", kercol_width, kernel_col_name);
    for iv in variant_ids.indices do
      channel.writef("%s%-*s", sepchr, varcol_width[iv], variant_ids[iv]);
    channel.writeln();

    //
    // Print row of data for variants of each kernel.
    //
    for kern in kernels {
      channel.writef("%-*s", kercol_width, kern.getName());
      for iv in variant_ids.indices {
        const vid = variant_ids[iv];
        channel.write(sepchr);
        if (mode == CSVRepMode.Speedup) &&
           (!kern.hasVariantDefined(reference_vid) ||
            !kern.hasVariantDefined(vid)) then
          channel.writef("%*s", varcol_width[iv], "Not run");
        else if (mode == CSVRepMode.Timing) &&
                !kern.hasVariantDefined(vid) then
          channel.writef("%*s", varcol_width[iv], "Not run");
        else
          cprintf(channel, "%*.*Lf",
              varcol_width[iv], prec, getReportDataEntry(mode, kern, vid));
      }
      channel.writeln();
    }

    channel.flush();
  } // note files and channels are closed when their variables go out of scope

  proc getReportTitle(mode:CSVRepMode) {
    var title:string;
    select mode {
      when CSVRepMode.Timing do title = "Mean Runtime Report (sec.) ";
      when CSVRepMode.Speedup do
        if haveReferenceVariant() then
          title = "Speedup Report (T_ref/T_var)" +
                  ": ref var = " + reference_vid:string + " ";
      otherwise do writeln("\n Unknown CSV report mode = " + mode:string);
    };
    return title;
  }

  proc getReportDataEntry(mode:CSVRepMode, kern, vid:VariantID): longdouble {
    var retval:longdouble = 0.0;
    select mode {
      when CSVRepMode.Timing do
        retval = kern.getTotTime(vid):longdouble /
                 RunParams.getNumPasses();
      when CSVRepMode.Speedup {
        if haveReferenceVariant() {
          if kern.hasVariantDefined(reference_vid) &&
             kern.hasVariantDefined(vid) then
            retval = kern.getTotTime(reference_vid):longdouble /
                     kern.getTotTime(vid);
          else
            retval = 0.0;
          if false {
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

  proc writeChecksumReport(const ref filename:string) throws {
    var channel;

    try {
      var f = open(filename, iomode.cw);
      channel = f.writer();
    } catch {
      writeln(" ERROR: Can't open output file " + filename);
      return;
    }

    //
    // Set basic table formatting parameters.
    //
    const equal_line = "===================================================================================================";
    const dash_line = "----------------------------------------------------------------------------------------";
    const dash_line_short = "-------------------------------------------------------";
    const dot_line = "........................................................";

    var prec = 20;
    var checksum_width = (prec + 8);

    var namecol_width = 0;
    for kern in kernels do
      namecol_width = max(namecol_width, kern.getName().size);
    for vid in variant_ids do
      namecol_width = max(namecol_width,   (vid:string).size);
    namecol_width += 1;

    //
    // Print title.
    //
    channel.writeln(equal_line);
    channel.writeln("Checksum Report ");
    channel.writeln(equal_line);

    //
    // Print column title line.
    //
    channel.writef("%-*s", namecol_width, "Kernel  ");
    channel.writeln();
    channel.writeln(dot_line);
    channel.writef("%-*s", namecol_width, "Variants  ");
    channel.writef("%-*s", checksum_width, "Checksum  ");
    channel.writef("%-*s", checksum_width, "Checksum Diff (vs. first variant listed)");
    channel.writeln();
    channel.writeln(dash_line);

    //
    // Print checksum and diff against baseline for each kernel variant.
    //
    for kern in kernels {
      channel.writef("%-*s\n", namecol_width, kern.getName());
      channel.writeln(dot_line);

      var cksum_ref:Checksum_type = 0.0;
      for vid in variant_ids do
        if kern.wasVariantRun(vid) {
          cksum_ref = kern.getChecksum(vid);
          break;
        }

      for vid in variant_ids {
        if kern.wasVariantRun(vid) {
          var vcheck_sum = kern.getChecksum(vid);
          var diff = cksum_ref - kern.getChecksum(vid);

          channel.writef("%-*s", namecol_width, vid);
          cprintf(channel, "%#-*.*Lg", checksum_width, prec, vcheck_sum);
          cprintf(channel, "%#-*.*Lg", checksum_width, prec, diff);
          channel.writeln();
        } else {
          channel.writef("%-*s", namecol_width, vid);
          channel.writef("%-*s", checksum_width, "Not Run");
          channel.writef("%-*s", checksum_width, "Not Run");
          channel.writeln();
        }
      }

      channel.writeln();
      channel.writeln(dash_line_short);
    }

    channel.flush();
  } // note files and channels are closed when their variables go out of scope

  proc writeFOMReport(const ref filename:string) throws {
    var fom_groups = getFOMGroups();
    if fom_groups.isEmpty() then return;

    var ik:int;
    var channel;

    try {
      channel = open(filename, iomode.cw).writer();
    } catch {
      writeln(" ERROR: Can't open output file " + filename);
      return;
    }

    //
    // Set basic table formatting parameters.
    //
    const kernel_col_name = "Kernel  ";
    const sepchr = " , ";
    const prec = 2;

    var kercol_width = kernel_col_name.size;
    for kern in kernels do
      kercol_width = max(kercol_width, kern.getName().size);
    kercol_width += 1;

    const fom_col_width = prec+14;

    var ncols = 0;
    for group in fom_groups do
      ncols += group.variants.size;  // num variants to compare
                                     // to each PM baseline
    var col_exec_count: [0..<ncols] int = 0;
    var col_min: [0..<ncols] real =  INFINITY;
    var col_max: [0..<ncols] real = -INFINITY;
    var col_avg: [0..<ncols] real = 0.0;
    var col_stddev: [0..<ncols] real = 0.0;
    var pct_diff: [0..<kernels.size, 0..<ncols] real = 0.0;

    //
    // Print title line.
    //
    channel.write("FOM Report : signed speedup(-)/slowdown(+) for each PM (base vs. RAJA) -> (T_RAJA - T_base) / T_base )");
    channel.write(for iv in 0..<ncols*2 do sepchr);
    channel.writeln();

    channel.write("'OVER_TOL' in column to right if RAJA speedup is over tolerance");
    channel.write(for iv in 0..<ncols*2 do sepchr);
    channel.writeln();

    const pass = ",        ";
    const fail = ",OVER_TOL";

    //
    // Print column title line.
    //
    channel.writef("%-*s", kercol_width, kernel_col_name);
    for group in fom_groups do
      for vid in group.variants do
        channel.writef("%s%-*s%s", sepchr, fom_col_width, vid:string, pass);
    channel.writeln();

    //
    // Write CSV file contents for FOM report.
    //

    //
    // Print row of FOM data for each kernel.
    //
    ik = 0; for kern in kernels {
      channel.writef("%-*s", kercol_width, kern.getName());

      var col = 0;
      for group in fom_groups {
        const base_vid = group.base;
        for comp_vid in group.variants {
          //
          // If kernel variant was run, generate data for it and
          // print (signed) percentage difference from baseline.
          //
          if kern.wasVariantRun(comp_vid) {
            col_exec_count[col] += 1;

            pct_diff[ik, col] =
              (kern.getTotTime(comp_vid) - kern.getTotTime(base_vid)) /
               kern.getTotTime(base_vid);

            var pfstring = if pct_diff[ik, col] <= RunParams.getPFTolerance()
                           then pass else fail;

            channel.writef("%s%-*.*r%s", sepchr,
                fom_col_width, prec, pct_diff[ik, col], pfstring);

            //
            // Gather data for column summaries (unsigned).
            //
            col_min[col] = min(col_min[col], pct_diff[ik, col]);
            col_max[col] = max(col_max[col], pct_diff[ik, col]);
            col_avg[col] += pct_diff[ik, col];
          } else  // variant was not run, print a big fat goose egg...
            channel.writef("%s%-*.*r%s", sepchr,
                fom_col_width, prec, 0.0, pass);

          col += 1;
        }  // loop over group variants
      }  // loop over fom_groups (i.e., columns)
      channel.writeln();
      ik += 1;
    }  // loop over kernels

    //
    // Compute column summary data.
    //

    // Column average...
    for col in 0..<ncols do
      if col_exec_count[col] > 0 then
        col_avg[col] /= col_exec_count[col];
      else
        col_avg[col] = 0.0;

    // Column standard deviaation...
    ik = 0; for kern in kernels {
      var col = 0;
      for group in fom_groups {
        for comp_vid in group.variants {
          if kern.wasVariantRun(comp_vid) then
            col_stddev[col] += ( pct_diff[ik, col] - col_avg[col] ) *
                               ( pct_diff[ik, col] - col_avg[col] );
          col += 1;
        }  // loop over group variants
      }  // loop over groups
      ik += 1;
    }  // loop over kernels

    for col in 0..<ncols do
      if col_exec_count[col] > 0 then
        col_stddev[col] /= col_exec_count[col];
      else
        col_stddev[col] = 0.0;

    //
    // Print column summaries.
    //
    channel.writef("%-*s", kercol_width, " ");
    for iv in 0..<ncols do
      channel.writef("%s%-*s%s", sepchr, fom_col_width, "  ", pass);
    channel.writeln();

    channel.writef("%-*s", kercol_width, "Col Min");
    for c_min in col_min do
      channel.writef("%s%-*.*r%s", sepchr, fom_col_width, prec, c_min, pass);
    channel.writeln();

    channel.writef("%-*s", kercol_width, "Col Max");
    for c_max in col_max do
      channel.writef("%s%-*.*r%s", sepchr, fom_col_width, prec, c_max, pass);
    channel.writeln();

    channel.writef("%-*s", kercol_width, "Col Avg");
    for c_avg in col_avg do
      channel.writef("%s%-*.*r%s", sepchr, fom_col_width, prec, c_avg, pass);
    channel.writeln();

    channel.writef("%-*s", kercol_width, "Col Std Dev");
    for c_std in col_stddev do
      channel.writef("%s%-*.*r%s", sepchr, fom_col_width, prec, c_std, pass);
    channel.writeln();

    channel.flush();
  } // note files and channels are closed when their variables go out of scope

  // This method needs to be adjusted to match our kernels with RAJAPerf's
  proc getFOMGroups() {
    var fom_groups:list(FOMGroup);

    for iv in variant_ids.indices {
      const vid = variant_ids[iv];
      const vname = vid:string;

      if vname.find("Base") != -1 {
        var group:FOMGroup;
        group.base = vid;

        const pos = vname.find("_");
        const pm = vname[pos+1..];

        for ivs in iv+1..<variant_ids.size {
          const vids = variant_ids[ivs];
          if (vids:string).find(pm) != -1 then
            group.variants.append(vids);
        }

        if !group.variants.isEmpty() then
          fom_groups.append(group);
      }  // if variant name contains 'Base'
    }  // iterate over variant ids to run

    if false {  //  RDH DEBUG   (leave this here, it's useful for debugging!)
      writeln("\nFOMGroups...");
      for group in fom_groups {
        writeln("\tBase : " + group.base:string);
        for vid in group.variants do
          writeln("\t\t " + vid:string);
      }
    }

    return fom_groups;
  }

  proc writeKernelInfoSummary(writer_or_filename, param to_file:bool) throws {
    var channel;

    if !to_file then
      channel = writer_or_filename;
    else
      try {
        channel = open(writer_or_filename, iomode.cw).writer();
      } catch {
        writeln(" ERROR: Can't open output file " + writer_or_filename);
        return;
      }

    //
    // Set up column headers and column widths for kernel summary output.
    //
    const kern_head = "Kernels";
    var kercol_width = kern_head.size;

    var psize_width:Index_type = 0;
    var reps_width:Index_type = 0;
    var itsrep_width:Index_type = 0;
    var bytesrep_width:Index_type = 0;
    var flopsrep_width:Index_type = 0;
    var dash_width:Index_type = 0;

    for kern in kernels {
      kercol_width   = max(kercol_width,   kern.getName().size);
      psize_width    = max(psize_width,    kern.getActualProblemSize());
      reps_width     = max(reps_width,     kern.getRunReps());
      itsrep_width   = max(itsrep_width,   kern.getItsPerRep());
      bytesrep_width = max(bytesrep_width, kern.getBytesPerRep());
      flopsrep_width = max(flopsrep_width, kern.getFLOPsPerRep());
    }

    const sepchr = " , ";

    kercol_width += 2;
    dash_width += kercol_width;

    const psize = log10(psize_width:real);
    const psize_head = "Problem size";
    psize_width = max(psize_head.size, psize:Index_type) + 3;
    dash_width += psize_width + sepchr.size;

    const rsize = log10(reps_width:real);
    const rsize_head = "Reps";
    reps_width = max(rsize_head.size, rsize:Index_type) + 3;
    dash_width += reps_width + sepchr.size;

    const irsize = log10(itsrep_width:real);
    const itsrep_head = "Iterations/rep";
    itsrep_width = max(itsrep_head.size, irsize:Index_type) + 3;
    dash_width += itsrep_width + sepchr.size;

    const kernsrep_head = "Kernels/rep";
    const kernsrep_width = max(kernsrep_head.size:Index_type, 4);
    dash_width += kernsrep_width + sepchr.size;

    const brsize = log10(bytesrep_width:real);
    const bytesrep_head = "Bytes/rep";
    bytesrep_width = max(bytesrep_head.size, brsize:Index_type) + 3;
    dash_width += bytesrep_width + sepchr.size;

    const frsize = log10(flopsrep_width:real);
    const flopsrep_head = "FLOPS/rep";
    flopsrep_width = max(flopsrep_head.size, frsize:Index_type) + 3;
    dash_width += flopsrep_width + sepchr.size;

    channel.writef("%-*s", kercol_width, kern_head);
    channel.writef("%s%*s", sepchr, psize_width,    psize_head);
    channel.writef("%s%*s", sepchr, reps_width,     rsize_head);
    channel.writef("%s%*s", sepchr, itsrep_width,   itsrep_head);
    channel.writef("%s%*s", sepchr, kernsrep_width, kernsrep_head);
    channel.writef("%s%*s", sepchr, bytesrep_width, bytesrep_head);
    channel.writef("%s%*s", sepchr, flopsrep_width, flopsrep_head);
    channel.writeln();

    if !to_file then channel.writef("%*s\n", dash_width, "-");

    for kern in kernels {
      channel.writef("%-*s", kercol_width, kern.getName());
      channel.writef("%s%*i", sepchr, psize_width,    kern.getActualProblemSize());
      channel.writef("%s%*i", sepchr, reps_width,     kern.getRunReps());
      channel.writef("%s%*i", sepchr, itsrep_width,   kern.getItsPerRep());
      channel.writef("%s%*i", sepchr, kernsrep_width, kern.getKernelsPerRep());
      channel.writef("%s%*i", sepchr, bytesrep_width, kern.getBytesPerRep());
      channel.writef("%s%*i", sepchr, flopsrep_width, kern.getFLOPsPerRep());
      channel.writeln();
    }

    channel.flush();
  }

  proc reportRunSummary(writer) throws {
    const in_state = RunParams.getInputState();

    if in_state == RunParams.InputOpt.BadInput {

      writer.write("\nRunParams state:\n");
      writer.write("----------------");
      RunParams.print(writer);

      writer.write("\n\nSuite will not be run now due to bad input.");
      writer.write("\n  See run parameters or option messages above.\n\n");
    } else if in_state == RunParams.InputOpt.PerfRun ||
              in_state == RunParams.InputOpt.DryRun ||
              in_state == RunParams.InputOpt.CheckRun {
      if in_state == RunParams.InputOpt.DryRun {
        writer.write("\n\nRAJA performance suite dry run summary....");
        writer.write("\n--------------------------------------\n");

        writer.write("\nInput state:");
        writer.write("\n------------");
        RunParams.print(writer);
      }

      if in_state == RunParams.InputOpt.PerfRun ||
         in_state == RunParams.InputOpt.CheckRun {
        writer.write("\n\nRAJA performance suite run summary....");
        writer.write("\n--------------------------------------\n");
      }

      var ofiles:string;
      if !RunParams.getOutputDirName().isEmpty() {
        ofiles = RunParams.getOutputDirName();
      } else {
        ofiles = ".";
      }
      ofiles += "/" + RunParams.getOutputFilePrefix() + "*";

      writer.writeln("\nHow suite will be run:");
      writer.writeln("\t # passes = " + RunParams.getNumPasses():string);
      if RunParams.getSizeMeaning() == RunParams.SizeMeaning.Factor then
        writer.writeln("\t Kernel size factor = " + RunParams.getSizeFactor():string);
      else if RunParams.getSizeMeaning() == RunParams.SizeMeaning.Direct then
        writer.writeln("\t Kernel size = " + RunParams.getSize():string);
      writer.writeln("\t Kernel rep factor = " + RunParams.getRepFactor():string);
      writer.writeln("\t Output files will be named " + ofiles);

      writer.writeln("\nThe following kernels and variants (when available for a kernel) will be run:");

      writer.write("\nVariants");
      writer.write("\n--------\n");
      for vid in variant_ids do
        writer.writeln(vid:string);
      writer.writeln();

      writeKernelInfoSummary(writer, false);
    }

    writer.flush();
  }
}
