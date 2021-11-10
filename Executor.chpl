module Executor {
  private use IO;
  private use List;
  private use OrderedSet;
  private import FileSystem;

  private use DataTypes;
  private use Enums;
  private use KernelBase;
  private use LongDouble;
  private use Utils;
  private import RunParams;

  private import algorithm;
  private import apps;
  private import basic;
  private import lcals;
  private import stream;

  record FOMGroup {
    var base: VariantID;
    var variants: list(VariantID);
  };

  var kernels: list(unmanaged KernelBase);
  var variant_ids: list(VariantID);

  // in RAJAPerf they have `reference_vid', and they mark it `invalid' with an
  // extra enum symbol `NumVariants'
  var reference_vid_int: int = VariantID.size;
  proc reference_vid          { return try! reference_vid_int:VariantID; }
  proc haveReferenceVariant() { return reference_vid_int < VariantID.size; }

  proc setupSuite() throws {
    const in_state = RunParams.getInputState();
    if in_state == InputOpt.InfoRequest || in_state == InputOpt.BadInput then
      return;

    writeln("\nSetting up suite based on input...");

    type Slist = list(string);
    type Svector = list(string);
    //type KIDset = orderedSet(KernelID);
    //type VIDset = orderedSet(VariantID);

    //
    // Determine which kernels to exclude from input.
    // exclude_kern will be non-duplicated ordered set of IDs of kernel to
    // exclude.
    //
    const exclude_kernel_input = RunParams.getExcludeKernelInput();
    const exclude_feature_input = RunParams.getExcludeFeatureInput();

    //var exclude_kern: KIDset;
    var exclude_kern = new orderedSet(KernelID);

    if !exclude_kernel_input.isEmpty() {
      // Make list copy of exclude kernel name input to manipulate for
      // processing potential group names and/or kernel names, next
      var exclude_kern_names = exclude_kernel_input;

      //
      // Search exclude_kern_names for matching group names.
      // groups2exclude will contain names of groups to exclude.
      //
      var groups2exclude:Svector;
      for kern_name in exclude_kern_names {
        for gid in GroupID {
          const group_name = getGroupName(gid);
          if group_name == kern_name then groups2exclude.append(group_name);
        }
      }

      //
      // If group name(s) found in exclude_kern_names, assemble kernels in
      // group(s) to run and remove those group name(s) from exclude_kern_names
      // list.
      //
      for gname in groups2exclude {
        for kid in KernelID do
          if getFullKernelName(kid).find(gname) != -1 then
            exclude_kern.add(kid);

        exclude_kern_names.remove(gname);
      }

      //
      // Look for matching names of individual kernels in remaining
      // exclude_kern_names.
      //
      // Assemble invalid input for warning message.
      //
      var invalid:Svector;

      for kern_name in exclude_kern_names {
        var found_it = false;

        for kid in KernelID {
          if found_it then break;

          if getKernelName(kid) == kern_name ||
             getFullKernelName(kid) == kern_name {
            exclude_kern.add(kid);
            found_it = true;
          }
        }

        if !found_it then invalid.append(kern_name);
      }

      RunParams.setInvalidExcludeKernelInput(invalid);
    }

    if !exclude_feature_input.isEmpty() {

      // First, check for invalid exclude_feature input.
      // Assemble invalid input for warning message.
      //
      var invalid:Svector;

      for feature in exclude_feature_input {
        var found_it = false;

        for tfid in FeatureID {
          if found_it then break;

          if getFeatureName(tfid) == feature then
            found_it = true;
        }

        if !found_it then invalid.append(feature);
      }

      RunParams.setInvalidExcludeFeatureInput(invalid);

      //
      // If feature input is valid, determine which kernels use
      // input-specified features and add to set of kernels to run.
      //
      if RunParams.getInvalidExcludeFeatureInput().isEmpty() {
        for feature in exclude_feature_input {
          var found_it = false;

          for tfid in FeatureID {
            if found_it then break;

            if getFeatureName(tfid) == feature {
              found_it = true;

              for tkid in KernelID {
                var kern = getKernelObject(tkid);
                if kern.usesFeature(tfid) then exclude_kern.add(tkid);
                delete kern;
              }  // loop over kernels
            }  // if input feature name matches feature id
          }  // loop over feature ids until name match is found
        }  // loop over feature name input
      }  // if feature name input is valid
    }

    //
    // Determine which kernels to execute from input.
    // run_kern will be non-duplicated ordered set of IDs of kernel to run.
    //
    const kernel_input = RunParams.getKernelInput();
    const feature_input = RunParams.getFeatureInput();

    //var run_kern: KIDset;
    var run_kern = new orderedSet(KernelID);

    if kernel_input.isEmpty() && feature_input.isEmpty() {

      //
      // No kernels or features specified in input, run them all...
      //
      for tkid in KernelID do
        if !exclude_kern.contains(tkid) then run_kern.add(tkid);

    } else {

      //
      // Need to parse input to determine which kernels to run
      //

      //
      // Look for kernels using features if such input provided
      //
      if !feature_input.isEmpty() {

        // First, check for invalid feature input.
        // Assemble invalid input for warning message.
        //
        var invalid:Svector;

        for feature in feature_input {
          var found_it = false;

          for tfid in FeatureID {
            if found_it then break;
            if getFeatureName(tfid) == feature then found_it = true;
          }

          if !found_it then invalid.append(feature);
        }

        RunParams.setInvalidFeatureInput(invalid);

        //
        // If feature input is valid, determine which kernels use
        // input-specified features and add to set of kernels to run.
        //
        if RunParams.getInvalidFeatureInput().isEmpty() {
          for feature in feature_input {
            var found_it = false;

            for tfid in FeatureID {
              if found_it then break;

              if getFeatureName(tfid) == feature {
                found_it = true;

                for tkid in KernelID {
                  var kern = getKernelObject(tkid);
                  if kern.usesFeature(tfid) &&
                     !exclude_kern.contains(tkid) then
                    run_kern.add(tkid);
                  delete kern;
                }  // loop over kernels
              }  // if input feature name matches feature id
            }  // loop over feature ids until name match is found
          }  // loop over feature name input
        }  // if feature name input is valid
      } // if !feature_input.isEmpty()

      // Make list copy of kernel name input to manipulate for
      // processing potential group names and/or kernel names, next
      var kern_names = kernel_input;

      //
      // Search kern_names for matching group names.
      // groups2run will contain names of groups to run.
      //
      var groups2run:Svector;
      for kern_name in kern_names {
        for gid in GroupID {
          const group_name = getGroupName(gid);
          if group_name == kern_name then groups2run.append(group_name);
        }
      }

      //
      // If group name(s) found in kern_names, assemble kernels in group(s)
      // to run and remove those group name(s) from kern_names list.
      //
      for gname in groups2run {
        for tkid in KernelID {
          if getFullKernelName(tkid).find(gname) != -1 &&
             !exclude_kern.contains(tkid) then
            run_kern.add(tkid);
        }

        kern_names.remove(gname);
      }

      //
      // Look for matching names of individual kernels in remaining kern_names.
      //
      // Assemble invalid input for warning message.
      //
      var invalid:Svector;

      for kern_name in kern_names {
        var found_it = false;

        for tkid in KernelID {
          if found_it then break;

          if getKernelName(tkid) == kern_name ||
             getFullKernelName(tkid) == kern_name {
            if !exclude_kern.contains(tkid) then run_kern.add(tkid);
            found_it = true;
          }
        }

        if !found_it then invalid.append(kern_name);
      }

      RunParams.setInvalidKernelInput(invalid);
    }

    //
    // Assemble set of available variants to run
    // (based on compile-time configuration).
    //
    //var available_var: VIDset;
    var available_var = new orderedSet(VariantID);
    for vid in VariantID do
      if isVariantAvailable(vid) then
        available_var.add(vid);

    //
    // Determine variants to execute from input.
    // run_var will be non-duplicated ordered set of IDs of variants to run.
    //
    const exclude_variant_names = RunParams.getExcludeVariantInput();

    //var exclude_var: VIDset;
    var exclude_var = new orderedSet(VariantID);

    if !exclude_variant_names.isEmpty() {
      //
      // Parse input to determine which variants to exclude.
      //
      // Assemble invalid input for warning message.
      //
      var invalid:Svector;

      for variant_name in exclude_variant_names {
        var found_it = false;

        for vid in available_var {
          if getVariantName(vid) == variant_name {
            exclude_var.add(vid);
            found_it = true;
          }
        }

        if !found_it then invalid.append(variant_name);
      }

      RunParams.setInvalidExcludeVariantInput(invalid);
    }

    //
    // Determine variants to execute from input.
    // run_var will be non-duplicated ordered set of IDs of variants to run.
    //
    const variant_names = RunParams.getVariantInput();

    //var run_var: VIDset;
    var run_var = new orderedSet(VariantID);

    if variant_names.isEmpty() {

      //
      // No variants specified in input options, run all available.
      // Also, set reference variant if specified.
      //
      for vid in available_var {
        if !exclude_var.contains(vid) {
          run_var.add(vid);
          if getVariantName(vid) == RunParams.getReferenceVariant() then
            reference_vid_int = vid:int;
        }
      }

      //
      // Set reference variant if not specified.
      //
      if RunParams.getReferenceVariant().isEmpty() && !run_var.isEmpty() then
        reference_vid_int = run_var.toArray()[0]:int;

    } else {
      //
      // Parse input to determine which variants to run:
      //   - variants to run will be the intersection of available variants
      //     and those specified in input
      //   - reference variant will be set to specified input if available
      //     and variant will be run; else first variant that will be run.
      //
      // Assemble invalid input for warning message.
      //

      var invalid:Svector;

      for variant_name in variant_names {
        var found_it = false;

        for vid in available_var {
          if getVariantName(vid) == variant_name {
            if !exclude_var.contains(vid) {
              run_var.add(vid);
              if getVariantName(vid) == RunParams.getReferenceVariant() then
                reference_vid_int = vid:int;
            }
            found_it = true;
          }
        }

        if !found_it then invalid.append(variant_name);
      }

      //
      // Set reference variant if not specified.
      //
      if RunParams.getReferenceVariant().isEmpty() && !run_var.isEmpty() then
        reference_vid_int = run_var.toArray()[0]:int;

      RunParams.setInvalidVariantInput(invalid);
    }

    //
    // Create kernel objects and variants to execute. If invalid input is not
    // empty for either case, then there were unmatched input items.
    //
    // A message will be emitted later so user can sort it out...
    //

    if !RunParams.getInvalidKernelInput().isEmpty() ||
       !RunParams.getInvalidExcludeKernelInput().isEmpty() {

      RunParams.setInputState(InputOpt.BadInput);

    } else if !RunParams.getInvalidFeatureInput().isEmpty() ||
              !RunParams.getInvalidExcludeFeatureInput().isEmpty() {

      RunParams.setInputState(InputOpt.BadInput);

    } else {  // kernel and feature input looks good

      for kid in run_kern do
        kernels.append(getKernelObject(kid));

      if !RunParams.getInvalidVariantInput().isEmpty() ||
         !RunParams.getInvalidExcludeVariantInput().isEmpty() {

        RunParams.setInputState(InputOpt.BadInput);

      } else {  // variant input lools good

        for vid in run_var do
          variant_ids.append(vid);

        //
        // If we've gotten to this point, we have good input to run.
        //
        if RunParams.getInputState() != InputOpt.DryRun &&
           RunParams.getInputState() != InputOpt.CheckRun then
          RunParams.setInputState(InputOpt.PerfRun);
      }  // kernel and variant input both look good
    }  // if kernel input looks good
  }

  proc reportRunSummary(writer) throws {
    const in_state = RunParams.getInputState();

    if in_state == InputOpt.BadInput {
      writer.write("\nRunParams state:\n");
      writer.write("----------------");
      RunParams.print(writer);

      writer.write("\n\nSuite will not be run now due to bad input.");
      writer.write("\n  See run parameters or option messages above.\n\n");
    } else if in_state == InputOpt.PerfRun ||
              in_state == InputOpt.DryRun ||
              in_state == InputOpt.CheckRun {
      if in_state == InputOpt.DryRun {
        writer.write("\n\nRAJA performance suite dry run summary....");
        writer.write("\n--------------------------------------\n");

        writer.write("\nInput state:");
        writer.write("\n------------");
        RunParams.print(writer);
      }

      if in_state == InputOpt.PerfRun || in_state == InputOpt.CheckRun {
        writer.write("\n\nRAJA performance suite run summary....");
        writer.write("\n--------------------------------------\n");
      }

      var ofiles:string;
      if !RunParams.getOutputDirName().isEmpty() then
        ofiles = RunParams.getOutputDirName();
      else
        ofiles = ".";
      ofiles += "/" + RunParams.getOutputFilePrefix() + "*";

      writer.writeln("\nHow suite will be run:");
      writer.writeln("\t # passes = " + RunParams.getNumPasses():string);
      if RunParams.getSizeMeaning() == SizeMeaning.Factor then
        writer.writeln("\t Kernel size factor = " + RunParams.getSizeFactor():string);
      else if RunParams.getSizeMeaning() == SizeMeaning.Direct then
        writer.writeln("\t Kernel size = " + RunParams.getSize():string);
      writer.writeln("\t Kernel rep factor = " + RunParams.getRepFactor():string);
      writer.writeln("\t Output files will be named " + ofiles);

      writer.writeln("\nThe following kernels and variants (when available for a kernel) will be run:");

      writer.write("\nVariants");
      writer.write("\n--------\n");
      for vid in variant_ids do
        writer.writeln(getVariantName(vid));
      writer.writeln();

      writeKernelInfoSummary(writer, false);
    }

    writer.flush();
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

    if !to_file then channel.writeln("-" * dash_width);

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

  proc runSuite() {
    const in_state = RunParams.getInputState();

    if in_state != InputOpt.PerfRun && in_state != InputOpt.CheckRun then
      return;

    writeln("\n\nRun warmup kernels...");

    var warmup_kernels = (
      new basic.DAXPY(),
      new basic.REDUCE3_INT(),
      new algorithm.SORT(),
    );

    for warmup_kernel in warmup_kernels {
      writeln("Kernel : " + warmup_kernel.getName());
      for vid in VariantID {
        if RunParams.showProgress() {
          if warmup_kernel.hasVariantDefined(vid) then
            write("   Running ");
          else
            write("   No ");
          writeln(getVariantName(vid) + " variant");
        }
        if warmup_kernel.hasVariantDefined(vid) then
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

        for vid in variant_ids {
          if RunParams.showProgress() {
            if kernel.hasVariantDefined(vid) then
              write("   Running ");
            else
              write("   No ");
            writeln(getVariantName(vid) + " variant");
          }
          if kernel.hasVariantDefined(vid) then
            kernel.execute(vid);
        }  // loop over variants
      }  // loop over kernels
    }  // loop over passes through suite
  }

  proc outputRunData() throws {
    const in_state = RunParams.getInputState();

    if in_state != InputOpt.PerfRun && in_state != InputOpt.CheckRun then
      return;

    writeln("\n\nGenerate run report files...");

    //
    // Generate output file prefix (including directory path).
    //
    const outdir = RunParams.getOutputDirName();
    const out_fprefix = "./" + RunParams.getOutputFilePrefix();

    if !outdir.isEmpty() {
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

    var varcol_width = for vid in variant_ids do
                         max(prec+2, getVariantName(vid).size);

    //
    // Print title line.
    //
    channel.write(getReportTitle(mode));

    //
    // Wrtie CSV file contents for report.
    //
    channel.writeln(sepchr * variant_ids.size);

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
          cprintf(channel, "%*.*Lf", varcol_width[iv], prec,
                  getReportDataEntry(mode, kern, vid));
      }
      channel.writeln();
    }

    channel.flush();
  } // note files and channels are closed when their variables go out of scope

  proc writeFOMReport(const ref filename:string) throws {
    var fom_groups = getFOMGroups();
    if fom_groups.isEmpty() then return;

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
    var col_min: [0..<ncols] real = max(real);
    var col_max: [0..<ncols] real = min(real);
    var col_avg: [0..<ncols] real = 0.0;
    var col_stddev: [0..<ncols] real = 0.0;
    var pct_diff: [0..<kernels.size, 0..<ncols] real = 0.0;

    //
    // Print title line.
    //
    channel.write("FOM Report : signed speedup(-)/slowdown(+) for each PM (base vs. RAJA) -> (T_RAJA - T_base) / T_base )");
    channel.writeln(sepchr * 2*ncols);

    channel.write("'OVER_TOL' in column to right if RAJA speedup is over tolerance");
    channel.writeln(sepchr * 2*ncols);

    const pass = ",        ";
    const fail = ",OVER_TOL";

    //
    // Print column title line.
    //
    channel.writef("%-*s", kercol_width, kernel_col_name);
    for group in fom_groups do
      for vid in group.variants do
        channel.writef("%s%-*s%s", sepchr, fom_col_width, getVariantName(vid), pass);
    channel.writeln();

    //
    // Write CSV file contents for FOM report.
    //

    //
    // Print row of FOM data for each kernel.
    //
    for (ik, kern) in zip(0..<kernels.size, kernels) {
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
    for (ik, kern) in zip(0..<kernels.size, kernels) {
      var col = 0;
      for group in fom_groups {
        for comp_vid in group.variants {
          if kern.wasVariantRun(comp_vid) then
            col_stddev[col] += ( pct_diff[ik, col] - col_avg[col] ) *
                               ( pct_diff[ik, col] - col_avg[col] );
          col += 1;
        }  // loop over group variants
      }  // loop over groups
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
    var checksum_width = prec + 8;

    var namecol_width = 0;
    for kern in kernels do
      namecol_width = max(namecol_width, kern.getName().size);
    for vid in variant_ids do
      namecol_width = max(namecol_width, getVariantName(vid).size);
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

  proc getReportTitle(mode:CSVRepMode) {
    var title:string;
    select mode {
      when CSVRepMode.Timing do title = "Mean Runtime Report (sec.) ";
      when CSVRepMode.Speedup do
        if haveReferenceVariant() then
          title = "Speedup Report (T_ref/T_var)" +
                  ": ref var = " + getVariantName(reference_vid) + " ";
      otherwise do writeln("\n Unknown CSV report mode = " + mode:string);
    };
    return title;
  }

  proc getReportDataEntry(mode:CSVRepMode, kern, vid:VariantID): longdouble {
    var retval:longdouble = 0.0;
    select mode {
      when CSVRepMode.Timing do
        retval = kern.getTotTime(vid):longdouble / RunParams.getNumPasses();
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

  // This method needs to be adjusted to match our kernels with RAJAPerf's
  proc getFOMGroups() {
    var fom_groups:list(FOMGroup);

    for iv in variant_ids.indices {
      const vid = variant_ids[iv];
      const vname = getVariantName(vid);

      if vname.find("Base") != -1 {
        var group:FOMGroup;
        group.base = vid;

        const pos = vname.find("_");
        const pm = vname[pos+1..];

        for ivs in iv+1..<variant_ids.size {
          const vids = variant_ids[ivs];
          if getVariantName(vids).find(pm) != -1 then
            group.variants.append(vids);
        }

        if !group.variants.isEmpty() then
          fom_groups.append(group);
      }  // if variant name contains 'Base'
    }  // iterate over variant ids to run

    if false {  //  RDH DEBUG   (leave this here, it's useful for debugging!)
      writeln("\nFOMGroups...");
      for group in fom_groups {
        writeln("\tBase : " + getVariantName(group.base));
        for vid in group.variants do
          writeln("\t\t " + getVariantName(vid));
      }
    }

    return fom_groups;
  }

  /* Construct and return *unmanaged* kernel object for given KernelID enum
     value. */
  proc getKernelObject(kid) {
    select kid {

      //
      // Basic kernels...
      //
      when KernelID.Basic_DAXPY               do return new unmanaged basic.DAXPY():KernelBase;
      when KernelID.Basic_IF_QUAD             do return new unmanaged basic.IF_QUAD():KernelBase;
      when KernelID.Basic_INIT3               do return new unmanaged basic.INIT3():KernelBase;
      when KernelID.Basic_INIT_VIEW1D         do return new unmanaged basic.INIT_VIEW1D():KernelBase;
      when KernelID.Basic_INIT_VIEW1D_OFFSET  do return new unmanaged basic.INIT_VIEW1D_OFFSET():KernelBase;
      when KernelID.Basic_MAT_MAT_SHARED      do return new unmanaged basic.MAT_MAT_SHARED():KernelBase;
      when KernelID.Basic_MULADDSUB           do return new unmanaged basic.MULADDSUB():KernelBase;
      when KernelID.Basic_NESTED_INIT         do return new unmanaged basic.NESTED_INIT():KernelBase;
      when KernelID.Basic_PI_ATOMIC           do return new unmanaged basic.PI_ATOMIC():KernelBase;
      when KernelID.Basic_PI_REDUCE           do return new unmanaged basic.PI_REDUCE():KernelBase;
      when KernelID.Basic_REDUCE3_INT         do return new unmanaged basic.REDUCE3_INT():KernelBase;
      when KernelID.Basic_TRAP_INT            do return new unmanaged basic.TRAP_INT():KernelBase;

      //
      // Lcals kernels...
      //
      when KernelID.Lcals_DIFF_PREDICT        do return new unmanaged lcals.DIFF_PREDICT():KernelBase;
      when KernelID.Lcals_EOS                 do return new unmanaged lcals.EOS():KernelBase;
      when KernelID.Lcals_FIRST_DIFF          do return new unmanaged lcals.FIRST_DIFF():KernelBase;
      when KernelID.Lcals_FIRST_MIN           do return new unmanaged lcals.FIRST_MIN():KernelBase;
      when KernelID.Lcals_FIRST_SUM           do return new unmanaged lcals.FIRST_SUM():KernelBase;
      when KernelID.Lcals_GEN_LIN_RECUR       do return new unmanaged lcals.GEN_LIN_RECUR():KernelBase;
      when KernelID.Lcals_HYDRO_1D            do return new unmanaged lcals.HYDRO_1D():KernelBase;
      when KernelID.Lcals_HYDRO_2D            do return new unmanaged lcals.HYDRO_2D():KernelBase;
      when KernelID.Lcals_INT_PREDICT         do return new unmanaged lcals.INT_PREDICT():KernelBase;
      when KernelID.Lcals_PLANCKIAN           do return new unmanaged lcals.PLANCKIAN():KernelBase;
      when KernelID.Lcals_TRIDIAG_ELIM        do return new unmanaged lcals.TRIDIAG_ELIM():KernelBase;

      //
      // Polybench kernels...
      //
      //when KernelID.Polybench_2MM do return new unmanaged polybench.2MM():KernelBase;
      //when KernelID.Polybench_3MM do return new unmanaged polybench.3MM():KernelBase;
      //when KernelID.Polybench_ADI do return new unmanaged polybench.ADI():KernelBase;
      //when KernelID.Polybench_ATAX do return new unmanaged polybench.ATAX():KernelBase;
      //when KernelID.Polybench_FDTD_2D do return new unmanaged polybench.FDTD_2D():KernelBase;
      //when KernelID.Polybench_FLOYD_WARSHALL do return new unmanaged polybench.FLOYD_WARSHALL():KernelBase;
      //when KernelID.Polybench_GEMM do return new unmanaged polybench.GEMM():KernelBase;
      //when KernelID.Polybench_GEMVER do return new unmanaged polybench.GEMVER():KernelBase;
      //when KernelID.Polybench_GESUMMV do return new unmanaged polybench.GESUMMV():KernelBase;
      //when KernelID.Polybench_HEAT_3D do return new unmanaged polybench.HEAT_3D():KernelBase;
      //when KernelID.Polybench_JACOBI_1D do return new unmanaged polybench.JACOBI_1D():KernelBase;
      //when KernelID.Polybench_JACOBI_2D do return new unmanaged polybench.JACOBI_2D():KernelBase;
      //when KernelID.Polybench_MVT do return new unmanaged polybench.MVT():KernelBase;

      //
      // Stream kernels...
      //
      when KernelID.Stream_ADD                do return new unmanaged stream.ADD():KernelBase;
      when KernelID.Stream_COPY               do return new unmanaged stream.COPY():KernelBase;
      when KernelID.Stream_DOT                do return new unmanaged stream.DOT():KernelBase;
      when KernelID.Stream_MUL                do return new unmanaged stream.MUL():KernelBase;
      when KernelID.Stream_TRIAD              do return new unmanaged stream.TRIAD():KernelBase;

      //
      // Apps kernels...
      //
      //when KernelID.Apps_COUPLE do return new unmanaged apps.COUPLE():KernelBase;
      when KernelID.Apps_DEL_DOT_VEC_2D       do return new unmanaged apps.DEL_DOT_VEC_2D():KernelBase;
      //when KernelID.Apps_DIFFUSION3DPA do return new unmanaged apps.DIFFUSION3DPA():KernelBase;
      //when KernelID.Apps_ENERGY do return new unmanaged apps.ENERGY():KernelBase;
      //when KernelID.Apps_FIR do return new unmanaged apps.FIR():KernelBase;
      //when KernelID.Apps_HALOEXCHANGE do return new unmanaged apps.HALOEXCHANGE():KernelBase;
      //when KernelID.Apps_HALOEXCHANGE_FUSED do return new unmanaged apps.HALOEXCHANGE_FUSED():KernelBase;
      //when KernelID.Apps_LTIMES do return new unmanaged apps.LTIMES():KernelBase;
      //when KernelID.Apps_LTIMES_NOVIEW do return new unmanaged apps.LTIMES_NOVIEW():KernelBase;
      //when KernelID.Apps_MASS3DPA do return new unmanaged apps.MASS3DPA():KernelBase;
      //when KernelID.Apps_PRESSURE do return new unmanaged apps.PRESSURE():KernelBase;
      //when KernelID.Apps_VOL3D do return new unmanaged apps.VOL3D():KernelBase;

      //
      // Algorithm kernels...
      //
      when KernelID.Algorithm_SORT            do return new unmanaged algorithm.SORT():KernelBase;
      when KernelID.Algorithm_SORTPAIRS       do return new unmanaged algorithm.SORTPAIRS():KernelBase;

      otherwise halt("\n Unknown Kernel ID = " + getFullKernelName(kid));
    }  // end switch on kernel id
  }
}
