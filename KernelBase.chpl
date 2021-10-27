module KernelBase {
  private use Set;
  private use Time;

  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private import RunParams;

  class KernelBase {
    //
    // Static properties of kernel, independent of run
    //
    var kernel_id: KernelID;
    var name = kernel_id:string;

    var default_prob_size: Index_type;
    var default_reps: Index_type;

    var actual_prob_size: Index_type;

    var uses_feature: set(FeatureID);
    var has_variant_defined: set(VariantID);

    //
    // Properties of kernel dependent on how kernel is run
    //
    var its_per_rep: Index_type;
    var kernels_per_rep: Index_type;
    var bytes_per_rep: Index_type;
    var flops_per_rep: Index_type;

    var running_variant:int;

    var num_exec: [0..<VariantID.size] int = 0;

    // Checksums
    var checksum: [0..<VariantID.size] Checksum_type = 0:Checksum_type;
    var checksum_scale_factor: Checksum_type;

    // Elapsed time in seconds
    var timer: Timer;
    var min_time: [0..<VariantID.size] Elapsed_type = 0;
    var max_time: [0..<VariantID.size] Elapsed_type = 0;
    var tot_time: [0..<VariantID.size] Elapsed_type = 0;

    proc init(kernel_id: KernelID) {
      this.kernel_id = kernel_id;
    }

    proc usesFeature(fid: FeatureID) { return uses_feature.contains(fid); };

    proc startTimer() { timer.start(); }

    proc stopTimer() {
      timer.stop();

      recordExecTime();
    }

    proc recordExecTime() {
      num_exec[running_variant] += 1;

      const exec_time = timer.elapsed():Elapsed_type;
      min_time[running_variant] = min(min_time[running_variant], exec_time);
      max_time[running_variant] = max(max_time[running_variant], exec_time);
      tot_time[running_variant] += exec_time;
    }

    proc resetTimer() { timer.clear(); }

    proc readTimer(unit:TimeUnits=TimeUnits.seconds) { return timer.elapsed(unit); }

    proc log(checksum) {
      //printf("%s: done in %f seconds (%f)\n", name.c_str(), readTimer(), checksum);
      writef("%s: done in %dr seconds (%dr)\n", kernel_id, readTimer(), checksum:real);
    }

    proc getKernelID()                          { return kernel_id; }
    proc getName() const ref                    { return name; }

    //
    // Methods called in kernel subclass constructors to set kernel properties
    // used to describe kernel and define how it will run
    //

    proc setDefaultProblemSize(size:Index_type) { default_prob_size = size; }
    proc setActualProblemSize(size:Index_type)  { actual_prob_size = size; }
    proc setDefaultReps(reps:Index_type)        { default_reps = reps; }
    proc setItsPerRep(its:Index_type)           { its_per_rep = its; };
    proc setKernelsPerRep(nkerns:Index_type)    { kernels_per_rep = nkerns; };
    proc setBytesPerRep(bytes_:Index_type)      { bytes_per_rep = bytes_;}
    proc setFLOPsPerRep(flops:Index_type)       { flops_per_rep = flops; }

    proc getTargetProblemSize(): Index_type
    {
      if RunParams.getSizeMeaning() == SizeMeaning.Factor then
        return (default_prob_size*RunParams.getSizeFactor()):Index_type;
      else if RunParams.getSizeMeaning() == SizeMeaning.Direct then
        return RunParams.getSize():Index_type;
      return 0;
    }

    proc getRunReps(): Index_type
    {
      if RunParams.getInputState() == InputOpt.CheckRun then
        return RunParams.getCheckRunReps():Index_type;
      return (default_reps*RunParams.getRepFactor()):Index_type;
    }

    proc setUsesFeature(fid:FeatureID)          { uses_feature.add(fid); }

    proc setVariantDefined(vid:VariantID)       { has_variant_defined.add(vid); }
    proc hasVariantDefined(vid:VariantID)       { return has_variant_defined.contains(vid); }

    //
    // Getter methods used to generate kernel execution summary
    // and kernel details report ouput.
    //

    proc getDefaultProblemSize(): Index_type    { return default_prob_size; }
    proc getActualProblemSize(): Index_type     { return actual_prob_size; }
    proc getDefaultReps(): Index_type           { return default_reps; }
    proc getItsPerRep(): Index_type             { return its_per_rep; }
    proc getKernelsPerRep(): Index_type         { return kernels_per_rep; }
    proc getBytesPerRep(): Index_type           { return bytes_per_rep; }
    proc getFLOPsPerRep(): Index_type           { return flops_per_rep; }

    proc getVariants() ref { return has_variant_defined; }

    //
    // Methods to get information about kernel execution for reports containing
    // kernel execution information
    //
    proc wasVariantRun(vid:VariantID)              { return num_exec[vid:int] > 0; }

    proc getMinTime(vid:VariantID): Elapsed_type   { return min_time[vid:int]; }
    proc getMaxTime(vid:VariantID): Elapsed_type   { return max_time[vid:int]; }
    proc getTotTime(vid:VariantID): Elapsed_type   { return tot_time[vid:int]; }
    proc getChecksum(vid:VariantID): Checksum_type { return checksum[vid:int]; }

    proc execute(vid:VariantID) {
      running_variant = vid:int;

      resetTimer();
      resetDataInitCount();

      run(vid);

      running_variant = VariantID.size;
    }

    proc run(vid:VariantID) { halt("Error: Called base method!"); }
  };
}
