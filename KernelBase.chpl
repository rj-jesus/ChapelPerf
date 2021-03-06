//----------------------------------------------------------------------------
// Copyright (C) 2022 Ricardo Jesus, EPCC, United Kingdom.
// All rights reserved.
//
// Redistribution and use of this software, with or without modification, is
// permitted provided that the following conditions are met:
//
// 1. Redistributions of this software must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
// THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
// EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//----------------------------------------------------------------------------

module KernelBase {
  private use Barriers;
  private use Set;
  private use Time;

  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private import RunParams;

  const barrier = new Barrier(numLocales);

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

    var running_variant_int:int;
    proc running_variant return try! running_variant_int:VariantID;

    var num_exec: [VariantID.first..VariantID.last] int = 0;

    // Checksums
    var checksum: [VariantID.first..VariantID.last] Checksum_type = 0:Checksum_type;
    var checksum_scale_factor: Checksum_type;

    // Elapsed time in seconds
    var timer: Timer;
    var min_time: [VariantID.first..VariantID.last] Elapsed_type = 0;
    var max_time: [VariantID.first..VariantID.last] Elapsed_type = 0;
    var tot_time: [VariantID.first..VariantID.last] Elapsed_type = 0;

    proc init(kernel_id: KernelID) {
      this.kernel_id = kernel_id;
    }

    proc usesFeature(fid: FeatureID) return uses_feature.contains(fid);

    proc startTimer() {
      if numLocales > 1 then barrier.barrier();
      timer.start();
    }

    proc stopTimer()  {
      if numLocales > 1 then barrier.barrier();
      timer.stop();
      recordExecTime();
    }

    proc resetTimer() { timer.clear(); }

    proc recordExecTime() {
      num_exec[running_variant] += 1;

      const exec_time = timer.elapsed():Elapsed_type;
      min_time[running_variant] = min(min_time[running_variant], exec_time);
      max_time[running_variant] = max(max_time[running_variant], exec_time);
      tot_time[running_variant] += exec_time;
    }

    proc readTimer(unit:TimeUnits=TimeUnits.seconds) { return timer.elapsed(unit); }

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

    //
    // Methods to get information about kernel execution for reports containing
    // kernel execution information
    //
    proc wasVariantRun(vid:VariantID)              { return num_exec[vid] > 0; }

    proc getMinTime(vid:VariantID): Elapsed_type   { return min_time[vid]; }
    proc getMaxTime(vid:VariantID): Elapsed_type   { return max_time[vid]; }
    proc getTotTime(vid:VariantID): Elapsed_type   { return tot_time[vid]; }
    proc getChecksum(vid:VariantID): Checksum_type { return checksum[vid]; }

    proc execute(vid:VariantID) {
      running_variant_int = vid:int;

      resetTimer();
      resetDataInitCount();

      runVariant(vid);

      running_variant_int = VariantID.size;
    }

    proc runVariant(vid:VariantID) { halt("Error: Called base method!"); }
  };
}
