module lcals {
  use Reflection;

  private use DataTypes;
  private use DataUtils;
  private use KernelBase;
  private use Utils;
  private use Enums;

  class DIFF_PREDICT : KernelBase {
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Lcals_DIFF_PREDICT);

      setDefaultProblemSize(1000000);
      setDefaultReps(200);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());

      setKernelsPerRep(1);
      setBytesPerRep((10*sizeof(Real_type) + 10*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(9 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      m_array_length = getActualProblemSize() * 14;
      const offset = getActualProblemSize();

      var px = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var cx = allocAndInitData(Real_type, m_array_length, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              var ar, br, cr: Real_type;

              ar                  =      cx[i + offset *  4];
              br                  = ar - px[i + offset *  4];
              px[i + offset *  4] = ar;
              cr                  = br - px[i + offset *  5];
              px[i + offset *  5] = br;
              ar                  = cr - px[i + offset *  6];
              px[i + offset *  6] = cr;
              br                  = ar - px[i + offset *  7];
              px[i + offset *  7] = ar;
              cr                  = br - px[i + offset *  8];
              px[i + offset *  8] = br;
              ar                  = cr - px[i + offset *  9];
              px[i + offset *  9] = cr;
              br                  = ar - px[i + offset * 10];
              px[i + offset * 10] = ar;
              cr                  = br - px[i + offset * 11];
              px[i + offset * 11] = br;
              px[i + offset * 13] = cr - px[i + offset * 12];
              px[i + offset * 12] = cr;
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend {
              var ar, br, cr: Real_type;

              ar                  =      cx[i + offset *  4];
              br                  = ar - px[i + offset *  4];
              px[i + offset *  4] = ar;
              cr                  = br - px[i + offset *  5];
              px[i + offset *  5] = br;
              ar                  = cr - px[i + offset *  6];
              px[i + offset *  6] = cr;
              br                  = ar - px[i + offset *  7];
              px[i + offset *  7] = ar;
              cr                  = br - px[i + offset *  8];
              px[i + offset *  8] = br;
              ar                  = cr - px[i + offset *  9];
              px[i + offset *  9] = cr;
              br                  = ar - px[i + offset * 10];
              px[i + offset * 10] = ar;
              cr                  = br - px[i + offset * 11];
              px[i + offset * 11] = br;
              px[i + offset * 13] = cr - px[i + offset * 12];
              px[i + offset * 12] = cr;
            }
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(px, m_array_length);
    }
  }

  class EOS : KernelBase {
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Lcals_EOS);

      setDefaultProblemSize(1000000);
      setDefaultReps(500);

      setActualProblemSize(getTargetProblemSize());

      m_array_length = getActualProblemSize() + 7;

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize() +
                      (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_array_length );
      setFLOPsPerRep(16 * getActualProblemSize());

      checksum_scale_factor = 0.0001:Checksum_type * (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    inline proc EOS_BODY(x, y, z, u, q, r, t, i) {
      x[i] = u[i] + r*( z[i] + r*y[i] ) +
                    t*( u[i+3] + r*( u[i+2] + r*u[i+1] ) +
                                 t*( u[i+6] + q*( u[i+5] + q*u[i+4] ) ) );
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var y = allocAndInitData(Real_type, m_array_length, vid);
      var z = allocAndInitData(Real_type, m_array_length, vid);
      var u = allocAndInitData(Real_type, m_array_length, vid);

      var q = initData(Real_type, vid);
      var r = initData(Real_type, vid);
      var t = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              EOS_BODY(x, y, z, u, q, r, t, i);
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              EOS_BODY(x, y, z, u, q, r, t, i);
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            EOS_BODY(x, y, z, u, q, r, t, I);
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize(), checksum_scale_factor:Real_type);
    }
  }

  class FIRST_DIFF: KernelBase {
    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Lcals_FIRST_DIFF);

      setDefaultProblemSize(1000000);
      setDefaultReps(2000);

      setActualProblemSize(getTargetProblemSize());

      m_N = getActualProblemSize()+1;

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type) + 0*sizeof(Real_type)) * getActualProblemSize() +
                      (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_N );
      setFLOPsPerRep(1 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      var y = allocAndInitData(Real_type, m_N, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              x[i] = y[i+1] - y[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              x[i] = y[i+1] - y[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            x[I] = y[I+1]-y[I];
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize());
    }
  }

  class FIRST_MIN: KernelBase {
    var m_N: Index_type;
    var m_minloc = -1;

    proc init() {
      super.init(KernelID.Lcals_FIRST_MIN);

      setDefaultProblemSize(1000000);
      setDefaultReps(100);

      setActualProblemSize(getTargetProblemSize());

      m_N = getActualProblemSize();

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type ) + 1*sizeof(Real_type )) +
                      (1*sizeof(Index_type) + 1*sizeof(Index_type)) +
                      (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_N );
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.Reduction);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Reduction_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      x[m_N/2] = -1.0e+10;
      const xmin_init = x[0];
      const initloc = 0;
      m_minloc = -1;

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {
        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            var mymin = (xmin_init, initloc);

            for i in ibegin..<iend do
              if x[i] < mymin(0) then
                mymin = (x[i], i);

            m_minloc = max(m_minloc, mymin(1));
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            var mymin = (xmin_init, initloc);

            forall i in ibegin..<iend with (minloc reduce mymin) do
              // if I comment the `if' it becomes ~10x slower
              if x[i] < mymin(0) then
                mymin reduce= (x[i], i);

            m_minloc = max(m_minloc, mymin(1));
          }

          stopTimer();
        }

        when VariantID.Reduction_Chpl {
          startTimer();

          for 0..#run_reps {
            var (myminval, myminloc) = minloc reduce zip(x, x.domain);
            m_minloc = max(m_minloc, myminloc);
          }

          stopTimer();
        }
      }

      // update checksum
      checksum[vid] += m_minloc:Checksum_type;
    }
  }

  class FIRST_SUM: KernelBase {
    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Lcals_FIRST_SUM);

      setDefaultProblemSize(1000000);
      setDefaultReps(2000);

      setActualProblemSize(getTargetProblemSize());

      m_N = getActualProblemSize();

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 0*sizeof(Real_type)) * (m_N-1) +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_N);
      setFLOPsPerRep(1 * (getActualProblemSize()-1));

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      var y = allocAndInitData(Real_type, m_N, vid);

      const run_reps = getRunReps();
      const ibegin = 1;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              x[i] = y[i-1] + y[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              x[i] = y[i-1] + y[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            x[I] = y[I-1] + y[I];
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize());
    }
  }

  class GEN_LIN_RECUR: KernelBase {
    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Lcals_GEN_LIN_RECUR);

      setDefaultProblemSize(1000000);
      setDefaultReps(500);

      setActualProblemSize(getTargetProblemSize());

      m_N = getActualProblemSize();

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(2);
      setBytesPerRep((2*sizeof(Real_type) + 3*sizeof(Real_type)) * m_N +
                     (2*sizeof(Real_type) + 3*sizeof(Real_type)) * m_N);
      setFLOPsPerRep((3 + 3) * getActualProblemSize());

      checksum_scale_factor = 0.01:Checksum_type *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var N = m_N;
      var kb5i = 0:Index_type;

      var b5 = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      var stb5 = allocAndInitData(Real_type, m_N, vid);
      var sa = allocAndInitData(Real_type, m_N, vid);
      var sb = allocAndInitData(Real_type, m_N, vid);

      const run_reps = getRunReps();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for k in 0..<N {
              b5[k+kb5i] = sa[k] + stb5[k]*sb[k];
              stb5[k] = b5[k+kb5i] - stb5[k];
            }

            for k in 0..<N by -1 {
              b5[k+kb5i] = sa[k] + stb5[k]*sb[k];
              stb5[k] = b5[k+kb5i] - stb5[k];
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall k in 0..<N {
              b5[k+kb5i] = sa[k] + stb5[k]*sb[k];
              stb5[k] = b5[k+kb5i] - stb5[k];
            }

            forall k in 0..<N by -1 {
              b5[k+kb5i] = sa[k] + stb5[k]*sb[k];
              stb5[k] = b5[k+kb5i] - stb5[k];
            }
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(b5, getActualProblemSize(), checksum_scale_factor:Real_type);
    }
  }

  class HYDRO_1D: KernelBase {
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Lcals_HYDRO_1D);

      setDefaultProblemSize(1000000);
      setDefaultReps(1000);

      setActualProblemSize(getTargetProblemSize());

      m_array_length = getActualProblemSize() + 12;

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) * getActualProblemSize() +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * (getActualProblemSize()+1));
      setFLOPsPerRep(5 * getActualProblemSize());

      checksum_scale_factor = 0.001:Checksum_type *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var y = allocAndInitData(Real_type, m_array_length, vid);
      var z = allocAndInitData(Real_type, m_array_length, vid);

      const q = initData(Real_type, vid);
      const r = initData(Real_type, vid);
      const t = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              x[i] = q + y[i]*(r*z[i+10] + t*z[i+11]);
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend {
              x[i] = q + y[i]*(r*z[i+10] + t*z[i+11]);
            }
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            x[I] = q + y[I]*(r*z[I+10] + t*z[I+11]);
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(x, getActualProblemSize(), checksum_scale_factor:Real_type);
    }
  }

  class HYDRO_2D: KernelBase {
    var m_jn: Index_type;
    var m_kn: Index_type;

    var m_s: Real_type;
    var m_t: Real_type;

    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Lcals_HYDRO_2D);

      m_jn = 1000;
      m_kn = 1000;

      m_s = 0.0041;
      m_t = 0.0037;

      setDefaultProblemSize(m_kn * m_jn);
      setDefaultReps(100);

      m_jn = sqrt(getTargetProblemSize()):Index_type;
      m_kn = m_jn;
      m_array_length = m_kn * m_jn;

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(3 * getActualProblemSize());
      setKernelsPerRep(3);
      setBytesPerRep((2*sizeof(Real_type) + 0*sizeof(Real_type)) * (m_kn-2) * (m_jn-2) +
                     (0*sizeof(Real_type) + 4*sizeof(Real_type)) * m_array_length +
                     (2*sizeof(Real_type) + 0*sizeof(Real_type)) * (m_kn-2) * (m_jn-2) +
                     (0*sizeof(Real_type) + 4*sizeof(Real_type)) * m_array_length +
                     (2*sizeof(Real_type) + 4*sizeof(Real_type)) * (m_kn-2) * (m_jn-2));
      setFLOPsPerRep((14 + 26 + 4) * (m_jn-2)*(m_kn-2));

      checksum_scale_factor = 0.001:Checksum_type *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var zroutdat = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var zzoutdat = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var zadat = allocAndInitData(Real_type, m_array_length, vid);
      var zbdat = allocAndInitData(Real_type, m_array_length, vid);
      var zmdat = allocAndInitData(Real_type, m_array_length, vid);
      var zpdat = allocAndInitData(Real_type, m_array_length, vid);
      var zqdat = allocAndInitData(Real_type, m_array_length, vid);
      var zrdat = allocAndInitData(Real_type, m_array_length, vid);
      var zudat = allocAndInitData(Real_type, m_array_length, vid);
      var zvdat = allocAndInitData(Real_type, m_array_length, vid);
      var zzdat = allocAndInitData(Real_type, m_array_length, vid);

      const run_reps = getRunReps();
      const kbeg = 1;
      const kend = m_kn - 1;
      const jbeg = 1;
      const jend = m_jn - 1;

      const s: Real_type = m_s;
      const t: Real_type = m_t;

      const kn: Index_type = m_kn;
      const jn: Index_type = m_jn;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for k in kbeg..<kend {
              for j in jbeg..<jend {
                zadat[j+k*jn] = ( zpdat[j-1+(k+1)*jn] + zqdat[j-1+(k+1)*jn] -
                                  zpdat[j-1+k*jn] - zqdat[j-1+k*jn] ) *
                                ( zrdat[j+k*jn] + zrdat[j-1+k*jn] ) /
                                ( zmdat[j-1+k*jn] + zmdat[j-1+(k+1)*jn] );
                zbdat[j+k*jn] = ( zpdat[j-1+k*jn] + zqdat[j-1+k*jn] -
                                  zpdat[j+k*jn] - zqdat[j+k*jn] ) *
                                ( zrdat[j+k*jn] + zrdat[j+(k-1)*jn] ) /
                                ( zmdat[j+k*jn] + zmdat[j-1+k*jn] );
              }
            }

            for k in kbeg..<kend {
              for j in jbeg..<jend {
                zudat[j+k*jn] += s*( zadat[j+k*jn] * ( zzdat[j+k*jn] - zzdat[j+1+k*jn] ) -
                                  zadat[j-1+k*jn] * ( zzdat[j+k*jn] - zzdat[j-1+k*jn] ) -
                                  zbdat[j+k*jn] * ( zzdat[j+k*jn] - zzdat[j+(k-1)*jn] ) +
                                  zbdat[j+(k+1)*jn] * ( zzdat[j+k*jn] - zzdat[j+(k+1)*jn] ) );
                zvdat[j+k*jn] += s*( zadat[j+k*jn] * ( zrdat[j+k*jn] - zrdat[j+1+k*jn] ) -
                                  zadat[j-1+k*jn] * ( zrdat[j+k*jn] - zrdat[j-1+k*jn] ) -
                                  zbdat[j+k*jn] * ( zrdat[j+k*jn] - zrdat[j+(k-1)*jn] ) +
                                  zbdat[j+(k+1)*jn] * ( zrdat[j+k*jn] - zrdat[j+(k+1)*jn] ) );
              }
            }

            for k in kbeg..<kend {
              for j in jbeg..<jend {
                zroutdat[j+k*jn] = zrdat[j+k*jn] + t*zudat[j+k*jn];
                zzoutdat[j+k*jn] = zzdat[j+k*jn] + t*zvdat[j+k*jn];
              }
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            forall k in kbeg..<kend do
              for j in jbeg..<jend {
                zadat[j+k*jn] = ( zpdat[j-1+(k+1)*jn] + zqdat[j-1+(k+1)*jn] -
                                  zpdat[j-1+k*jn] - zqdat[j-1+k*jn] ) *
                                ( zrdat[j+k*jn] + zrdat[j-1+k*jn] ) /
                                ( zmdat[j-1+k*jn] + zmdat[j-1+(k+1)*jn] );
                zbdat[j+k*jn] = ( zpdat[j-1+k*jn] + zqdat[j-1+k*jn] -
                                  zpdat[j+k*jn] - zqdat[j+k*jn] ) *
                                ( zrdat[j+k*jn] + zrdat[j+(k-1)*jn] ) /
                                ( zmdat[j+k*jn] + zmdat[j-1+k*jn] );
              }

            forall k in kbeg..<kend do
              for j in jbeg..<jend {
                zudat[j+k*jn] += s*( zadat[j+k*jn] * ( zzdat[j+k*jn] - zzdat[j+1+k*jn] ) -
                                  zadat[j-1+k*jn] * ( zzdat[j+k*jn] - zzdat[j-1+k*jn] ) -
                                  zbdat[j+k*jn] * ( zzdat[j+k*jn] - zzdat[j+(k-1)*jn] ) +
                                  zbdat[j+(k+1)*jn] * ( zzdat[j+k*jn] - zzdat[j+(k+1)*jn] ) );
                zvdat[j+k*jn] += s*( zadat[j+k*jn] * ( zrdat[j+k*jn] - zrdat[j+1+k*jn] ) -
                                  zadat[j-1+k*jn] * ( zrdat[j+k*jn] - zrdat[j-1+k*jn] ) -
                                  zbdat[j+k*jn] * ( zrdat[j+k*jn] - zrdat[j+(k-1)*jn] ) +
                                  zbdat[j+(k+1)*jn] * ( zrdat[j+k*jn] - zrdat[j+(k+1)*jn] ) );
              }

            forall k in kbeg..<kend do
              for j in jbeg..<jend {
                zroutdat[j+k*jn] = zrdat[j+k*jn] + t*zudat[j+k*jn];
                zzoutdat[j+k*jn] = zzdat[j+k*jn] + t*zvdat[j+k*jn];
              }

          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(zzoutdat, m_array_length, checksum_scale_factor:Real_type);
      checksum[vid] += calcChecksum(zroutdat, m_array_length, checksum_scale_factor:Real_type);
    }
  }

  class INT_PREDICT: KernelBase {
    var m_array_length: Index_type;
    var m_offset: Index_type;

    var m_px_initval: Real_type;

    proc init() {
      super.init(KernelID.Lcals_INT_PREDICT);

      setDefaultProblemSize(1000000);
      setDefaultReps(400);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 10*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(17 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      m_array_length = getActualProblemSize() * 13;
      const offset = getActualProblemSize();

      m_px_initval = 1.0;
      var px = allocAndInitDataConst(Real_type, m_array_length, m_px_initval, vid);

      var dm22 = initData(Real_type, vid);
      var dm23 = initData(Real_type, vid);
      var dm24 = initData(Real_type, vid);
      var dm25 = initData(Real_type, vid);
      var dm26 = initData(Real_type, vid);
      var dm27 = initData(Real_type, vid);
      var dm28 = initData(Real_type, vid);
      var c0   = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              px[i] = dm28*px[i + offset * 12] + dm27*px[i + offset * 11] +
                      dm26*px[i + offset * 10] + dm25*px[i + offset *  9] +
                      dm24*px[i + offset *  8] + dm23*px[i + offset *  7] +
                      dm22*px[i + offset *  6] +
                      c0*( px[i + offset *  4] + px[i + offset *  5] ) +
                      px[i + offset *  2];
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              px[i] = dm28*px[i + offset * 12] + dm27*px[i + offset * 11] +
                      dm26*px[i + offset * 10] + dm25*px[i + offset *  9] +
                      dm24*px[i + offset *  8] + dm23*px[i + offset *  7] +
                      dm22*px[i + offset *  6] +
                      c0*( px[i + offset *  4] + px[i + offset *  5] ) +
                      px[i + offset *  2];
          }

          stopTimer();
        }

      }

      // update checksum
      px -= m_px_initval;
      checksum[vid] += calcChecksum(px, getActualProblemSize());
    }
  }

  class PLANCKIAN: KernelBase {
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Lcals_PLANCKIAN);

      setDefaultProblemSize(1000000);
      setDefaultReps(50);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((2*sizeof(Real_type) + 3*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(4 * getActualProblemSize()); // 1 exp

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var y = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var u = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var v = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var w = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              y[i] = u[i] / v[i];
              w[i] = x[i] / (exp(y[i]) - 1.0);
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend {
              y[i] = u[i] / v[i];
              w[i] = x[i] / (exp(y[i]) - 1.0);
            }
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(w, getActualProblemSize());
    }
  }

  class TRIDIAG_ELIM: KernelBase {
    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Lcals_TRIDIAG_ELIM);

      setDefaultProblemSize(1000000);
      setDefaultReps(1000);

      setActualProblemSize(getTargetProblemSize());

      m_N = getActualProblemSize();

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 3*sizeof(Real_type)) * (m_N-1));
      setFLOPsPerRep(2 * (getActualProblemSize()-1));

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var xout = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      var xin = allocAndInitData(Real_type, m_N, vid);
      var y = allocAndInitData(Real_type, m_N, vid);
      var z = allocAndInitData(Real_type, m_N, vid);

      const run_reps = getRunReps();
      const ibegin = 1;
      const iend = m_N;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              xout[i] = z[i] * (y[i] - xin[i-1]);
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              xout[i] = z[i] * (y[i] - xin[i-1]);
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            xout[I] = z[I] * (y[I] - xin[I-1]);
          }

          stopTimer();
        }

      }

      // update checksum
      checksum[vid] += calcChecksum(xout, getActualProblemSize());
    }
  }
}
