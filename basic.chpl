module basic {
  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private use KernelBase;
  private use Utils;

  class DAXPY: KernelBase {

    proc init() {
      super.init(KernelID.Basic_DAXPY);

      setDefaultProblemSize(1000000);
      setDefaultReps(500);

      setActualProblemSize( getTargetProblemSize() );

      setItsPerRep( getActualProblemSize() );
      setKernelsPerRep(1);
      setBytesPerRep( (1*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize() );
      setFLOPsPerRep(2 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var y = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var x = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var a = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend do
              y[i] += a * x[i];
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend do
              y[i] += a * x[i];
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps do
            y += a * x;

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(y, getActualProblemSize());
    }
  }

  class IF_QUAD: KernelBase {

    proc init() {
      super.init(KernelID.Basic_IF_QUAD);

      setDefaultProblemSize(1000000);
      setDefaultReps(180);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((2*sizeof(Real_type) + 3*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(11 * getActualProblemSize());  // 1 sqrt

      checksum_scale_factor = 0.0001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var a = allocAndInitDataRandSign(Real_type, getActualProblemSize(), vid);
      var b = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var c = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var x1 = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var x2 = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              var s: Real_type = b[i]*b[i] - 4.0*a[i]*c[i];
              if s >= 0 {
                s = sqrt(s);
                x2[i] = (-b[i]+s)/(2.0*a[i]);
                x1[i] = (-b[i]-s)/(2.0*a[i]);
              } else {
                x2[i] = 0.0;
                x1[i] = 0.0;
              }
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend {
              var s: Real_type = b[i]*b[i] - 4.0*a[i]*c[i];
              if s >= 0 {
                s = sqrt(s);
                x2[i] = (-b[i]+s)/(2.0*a[i]);
                x1[i] = (-b[i]-s)/(2.0*a[i]);
              } else {
                x2[i] = 0.0;
                x1[i] = 0.0;
              }
            }
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(x1, getActualProblemSize(), checksum_scale_factor:Real_type);
      checksum[vid] += calcChecksum(x2, getActualProblemSize(), checksum_scale_factor:Real_type);
    }
  }

  class INIT_VIEW1D_OFFSET: KernelBase {

    proc init() {
      super.init(KernelID.Basic_INIT_VIEW1D_OFFSET);

      setDefaultProblemSize(1000000);
      setDefaultReps(2500);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 0*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(1 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.View);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var a = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      const v: Real_type = 0.00000123;

      const run_reps = getRunReps();
      const ibegin = 1;
      const iend = getActualProblemSize()+1;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              a[i-ibegin] = i * v;
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend {
              a[i-ibegin] = i * v;
            }
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            a[I-ibegin] = I * v;
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(a, getActualProblemSize());
    }
  }

  class INIT_VIEW1D: KernelBase {

    proc init() {
      super.init(KernelID.Basic_INIT_VIEW1D);

      setDefaultProblemSize(1000000);
      setDefaultReps(2500);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 0*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(1 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.View);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var a = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      const v: Real_type = 0.00000123;

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              a[i] = (i+1) * v;
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend {
              a[i] = (i+1) * v;
            }
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            a[I] = (I+1) * v;
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(a, getActualProblemSize());
    }
  }

  class INIT3: KernelBase {

    proc init() {
      super.init(KernelID.Basic_INIT3);

      setDefaultProblemSize(1000000);
      setDefaultReps(500);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((3*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep(1 * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var out1 = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var out2 = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var out3 = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var in1 = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var in2 = allocAndInitData(Real_type, getActualProblemSize(), vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for i in ibegin..<iend {
              out1[i] = -in1[i] - in2[i];
              out2[i] = -in1[i] - in2[i];
              out3[i] = -in1[i] - in2[i];
            }
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {
            forall i in ibegin..<iend {
              out1[i] = -in1[i] - in2[i];
              out2[i] = -in1[i] - in2[i];
              out3[i] = -in1[i] - in2[i];
            }
          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {
            const I = ibegin..<iend;
            out1[I] = -in1[I] - in2[I];
            out2[I] = -in1[I] - in2[I];
            out3[I] = -in1[I] - in2[I];
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(out1, getActualProblemSize());
      checksum[vid] += calcChecksum(out2, getActualProblemSize());
      checksum[vid] += calcChecksum(out3, getActualProblemSize());
    }
  }

  class MAT_MAT_SHARED: KernelBase {

    const TL_SZ: Index_type = 16;

    var m_N_default: Index_type;
    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Basic_MAT_MAT_SHARED);

      m_N_default = 1000;
      setDefaultProblemSize(m_N_default*m_N_default);
      setDefaultReps(5);

      m_N = max(sqrt(getTargetProblemSize()):Index_type, 1:Index_type);

      setActualProblemSize(m_N * m_N);

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);

      setBytesPerRep(m_N*m_N*sizeof(Real_type) + m_N*m_N*sizeof(Real_type));

      const no_tiles = (TL_SZ + m_N - 1) / TL_SZ;
      const no_blocks = divceil(m_N, TL_SZ);
      setFLOPsPerRep(2 * TL_SZ * TL_SZ * TL_SZ * no_tiles * no_blocks * no_blocks);

      checksum_scale_factor = 1e-6 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Teams);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      const NN = m_N * m_N;

      var A = allocAndInitDataConst(Real_type, NN, 1.0, vid);
      var B = allocAndInitDataConst(Real_type, NN, 1.0, vid);
      var C = allocAndInitDataConst(Real_type, NN, 0.0, vid);

      const run_reps = getRunReps();
      const N = m_N;

      const Nx = divceil(N, TL_SZ);
      const Ny = divceil(N, TL_SZ);

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {
            for cy in 0..<Ny {
              for cx in 0..<Nx {
                var As: [0..<TL_SZ, 0..<TL_SZ] Real_type;
                var Bs: [0..<TL_SZ, 0..<TL_SZ] Real_type;
                var Cs: [0..<TL_SZ, 0..<TL_SZ] Real_type;

                for ty in 0..<TL_SZ {
                  for tx in 0..<TL_SZ {
                    Cs[ty, tx] = 0;
                  }
                }

                for k in 0..<(TL_SZ+N-1)/TL_SZ {

                  for ty in 0..<TL_SZ {
                    for tx in 0..<TL_SZ {
                      const Row = cy * TL_SZ + ty;
                      const Col = cx * TL_SZ + tx;
                      if k * TL_SZ + tx < N && Row < N then
                        As[ty, tx] = A[Row * N + k * TL_SZ + tx];
                      else
                        As[ty, tx] = 0.0;
                      if k * TL_SZ + ty < N && Col < N then
                        Bs[ty, tx] = B[(k * TL_SZ + ty) * N + Col];
                      else
                        Bs[ty, tx] = 0.0;
                    }
                  }

                  for ty in 0..<TL_SZ {
                    for tx in 0..<TL_SZ {
                      for n in 0..<TL_SZ {
                        Cs[ty, tx] += As[ty, n] * Bs[n, tx];
                      }
                    }
                  }

                }  // Sequential loop

                for ty in 0..<TL_SZ {
                  for tx in 0..<TL_SZ {
                    const Row = cy * TL_SZ + ty;
                    const Col = cx * TL_SZ + tx;
                    if Row < N && Col < N then
                      C[Col + N * Row] = Cs[ty, tx];
                  }
                }
              }
            }
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(C, m_N*m_N, checksum_scale_factor:Real_type);
    }
  }

  class REDUCE3_INT: KernelBase {

    var m_vsum: Int_type;
    var m_vsum_init: Int_type;
    var m_vmin: Int_type;
    var m_vmin_init: Int_type;
    var m_vmax: Int_type;
    var m_vmax_init: Int_type;

    proc init() {
      super.init(KernelID.Basic_REDUCE3_INT);

      setDefaultProblemSize(1000000);
      setDefaultReps(50);

      setActualProblemSize( getTargetProblemSize() );

      setItsPerRep( getActualProblemSize() );
      setKernelsPerRep(1);
      setBytesPerRep( (3*sizeof(Int_type) + 3*sizeof(Int_type)) +
                      (0*sizeof(Int_type) + 1*sizeof(Int_type)) * getActualProblemSize() );
      setFLOPsPerRep(1 * getActualProblemSize() + 1);

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.Reduction);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Reduction_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var vec = allocAndInitData(Int_type, getActualProblemSize(), vid);

      m_vsum = 0;
      m_vsum_init = 0;
      m_vmin = max(Int_type);
      m_vmin_init = max(Int_type);
      m_vmax = min(Int_type);
      m_vmax_init = min(Int_type);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            var vsum = m_vsum_init;
            var vmin = m_vmin_init;
            var vmax = m_vmax_init;

            for i in ibegin..<iend {
              vsum += vec[i];
              vmin = min(vmin, vec[i]);
              vmax = max(vmax, vec[i]);
            }

            m_vsum += vsum;
            m_vmin  = min(m_vmin, vmin);
            m_vmax  = max(m_vmax, vmax);

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            var vsum = m_vsum_init;
            var vmin = m_vmin_init;
            var vmax = m_vmax_init;

            forall i in ibegin..<iend with (+   reduce vsum,
                                            min reduce vmin,
                                            max reduce vmax) {
              vsum += vec[i];
              vmin reduce= vec[i];
              vmax reduce= vec[i];
            }

            m_vsum += vsum;
            m_vmin  = min(m_vmin, vmin);
            m_vmax  = max(m_vmax, vmax);

          }

          stopTimer();
        }

        when VariantID.Reduction_Chpl {
          startTimer();

          for 0..#run_reps {

            var vsum = +   reduce vec;
            var vmin = min reduce vec;
            var vmax = max reduce vec;

            m_vsum += vsum;
            m_vmin  = min(m_vmin, vmin);
            m_vmax  = max(m_vmax, vmax);

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += m_vsum;
      checksum[vid] += m_vmin;
      checksum[vid] += m_vmax;
    }
  }
}
