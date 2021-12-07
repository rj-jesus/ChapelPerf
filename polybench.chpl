module polybench {
  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private use KernelBase;
  private use Utils;

  class P2MM: KernelBase {

    var m_ni: Index_type;
    var m_nj: Index_type;
    var m_nk: Index_type;
    var m_nl: Index_type;

    var m_alpha: Real_type;
    var m_beta: Real_type;

    proc init() {
      super.init(KernelID.Polybench_2MM);

      var ni_default: Index_type = 1000;
      var nj_default: Index_type = 1000;
      var nk_default: Index_type = 1120;
      var nl_default: Index_type = 1000;

      setDefaultProblemSize(max(ni_default*nj_default, ni_default*nl_default));
      setDefaultReps(2);

      m_ni = sqrt(getTargetProblemSize()):Index_type + 1;
      m_nj = m_ni;
      m_nk = nk_default;
      m_nl = m_ni;

      m_alpha = 1.5;
      m_beta = 1.2;

      setActualProblemSize(max(m_ni*m_nj, m_ni*m_nl));

      setItsPerRep(m_ni*m_nj + m_ni*m_nl);
      setKernelsPerRep(2);
      setBytesPerRep((1*sizeof(Real_type ) + 0*sizeof(Real_type )) * m_ni * m_nj +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_ni * m_nk +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_nj * m_nk +

                     (1*sizeof(Real_type ) + 0*sizeof(Real_type )) * m_ni * m_nl +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_ni * m_nj +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_nj * m_nl);
      setFLOPsPerRep(3 * m_ni*m_nj*m_nk + 2 * m_ni*m_nj*m_nl);

      checksum_scale_factor = 0.000001 *
        getDefaultProblemSize():Checksum_type/getActualProblemSize();

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var tmp = allocAndInitData(Real_type, m_ni * m_nj, vid);
      var A = allocAndInitData(Real_type, m_ni * m_nk, vid);
      var B = allocAndInitData(Real_type, m_nk * m_nj, vid);
      var C = allocAndInitData(Real_type, m_nj * m_nl, vid);
      var D = allocAndInitDataConst(Real_type, m_ni * m_nl, 0.0, vid);

      const run_reps = getRunReps();

      var alpha: Real_type = m_alpha;
      var beta: Real_type = m_beta;

      const ni = m_ni;
      const nj = m_nj;
      const nk = m_nk;
      const nl = m_nl;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for i in 0..<ni {
              for j in 0..<nj {
                var dot: Real_type = 0.0;

                for k in 0..<nk do
                  dot += alpha * A[k + i*nk] * B[j + k*nj];

                tmp[j + i*nj] = dot;
              }
            }

            for i in 0..<ni {
              for l in 0..<nl {
                var dot: Real_type = beta;

                for j in 0..<nj do
                  dot += tmp[j + i*nj] * C[l + j*nl];

                D[l + i*nl] = dot;
              }
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(D, m_ni * m_nl, checksum_scale_factor:Real_type);
    }
  }

  class P3MM: KernelBase {

    var m_ni: Index_type;
    var m_nj: Index_type;
    var m_nk: Index_type;
    var m_nl: Index_type;
    var m_nm: Index_type;

    var m_alpha: Real_type;
    var m_beta: Real_type;

    proc init() {
      super.init(KernelID.Polybench_3MM);

      var ni_default: Index_type = 1000;
      var nj_default: Index_type = 1000;
      var nk_default: Index_type = 1010;
      var nl_default: Index_type = 1000;
      var nm_default: Index_type = 1200;

      setDefaultProblemSize(max(max(ni_default*nj_default,
                                    nj_default*nl_default),
                                ni_default*nl_default));
      setDefaultProblemSize(ni_default * nj_default);
      setDefaultReps(2);

      m_ni = sqrt(getTargetProblemSize()):Index_type + 1;
      m_nj = m_ni;
      m_nk = nk_default;
      m_nl = m_ni;
      m_nm = nm_default;

      setActualProblemSize(max(max(m_ni*m_nj, m_nj*m_nl),
                               m_ni*m_nl));

      setItsPerRep(m_ni*m_nj + m_nj*m_nl + m_ni*m_nl);
      setKernelsPerRep(3);
      setBytesPerRep((1*sizeof(Real_type ) + 0*sizeof(Real_type )) * m_ni * m_nj +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_ni * m_nk +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_nj * m_nk +

                     (1*sizeof(Real_type ) + 0*sizeof(Real_type )) * m_nj * m_nl +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_nj * m_nm +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_nl * m_nm +

                     (1*sizeof(Real_type ) + 0*sizeof(Real_type )) * m_ni * m_nl +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_ni * m_nj +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_nj * m_nl);
      setFLOPsPerRep(2 * m_ni*m_nj*m_nk +
                     2 * m_nj*m_nl*m_nm +
                     2 * m_ni*m_nj*m_nl);

      checksum_scale_factor = 0.000000001 *
        (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var A = allocAndInitData(Real_type, m_ni * m_nk, vid);
      var B = allocAndInitData(Real_type, m_nk * m_nj, vid);
      var C = allocAndInitData(Real_type, m_nj * m_nm, vid);
      var D = allocAndInitData(Real_type, m_nm * m_nl, vid);
      var E = allocAndInitDataConst(Real_type, m_ni * m_nj, 0.0, vid);
      var F = allocAndInitDataConst(Real_type, m_nj * m_nl, 0.0, vid);
      var G = allocAndInitDataConst(Real_type, m_ni * m_nl, 0.0, vid);

      const run_reps = getRunReps();

      const ni = m_ni;
      const nj = m_nj;
      const nk = m_nk;
      const nl = m_nl;
      const nm = m_nm;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for i in 0..<ni {
              for j in 0..<nj {
                var dot = 0.0:Real_type;

                for k in 0..<nk do
                  dot += A[k + i*nk] * B[j + k*nj];

                E[j + i*nj] = dot;
              }
            }

            for j in 0..<nj {
              for l in 0..<nl {
                var dot = 0.0:Real_type;

                for m in 0..<nm do
                  dot += C[m + j*nm] * D[l + m*nl];

                F[l + j*nl] = dot;
              }
            }

            for i in 0..<ni {
              for l in 0..<nl {
                var dot = 0.0:Real_type;

                for j in 0..<nj do
                  dot += E[j + i*nj] * F[l + j*nl];

                G[l + i*nl] = dot;
              }
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(G, m_ni * m_nl, checksum_scale_factor:Real_type);
    }
  }

  class ADI: KernelBase {

    var m_n: Index_type;
    var m_tsteps: Index_type;

    proc init() {
      super.init(KernelID.Polybench_ADI);

      var n_default: Index_type = 1000;

      setDefaultProblemSize((n_default-2) * (n_default-2));
      setDefaultReps(4);

      m_n = sqrt(getTargetProblemSize()):Index_type + 1;
      m_tsteps = 4;

      setItsPerRep(m_tsteps * ((m_n-2) + (m_n-2)));


      setActualProblemSize((m_n-2) * (m_n-2));

      setKernelsPerRep(m_tsteps * 2);
      setBytesPerRep(m_tsteps * ((3*sizeof(Real_type) + 3*sizeof(Real_type)) * m_n * (m_n-2) +
                                 (3*sizeof(Real_type) + 3*sizeof(Real_type)) * m_n * (m_n-2)));
      setFLOPsPerRep(m_tsteps * ((15 + 2) * (m_n-2)*(m_n-2) +
                                 (15 + 2) * (m_n-2)*(m_n-2)));

      checksum_scale_factor = 0.0000001 *
        (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var U = allocAndInitDataConst(Real_type, m_n * m_n, 0.0, vid);
      var V = allocAndInitData(Real_type, m_n * m_n, vid);
      var P = allocAndInitData(Real_type, m_n * m_n, vid);
      var Q = allocAndInitData(Real_type, m_n * m_n, vid);

      const run_reps = getRunReps();

      const n = m_n;
      const tsteps = m_tsteps;

      var DX: Real_type = 1.0/n:Real_type;
      var DY: Real_type = 1.0/n:Real_type;
      var DT: Real_type = 1.0/tsteps:Real_type;
      var B1: Real_type = 2.0;
      var B2: Real_type = 1.0;
      var mul1: Real_type = B1 * DT / (DX * DX);
      var mul2: Real_type = B2 * DT / (DY * DY);
      var a: Real_type = -mul1 / 2.0;
      var b: Real_type = 1.0 + mul1;
      var c: Real_type = a;
      var d: Real_type = -mul2 /2.0;
      var e: Real_type = 1.0 + mul2;
      var f: Real_type = d;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 1..tsteps {

              for i in 1..<n-1 {
                V[0 * n + i] = 1.0;
                P[i * n + 0] = 0.0;
                Q[i * n + 0] = V[0 * n + i];

                for j in 1..<n-1 {
                  P[i * n + j] = -c / (a * P[i * n + j-1] + b);
                  Q[i * n + j] = (-d * U[j * n + i-1] + (1.0 + 2.0*d) * U[j * n + i] -
                                   f * U[j * n + i + 1] - a * Q[i * n + j-1]) /
                                 (a * P[i * n + j-1] + b);
                }

                V[(n-1) * n + i] = 1.0;

                for k in 1..n-2 by -1 do
                  V[k * n + i] = P[i * n + k] * V[(k+1) * n + i] + Q[i * n + k];
              }

              for i in 1..<n-1 {
                U[i * n + 0] = 1.0;
                P[i * n + 0] = 0.0;
                Q[i * n + 0] = U[i * n + 0];

                for j in 1..<n-1 {
                  P[i * n + j] = -f / (d * P[i * n + j-1] + e);
                  Q[i * n + j] = (-a * V[(i-1) * n + j] + (1.0 + 2.0*a) * V[i * n + j] -
                                   c * V[(i + 1) * n + j] - d * Q[i * n + j-1]) /
                                 (d * P[i * n + j-1] + e);
                }

                U[i * n + n-1] = 1.0;

                for k in 1..n-2 by -1 do
                  U[i * n + k] = P[i * n + k] * U[i * n + k +1] + Q[i * n + k];
              }

            }  // tstep loop

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(U, m_n * m_n, checksum_scale_factor:Real_type);
    }
  }

  class ATAX: KernelBase {

    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Polybench_ATAX);

      var N_default: Index_type = 1000;

      setDefaultProblemSize(N_default * N_default);
      setDefaultReps(100);

      m_N = sqrt(getTargetProblemSize()):Index_type+1;

      setActualProblemSize(m_N * m_N);

      setItsPerRep(m_N + m_N);
      setKernelsPerRep(2);
      setBytesPerRep((2*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_N +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_N * m_N +

                     (1*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_N +
                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_N * m_N);
      setFLOPsPerRep(2 * m_N*m_N + 2 * m_N*m_N);

      checksum_scale_factor = 0.001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var tmp = allocAndInitData(Real_type, m_N, vid);
      var x = allocAndInitData(Real_type, m_N, vid);
      var A = allocAndInitData(Real_type, m_N * m_N, vid);
      var y = allocAndInitDataConst(Real_type, m_N, 0.0, vid);

      const run_reps = getRunReps();

      const N = m_N;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for i in 0..<N {
              y[i] = 0.0;
              var dot = 0.0:Real_type;

              for j in 0..<N do
                dot += A[j + i*N] * x[j];

              tmp[i] = dot;
            }

            for j in 0..<N {
              var dot: Real_type = y[j];

              for i in 0..<N do
                dot += A[j + i*N] * tmp[i];

              y[j] = dot;
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(y, m_N, checksum_scale_factor:Real_type);
    }
  }

  class FDTD_2D: KernelBase {

    var m_nx: Index_type;
    var m_ny: Index_type;
    var m_tsteps: Index_type;

    proc init() {
      super.init(KernelID.Polybench_FDTD_2D);

      var nx_default: Index_type = 1000;
      var ny_default: Index_type = 1000;

      setDefaultProblemSize(max((nx_default-1) * ny_default,
                                nx_default * (ny_default-1)));
      setDefaultReps(8);

      m_nx = sqrt(getTargetProblemSize()):Index_type + 1;
      m_ny = m_nx;
      m_tsteps = 40;

      setActualProblemSize(max((m_nx-1)*m_ny, m_nx*(m_ny-1)));

      setItsPerRep(m_tsteps * (m_ny + (m_nx-1)*m_ny +
                               m_nx*(m_ny-1) + (m_nx-1)*(m_ny-1)));
      setKernelsPerRep(m_tsteps * 4);
      setBytesPerRep(m_tsteps * ((0*sizeof(Real_type) + 1*sizeof(Real_type)) +
                                 (1*sizeof(Real_type) + 0*sizeof(Real_type)) * m_ny +

                                 (1*sizeof(Real_type) + 1*sizeof(Real_type)) * (m_nx-1) * m_ny +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_nx * m_ny +

                                 (1*sizeof(Real_type) + 1*sizeof(Real_type)) * m_nx * (m_ny-1) +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_nx * m_ny +

                                 (1*sizeof(Real_type) + 1*sizeof(Real_type)) * (m_nx-1) * (m_ny-1) +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) * (m_nx-1) * m_ny +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_nx * (m_ny-1)));
      setFLOPsPerRep(m_tsteps * ( 0 * m_ny +
                                  3 * (m_nx-1)*m_ny +
                                  3 * m_nx*(m_ny-1) +
                                  5 * (m_nx-1)*(m_ny-1)));

      checksum_scale_factor = 0.001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var hz = allocAndInitDataConst(Real_type, m_nx * m_ny, 0.0, vid);
      var ex = allocAndInitData(Real_type, m_nx * m_ny, vid);
      var ey = allocAndInitData(Real_type, m_nx * m_ny, vid);
      var fict = allocAndInitData(Real_type, m_tsteps, vid);

      const run_reps = getRunReps();

      const nx = m_nx;
      const ny = m_ny;
      const tsteps = m_tsteps;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              for j in 0..<ny do
                ey[j + 0*ny] = fict[t];

              for i in 1..<nx do
                for j in 0..<ny do
                  ey[j + i*ny] = ey[j + i*ny] - 0.5*(hz[j + i*ny] - hz[j + (i-1)*ny]);

              for i in 0..<nx do
                for j in 1..<ny do
                  ex[j + i*ny] = ex[j + i*ny] - 0.5*(hz[j + i*ny] - hz[j-1 + i*ny]);

              for i in 0..<nx - 1 do
                for j in 0..<ny - 1 do
                  hz[j + i*ny] = hz[j + i*ny] - 0.7*(ex[j+1 + i*ny] - ex[j + i*ny] +
                                                     ey[j + (i+1)*ny] - ey[j + i*ny]);

            }  // tstep loop

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(hz, m_nx * m_ny, checksum_scale_factor:Real_type);
    }
  }

  class FLOYD_WARSHALL: KernelBase {

    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Polybench_FLOYD_WARSHALL);

      var N_default: Index_type = 1000;

      setDefaultProblemSize(N_default * N_default);
      setDefaultReps(8);

      m_N = sqrt(getTargetProblemSize()):Index_type + 1;

      setActualProblemSize(m_N * m_N);

      setItsPerRep(m_N*m_N);
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) * m_N * m_N);
      setFLOPsPerRep(1 * m_N*m_N*m_N);

      checksum_scale_factor = 1.0 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var pin = allocAndInitDataRandSign(Real_type, m_N*m_N, vid);
      var pout = allocAndInitDataConst(Real_type, m_N*m_N, 0.0, vid);

      const run_reps = getRunReps();

      const N = m_N;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for k in 0..<N do
              for i in 0..<N do
                for j in 0..<N do
                  pout[j + i*N] = if pin[j + i*N] < pin[k + i*N] + pin[j + k*N]
                                  then pin[j + i*N]
                                  else pin[k + i*N] + pin[j + k*N];

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            for k in 0..<N do
              forall (i, j) in {0..<N, 0..<N} do
                pout[j + i*N] = if pin[j + i*N] < pin[k + i*N] + pin[j + k*N]
                                then pin[j + i*N]
                                else pin[k + i*N] + pin[j + k*N];

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(pout, m_N*m_N, checksum_scale_factor:Real_type);
    }
  }

  class GEMM: KernelBase {

    var m_ni: Index_type;
    var m_nj: Index_type;
    var m_nk: Index_type;

    var m_alpha: Real_type;
    var m_beta: Real_type;

    proc init() {
      super.init(KernelID.Polybench_GEMM);

      var ni_default: Index_type = 1000;
      var nj_default: Index_type = 1000;
      var nk_default: Index_type = 1200;

      setDefaultProblemSize(ni_default * nj_default);
      setDefaultReps(4);

      m_ni = sqrt(getTargetProblemSize()):Index_type + 1;
      m_nj = m_ni;
      m_nk = nk_default;

      m_alpha = 0.62;
      m_beta = 1.002;

      setActualProblemSize(m_ni * m_nj);

      setItsPerRep(m_ni * m_nj);
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 0*sizeof(Real_type)) * m_ni * m_nj +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_ni * m_nk +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_nj * m_nk);
      setFLOPsPerRep((1 + 3 * m_nk) * m_ni*m_nj);

      checksum_scale_factor = 0.001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var A = allocAndInitData(Real_type, m_ni * m_nk, vid);
      var B = allocAndInitData(Real_type, m_nk * m_nj, vid);
      var C = allocAndInitDataConst(Real_type, m_ni * m_nj, 0.0, vid);

      const run_reps = getRunReps();

      const ni = m_ni;
      const nj = m_nj;
      const nk = m_nk;

      var alpha: Real_type = m_alpha;
      var beta: Real_type = m_beta;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for i in 0..<ni {
              for j in 0..<nj {
                var dot: Real_type = 0.0;
                C[j + i*nj] *= beta;

                for k in 0..<nk do
                  dot += alpha * A[k + i*nk] * B[j + k*nj];

                C[j + i*nj] = dot;
              }
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            forall (i, j) in {0..<ni, 0..<nj} {
              var dot: Real_type = 0.0;
              C[j + i*nj] *= beta;

              for k in 0..<nk do
                dot += alpha * A[k + i*nk] * B[j + k*nj];

              C[j + i*nj] = dot;
            }

          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          var Aview = makeArrayFromArray(A, (ni, nk));
          var Bview = makeArrayFromArray(B, (nk, nj));
          var Cview = makeArrayFromArray(C, (ni, nj));

          startTimer();

          for 0..#run_reps {

            const K = 0..<nk;

            forall (i, j) in {0..<ni, 0..<nj} do
              Cview[i, j] = +reduce (alpha*Aview[i, K]*Bview[K, j]);

          }

          stopTimer();

          vcopy(A, Aview);
          vcopy(B, Bview);
          vcopy(C, Cview);
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(C, m_ni * m_nj, checksum_scale_factor:Real_type);
    }
  }

  class GEMVER: KernelBase {

    var m_n: Index_type;
    var m_alpha: Real_type;
    var m_beta: Real_type;

    proc init() {
      super.init(KernelID.Polybench_GEMVER);

      var n_default: Index_type = 1000;

      setDefaultProblemSize(n_default * n_default);
      setDefaultReps(20);

      m_n =  sqrt(getTargetProblemSize()):Index_type + 1;

      m_alpha = 1.5;
      m_beta = 1.2;

      setActualProblemSize(m_n * m_n);

      setItsPerRep(m_n*m_n + m_n*m_n + m_n + m_n*m_n);
      setKernelsPerRep(4);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) * m_n * m_n +
                     (0*sizeof(Real_type ) + 4*sizeof(Real_type )) * m_n +

                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_n * m_n +
                     (1*sizeof(Real_type ) + 2*sizeof(Real_type )) * m_n +

                     (1*sizeof(Real_type ) + 2*sizeof(Real_type )) * m_n +

                     (0*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_n * m_n +
                     (1*sizeof(Real_type ) + 2*sizeof(Real_type )) * m_n);
      setFLOPsPerRep(4 * m_n*m_n +
                     3 * m_n*m_n +
                     1 * m_n +
                     3 * m_n*m_n);

      checksum_scale_factor = 0.001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Forall);
      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Reduction_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var A  = allocAndInitData(Real_type, m_n * m_n, vid);
      var u1 = allocAndInitData(Real_type, m_n, vid);
      var v1 = allocAndInitData(Real_type, m_n, vid);
      var u2 = allocAndInitData(Real_type, m_n, vid);
      var v2 = allocAndInitData(Real_type, m_n, vid);
      var w = allocAndInitDataConst(Real_type, m_n, 0.0, vid);
      var x = allocAndInitData(Real_type, m_n, vid);
      var y = allocAndInitData(Real_type, m_n, vid);
      var z = allocAndInitData(Real_type, m_n, vid);

      const run_reps = getRunReps();

      var alpha: Real_type = m_alpha;
      var beta: Real_type = m_beta;

      const n = m_n;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for i in 0..<n do
              for j in 0..<n do
                A[j + i*n] += u1[i] * v1[j] + u2[i] * v2[j];

            for i in 0..<n {
              var dot: Real_type = 0.0;

              for j in 0..<n do
                dot +=  beta * A[i + j*n] * y[j];

              x[i] += dot;
            }

            for i in 0..<n do
              x[i] += z[i];

            for i in 0..<n {
              var dot: Real_type = w[i];

              for j in 0..<n do
                dot +=  alpha * A[j + i*n] * x[j];

              w[i] = dot;
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            forall i in 0..<n do
              for j in 0..<n do
                A[j + i*n] += u1[i] * v1[j] + u2[i] * v2[j];

            forall i in 0..<n {
              var dot: Real_type = 0.0;

              for j in 0..<n do
                dot +=  beta * A[i + j*n] * y[j];

              x[i] += dot;
            }

            forall i in 0..<n do
              x[i] += z[i];

            forall i in 0..<n {
              var dot: Real_type = w[i];

              for j in 0..<n do
                dot +=  alpha * A[j + i*n] * x[j];

              w[i] = dot;
            }

          }

          stopTimer();
        }

        when VariantID.Reduction_Chpl {
          var Aview = makeArrayFromArray(A, (n, n));

          startTimer();

          for 0..#run_reps {

            forall (i, j) in {0..<n, 0..<n} do
              Aview[i, j] += u1[i] * v1[j] + u2[i] * v2[j];

            forall j in 0..<n do
              x[j] += +reduce (beta*Aview[.., j]*y);

            forall i in 0..<n do
              x[i] += z[i];

            forall i in 0..<n do
              w[i] += +reduce (alpha*Aview[i, ..]*x);

          }

          stopTimer();

          vcopy(A, Aview);
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(w, m_n, checksum_scale_factor:Real_type);
    }
  }

  class GESUMMV: KernelBase {

    var m_N: Index_type;

    var m_alpha: Real_type;
    var m_beta: Real_type;

    proc init() {
      super.init(KernelID.Polybench_GESUMMV);

      var N_default: Index_type = 1000;

      setDefaultProblemSize(N_default * N_default);
      setDefaultReps(120);

      m_N = sqrt(getTargetProblemSize()):Index_type + 1;

      m_alpha = 0.62;
      m_beta = 1.002;

      setActualProblemSize(m_N * m_N);

      setItsPerRep(m_N);
      setKernelsPerRep(1);
      setBytesPerRep((2*sizeof(Real_type ) + 1*sizeof(Real_type )) * m_N +
                     (0*sizeof(Real_type ) + 2*sizeof(Real_type )) * m_N * m_N);
      setFLOPsPerRep((4*m_N + 3) * m_N);

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Reduction_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitData(Real_type, m_N, vid);
      var y = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      var A = allocAndInitData(Real_type, m_N * m_N, vid);
      var B = allocAndInitData(Real_type, m_N * m_N, vid);

      const run_reps = getRunReps();

      const N = m_N;

      var alpha: Real_type = m_alpha;
      var beta: Real_type = m_beta;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for i in 0..<N {
              var tmpdot: Real_type = 0.0;
              var ydot: Real_type = 0.0;

              for j in 0..<N {
                tmpdot += A[j + i*N] * x[j];
                ydot   += B[j + i*N] * x[j];
              }

              y[i] = alpha * tmpdot + beta * ydot;
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            forall i in 0..<N {
              var tmpdot: Real_type = 0.0;
              var ydot: Real_type = 0.0;

              for j in 0..<N {
                tmpdot += A[j + i*N] * x[j];
                ydot   += B[j + i*N] * x[j];
              }

              y[i] = alpha * tmpdot + beta * ydot;
            }

          }

          stopTimer();
        }

        when VariantID.Reduction_Chpl {
          startTimer();

          for 0..#run_reps {

            forall i in 0..<N {
              const J = 0..<N;

              var tmpdot = +reduce (A[i*N+J]*x);
              var ydot   = +reduce (B[i*N+J]*x);

              y[i] = alpha * tmpdot + beta * ydot;
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(y, m_N);
    }
  }

  class HEAT_3D: KernelBase {

    var m_N: Index_type;
    var m_tsteps: Index_type;

    proc init() {
      super.init(KernelID.Polybench_HEAT_3D);

      var N_default: Index_type = 100;

      setDefaultProblemSize((N_default-2)*(N_default-2)*(N_default-2));
      setDefaultReps(20);

      m_N = cbrt(getTargetProblemSize()):Index_type + 1;
      m_tsteps = 20;

      setActualProblemSize((m_N-2) * (m_N-2) * (m_N-2));

      setItsPerRep(m_tsteps * (2 * getActualProblemSize()));
      setKernelsPerRep(m_tsteps * 2);
      setBytesPerRep(m_tsteps * ((1*sizeof(Real_type) + 0*sizeof(Real_type)) *
                                 (m_N-2) * (m_N-2) * (m_N-2) +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) *
                                 (m_N * m_N * m_N - 12*(m_N-2) - 8) +
                                 (1*sizeof(Real_type) + 0*sizeof(Real_type)) *
                                 (m_N-2) * (m_N-2) * (m_N-2) +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) *
                                 (m_N * m_N * m_N - 12*(m_N-2) - 8)));
      setFLOPsPerRep(m_tsteps * (15 * (m_N-2) * (m_N-2) * (m_N-2) +
                                 15 * (m_N-2) * (m_N-2) * (m_N-2)));

      checksum_scale_factor = 0.0001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var A = allocAndInitData(Real_type, m_N*m_N*m_N, vid);
      var B = allocAndInitData(Real_type, m_N*m_N*m_N, vid);
      allocAndInitDataConst(Real_type, m_N*m_N*m_N, 0.0, vid);
      allocAndInitDataConst(Real_type, m_N*m_N*m_N, 0.0, vid);

      const run_reps = getRunReps();

      const N = m_N;
      const tsteps = m_tsteps;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              for i in 1..<N-1 do
                for j in 1..<N-1 do
                  for k in 1..<N-1 do
                    B[k + N*(j + N*i)] = 0.125*(A[k + N*(j + N*(i+1))] - 2.0*A[k + N*(j + N*i)] +
                                                A[k + N*(j + N*(i-1))]) +
                                         0.125*(A[k + N*(j+1 + N*i)]   - 2.0*A[k + N*(j + N*i)] +
                                                A[k + N*(j-1 + N*i)]) +
                                         0.125*(A[k+1 + N*(j + N*i)]   - 2.0*A[k + N*(j + N*i)] +
                                                A[k-1 + N*(j + N*i)]) +
                                         A[k + N*(j + N*i)];

              for i in 1..<N-1 do
                for j in 1..<N-1 do
                  for k in 1..<N-1 do
                    A[k + N*(j + N*i)] = 0.125*(B[k + N*(j + N*(i+1))] - 2.0*B[k + N*(j + N*i)] +
                                                B[k + N*(j + N*(i-1))]) +
                                         0.125*(B[k + N*(j+1 + N*i)]   - 2.0*B[k + N*(j + N*i)] +
                                                B[k + N*(j-1 + N*i)]) +
                                         0.125*(B[k+1 + N*(j + N*i)]   - 2.0*B[k + N*(j + N*i)] +
                                                B[k-1 + N*(j + N*i)]) +
                                         B[k + N*(j + N*i)];

            }  // tstep loop

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              forall (i, j) in zip(1..<N-1, 1..<N-1) do
                for k in 1..<N-1 do
                  B[k + N*(j + N*i)] = 0.125*(A[k + N*(j + N*(i+1))] - 2.0*A[k + N*(j + N*i)] +
                                              A[k + N*(j + N*(i-1))]) +
                                       0.125*(A[k + N*(j+1 + N*i)]   - 2.0*A[k + N*(j + N*i)] +
                                              A[k + N*(j-1 + N*i)]) +
                                       0.125*(A[k+1 + N*(j + N*i)]   - 2.0*A[k + N*(j + N*i)] +
                                              A[k-1 + N*(j + N*i)]) +
                                       A[k + N*(j + N*i)];

              for (i, j) in zip(1..<N-1, 1..<N-1) do
                for k in 1..<N-1 do
                  A[k + N*(j + N*i)] = 0.125*(B[k + N*(j + N*(i+1))] - 2.0*B[k + N*(j + N*i)] +
                                              B[k + N*(j + N*(i-1))]) +
                                       0.125*(B[k + N*(j+1 + N*i)]   - 2.0*B[k + N*(j + N*i)] +
                                              B[k + N*(j-1 + N*i)]) +
                                       0.125*(B[k+1 + N*(j + N*i)]   - 2.0*B[k + N*(j + N*i)] +
                                              B[k-1 + N*(j + N*i)]) +
                                       B[k + N*(j + N*i)];

            }  // tstep loop

          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          var Aview = makeArrayFromArray(A, (N, N, N)),
              Bview = makeArrayFromArray(B, (N, N, N));

          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              const I = 1..<N-1, J = 1..<N-1, K = 1..<N-1;

              Bview[I, J, K] = 0.125*(Aview[I+1, J, K] - 2.0*Aview[I, J, K] + Aview[I-1, J, K]) +
                               0.125*(Aview[I, J+1, K] - 2.0*Aview[I, J, K] + Aview[I, J-1, K]) +
                               0.125*(Aview[I, J, K+1] - 2.0*Aview[I, J, K] + Aview[I, J, K-1]) +
                               Aview[I, J, K];

              Aview[I, J, K] = 0.125*(Bview[I+1, J, K] - 2.0*Bview[I, J, K] + Bview[I-1, J, K]) +
                               0.125*(Bview[I, J+1, K] - 2.0*Bview[I, J, K] + Bview[I, J-1, K]) +
                               0.125*(Bview[I, J, K+1] - 2.0*Bview[I, J, K] + Bview[I, J, K-1]) +
                               Bview[I, J, K];

            }  // tstep loop

          }

          stopTimer();

          vcopy(A, Aview);
          vcopy(B, Bview);
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(A, m_N*m_N*m_N, checksum_scale_factor:Real_type);
      checksum[vid] += calcChecksum(B, m_N*m_N*m_N, checksum_scale_factor:Real_type);
    }
  }

  class JACOBI_1D: KernelBase {

    var m_N: Index_type;
    var m_tsteps: Index_type;

    proc init() {
      super.init(KernelID.Polybench_JACOBI_1D);

      var N_default: Index_type = 1000000;

      setDefaultProblemSize(N_default-2);
      setDefaultReps(100);

      m_N = getTargetProblemSize();
      m_tsteps = 16;

      setActualProblemSize(m_N-2);

      setItsPerRep(m_tsteps * (2 * getActualProblemSize()));
      setKernelsPerRep(m_tsteps * 2);
      setBytesPerRep(m_tsteps * ((1*sizeof(Real_type ) + 0*sizeof(Real_type )) *
                                 (m_N-2) +
                                 (0*sizeof(Real_type ) + 1*sizeof(Real_type )) *
                                 m_N +
                                 (1*sizeof(Real_type ) + 0*sizeof(Real_type )) *
                                 (m_N-2) +
                                 (0*sizeof(Real_type ) + 1*sizeof(Real_type )) *
                                 m_N));
      setFLOPsPerRep(m_tsteps * (3 * (m_N-2) + 3 * (m_N-2)));

      checksum_scale_factor = 0.0001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var A = allocAndInitData(Real_type, m_N, vid);
      var B = allocAndInitData(Real_type, m_N, vid);
      allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      allocAndInitDataConst(Real_type, m_N, 0.0, vid);

      const run_reps = getRunReps();

      const N = m_N;
      const tsteps = m_tsteps;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              for i in 1..<N-1 do
                B[i] = 0.33333 * (A[i-1] + A[i] + A[i+1]);
              for i in 1..<N-1 do
                A[i] = 0.33333 * (B[i-1] + B[i] + B[i+1]);

            }  // tstep loop

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              forall i in 1..<N-1 do
                B[i] = 0.33333 * (A[i-1] + A[i] + A[i+1]);
              forall i in 1..<N-1 do
                A[i] = 0.33333 * (B[i-1] + B[i] + B[i+1]);

            }  // tstep loop

          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              const I = 1..<N-1;

              B[I] = 0.33333 * (A[I-1] + A[I] + A[I+1]);
              A[I] = 0.33333 * (B[I-1] + B[I] + B[I+1]);

            }  // tstep loop

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(A, m_N, checksum_scale_factor:Real_type);
      checksum[vid] += calcChecksum(B, m_N, checksum_scale_factor:Real_type);
    }
  }

  class JACOBI_2D: KernelBase {

    var m_N: Index_type;
    var m_tsteps: Index_type;

    proc init() {
      super.init(KernelID.Polybench_JACOBI_2D);

      var N_default: Index_type = 1000;

      setDefaultProblemSize(N_default * N_default);
      setDefaultReps(50);

      m_N = sqrt(getTargetProblemSize()):Index_type + 1;
      m_tsteps = 40;

      setActualProblemSize((m_N-2) * (m_N-2));

      setItsPerRep(m_tsteps * (2 * (m_N-2) * (m_N-2)));
      setKernelsPerRep(2);
      setBytesPerRep(m_tsteps * ((1*sizeof(Real_type) + 0*sizeof(Real_type)) *
                                 (m_N-2) * (m_N-2) +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) *
                                 (m_N * m_N - 4) +
                                 (1*sizeof(Real_type) + 0*sizeof(Real_type)) *
                                 (m_N-2) * (m_N-2) +
                                 (0*sizeof(Real_type) + 1*sizeof(Real_type)) *
                                 (m_N * m_N  - 4)));
      setFLOPsPerRep(m_tsteps * (5 * (m_N-2)*(m_N-2) +
                                 5 * (m_N-2)*(m_N-2)));

      checksum_scale_factor = 0.0001 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var A = allocAndInitData(Real_type, m_N*m_N, vid);
      var B = allocAndInitData(Real_type, m_N*m_N, vid);
      allocAndInitDataConst(Real_type, m_N*m_N, 0.0, vid);
      allocAndInitDataConst(Real_type, m_N*m_N, 0.0, vid);

      const run_reps = getRunReps();

      const N = m_N;
      const tsteps = m_tsteps;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              for i in 1..<N-1 do
                for j in 1..<N-1 do
                  B[j + i*N] = 0.2 * (A[j + i*N] + A[j-1 + i*N] + A[j+1 + i*N] +
                                                   A[j + (i+1)*N] + A[j + (i-1)*N]);

              for i in 1..<N-1 do
                for j in 1..<N-1 do
                  A[j + i*N] = 0.2 * (B[j + i*N] + B[j-1 + i*N] + B[j+1 + i*N] +
                                                   B[j + (i+1)*N] + B[j + (i-1)*N]);

            }  // tstep loop

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              forall i in 1..<N-1 do
                for j in 1..<N-1 do
                  B[j + i*N] = 0.2 * (A[j + i*N] + A[j-1 + i*N] + A[j+1 + i*N] +
                                                   A[j + (i+1)*N] + A[j + (i-1)*N]);

              forall i in 1..<N-1 do
                for j in 1..<N-1 do
                  A[j + i*N] = 0.2 * (B[j + i*N] + B[j-1 + i*N] + B[j+1 + i*N] +
                                                   B[j + (i+1)*N] + B[j + (i-1)*N]);

            }  // tstep loop

          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          var Aview = makeArrayFromArray(A, (N, N)),
              Bview = makeArrayFromArray(B, (N, N));

          startTimer();

          for 0..#run_reps {

            for t in 0..<tsteps {

              const I = 1..<N-1, J = 1..<N-1;

              Bview[I, J] = 0.2 * (Aview[I, J] + Aview[I, J-1] + Aview[I, J+1] +
                                                 Aview[I+1, J] + Aview[I-1, J]);

              Aview[I, J] = 0.2 * (Bview[I, J] + Bview[I, J-1] + Bview[I, J+1] +
                                                 Bview[I+1, J] + Bview[I-1, J]);

            }  // tstep loop

          }

          stopTimer();

          vcopy(A, Aview);
          vcopy(B, Bview);
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(A, m_N*m_N, checksum_scale_factor:Real_type);
      checksum[vid] += calcChecksum(B, m_N*m_N, checksum_scale_factor:Real_type);
    }
  }

  class MVT: KernelBase {

    var m_N: Index_type;

    proc init() {
      super.init(KernelID.Polybench_MVT);

      var N_default: Index_type = 1000;

      setDefaultProblemSize(N_default * N_default);
      setDefaultReps(100);

      m_N = sqrt(getTargetProblemSize()):Index_type + 1;

      setActualProblemSize(m_N * m_N);

      setItsPerRep(2 * m_N);
      setKernelsPerRep(2);
      setBytesPerRep((1*sizeof(Real_type) + 2*sizeof(Real_type)) * m_N +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_N * m_N +
                     (1*sizeof(Real_type) + 2*sizeof(Real_type)) * m_N +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_N * m_N);
      setFLOPsPerRep(2 * m_N*m_N + 2 * m_N*m_N);

      checksum_scale_factor = 1.0 *
        (getDefaultProblemSize():Checksum_type / getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
      setVariantDefined(VariantID.Promotion_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var y1 = allocAndInitData(Real_type, m_N, vid);
      var y2 = allocAndInitData(Real_type, m_N, vid);
      var  A = allocAndInitData(Real_type, m_N * m_N, vid);
      var x1 = allocAndInitDataConst(Real_type, m_N, 0.0, vid);
      var x2 = allocAndInitDataConst(Real_type, m_N, 0.0, vid);

      const run_reps = getRunReps();

      const N = m_N;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..#run_reps {

            for i in 0..<N {
              var dot: Real_type = 0.0;

              for j in 0..<N do
                dot += A[j + i*N] * y1[j];

              x1[i] += dot;
            }

            for i in 0..<N {
              var dot: Real_type = 0.0;

              for j in 0..<N do
                dot += A[i + j*N] * y2[i];

              x2[i] += dot;
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..#run_reps {

            forall i in 0..<N {
              var dot: Real_type = 0.0;

              for j in 0..<N do
                dot += A[j + i*N] * y1[j];

              x1[i] += dot;
            }

            forall i in 0..<N {
              var dot: Real_type = 0.0;

              for j in 0..<N do
                dot += A[i + j*N] * y2[i];

              x2[i] += dot;
            }

          }

          stopTimer();
        }

        when VariantID.Promotion_Chpl {
          startTimer();

          for 0..#run_reps {

            const I = 0..<N, J = 0..<N;

            forall i in I do
              x1[i] += +reduce (A[i*N+J]*y1);

            forall j in J do
              x2[j] += +reduce (A[I*N+j]*y2);

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(x1, m_N, checksum_scale_factor:Real_type);
      checksum[vid] += calcChecksum(x2, m_N, checksum_scale_factor:Real_type);
    }
  }

}
