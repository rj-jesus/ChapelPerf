module apps {
  private use DataTypes;
  private use DataUtils;
  private use Enums;
  private use KernelBase;
  private use Utils;

  class ADomain {
    var ndims: Index_type;
    var NPNL: Index_type;
    var NPNR: Index_type;

    var imin: Index_type;
    var jmin: Index_type;
    var kmin: Index_type;
    var imax: Index_type;
    var jmax: Index_type;
    var kmax: Index_type;

    var jp: Index_type;
    var kp: Index_type;
    var nnalls: Index_type;

    var fpn: Index_type;
    var lpn: Index_type;
    var frn: Index_type;
    var lrn: Index_type;

    var fpz: Index_type;
    var lpz: Index_type;

    var real_zones: [0..<nnalls] Index_type;
    var n_real_zones: Index_type;

    proc init(rzmax: Index_type, ndims: Index_type) {
      this.ndims = ndims;
      this.NPNL = 2;
      this.NPNR = 1;

      imin = NPNL;
      jmin = NPNL;
      imax = rzmax + NPNR;
      jmax = rzmax + NPNR;
      jp = imax - imin + 1 + NPNL + NPNR;

      if ndims == 2 {
        kmin = 0;
        kmax = 0;
        kp = 0;
        nnalls = jp * (jmax - jmin + 1 + NPNL + NPNR) ;
      } else if ndims == 3 {
        kmin = NPNL;
        kmax = rzmax + NPNR;
        kp = jp * (jmax - jmin + 1 + NPNL + NPNR);
        nnalls = kp * (kmax - kmin + 1 + NPNL + NPNR) ;
      }

      fpn = 0;
      lpn = nnalls - 1;
      frn = fpn + NPNL * (kp + jp) + NPNL;
      lrn = lpn - NPNR * (kp + jp) - NPNR;

      fpz = frn - jp - kp - 1;
      lpz = lrn;

      real_zones = -1;

      n_real_zones = 0;

      if ndims == 2 {

        for j in jmin..<jmax {
          for i in imin..<imax {
            var ip: Index_type = i + j*jp ;
            var id: Index_type = n_real_zones;
            real_zones[id] = ip;
            n_real_zones += 1;
          }
        }

      } else if ndims == 3 {

        for k in kmin..<kmax {
          for j in jmin..<jmax {
            for i in imin..<imax {
              var ip: Index_type = i + j*jp + kp*k ;
              var id: Index_type = n_real_zones;
              real_zones[id] = ip;
              n_real_zones += 1;
            }
          }
        }

      }
    }
  }

  //
  // Set mesh positions for 2d mesh.
  //
  proc setMeshPositions_2d(ref x: [] Real_type, dx: Real_type,
                           ref y: [] Real_type, dy: Real_type,
                           const ref dom: ADomain)
  {
    if dom.ndims != 2 then
      halt("\n******* ERROR!!! domain is not 2d *******");

    var imin: Index_type = dom.imin;
    var imax: Index_type = dom.imax;
    var jmin: Index_type = dom.jmin;
    var jmax: Index_type = dom.jmax;

    var jp: Index_type = dom.jp;

    var npnl: Index_type = dom.NPNL;
    var npnr: Index_type = dom.NPNR;

    // NDSET2D(domain.jp, x, x1,x2,x3,x4);
    ref x4 = reindex(x,            0..);
    ref x1 = reindex(x4[1..],      0..);
    ref x2 = reindex(x1[dom.jp..], 0..);
    ref x3 = reindex(x4[dom.jp..], 0..);

    // NDSET2D(domain.jp, y, y1,y2,y3,y4);
    ref y4 = reindex(y,            0..);
    ref y1 = reindex(y4[1..],      0..);
    ref y2 = reindex(y1[dom.jp..], 0..);
    ref y3 = reindex(y4[dom.jp..], 0..);

    for j in jmin-npnl..<jmax+npnr {
      for i in imin-npnl..<imax+npnr {
        const iz: Index_type = i + j*jp;

        setv(x3[iz], x4[iz], i*dx);
        setv(x1[iz], x2[iz], (i+1)*dx);

        setv(y1[iz], y4[iz], j*dy);
        setv(y2[iz], y3[iz], (j+1)*dy);

      }
    }
  }

  //
  // Set mesh positions for 3d mesh.
  //
  proc setMeshPositions_3d(ref x: [] Real_type, dx: Real_type,
                           ref y: [] Real_type, dy: Real_type,
                           ref z: [] Real_type, dz: Real_type,
                           const ref dom: ADomain)
  {
    if dom.ndims != 3 then
      halt("\n******* ERROR!!! dom is not 3d *******");

    var imin: Index_type = dom.imin;
    var imax: Index_type = dom.imax;
    var jmin: Index_type = dom.jmin;
    var jmax: Index_type = dom.jmax;
    var kmin: Index_type = dom.kmin;
    var kmax: Index_type = dom.kmax;

    var jp: Index_type = dom.jp;
    var kp: Index_type = dom.kp;

    var npnl: Index_type = dom.NPNL;
    var npnr: Index_type = dom.NPNR;

    // NDPTRSET(dom.jp, dom.kp, x,x0,x1,x2,x3,x4,x5,x6,x7);
    ref x0 = reindex(x,            0..);
    ref x1 = reindex(x0[1..],      0..);
    ref x2 = reindex(x0[dom.jp..], 0..);
    ref x3 = reindex(x1[dom.jp..], 0..);
    ref x4 = reindex(x0[dom.kp..], 0..);
    ref x5 = reindex(x1[dom.kp..], 0..);
    ref x6 = reindex(x2[dom.kp..], 0..);
    ref x7 = reindex(x3[dom.kp..], 0..);

    // NDPTRSET(dom.jp, dom.kp, y,y0,y1,y2,y3,y4,y5,y6,y7);
    ref y0 = reindex(y,            0..);
    ref y1 = reindex(y0[1..],      0..);
    ref y2 = reindex(y0[dom.jp..], 0..);
    ref y3 = reindex(y1[dom.jp..], 0..);
    ref y4 = reindex(y0[dom.kp..], 0..);
    ref y5 = reindex(y1[dom.kp..], 0..);
    ref y6 = reindex(y2[dom.kp..], 0..);
    ref y7 = reindex(y3[dom.kp..], 0..);

    // NDPTRSET(dom.jp, dom.kp, z,z0,z1,z2,z3,z4,z5,z6,z7);
    ref z0 = reindex(z,            0..);
    ref z1 = reindex(z0[1..],      0..);
    ref z2 = reindex(z0[dom.jp..], 0..);
    ref z3 = reindex(z1[dom.jp..], 0..);
    ref z4 = reindex(z0[dom.kp..], 0..);
    ref z5 = reindex(z1[dom.kp..], 0..);
    ref z6 = reindex(z2[dom.kp..], 0..);
    ref z7 = reindex(z3[dom.kp..], 0..);

    for k in kmin-npnl..<kmax+npnr {
      for j in jmin-npnl..<jmax+npnr {
        for i in imin-npnl..<imax+npnr {
          const iz: Index_type = i + j*jp + kp*k;

          setv(x0[iz], x2[iz], x4[iz], x6[iz], i*dx);
          setv(x1[iz], x3[iz], x5[iz], x7[iz], (i+1)*dx);

          setv(y0[iz], y1[iz], y4[iz], y5[iz], j*dy);
          setv(y2[iz], y3[iz], y6[iz], y7[iz], (j+1)*dy);

          setv(z0[iz], z1[iz], z2[iz], z3[iz], k*dz);
          setv(z4[iz], z5[iz], z6[iz], z7[iz], (k+1)*dz);

        }
      }
    }
  }

  class DEL_DOT_VEC_2D: KernelBase {

    var m_domain: ADomain;
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Apps_DEL_DOT_VEC_2D);

      setDefaultProblemSize(1000*1000);  // See rzmax in ADomain struct
      setDefaultReps(100);

      var rzmax = sqrt(getTargetProblemSize()):Index_type+1;
      m_domain = new ADomain(rzmax, /* ndims = */ 2);

      m_array_length = m_domain.nnalls;

      setActualProblemSize(m_domain.n_real_zones);

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((0*sizeof(Index_type) + 1*sizeof(Index_type)) * getItsPerRep() +
                     (1*sizeof(Real_type)  + 0*sizeof(Real_type) ) * getItsPerRep() +
                     (0*sizeof(Real_type)  + 4*sizeof(Real_type) ) * (m_domain.imax+1-m_domain.imin)*(m_domain.jmax+1-m_domain.jmin));  // touched data size, not actual number of stores and loads
      setFLOPsPerRep(54 * m_domain.n_real_zones);

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var y = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);

      var dx: Real_type = 0.2;
      var dy: Real_type = 0.1;
      setMeshPositions_2d(x, dx, y, dy, m_domain);

      var xdot = allocAndInitData(Real_type, m_array_length, vid);
      var ydot = allocAndInitData(Real_type, m_array_length, vid);

      var div = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);

      var ptiny: Real_type = 1.0e-20;
      var half: Real_type = 0.5;

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = m_domain.n_real_zones;

      ref real_zones = m_domain.real_zones;

      // NDSET2D(m_domain.jp, x,x1,x2,x3,x4);
      ref x4 = reindex(x,                 0..);
      ref x1 = reindex(x4[1..],           0..);
      ref x2 = reindex(x1[m_domain.jp..], 0..);
      ref x3 = reindex(x4[m_domain.jp..], 0..);

      // NDSET2D(m_domain.jp, y,y1,y2,y3,y4);
      ref y4 = reindex(y,                 0..);
      ref y1 = reindex(y4[1..],           0..);
      ref y2 = reindex(y1[m_domain.jp..], 0..);
      ref y3 = reindex(y4[m_domain.jp..], 0..);

      // NDSET2D(m_domain.jp, xdot,fx1,fx2,fx3,fx4);
      ref fx4 = reindex(xdot,               0..);
      ref fx1 = reindex(fx4[1..],           0..);
      ref fx2 = reindex(fx1[m_domain.jp..], 0..);
      ref fx3 = reindex(fx4[m_domain.jp..], 0..);

      // NDSET2D(m_domain.jp, ydot,fy1,fy2,fy3,fy4);
      ref fy4 = reindex(ydot,               0..);
      ref fy1 = reindex(fx4[1..],           0..);
      ref fy2 = reindex(fx1[m_domain.jp..], 0..);
      ref fy3 = reindex(fx4[m_domain.jp..], 0..);

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..<run_reps do
            for ii in ibegin..<iend {
              const i: Index_type = real_zones[ii];

              const  xi: Real_type = half*( x1[i] +  x2[i] -  x3[i] -  x4[i]);
              const  xj: Real_type = half*( x2[i] +  x3[i] -  x4[i] -  x1[i]);

              const  yi: Real_type = half*( y1[i] +  y2[i] -  y3[i] -  y4[i]);
              const  yj: Real_type = half*( y2[i] +  y3[i] -  y4[i] -  y1[i]);

              const fxi: Real_type = half*(fx1[i] + fx2[i] - fx3[i] - fx4[i]);
              const fxj: Real_type = half*(fx2[i] + fx3[i] - fx4[i] - fx1[i]);

              const fyi: Real_type = half*(fy1[i] + fy2[i] - fy3[i] - fy4[i]);
              const fyj: Real_type = half*(fy2[i] + fy3[i] - fy4[i] - fy1[i]);

              const rarea: Real_type  = 1.0/(xi*yj - xj*yi + ptiny);

              const dfxdx: Real_type  = rarea*(fxi*yj - fxj*yi);

              const dfydy: Real_type  = rarea*(fyj*xi - fyi*xj);

              const affine: Real_type = (fy1[i] + fy2[i] + fy3[i] + fy4[i]) /
                                        ( y1[i] +  y2[i] +  y3[i] +  y4[i]);

              div[i] = dfxdx + dfydy + affine;
            }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..<run_reps do
            forall ii in ibegin..<iend {
              const i: Index_type = real_zones[ii];

              const  xi: Real_type = half*( x1[i] +  x2[i] -  x3[i] -  x4[i]);
              const  xj: Real_type = half*( x2[i] +  x3[i] -  x4[i] -  x1[i]);

              const  yi: Real_type = half*( y1[i] +  y2[i] -  y3[i] -  y4[i]);
              const  yj: Real_type = half*( y2[i] +  y3[i] -  y4[i] -  y1[i]);

              const fxi: Real_type = half*(fx1[i] + fx2[i] - fx3[i] - fx4[i]);
              const fxj: Real_type = half*(fx2[i] + fx3[i] - fx4[i] - fx1[i]);

              const fyi: Real_type = half*(fy1[i] + fy2[i] - fy3[i] - fy4[i]);
              const fyj: Real_type = half*(fy2[i] + fy3[i] - fy4[i] - fy1[i]);

              const rarea: Real_type  = 1.0/(xi*yj - xj*yi + ptiny);

              const dfxdx: Real_type  = rarea*(fxi*yj - fxj*yi);

              const dfydy: Real_type  = rarea*(fyj*xi - fyi*xj);

              const affine: Real_type = (fy1[i] + fy2[i] + fy3[i] + fy4[i]) /
                                        ( y1[i] +  y2[i] +  y3[i] +  y4[i]);

              div[i] = dfxdx + dfydy + affine;
            }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(div, m_array_length);
    }
  }

  class DIFFUSION3DPA: KernelBase {

    param DPA_D1D = 3;
    param DPA_Q1D = 4;
    param SYM = 6;

    var m_NE_default: Index_type;
    var m_NE: Index_type;

    proc init() {
      super.init(KernelID.Apps_DIFFUSION3DPA);

      m_NE_default = 15625;

      setDefaultProblemSize(m_NE_default*DPA_Q1D*DPA_Q1D*DPA_Q1D);
      setDefaultReps(50);

      m_NE = max(getTargetProblemSize()/(DPA_Q1D*DPA_Q1D*DPA_Q1D), 1:Index_type);

      setActualProblemSize(m_NE*DPA_Q1D*DPA_Q1D*DPA_Q1D);

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);

      setBytesPerRep(2*DPA_Q1D*DPA_D1D*sizeof(Real_type)  +
                     DPA_Q1D*DPA_Q1D*DPA_Q1D*SYM*m_NE*sizeof(Real_type) +
                     DPA_D1D*DPA_D1D*DPA_D1D*m_NE*sizeof(Real_type) +
                     DPA_D1D*DPA_D1D*DPA_D1D*m_NE*sizeof(Real_type));

      setFLOPsPerRep(m_NE * (DPA_Q1D * DPA_D1D +
                                   5 * DPA_D1D * DPA_D1D * DPA_Q1D * DPA_D1D +
                                   7 * DPA_D1D * DPA_D1D * DPA_Q1D * DPA_Q1D +
                                   7 * DPA_Q1D * DPA_D1D * DPA_Q1D * DPA_Q1D +
                                  15 * DPA_Q1D * DPA_Q1D * DPA_Q1D +
                             DPA_Q1D * DPA_D1D +
                                   7 * DPA_Q1D * DPA_Q1D * DPA_D1D * DPA_Q1D +
                                   7 * DPA_Q1D * DPA_Q1D * DPA_D1D * DPA_D1D +
                                   7 * DPA_D1D * DPA_Q1D * DPA_D1D * DPA_D1D +
                                   3 * DPA_D1D * DPA_D1D * DPA_D1D));

      setUsesFeature(FeatureID.Teams);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var  Basis = allocAndInitDataConst(Real_type, DPA_Q1D*DPA_D1D,                  1.0:Real_type, vid);
      var dBasis = allocAndInitDataConst(Real_type, DPA_Q1D*DPA_D1D,                  1.0:Real_type, vid);
      var      D = allocAndInitDataConst(Real_type, DPA_Q1D*DPA_Q1D*DPA_Q1D*SYM*m_NE, 1.0:Real_type, vid);
      var      X = allocAndInitDataConst(Real_type, DPA_D1D*DPA_D1D*DPA_D1D*m_NE,     1.0:Real_type, vid);
      var      Y = allocAndInitDataConst(Real_type, DPA_D1D*DPA_D1D*DPA_D1D*m_NE,     0.0:Real_type, vid);

      const run_reps = getRunReps();

      var NE: Index_type = m_NE;
      const symmetric = true;

      inline proc b(x, y) ref return  Basis[x + DPA_Q1D * y];
      inline proc g(x, y) ref return dBasis[x + DPA_Q1D * y];
      inline proc dpaX_(dx, dy, dz, e) ref return
        X[dx + DPA_D1D * dy + DPA_D1D * DPA_D1D * dz + DPA_D1D * DPA_D1D * DPA_D1D * e];
      inline proc dpaY_(dx, dy, dz, e) ref return
        Y[dx + DPA_D1D * dy + DPA_D1D * DPA_D1D * dz + DPA_D1D * DPA_D1D * DPA_D1D * e];
      inline proc d(qx, qy, qz, s, e) ref return
        D[qx + DPA_Q1D * qy + DPA_Q1D * DPA_Q1D * qz + DPA_Q1D * DPA_Q1D * DPA_Q1D * s +  DPA_Q1D * DPA_Q1D * DPA_Q1D * SYM * e];

      // Half of B and G are stored in shared to get B, Bt, G and Gt.
      // Indices computation for SmemPADiffusionApply3D.
      inline proc   qi(const q, const d, const Q) return if q <= d then     q else Q-1-q;
      inline proc   dj(const q, const d, const D) return if q <= d then     d else D-1-d;
      inline proc   qk(const q, const d, const Q) return if q <= d then Q-1-q else     q;
      inline proc   dl(const q, const d, const D) return if q <= d then D-1-d else     d;
      inline proc sign(const q, const d)          return if q <= d then  -1.0 else   1.0;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..<run_reps {
            for e in 0..<NE {

              // note: we use procedures to create ``fake views'' over other
              // arrays since Chapel does not seem to support array reindexing
              // with rank change, see
              // https://matrix.to/#/!PBYDSerrfYujeStENM:gitter.im/$t-4bYrNGMVjXNSRlTiYQ1Pjp_Q7dTHWo2K4E5wERAIU
              // and subsequent messages
              const MQ1: int = DPA_Q1D;
              const MD1: int = DPA_D1D;
              const MDQ: int = if MQ1 > MD1 then MQ1 else MD1;
              var sBG: [0..<MQ1*MD1] real;
              inline proc  B(i, j) ref return sBG[i*MD1+j];
              inline proc  G(i, j) ref return sBG[i*MD1+j];
              inline proc Bt(i, j) ref return sBG[i*MQ1+j];
              inline proc Gt(i, j) ref return sBG[i*MQ1+j];
              var sm0: [0..<3][0..<MDQ*MDQ*MDQ] real;
              var sm1: [0..<3][0..<MDQ*MDQ*MDQ] real;
              inline proc  s_X(i, j, k) ref return sm0[2][(i*MD1+j)*MD1+k];
              inline proc DDQ0(i, j, k) ref return sm0[0][(i*MD1+j)*MQ1+k];
              inline proc DDQ1(i, j, k) ref return sm0[1][(i*MD1+j)*MQ1+k];
              inline proc DQQ0(i, j, k) ref return sm0[0][(i*MQ1+j)*MQ1+k];
              inline proc DQQ1(i, j, k) ref return sm0[1][(i*MQ1+j)*MQ1+k];
              inline proc DQQ2(i, j, k) ref return sm0[2][(i*MQ1+j)*MQ1+k];
              inline proc QQQ0(i, j, k) ref return sm0[0][(i*MQ1+j)*MQ1+k];
              inline proc QQQ1(i, j, k) ref return sm0[1][(i*MQ1+j)*MQ1+k];
              inline proc QQQ2(i, j, k) ref return sm0[2][(i*MQ1+j)*MQ1+k];
              inline proc QQD0(i, j, k) ref return sm0[0][(i*MQ1+j)*MD1+k];
              inline proc QQD1(i, j, k) ref return sm0[1][(i*MQ1+j)*MD1+k];
              inline proc QQD2(i, j, k) ref return sm0[2][(i*MQ1+j)*MD1+k];
              inline proc QDD0(i, j, k) ref return sm0[0][(i*MD1+j)*MD1+k];
              inline proc QDD1(i, j, k) ref return sm0[1][(i*MD1+j)*MD1+k];
              inline proc QDD2(i, j, k) ref return sm0[2][(i*MD1+j)*MD1+k];

              for dy in 0..<DPA_D1D {
                for dx in 0..<DPA_D1D {
                  // DIFFUSION3DPA_1
                  for dz in 0..<DPA_D1D do
                    s_X[dz, dy, dx] = dpaX_(dx, dy, dz, e);
                }
                for qx in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_2
                  const i = qi(qx, dy, DPA_Q1D);
                  const j = dj(qx, dy, DPA_D1D);
                  const k = qk(qx, dy, DPA_Q1D);
                  const l = dl(qx, dy, DPA_D1D);
                  B[i, j] = b(qx, dy);
                  G[k, l] = g(qx, dy) * sign(qx, dy);
                }
              }

              for dy in 0..<DPA_D1D {
                for qx in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_3;
                  var u: [0..<DPA_D1D] real = 0.0;
                  var v: [0..<DPA_D1D] real = 0.0;
                  for dx in 0..<DPA_D1D {
                    const i = qi(qx, dx, DPA_Q1D);
                    const j = dj(qx, dx, DPA_D1D);
                    const k = qk(qx, dx, DPA_Q1D);
                    const l = dl(qx, dx, DPA_D1D);
                    const s = sign(qx, dx);
                    for dz in 0..<DPA_D1D {
                      const coords = s_X[dz, dy, dx];
                      u[dz] += coords * B[i, j];
                      v[dz] += coords * G[k, l] * s;
                    }
                  }
                  for dz in 0..<DPA_D1D {
                    DDQ0[dz, dy, qx] = u[dz];
                    DDQ1[dz, dy, qx] = v[dz];
                  }
                }
              }

              for qy in 0..<DPA_Q1D {
                for qx in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_4;
                  var u: [0..<DPA_D1D] real = 0.0;
                  var v: [0..<DPA_D1D] real = 0.0;
                  var w: [0..<DPA_D1D] real = 0.0;
                  for dy in 0..<DPA_D1D {
                    const i = qi(qy, dy, DPA_Q1D);
                    const j = dj(qy, dy, DPA_D1D);
                    const k = qk(qy, dy, DPA_Q1D);
                    const l = dl(qy, dy, DPA_D1D);
                    const s = sign(qy, dy);
                    for dz in 0..<DPA_D1D {
                      u[dz] += DDQ1[dz, dy, qx] * B[i, j];
                      v[dz] += DDQ0[dz, dy, qx] * G[k, l] * s;
                      w[dz] += DDQ0[dz, dy, qx] * B[i, j];
                    }
                  }
                  for dz in 0..<DPA_D1D {
                    DQQ0[dz, qy, qx] = u[dz];
                    DQQ1[dz, qy, qx] = v[dz];
                    DQQ2[dz, qy, qx] = w[dz];
                  }
                }
              }

              for qy in 0..<DPA_Q1D {
                for qx in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_5;
                  var u: [0..<DPA_Q1D] real = 0.0;
                  var v: [0..<DPA_Q1D] real = 0.0;
                  var w: [0..<DPA_Q1D] real = 0.0;
                  for dz in 0..<DPA_D1D {
                    for qz in 0..<DPA_Q1D {
                      const i = qi(qz, dz, DPA_Q1D);
                      const j = dj(qz, dz, DPA_D1D);
                      const k = qk(qz, dz, DPA_Q1D);
                      const l = dl(qz, dz, DPA_D1D);
                      const s = sign(qz, dz);
                      u[qz] += DQQ0[dz, qy, qx] * B[i, j];
                      v[qz] += DQQ1[dz, qy, qx] * B[i, j];
                      w[qz] += DQQ2[dz, qy, qx] * G[k, l] * s;
                    }
                  }
                  for qz in 0..<DPA_Q1D {
                    const O11 = d(qx, qy, qz, 0, e): real;
                    const O12 = d(qx, qy, qz, 1, e): real;
                    const O13 = d(qx, qy, qz, 2, e): real;
                    const O21 = (if symmetric then                O12 else d(qx, qy, qz, 3, e)): real;
                    const O22 = (if symmetric then d(qx, qy, qz, 3,e) else d(qx, qy, qz, 4, e)): real;
                    const O23 = (if symmetric then d(qx, qy, qz, 4,e) else d(qx, qy, qz, 5, e)): real;
                    const O31 = (if symmetric then                O13 else d(qx, qy, qz, 6, e)): real;
                    const O32 = (if symmetric then                O23 else d(qx, qy, qz, 7, e)): real;
                    const O33 = (if symmetric then d(qx, qy, qz, 5,e) else d(qx, qy, qz, 8, e)): real;
                    const  gX = u[qz];
                    const  gY = v[qz];
                    const  gZ = w[qz];
                    QQQ0[qz, qy, qx] = (O11*gX) + (O12*gY) + (O13*gZ);
                    QQQ1[qz, qy, qx] = (O21*gX) + (O22*gY) + (O23*gZ);
                    QQQ2[qz, qy, qx] = (O31*gX) + (O32*gY) + (O33*gZ);
                  }
                }
              }

              for d in 0..<DPA_D1D {
                for q in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_6;
                  const i = qi(q, d, DPA_Q1D);
                  const j = dj(q, d, DPA_D1D);
                  const k = qk(q, d, DPA_Q1D);
                  const l = dl(q, d, DPA_D1D);
                  Bt[j, i] = b(q, d);
                  Gt[l, k] = g(q, d) * sign(q, d);
                }
              }

              for qy in 0..<DPA_Q1D {
                for dx in 0..<DPA_D1D {
                  // DIFFUSION3DPA_7;
                  var u: [0..<DPA_Q1D] real = 0.0;
                  var v: [0..<DPA_Q1D] real = 0.0;
                  var w: [0..<DPA_Q1D] real = 0.0;
                  for qx in 0..<DPA_Q1D {
                    const i = qi(qx, dx, DPA_Q1D);
                    const j = dj(qx, dx, DPA_D1D);
                    const k = qk(qx, dx, DPA_Q1D);
                    const l = dl(qx, dx, DPA_D1D);
                    const s = sign(qx, dx);
                    for qz in 0..<DPA_Q1D {
                      u[qz] += QQQ0[qz, qy, qx] * Gt[l, k] * s;
                      v[qz] += QQQ1[qz, qy, qx] * Bt[j, i];
                      w[qz] += QQQ2[qz, qy, qx] * Bt[j, i];
                    }
                  }
                  for qz in 0..<DPA_Q1D {
                    QQD0[qz, qy, dx] = u[qz];
                    QQD1[qz, qy, dx] = v[qz];
                    QQD2[qz, qy, dx] = w[qz];
                  }
                }
              }

              for dy in 0..<DPA_D1D {
                for dx in 0..<DPA_D1D {
                  // DIFFUSION3DPA_8;
                  var u: [0..<DPA_Q1D] real = 0.0;
                  var v: [0..<DPA_Q1D] real = 0.0;
                  var w: [0..<DPA_Q1D] real = 0.0;
                  for qy in 0..<DPA_Q1D {
                    const i = qi(qy, dy, DPA_Q1D);
                    const j = dj(qy, dy, DPA_D1D);
                    const k = qk(qy, dy, DPA_Q1D);
                    const l = dl(qy, dy, DPA_D1D);
                    const s = sign(qy, dy);
                    for qz in 0..<DPA_Q1D {
                      u[qz] += QQD0[qz, qy, dx] * Bt[j, i];
                      v[qz] += QQD1[qz, qy, dx] * Gt[l, k] * s;
                      w[qz] += QQD2[qz, qy, dx] * Bt[j, i];
                    }
                  }
                  for qz in 0..<DPA_Q1D {
                    QDD0[qz, dy, dx] = u[qz];
                    QDD1[qz, dy, dx] = v[qz];
                    QDD2[qz, dy, dx] = w[qz];
                  }
                }
              }

              for dy in 0..<DPA_D1D {
                for dx in 0..<DPA_D1D {
                  // DIFFUSION3DPA_9;
                  var u: [0..<DPA_D1D] real = 0.0;
                  var v: [0..<DPA_D1D] real = 0.0;
                  var w: [0..<DPA_D1D] real = 0.0;
                  for qz in 0..<DPA_Q1D {
                    for dz in 0..<DPA_D1D {
                      const i = qi(qz, dz, DPA_Q1D);
                      const j = dj(qz, dz, DPA_D1D);
                      const k = qk(qz, dz, DPA_Q1D);
                      const l = dl(qz, dz, DPA_D1D);
                      const s = sign(qz, dz);
                      u[dz] += QDD0[qz, dy, dx] * Bt[j, i];
                      v[dz] += QDD1[qz, dy, dx] * Bt[j, i];
                      w[dz] += QDD2[qz, dy, dx] * Gt[l, k] * s;
                    }
                  }
                  for dz in 0..<DPA_D1D do
                    dpaY_(dx, dy, dz, e) += (u[dz] + v[dz] + w[dz]);
                }
              }

            }  // element loop
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(Y, DPA_D1D*DPA_D1D*DPA_D1D*m_NE);
    }
  }

  class ENERGY: KernelBase {

    proc init() {
      super.init(KernelID.Apps_ENERGY);

      setDefaultProblemSize(1000000);
      setDefaultReps(130);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(6 * getActualProblemSize());
      setKernelsPerRep(6);
      // some branches are never taken due to the nature of the initialization of delvc
      // the additional reads and writes that would be done if those branches were taken are noted in the comments
      setBytesPerRep((1*sizeof(Real_type) + 5*sizeof(Real_type)) * getActualProblemSize() +
                     (1*sizeof(Real_type) + 1*sizeof(Real_type)) * getActualProblemSize() + /* 1 +  8 */
                     (1*sizeof(Real_type) + 6*sizeof(Real_type)) * getActualProblemSize() +
                     (1*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize() +
                     (1*sizeof(Real_type) + 7*sizeof(Real_type)) * getActualProblemSize() + /* 1 + 12 */
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * getActualProblemSize()); /* 1 +  8 */
      setFLOPsPerRep(( 6 +
                      11 + // 1 sqrt
                       8 +
                       2 +
                      19 + // 1 sqrt
                       9   // 1 sqrt
                     ) * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var e_new        = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var e_old        = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var delvc        = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var p_new        = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var p_old        = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var q_new        = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var q_old        = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var work         = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var compHalfStep = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var pHalfStep    = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var bvc          = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var pbvc         = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var ql_old       = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var qq_old       = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var vnewc        = allocAndInitData(Real_type, getActualProblemSize(), vid);

      const rho0  = initData(Real_type, vid);
      const e_cut = initData(Real_type, vid);
      const emin  = initData(Real_type, vid);
      const q_cut = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..<run_reps {

            for i in ibegin..<iend {
              e_new[i] = e_old[i] - 0.5 * delvc[i] * (p_old[i] + q_old[i])
                                  + 0.5 * work[i];
            }

            for i in ibegin..<iend {
              if delvc[i] > 0.0 then
                q_new[i] = 0.0;
              else {
                var vhalf: Real_type = 1.0/(1.0 + compHalfStep[i]);
                var ssc: Real_type = (pbvc[i] * e_new[i] +
                                      vhalf * vhalf * bvc[i] * pHalfStep[i])/rho0;
                  if ssc <= 0.1111111e-36 then ssc = 0.3333333e-18;
                                          else ssc = sqrt(ssc);
                q_new[i] = (ssc*ql_old[i] + qq_old[i]);
              }
            }

            for i in ibegin..<iend {
              e_new[i] = e_new[i] + 0.5 * delvc[i] * (3.0*(p_old[i] + q_old[i]) -
                                                      4.0*(pHalfStep[i] + q_new[i]));
            }

            for i in ibegin..<iend {
              e_new[i] += 0.5 * work[i];
              if abs(e_new[i]) < e_cut then e_new[i] = 0.0;
              if e_new[i] < emin then e_new[i] = emin;
            }

            for i in ibegin..<iend {
              var q_tilde: Real_type;
              if delvc[i] > 0.0 then
                q_tilde = 0.0;
              else {
                var ssc: Real_type = (pbvc[i] * e_new[i] +
                                      vnewc[i] * vnewc[i] * bvc[i] * p_new[i])/rho0;
                if ssc <= 0.1111111e-36 then ssc = 0.3333333e-18;
                                        else ssc = sqrt(ssc);
                q_tilde = (ssc*ql_old[i] + qq_old[i]);
              }
              e_new[i] = e_new[i] - (7.0*(p_old[i] + q_old[i]) -
                                     8.0*(pHalfStep[i] + q_new[i]) +
                                     (p_new[i] + q_tilde)) * delvc[i] / 6.0;
              if abs(e_new[i]) < e_cut then e_new[i] = 0.0;
              if e_new[i] < emin then e_new[i] = emin;
            }

            for i in ibegin..<iend {
              if delvc[i] <= 0.0 {
                var ssc: Real_type = (pbvc[i] * e_new[i] +
                                      vnewc[i] * vnewc[i] * bvc[i] * p_new[i])/rho0;
                if ssc <= 0.1111111e-36 then ssc = 0.3333333e-18;
                                        else ssc = sqrt(ssc);
                q_new[i] = (ssc*ql_old[i] + qq_old[i]);
                if abs(q_new[i]) < q_cut then q_new[i] = 0.0;
              }
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..<run_reps {

            forall i in ibegin..<iend {
              e_new[i] = e_old[i] - 0.5 * delvc[i] * (p_old[i] + q_old[i])
                                  + 0.5 * work[i];
            }

            forall i in ibegin..<iend {
              if delvc[i] > 0.0 then
                q_new[i] = 0.0;
              else {
                var vhalf: Real_type = 1.0/(1.0 + compHalfStep[i]);
                var ssc: Real_type = (pbvc[i] * e_new[i] +
                                      vhalf * vhalf * bvc[i] * pHalfStep[i])/rho0;
                  if ssc <= 0.1111111e-36 then ssc = 0.3333333e-18;
                                          else ssc = sqrt(ssc);
                q_new[i] = (ssc*ql_old[i] + qq_old[i]);
              }
            }

            forall i in ibegin..<iend {
              e_new[i] = e_new[i] + 0.5 * delvc[i] * (3.0*(p_old[i] + q_old[i]) -
                                                      4.0*(pHalfStep[i] + q_new[i]));
            }

            forall i in ibegin..<iend {
              e_new[i] += 0.5 * work[i];
              if abs(e_new[i]) < e_cut then e_new[i] = 0.0;
              if e_new[i] < emin then e_new[i] = emin;
            }

            forall i in ibegin..<iend {
              var q_tilde: Real_type;
              if delvc[i] > 0.0 then
                q_tilde = 0.0;
              else {
                var ssc: Real_type = (pbvc[i] * e_new[i] +
                                      vnewc[i] * vnewc[i] * bvc[i] * p_new[i])/rho0;
                if ssc <= 0.1111111e-36 then ssc = 0.3333333e-18;
                                        else ssc = sqrt(ssc);
                q_tilde = (ssc*ql_old[i] + qq_old[i]);
              }
              e_new[i] = e_new[i] - (7.0*(p_old[i] + q_old[i]) -
                                     8.0*(pHalfStep[i] + q_new[i]) +
                                     (p_new[i] + q_tilde)) * delvc[i] / 6.0;
              if abs(e_new[i]) < e_cut then e_new[i] = 0.0;
              if e_new[i] < emin then e_new[i] = emin;
            }

            forall i in ibegin..<iend {
              if delvc[i] <= 0.0 {
                var ssc: Real_type = (pbvc[i] * e_new[i] +
                                      vnewc[i] * vnewc[i] * bvc[i] * p_new[i])/rho0;
                if ssc <= 0.1111111e-36 then ssc = 0.3333333e-18;
                                        else ssc = sqrt(ssc);
                q_new[i] = (ssc*ql_old[i] + qq_old[i]);
                if abs(q_new[i]) < q_cut then q_new[i] = 0.0;
              }
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(e_new, getActualProblemSize());
      checksum[vid] += calcChecksum(q_new, getActualProblemSize());
    }
  }

  class FIR: KernelBase {

    param FIR_COEFFLEN = 16;

    var m_coefflen: Index_type;

    proc init() {
      super.init(KernelID.Apps_FIR);

      setDefaultProblemSize(1000000);
      setDefaultReps(160);

      m_coefflen = FIR_COEFFLEN;

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(getActualProblemSize() - m_coefflen);
      setKernelsPerRep(1);
      setBytesPerRep((1*sizeof(Real_type) + 0*sizeof(Real_type)) * getItsPerRep() +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep((2 * m_coefflen) * (getActualProblemSize() - m_coefflen));

      checksum_scale_factor = 0.0001 *
        (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var in_  = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var out_ = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize() - m_coefflen;

      var coeff_array: [0..<FIR_COEFFLEN] Real_type = [ 3.0, -1.0, -1.0, -1.0,
                                                       -1.0,  3.0, -1.0, -1.0,
                                                       -1.0, -1.0,  3.0, -1.0,
                                                       -1.0, -1.0, -1.0,  3.0, ];

      const coefflen = m_coefflen;

      var coeff: [0..<FIR_COEFFLEN] Real_type = coeff_array;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..<run_reps {

            for i in ibegin..<iend {
              var sum: Real_type = 0.0;
              for j in 0..<coefflen {
                sum += coeff[j]*in_[i+j];
              }
              out_[i] = sum;
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..<run_reps {

            forall i in ibegin..<iend {
              var sum: Real_type = 0.0;
              for j in 0..<coefflen {
                sum += coeff[j]*in_[i+j];
              }
              out_[i] = sum;
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(out_, getActualProblemSize(), checksum_scale_factor:Real_type);
    }
  }

  record Extent {
    var i_min: Index_type;
    var i_max: Index_type;
    var j_min: Index_type;
    var j_max: Index_type;
    var k_min: Index_type;
    var k_max: Index_type;
  };

  //
  // Function to generate index lists for packing.
  //
  proc create_pack_lists(const halo_width: Index_type, const ref grid_dims: [] Index_type,
                         const num_neighbors: Index_type, vid: VariantID)
  {
    var pack_index_list_extents: [0..<num_neighbors] Extent;

    // faces
    pack_index_list_extents[0]  = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[1]  = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[2]  = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[3]  = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[4]  = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[5]  = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);

    // edges
    pack_index_list_extents[6]  = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[7]  = new Extent(halo_width  , halo_width   + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[8]  = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[9]  = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[10] = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[11] = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);
    pack_index_list_extents[12] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[13] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);
    pack_index_list_extents[14] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[15] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);
    pack_index_list_extents[16] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[17] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);

    // corners
    pack_index_list_extents[18] = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[19] = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);
    pack_index_list_extents[20] = new Extent(halo_width  , halo_width   + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[21] = new Extent(halo_width  , halo_width   + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);
    pack_index_list_extents[22] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[23] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);
    pack_index_list_extents[24] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[25] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);

    const grid_i_stride: Index_type = 1;
    const grid_j_stride: Index_type = grid_dims[0] + 2*halo_width;
    const grid_k_stride: Index_type = grid_j_stride * (grid_dims[1] + 2*halo_width);

    var pack_index_list_lengths: [0..<num_neighbors] Index_type;

    for (l, extent) in zip(0..<num_neighbors, pack_index_list_extents) do
      pack_index_list_lengths[l] = (extent.i_max - extent.i_min) *
                                   (extent.j_max - extent.j_min) *
                                   (extent.k_max - extent.k_min);

    var pack_index_lists = for l in 0..<num_neighbors do
      allocAndInitData(Int_type, pack_index_list_lengths[l], vid);

    for (l, extent, pack_list) in zip(0..<num_neighbors, pack_index_list_extents, pack_index_lists) {
      var list_idx: Index_type = 0;

      for kk in extent.k_min..<extent.k_max {
        for jj in extent.j_min..<extent.j_max {
          for ii in extent.i_min..<extent.i_max {
            var pack_idx: Index_type = ii * grid_i_stride +
                                       jj * grid_j_stride +
                                       kk * grid_k_stride;

            pack_list[list_idx] = pack_idx:Int_type;

            list_idx += 1;
          }
        }
      }
    }

    return (pack_index_lists, pack_index_list_lengths);
  }

  //
  // Function to generate index lists for unpacking.
  //
  proc create_unpack_lists(const halo_width: Index_type, const ref grid_dims: [] Index_type,
                             const num_neighbors: Index_type, vid: VariantID)
    {
      var unpack_index_list_extents: [0..<num_neighbors] Extent;

      // faces
      unpack_index_list_extents[ 0] = new Extent(0                        ,                  halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[ 1] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[ 2] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 0                        ,                  halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[ 3] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[ 4] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[ 5] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);

      // edges
      unpack_index_list_extents[ 6] = new Extent(0                        ,                  halo_width,
                                                 0                        ,                  halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[ 7] = new Extent(0                        ,                  halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[ 8] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 0                        ,                  halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[ 9] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 halo_width               , grid_dims[2] +   halo_width);
      unpack_index_list_extents[10] = new Extent(0                        ,                  halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[11] = new Extent(0                        ,                  halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);
      unpack_index_list_extents[12] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[13] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 halo_width               , grid_dims[1] +   halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);
      unpack_index_list_extents[14] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 0                        ,                  halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[15] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 0                        ,                  halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);
      unpack_index_list_extents[16] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[17] = new Extent(halo_width               , grid_dims[0] +   halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);

      // corners
      unpack_index_list_extents[18] = new Extent(0                        ,                  halo_width,
                                                 0                        ,                  halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[19] = new Extent(0                        ,                  halo_width,
                                                 0                        ,                  halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);
      unpack_index_list_extents[20] = new Extent(0                        ,                  halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[21] = new Extent(0                        ,                  halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);
      unpack_index_list_extents[22] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 0                        ,                  halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[23] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 0                        ,                  halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);
      unpack_index_list_extents[24] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 0                        ,                  halo_width);
      unpack_index_list_extents[25] = new Extent(grid_dims[0] + halo_width, grid_dims[0] + 2*halo_width,
                                                 grid_dims[1] + halo_width, grid_dims[1] + 2*halo_width,
                                                 grid_dims[2] + halo_width, grid_dims[2] + 2*halo_width);

      const grid_i_stride: Index_type = 1;
      const grid_j_stride: Index_type = grid_dims[0] + 2*halo_width;
      const grid_k_stride: Index_type = grid_j_stride * (grid_dims[1] + 2*halo_width);

      var unpack_index_list_lengths: [0..<num_neighbors] Index_type;

      for (l, extent) in zip(0..<num_neighbors, unpack_index_list_extents) do
        unpack_index_list_lengths[l] = (extent.i_max - extent.i_min) *
                                       (extent.j_max - extent.j_min) *
                                       (extent.k_max - extent.k_min);

      var unpack_index_lists = for l in 0..<num_neighbors do
        allocAndInitData(Int_type, unpack_index_list_lengths[l], vid);

      for (l, extent, unpack_list) in zip(0..<num_neighbors, unpack_index_list_extents, unpack_index_lists) {
        var list_idx: Index_type = 0;

        for kk in extent.k_min..<extent.k_max {
          for jj in extent.j_min..<extent.j_max {
            for ii in extent.i_min..<extent.i_max {
              var unpack_idx: Index_type = ii * grid_i_stride +
                                           jj * grid_j_stride +
                                           kk * grid_k_stride;

              unpack_list[list_idx] = unpack_idx:Int_type;

              list_idx += 1;
            }
          }
        }
      }

      return (unpack_index_lists, unpack_index_list_lengths);
    }

  class HALOEXCHANGE_FUSED: KernelBase {

    record ptr_holder {
      var _buffer_b: Index_type;
      var _buffer_off: Index_type;
      var _list: Index_type;
      var _var: Index_type;
    };

    param s_num_neighbors = 26;

    var m_grid_dims_default: [0..<3] Index_type;
    var m_halo_width_default: Index_type;
    var m_num_vars_default: Index_type;

    var m_grid_dims: [0..<3] Index_type;
    var m_halo_width: Index_type;
    var m_num_vars: Index_type;

    var m_grid_plus_halo_dims: [0..<3] Index_type;
    var m_var_size: Index_type;
    var m_var_halo_size: Index_type;

    proc init() {
      super.init(KernelID.Apps_HALOEXCHANGE_FUSED);

      this.complete();

      m_grid_dims_default[0] = 100;
      m_grid_dims_default[1] = 100;
      m_grid_dims_default[2] = 100;
      m_halo_width_default   = 1;
      m_num_vars_default     = 3;

      setDefaultProblemSize(m_grid_dims_default[0] *
                            m_grid_dims_default[1] *
                            m_grid_dims_default[2]);
      setDefaultReps(50);

      var cbrt_run_size: real = cbrt(getTargetProblemSize());

      m_grid_dims[0] = cbrt_run_size:Index_type;
      m_grid_dims[1] = cbrt_run_size:Index_type;
      m_grid_dims[2] = cbrt_run_size:Index_type;
      m_halo_width   = m_halo_width_default;
      m_num_vars     = m_num_vars_default;

      m_grid_plus_halo_dims[0] = m_grid_dims[0] + 2*m_halo_width;
      m_grid_plus_halo_dims[1] = m_grid_dims[1] + 2*m_halo_width;
      m_grid_plus_halo_dims[2] = m_grid_dims[2] + 2*m_halo_width;
      m_var_size = m_grid_plus_halo_dims[0] *
                   m_grid_plus_halo_dims[1] *
                   m_grid_plus_halo_dims[2];

      setActualProblemSize(m_grid_dims[0] * m_grid_dims[1] * m_grid_dims[1]);

      setItsPerRep(m_num_vars * (m_var_size - getActualProblemSize()));
      setKernelsPerRep(2);
      setBytesPerRep((0*sizeof(Int_type)  + 1*sizeof(Int_type) ) * getItsPerRep() +
                     (1*sizeof(Real_type) + 1*sizeof(Real_type)) * getItsPerRep() +
                     (0*sizeof(Int_type)  + 1*sizeof(Int_type) ) * getItsPerRep() +
                     (1*sizeof(Real_type) + 1*sizeof(Real_type)) * getItsPerRep());
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Workgroup);

      setVariantDefined(VariantID.Base_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var vars = for v in 0..<m_num_vars do
        allocAndInitData(Real_type, m_var_size, vid);

      for v in 0..<m_num_vars do
        for i in 0..<m_var_size do
          vars[v][i] = i + v;

      var (pack_index_lists, pack_index_list_lengths) =
        create_pack_lists(m_halo_width, m_grid_dims, s_num_neighbors, vid);

      var (unpack_index_lists, unpack_index_list_lengths) =
        create_unpack_lists(m_halo_width, m_grid_dims, s_num_neighbors, vid);

      var buffers = for l in 0..<s_num_neighbors do
        allocAndInitData(Real_type, m_num_vars*pack_index_list_lengths[l], vid);

      const run_reps = getRunReps();

      var num_neighbors: Index_type = s_num_neighbors;
      var num_vars: Index_type = m_num_vars;

      // run
      select vid {

        when VariantID.Base_Chpl {

          var pack_ptr_holders: [0..<num_neighbors*num_vars] ptr_holder;
          var pack_lens: [0..<num_neighbors*num_vars] Index_type;
          var unpack_ptr_holders: [0..<num_neighbors*num_vars] ptr_holder;
          var unpack_lens: [0..<num_neighbors*num_vars] Index_type;

          startTimer();

          for irep in 0..<run_reps {

            var pack_index: Index_type = 0;

            for l in 0..<num_neighbors {
              var _buffer_b: Index_type = l;  // buffers[l];
              var _buffer_off: Index_type = 0;
              var _list: Index_type = l;  // pack_index_lists[l];
              var len: Index_type = pack_index_list_lengths[l];
              for v in 0..<num_vars {
                var _var = v;  // vars[v];
                pack_ptr_holders[pack_index] = new ptr_holder(_buffer_b, _buffer_off, _list, _var);
                pack_lens[pack_index]        = len;
                pack_index += 1;
                _buffer_off += len;
              }
            }
            for j in 0..<pack_index {
              ref buffer_ = reindex(buffers[pack_ptr_holders[j]._buffer_b][pack_ptr_holders[j]._buffer_off..], 0..);
              ref list_   = pack_index_lists[pack_ptr_holders[j]._list];
              ref var_    = vars[pack_ptr_holders[j]._var];
              var len: Index_type = pack_lens[j];
              for i in 0..<len do buffer_[i] = var_[list_[i]];
            }

            var unpack_index: Index_type = 0;

            for l in 0..<num_neighbors {
              var _buffer_b: Index_type = l;  // buffers[l];
              var _buffer_off: Index_type = 0;
              var _list: Index_type = l;  // unpack_index_lists[l];
              var len: Index_type = unpack_index_list_lengths[l];
              for v in 0..<num_vars {
                var _var = v;  // vars[v];
                unpack_ptr_holders[unpack_index] = new ptr_holder(_buffer_b, _buffer_off, _list, _var);
                unpack_lens[unpack_index]        = len;
                unpack_index += 1;
                _buffer_off += len;
              }
            }
            for j in 0..<unpack_index {
              ref buffer_ = reindex(buffers[unpack_ptr_holders[j]._buffer_b][unpack_ptr_holders[j]._buffer_off..], 0..);
              ref list_   = unpack_index_lists[unpack_ptr_holders[j]._list];
              ref var_    = vars[unpack_ptr_holders[j]._var];
              var len: Index_type = unpack_lens[j];
              for i in 0..<len do var_[list_[i]] = buffer_[i];
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      for v in vars do checksum[vid] += calcChecksum(v, m_var_size);
    }
  }

}
