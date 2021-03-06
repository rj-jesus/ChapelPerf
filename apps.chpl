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

module apps {
  private use IO;

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
      NPNL = 2;
      NPNR = 1;

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
      }
      else if ndims == 3 {
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
      ref fy1 = reindex(fy4[1..],           0..);
      ref fy2 = reindex(fy1[m_domain.jp..], 0..);
      ref fy3 = reindex(fy4[m_domain.jp..], 0..);

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
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var  Basis = allocAndInitDataConst(Real_type, (DPA_Q1D*DPA_D1D):Int_type,                  1.0:Real_type, vid);
      var dBasis = allocAndInitDataConst(Real_type, (DPA_Q1D*DPA_D1D):Int_type,                  1.0:Real_type, vid);
      var      D = allocAndInitDataConst(Real_type, (DPA_Q1D*DPA_Q1D*DPA_Q1D*SYM*m_NE):Int_type, 1.0:Real_type, vid);
      var      X = allocAndInitDataConst(Real_type, (DPA_D1D*DPA_D1D*DPA_D1D*m_NE):Int_type,     1.0:Real_type, vid);
      var      Y = allocAndInitDataConst(Real_type, (DPA_D1D*DPA_D1D*DPA_D1D*m_NE):Int_type,     0.0:Real_type, vid);

      const run_reps = getRunReps();

      var NE = m_NE: Index_type;
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

              // DIFFUSION3DPA_0_CPU
              const MQ1: int = DPA_Q1D;
              const MD1: int = DPA_D1D;
              const MDQ: int = if MQ1 > MD1 then MQ1 else MD1;
              var sBG: [0..<MDQ, 0..<MDQ] real;
              ref  B = sBG;  // (double (*)[MD1]) sBG
              ref  G = sBG;  // (double (*)[MD1]) sBG
              ref Bt = sBG;  // (double (*)[MQ1]) sBG;
              ref Gt = sBG;  // (double (*)[MQ1]) sBG;
              var sm0: [0..<3][0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              var sm1: [0..<3][0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              ref  s_X = sm0[2];  // (double (*)[MD1][MD1]) (sm0+2)
              ref DDQ0 = sm0[0];  // (double (*)[MD1][MQ1]) (sm0+0)
              ref DDQ1 = sm0[1];  // (double (*)[MD1][MQ1]) (sm0+1)
              ref DQQ0 = sm1[0];  // (double (*)[MQ1][MQ1]) (sm1+0)
              ref DQQ1 = sm1[1];  // (double (*)[MQ1][MQ1]) (sm1+1)
              ref DQQ2 = sm1[2];  // (double (*)[MQ1][MQ1]) (sm1+2)
              ref QQQ0 = sm0[0];  // (double (*)[MQ1][MQ1]) (sm0+0)
              ref QQQ1 = sm0[1];  // (double (*)[MQ1][MQ1]) (sm0+1)
              ref QQQ2 = sm0[2];  // (double (*)[MQ1][MQ1]) (sm0+2)
              ref QQD0 = sm1[0];  // (double (*)[MQ1][MD1]) (sm1+0)
              ref QQD1 = sm1[1];  // (double (*)[MQ1][MD1]) (sm1+1)
              ref QQD2 = sm1[2];  // (double (*)[MQ1][MD1]) (sm1+2)
              ref QDD0 = sm0[0];  // (double (*)[MD1][MD1]) (sm0+0)
              ref QDD1 = sm0[1];  // (double (*)[MD1][MD1]) (sm0+1)
              ref QDD2 = sm0[2];  // (double (*)[MD1][MD1]) (sm0+2)

              for dz in 0..<DPA_D1D {
                for dy in 0..<DPA_D1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_1
                    s_X[dz, dy, dx] = dpaX_(dx, dy, dz, e);
                  }
                }
              }

              for dy in 0..<DPA_D1D {
                for qx in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_2
                  const i = qi(qx, dy, DPA_Q1D):int;
                  const j = dj(qx, dy, DPA_D1D):int;
                  const k = qk(qx, dy, DPA_Q1D):int;
                  const l = dl(qx, dy, DPA_D1D):int;
                  B[i, j] = b(qx, dy);
                  G[k, l] = g(qx, dy) * sign(qx, dy);
                }
              }

              for dz in 0..<DPA_D1D {
                for dy in 0..<DPA_D1D {
                  for qx in 0..<DPA_Q1D {
                    // DIFFUSION3DPA_3
                    var u = 0.0:real, v = 0.0:real;
                    for dx in 0..<DPA_D1D {
                      const i = qi(qx, dx, DPA_Q1D):int;
                      const j = dj(qx, dx, DPA_D1D):int;
                      const k = qk(qx, dx, DPA_Q1D):int;
                      const l = dl(qx, dx, DPA_D1D):int;
                      const s = sign(qx, dx):real;
                      const coords = s_X[dz, dy, dx]:real;
                      u += coords * B[i, j];
                      v += coords * G[k, l] * s;
                    }
                    DDQ0[dz, dy, qx] = u;
                    DDQ1[dz, dy, qx] = v;
                  }
                }
              }

              for dz in 0..<DPA_D1D {
                for qy in 0..<DPA_Q1D {
                  for qx in 0..<DPA_Q1D {
                    // DIFFUSION3DPA_4
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for dy in 0..<DPA_D1D {
                      const i = qi(qy, dy, DPA_Q1D):int;
                      const j = dj(qy, dy, DPA_D1D):int;
                      const k = qk(qy, dy, DPA_Q1D):int;
                      const l = dl(qy, dy, DPA_D1D):int;
                      const s = sign(qy, dy):real;
                      u += DDQ1[dz, dy, qx] * B[i, j];
                      v += DDQ0[dz, dy, qx] * G[k, l] * s;
                      w += DDQ0[dz, dy, qx] * B[i, j];
                    }
                    DQQ0[dz, qy, qx] = u;
                    DQQ1[dz, qy, qx] = v;
                    DQQ2[dz, qy, qx] = w;
                  }
                }
              }

              for qz in 0..<DPA_Q1D {
                for qy in 0..<DPA_Q1D {
                  for qx in 0..<DPA_Q1D {
                    // DIFFUSION3DPA_5
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for dz in 0..<DPA_D1D {
                      const i = qi(qz, dz, DPA_Q1D):int;
                      const j = dj(qz, dz, DPA_D1D):int;
                      const k = qk(qz, dz, DPA_Q1D):int;
                      const l = dl(qz, dz, DPA_D1D):int;
                      const s = sign(qz, dz):real;
                      u += DQQ0[dz, qy, qx] * B[i, j];
                      v += DQQ1[dz, qy, qx] * B[i, j];
                      w += DQQ2[dz, qy, qx] * G[k, l] * s;
                    }
                    const O11 = d(qx, qy, qz, 0, e):real;
                    const O12 = d(qx, qy, qz, 1, e):real;
                    const O13 = d(qx, qy, qz, 2, e):real;
                    const O21 = (if symmetric then             O12 else d(qx,qy,qz,3,e)):real;
                    const O22 = (if symmetric then d(qx,qy,qz,3,e) else d(qx,qy,qz,4,e)):real;
                    const O23 = (if symmetric then d(qx,qy,qz,4,e) else d(qx,qy,qz,5,e)):real;
                    const O31 = (if symmetric then             O13 else d(qx,qy,qz,6,e)):real;
                    const O32 = (if symmetric then             O23 else d(qx,qy,qz,7,e)):real;
                    const O33 = (if symmetric then d(qx,qy,qz,5,e) else d(qx,qy,qz,8,e)):real;
                    const gX = u:real;
                    const gY = v:real;
                    const gZ = w:real;
                    QQQ0[qz, qy, qx] = (O11*gX) + (O12*gY) + (O13*gZ);
                    QQQ1[qz, qy, qx] = (O21*gX) + (O22*gY) + (O23*gZ);
                    QQQ2[qz, qy, qx] = (O31*gX) + (O32*gY) + (O33*gZ);
                  }
                }
              }

              for d in 0..<DPA_D1D {
                for q in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_6
                  const i = qi(q,d,DPA_Q1D):int;
                  const j = dj(q,d,DPA_D1D):int;
                  const k = qk(q,d,DPA_Q1D):int;
                  const l = dl(q,d,DPA_D1D):int;
                  Bt[j, i] = b(q,d);
                  Gt[l, k] = g(q,d) * sign(q,d);
                }
              }

              for qz in 0..<DPA_Q1D {
                for qy in 0..<DPA_Q1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_7
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for qx in 0..<DPA_Q1D {
                      const i = qi(qx,dx,DPA_Q1D):int;
                      const j = dj(qx,dx,DPA_D1D):int;
                      const k = qk(qx,dx,DPA_Q1D):int;
                      const l = dl(qx,dx,DPA_D1D):int;
                      const s = sign(qx,dx):real;
                      u += QQQ0[qz, qy, qx] * Gt[l, k] * s;
                      v += QQQ1[qz, qy, qx] * Bt[j, i];
                      w += QQQ2[qz, qy, qx] * Bt[j, i];
                    }
                    QQD0[qz, qy, dx] = u;
                    QQD1[qz, qy, dx] = v;
                    QQD2[qz, qy, dx] = w;
                  }
                }
              }

              for qz in 0..<DPA_Q1D {
                for dy in 0..<DPA_D1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_8
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for qy in 0..<DPA_Q1D {
                      const i = qi(qy,dy,DPA_Q1D):int;
                      const j = dj(qy,dy,DPA_D1D):int;
                      const k = qk(qy,dy,DPA_Q1D):int;
                      const l = dl(qy,dy,DPA_D1D):int;
                      const s = sign(qy,dy):real;
                      u += QQD0[qz, qy, dx] * Bt[j, i];
                      v += QQD1[qz, qy, dx] * Gt[l, k] * s;
                      w += QQD2[qz, qy, dx] * Bt[j, i];
                    }
                    QDD0[qz, dy, dx] = u;
                    QDD1[qz, dy, dx] = v;
                    QDD2[qz, dy, dx] = w;
                  }
                }
              }

              for dz in 0..<DPA_D1D {
                for dy in 0..<DPA_D1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_9
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for qz in 0..<DPA_Q1D {
                      const i = qi(qz,dz,DPA_Q1D):int;
                      const j = dj(qz,dz,DPA_D1D):int;
                      const k = qk(qz,dz,DPA_Q1D):int;
                      const l = dl(qz,dz,DPA_D1D):int;
                      const s = sign(qz,dz):real;
                      u += QDD0[qz, dy, dx] * Bt[j, i];
                      v += QDD1[qz, dy, dx] * Bt[j, i];
                      w += QDD2[qz, dy, dx] * Gt[l, k] * s;
                    }
                    dpaY_(dx,dy,dz,e) += (u + v + w);
                  }
                }
              }

            }  // element loop
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..<run_reps {
            forall e in 0..<NE {

              // DIFFUSION3DPA_0_CPU
              const MQ1: int = DPA_Q1D;
              const MD1: int = DPA_D1D;
              const MDQ: int = if MQ1 > MD1 then MQ1 else MD1;
              var sBG: [0..<MDQ, 0..<MDQ] real;
              ref  B = sBG;  // (double (*)[MD1]) sBG
              ref  G = sBG;  // (double (*)[MD1]) sBG
              ref Bt = sBG;  // (double (*)[MQ1]) sBG;
              ref Gt = sBG;  // (double (*)[MQ1]) sBG;
              var sm0: [0..<3][0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              var sm1: [0..<3][0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              ref  s_X = sm0[2];  // (double (*)[MD1][MD1]) (sm0+2)
              ref DDQ0 = sm0[0];  // (double (*)[MD1][MQ1]) (sm0+0)
              ref DDQ1 = sm0[1];  // (double (*)[MD1][MQ1]) (sm0+1)
              ref DQQ0 = sm1[0];  // (double (*)[MQ1][MQ1]) (sm1+0)
              ref DQQ1 = sm1[1];  // (double (*)[MQ1][MQ1]) (sm1+1)
              ref DQQ2 = sm1[2];  // (double (*)[MQ1][MQ1]) (sm1+2)
              ref QQQ0 = sm0[0];  // (double (*)[MQ1][MQ1]) (sm0+0)
              ref QQQ1 = sm0[1];  // (double (*)[MQ1][MQ1]) (sm0+1)
              ref QQQ2 = sm0[2];  // (double (*)[MQ1][MQ1]) (sm0+2)
              ref QQD0 = sm1[0];  // (double (*)[MQ1][MD1]) (sm1+0)
              ref QQD1 = sm1[1];  // (double (*)[MQ1][MD1]) (sm1+1)
              ref QQD2 = sm1[2];  // (double (*)[MQ1][MD1]) (sm1+2)
              ref QDD0 = sm0[0];  // (double (*)[MD1][MD1]) (sm0+0)
              ref QDD1 = sm0[1];  // (double (*)[MD1][MD1]) (sm0+1)
              ref QDD2 = sm0[2];  // (double (*)[MD1][MD1]) (sm0+2)

              for dz in 0..<DPA_D1D {
                for dy in 0..<DPA_D1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_1
                    s_X[dz, dy, dx] = dpaX_(dx, dy, dz, e);
                  }
                }
              }

              for dy in 0..<DPA_D1D {
                for qx in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_2
                  const i = qi(qx, dy, DPA_Q1D):int;
                  const j = dj(qx, dy, DPA_D1D):int;
                  const k = qk(qx, dy, DPA_Q1D):int;
                  const l = dl(qx, dy, DPA_D1D):int;
                  B[i, j] = b(qx, dy);
                  G[k, l] = g(qx, dy) * sign(qx, dy);
                }
              }

              for dz in 0..<DPA_D1D {
                for dy in 0..<DPA_D1D {
                  for qx in 0..<DPA_Q1D {
                    // DIFFUSION3DPA_3
                    var u = 0.0:real, v = 0.0:real;
                    for dx in 0..<DPA_D1D {
                      const i = qi(qx, dx, DPA_Q1D):int;
                      const j = dj(qx, dx, DPA_D1D):int;
                      const k = qk(qx, dx, DPA_Q1D):int;
                      const l = dl(qx, dx, DPA_D1D):int;
                      const s = sign(qx, dx):real;
                      const coords = s_X[dz, dy, dx]:real;
                      u += coords * B[i, j];
                      v += coords * G[k, l] * s;
                    }
                    DDQ0[dz, dy, qx] = u;
                    DDQ1[dz, dy, qx] = v;
                  }
                }
              }

              for dz in 0..<DPA_D1D {
                for qy in 0..<DPA_Q1D {
                  for qx in 0..<DPA_Q1D {
                    // DIFFUSION3DPA_4
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for dy in 0..<DPA_D1D {
                      const i = qi(qy, dy, DPA_Q1D):int;
                      const j = dj(qy, dy, DPA_D1D):int;
                      const k = qk(qy, dy, DPA_Q1D):int;
                      const l = dl(qy, dy, DPA_D1D):int;
                      const s = sign(qy, dy):real;
                      u += DDQ1[dz, dy, qx] * B[i, j];
                      v += DDQ0[dz, dy, qx] * G[k, l] * s;
                      w += DDQ0[dz, dy, qx] * B[i, j];
                    }
                    DQQ0[dz, qy, qx] = u;
                    DQQ1[dz, qy, qx] = v;
                    DQQ2[dz, qy, qx] = w;
                  }
                }
              }

              for qz in 0..<DPA_Q1D {
                for qy in 0..<DPA_Q1D {
                  for qx in 0..<DPA_Q1D {
                    // DIFFUSION3DPA_5
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for dz in 0..<DPA_D1D {
                      const i = qi(qz, dz, DPA_Q1D):int;
                      const j = dj(qz, dz, DPA_D1D):int;
                      const k = qk(qz, dz, DPA_Q1D):int;
                      const l = dl(qz, dz, DPA_D1D):int;
                      const s = sign(qz, dz):real;
                      u += DQQ0[dz, qy, qx] * B[i, j];
                      v += DQQ1[dz, qy, qx] * B[i, j];
                      w += DQQ2[dz, qy, qx] * G[k, l] * s;
                    }
                    const O11 = d(qx, qy, qz, 0, e):real;
                    const O12 = d(qx, qy, qz, 1, e):real;
                    const O13 = d(qx, qy, qz, 2, e):real;
                    const O21 = (if symmetric then             O12 else d(qx,qy,qz,3,e)):real;
                    const O22 = (if symmetric then d(qx,qy,qz,3,e) else d(qx,qy,qz,4,e)):real;
                    const O23 = (if symmetric then d(qx,qy,qz,4,e) else d(qx,qy,qz,5,e)):real;
                    const O31 = (if symmetric then             O13 else d(qx,qy,qz,6,e)):real;
                    const O32 = (if symmetric then             O23 else d(qx,qy,qz,7,e)):real;
                    const O33 = (if symmetric then d(qx,qy,qz,5,e) else d(qx,qy,qz,8,e)):real;
                    const gX = u:real;
                    const gY = v:real;
                    const gZ = w:real;
                    QQQ0[qz, qy, qx] = (O11*gX) + (O12*gY) + (O13*gZ);
                    QQQ1[qz, qy, qx] = (O21*gX) + (O22*gY) + (O23*gZ);
                    QQQ2[qz, qy, qx] = (O31*gX) + (O32*gY) + (O33*gZ);
                  }
                }
              }

              for d in 0..<DPA_D1D {
                for q in 0..<DPA_Q1D {
                  // DIFFUSION3DPA_6
                  const i = qi(q,d,DPA_Q1D):int;
                  const j = dj(q,d,DPA_D1D):int;
                  const k = qk(q,d,DPA_Q1D):int;
                  const l = dl(q,d,DPA_D1D):int;
                  Bt[j, i] = b(q,d);
                  Gt[l, k] = g(q,d) * sign(q,d);
                }
              }

              for qz in 0..<DPA_Q1D {
                for qy in 0..<DPA_Q1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_7
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for qx in 0..<DPA_Q1D {
                      const i = qi(qx,dx,DPA_Q1D):int;
                      const j = dj(qx,dx,DPA_D1D):int;
                      const k = qk(qx,dx,DPA_Q1D):int;
                      const l = dl(qx,dx,DPA_D1D):int;
                      const s = sign(qx,dx):real;
                      u += QQQ0[qz, qy, qx] * Gt[l, k] * s;
                      v += QQQ1[qz, qy, qx] * Bt[j, i];
                      w += QQQ2[qz, qy, qx] * Bt[j, i];
                    }
                    QQD0[qz, qy, dx] = u;
                    QQD1[qz, qy, dx] = v;
                    QQD2[qz, qy, dx] = w;
                  }
                }
              }

              for qz in 0..<DPA_Q1D {
                for dy in 0..<DPA_D1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_8
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for qy in 0..<DPA_Q1D {
                      const i = qi(qy,dy,DPA_Q1D):int;
                      const j = dj(qy,dy,DPA_D1D):int;
                      const k = qk(qy,dy,DPA_Q1D):int;
                      const l = dl(qy,dy,DPA_D1D):int;
                      const s = sign(qy,dy):real;
                      u += QQD0[qz, qy, dx] * Bt[j, i];
                      v += QQD1[qz, qy, dx] * Gt[l, k] * s;
                      w += QQD2[qz, qy, dx] * Bt[j, i];
                    }
                    QDD0[qz, dy, dx] = u;
                    QDD1[qz, dy, dx] = v;
                    QDD2[qz, dy, dx] = w;
                  }
                }
              }

              for dz in 0..<DPA_D1D {
                for dy in 0..<DPA_D1D {
                  for dx in 0..<DPA_D1D {
                    // DIFFUSION3DPA_9
                    var u = 0.0:real, v = 0.0:real, w = 0.0:real;
                    for qz in 0..<DPA_Q1D {
                      const i = qi(qz,dz,DPA_Q1D):int;
                      const j = dj(qz,dz,DPA_D1D):int;
                      const k = qk(qz,dz,DPA_Q1D):int;
                      const l = dl(qz,dz,DPA_D1D):int;
                      const s = sign(qz,dz):real;
                      u += QDD0[qz, dy, dx] * Bt[j, i];
                      v += QDD1[qz, dy, dx] * Bt[j, i];
                      w += QDD2[qz, dy, dx] * Gt[l, k] * s;
                    }
                    dpaY_(dx,dy,dz,e) += (u + v + w);
                  }
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

      var coeff: [0..<FIR_COEFFLEN] Real_type = [  3.0, -1.0, -1.0, -1.0,
                                                  -1.0,  3.0, -1.0, -1.0,
                                                  -1.0, -1.0,  3.0, -1.0,
                                                  -1.0, -1.0, -1.0,  3.0, ];

      const coefflen = m_coefflen;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for 0..<run_reps {

            for i in ibegin..<iend {
              var sum: Real_type = 0.0;
              for j in 0..<coefflen do
                sum += coeff[j]*in_[i+j];
              out_[i] = sum;
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for 0..<run_reps {

            forall i in ibegin..<iend do
              out_[i] = + reduce (coeff*in_[i..]);

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
    pack_index_list_extents[ 0] = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[ 1] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[ 2] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[ 3] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[ 4] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             halo_width  , halo_width   + halo_width);
    pack_index_list_extents[ 5] = new Extent(halo_width  , grid_dims[0] + halo_width,
                                             halo_width  , grid_dims[1] + halo_width,
                                             grid_dims[2], grid_dims[2] + halo_width);

    // edges
    pack_index_list_extents[ 6] = new Extent(halo_width  , halo_width   + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[ 7] = new Extent(halo_width  , halo_width   + halo_width,
                                             grid_dims[1], grid_dims[1] + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[ 8] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
                                             halo_width  , halo_width   + halo_width,
                                             halo_width  , grid_dims[2] + halo_width);
    pack_index_list_extents[ 9] = new Extent(grid_dims[0], grid_dims[0] + halo_width,
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

    var pack_index_list_lengths = for extent in pack_index_list_extents do
                                    (extent.i_max - extent.i_min) *
                                    (extent.j_max - extent.j_min) *
                                    (extent.k_max - extent.k_min);

    var pack_index_lists = for (l, len) in zip(0..<num_neighbors, pack_index_list_lengths) do
                             allocAndInitData(Int_type, len, vid);

    for (extent, pack_list) in zip(pack_index_list_extents, pack_index_lists) {
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

    return pack_index_lists;
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

    var unpack_index_list_lengths = for extent in unpack_index_list_extents do
                                      (extent.i_max - extent.i_min) *
                                      (extent.j_max - extent.j_min) *
                                      (extent.k_max - extent.k_min);

    var unpack_index_lists = for (l, len) in zip(0..<num_neighbors, unpack_index_list_lengths) do
                               allocAndInitData(Int_type, len, vid);

    for (extent, unpack_list) in zip(unpack_index_list_extents, unpack_index_lists) {
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

    return unpack_index_lists;
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
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var vars = for v in 0..<m_num_vars do
        allocAndInitData(Real_type, m_var_size, vid);

      for (v, i) in {0..<m_num_vars, 0..<m_var_size} do
        vars[v][i] = i + v;

      var pack_index_lists = create_pack_lists(m_halo_width, m_grid_dims, s_num_neighbors, vid);

      var unpack_index_lists = create_unpack_lists(m_halo_width, m_grid_dims, s_num_neighbors, vid);

      var buffers = for (l, pack_list) in zip(0..<s_num_neighbors, pack_index_lists) do
                      allocAndInitData(Real_type, m_num_vars*pack_list.size, vid);

      const run_reps = getRunReps();

      var num_neighbors: Index_type = s_num_neighbors;
      var num_vars: Index_type = m_num_vars;

      // run
      select vid {

        when VariantID.Base_Chpl {

          var pack_ptr_holders: [0..<num_neighbors*num_vars] ptr_holder;
          var unpack_ptr_holders: [0..<num_neighbors*num_vars] ptr_holder;

          startTimer();

          for irep in 0..<run_reps {

            var pack_index: Index_type = 0;

            for l in 0..<num_neighbors {
              var _buffer_b: Index_type = l;  // buffers[l];
              var _buffer_off: Index_type = 0;
              var _list: Index_type = l;  // pack_index_lists[l];
              var len: Index_type = pack_index_lists[l].size;
              for v in 0..<num_vars {
                var _var = v;  // vars[v];
                pack_ptr_holders[pack_index] = new ptr_holder(_buffer_b, _buffer_off, _list, _var);
                pack_index += 1;
                _buffer_off += len;
              }
            }

            for pack_ptr_holder in pack_ptr_holders {
              ref buffer_ = reindex(buffers[pack_ptr_holder._buffer_b][pack_ptr_holder._buffer_off..], 0..);
              ref list_ = pack_index_lists[pack_ptr_holder._list];
              ref var_ = vars[pack_ptr_holder._var];
              for (i, li_) in zip(0.., list_) do buffer_[i] = var_[li_];
            }

            var unpack_index: Index_type = 0;

            for l in 0..<num_neighbors {
              var _buffer_b: Index_type = l;  // buffers[l];
              var _buffer_off: Index_type = 0;
              var _list: Index_type = l;  // unpack_index_lists[l];
              var len: Index_type = unpack_index_lists[l].size;
              for v in 0..<num_vars {
                var _var = v;  // vars[v];
                unpack_ptr_holders[unpack_index] = new ptr_holder(_buffer_b, _buffer_off, _list, _var);
                unpack_index += 1;
                _buffer_off += len;
              }
            }

            for unpack_ptr_holder in unpack_ptr_holders {
              ref buffer_ = reindex(buffers[unpack_ptr_holder._buffer_b][unpack_ptr_holder._buffer_off..], 0..);
              ref list_ = unpack_index_lists[unpack_ptr_holder._list];
              ref var_ = vars[unpack_ptr_holder._var];
              for (i, li_) in zip(0.., list_) do var_[li_] = buffer_[i];
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {

          var pack_ptr_holders: [0..<num_neighbors*num_vars] ptr_holder;
          var unpack_ptr_holders: [0..<num_neighbors*num_vars] ptr_holder;

          startTimer();

          for irep in 0..<run_reps {

            var pack_index: Index_type = 0;

            for l in 0..<num_neighbors {
              var _buffer_b: Index_type = l;  // buffers[l];
              var _buffer_off: Index_type = 0;
              var _list: Index_type = l;  // pack_index_lists[l];
              var len: Index_type = pack_index_lists[l].size;
              for v in 0..<num_vars {
                var _var = v;  // vars[v];
                pack_ptr_holders[pack_index] = new ptr_holder(_buffer_b, _buffer_off, _list, _var);
                pack_index += 1;
                _buffer_off += len;
              }
            }

            forall pack_ptr_holder in pack_ptr_holders {
              ref buffer_ = reindex(buffers[pack_ptr_holder._buffer_b][pack_ptr_holder._buffer_off..], 0..);
              ref list_ = pack_index_lists[pack_ptr_holder._list];
              ref var_ = vars[pack_ptr_holder._var];
              for (i, li_) in zip(0.., list_) do buffer_[i] = var_[li_];
            }

            var unpack_index: Index_type = 0;

            for l in 0..<num_neighbors {
              var _buffer_b: Index_type = l;  // buffers[l];
              var _buffer_off: Index_type = 0;
              var _list: Index_type = l;  // unpack_index_lists[l];
              var len: Index_type = unpack_index_lists[l].size;
              for v in 0..<num_vars {
                var _var = v;  // vars[v];
                unpack_ptr_holders[unpack_index] = new ptr_holder(_buffer_b, _buffer_off, _list, _var);
                unpack_index += 1;
                _buffer_off += len;
              }
            }

            forall unpack_ptr_holder in unpack_ptr_holders {
              ref buffer_ = reindex(buffers[unpack_ptr_holder._buffer_b][unpack_ptr_holder._buffer_off..], 0..);
              ref list_ = unpack_index_lists[unpack_ptr_holder._list];
              ref var_ = vars[unpack_ptr_holder._var];
              for (i, li_) in zip(0.., list_) do var_[li_] = buffer_[i];
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

  class HALOEXCHANGE: KernelBase {

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
      super.init(KernelID.Apps_HALOEXCHANGE);

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
      setKernelsPerRep(2 * s_num_neighbors * m_num_vars);
      setBytesPerRep((0*sizeof(Int_type)  + 1*sizeof(Int_type) ) * getItsPerRep() +
                     (1*sizeof(Real_type) + 1*sizeof(Real_type)) * getItsPerRep() +
                     (0*sizeof(Int_type)  + 1*sizeof(Int_type) ) * getItsPerRep() +
                     (1*sizeof(Real_type) + 1*sizeof(Real_type)) * getItsPerRep());
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var vars = for v in 0..<m_num_vars do
        allocAndInitData(Real_type, m_var_size, vid);

      for v in 0..<m_num_vars do
        for i in 0..<m_var_size do
          vars[v][i] = i + v;

      var pack_index_lists = create_pack_lists(m_halo_width, m_grid_dims, s_num_neighbors, vid);

      var unpack_index_lists = create_unpack_lists(m_halo_width, m_grid_dims, s_num_neighbors, vid);

      var buffers = for (l, list) in zip(0..<s_num_neighbors, pack_index_lists) do
                      allocAndInitData(Real_type, m_num_vars*list.size, vid);

      const run_reps = getRunReps();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            for (buffer, list) in zip(buffers, pack_index_lists) {
              var k = 0;

              for vi in vars {
                for li in list {
                  buffer[k] = vi[li];
                  k += 1;
                }
              }
            }

            for (buffer, list) in zip(buffers, unpack_index_lists) {
              var k = 0;

              for vi in vars {
                for li in list {
                  vi[li] = buffer[k];
                  k += 1;
                }
              }
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            for (buffer, list) in zip(buffers, pack_index_lists) {
              var k = 0;

              for vi in vars {
                forall (li, bk) in zip(list, buffer[k..]) do
                  bk = vi[li];
                k += list.size;
              }
            }

            for (buffer, list) in zip(buffers, unpack_index_lists) {
              var k = 0;

              for vi in vars {
                forall (li, bk) in zip(list, buffer[k..]) do
                  vi[li] = bk;
                k += list.size;
              }
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

  class LTIMES_NOVIEW: KernelBase {

    var m_num_d_default: Index_type;
    var m_num_z_default: Index_type;
    var m_num_g_default: Index_type;
    var m_num_m_default: Index_type;

    var m_num_z: Index_type;
    var m_num_g: Index_type;
    var m_num_m: Index_type;
    var m_num_d: Index_type;

    var m_philen: Index_type;
    var m_elllen: Index_type;
    var m_psilen: Index_type;

    proc init() {
      super.init(KernelID.Apps_LTIMES_NOVIEW);

      m_num_d_default =  64;
      m_num_z_default = 488;
      m_num_g_default =  32;
      m_num_m_default =  25;

      setDefaultProblemSize(m_num_d_default * m_num_g_default * m_num_z_default);
      setDefaultReps(50);

      m_num_z = max(getTargetProblemSize()/(m_num_d_default * m_num_g_default),
                    1:Index_type);
      m_num_g = m_num_g_default;
      m_num_m = m_num_m_default;
      m_num_d = m_num_d_default;

      m_philen = m_num_m * m_num_g * m_num_z;
      m_elllen = m_num_d * m_num_m;
      m_psilen = m_num_d * m_num_g * m_num_z;

      setActualProblemSize(m_psilen);

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      // using total data size instead of writes and reads
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) * m_philen +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_elllen +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_psilen);
      setFLOPsPerRep(2 * m_num_z * m_num_g * m_num_m * m_num_d);

      checksum_scale_factor = 0.001 *
        (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var phidat = allocAndInitDataConst(Real_type, m_philen:Int_type, 0.0:Real_type, vid);
      var elldat = allocAndInitData(Real_type, m_elllen:Int_type, vid);
      var psidat = allocAndInitData(Real_type, m_psilen:Int_type, vid);

      const run_reps = getRunReps();

      var num_d: Index_type = m_num_d;
      var num_z: Index_type = m_num_z;
      var num_g: Index_type = m_num_g;
      var num_m: Index_type = m_num_m;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            for z in 0..< num_z do
              for g in 0..<num_g do
                for m in 0..<num_m do
                  for d in 0..<num_d do
                    phidat[m+ (g * num_m) + (z * num_m * num_g)] +=
                      elldat[d+ (m * num_d)] *
                      psidat[d+ (g * num_d) + (z * num_d * num_g)];

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            forall z in 0..< num_z do
              for g in 0..<num_g do
                for m in 0..<num_m do
                  for d in 0..<num_d do
                    phidat[m+ (g * num_m) + (z * num_m * num_g)] +=
                      elldat[d+ (m * num_d)] *
                      psidat[d+ (g * num_d) + (z * num_d * num_g)];

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(phidat, m_philen, checksum_scale_factor:Real_type);
    }
  }

  class LTIMES: KernelBase {

    var m_num_d_default: Index_type;
    var m_num_z_default: Index_type;
    var m_num_g_default: Index_type;
    var m_num_m_default: Index_type;

    var m_num_z: Index_type;
    var m_num_g: Index_type;
    var m_num_m: Index_type;
    var m_num_d: Index_type;

    var m_philen: Index_type;
    var m_elllen: Index_type;
    var m_psilen: Index_type;

    proc init() {
      super.init(KernelID.Apps_LTIMES);

      m_num_d_default =  64;
      m_num_z_default = 488;
      m_num_g_default =  32;
      m_num_m_default =  25;

      setDefaultProblemSize(m_num_d_default * m_num_g_default * m_num_z_default);
      setDefaultReps(50);

      m_num_z = max(getTargetProblemSize()/(m_num_d_default * m_num_g_default),
                    1:Index_type);
      m_num_g = m_num_g_default;
      m_num_m = m_num_m_default;
      m_num_d = m_num_d_default;

      m_philen = m_num_m * m_num_g * m_num_z;
      m_elllen = m_num_d * m_num_m;
      m_psilen = m_num_d * m_num_g * m_num_z;

      setActualProblemSize(m_psilen);

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      // using total data size instead of writes and reads
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) * m_philen +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_elllen +
                     (0*sizeof(Real_type) + 1*sizeof(Real_type)) * m_psilen);
      setFLOPsPerRep(2 * m_num_z * m_num_g * m_num_m * m_num_d);

      checksum_scale_factor = 0.001 *
        (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Kernel);
      setUsesFeature(FeatureID.View);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var phidat = allocAndInitDataConst(Real_type, m_philen:Int_type, 0.0:Real_type, vid);
      var elldat = allocAndInitData(Real_type, m_elllen:Int_type, vid);
      var psidat = allocAndInitData(Real_type, m_psilen:Int_type, vid);

      const run_reps = getRunReps();

      var num_d: Index_type = m_num_d;
      var num_z: Index_type = m_num_z;
      var num_g: Index_type = m_num_g;
      var num_m: Index_type = m_num_m;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            for z in 0..<num_z do
              for g in 0..<num_g do
                for m in 0..<num_m do
                  for d in 0..<num_d do
                    phidat[m+ (g * num_m) + (z * num_m * num_g)] +=
                      elldat[d+ (m * num_d)] *
                      psidat[d+ (g * num_d) + (z * num_d * num_g)];

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            forall z in 0..<num_z do
              for g in 0..<num_g do
                for m in 0..<num_m do
                  for d in 0..<num_d do
                    phidat[m+ (g * num_m) + (z * num_m * num_g)] +=
                      elldat[d+ (m * num_d)] *
                      psidat[d+ (g * num_d) + (z * num_d * num_g)];

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(phidat, m_philen, checksum_scale_factor:Real_type);
    }
  }

  class MASS3DPA: KernelBase {

    param MPA_D1D = 4;
    param MPA_Q1D = 5;

    var m_NE_default: Index_type;
    var m_NE: Index_type;

    proc init() {
      super.init(KernelID.Apps_MASS3DPA);

      m_NE_default = 8000;

      setDefaultProblemSize(m_NE_default*MPA_Q1D*MPA_Q1D*MPA_Q1D);
      setDefaultReps(50);

      m_NE = max(getTargetProblemSize()/(MPA_Q1D*MPA_Q1D*MPA_Q1D), 1:Index_type);

      setActualProblemSize(m_NE*MPA_Q1D*MPA_Q1D*MPA_Q1D);

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);

      setBytesPerRep(MPA_Q1D*MPA_D1D*sizeof(Real_type)  +
                     MPA_Q1D*MPA_D1D*sizeof(Real_type)  +
                     MPA_Q1D*MPA_Q1D*MPA_Q1D*m_NE*sizeof(Real_type) +
                     MPA_D1D*MPA_D1D*MPA_D1D*m_NE*sizeof(Real_type) +
                     MPA_D1D*MPA_D1D*MPA_D1D*m_NE*sizeof(Real_type));

      setFLOPsPerRep(m_NE * (2 * MPA_D1D * MPA_D1D * MPA_D1D * MPA_Q1D +
                             2 * MPA_D1D * MPA_D1D * MPA_Q1D * MPA_Q1D +
                             2 * MPA_D1D * MPA_Q1D * MPA_Q1D * MPA_Q1D + MPA_Q1D * MPA_Q1D * MPA_Q1D +
                             2 * MPA_Q1D * MPA_Q1D * MPA_Q1D * MPA_D1D +
                             2 * MPA_Q1D * MPA_Q1D * MPA_D1D * MPA_D1D +
                             2 * MPA_Q1D * MPA_D1D * MPA_D1D * MPA_D1D + MPA_D1D * MPA_D1D * MPA_D1D));
      setUsesFeature(FeatureID.Teams);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var B  = allocAndInitDataConst(Real_type, {0..<MPA_D1D, 0..<MPA_Q1D}, 1.0:Real_type, vid);
      var Bt = allocAndInitDataConst(Real_type, {0..<MPA_D1D, 0..<MPA_Q1D}, 1.0:Real_type, vid);
      var D  = allocAndInitDataConst(Real_type, {0..<m_NE, 0..<MPA_Q1D, 0..<MPA_Q1D, 0..<MPA_Q1D}, 1.0:Real_type, vid);
      var X  = allocAndInitDataConst(Real_type, {0..<m_NE, 0..<MPA_D1D, 0..<MPA_D1D, 0..<MPA_D1D}, 1.0:Real_type, vid);
      var Y  = allocAndInitDataConst(Real_type, {0..<m_NE, 0..<MPA_D1D, 0..<MPA_D1D, 0..<MPA_D1D}, 0.0:Real_type, vid);

      const run_reps = getRunReps();

      var NE: Index_type = m_NE;

      inline proc  B_(x, y) ref return  B[y, x];
      inline proc Bt_(x, y) ref return Bt[y, x];
      inline proc  X_(dx, dy, dz, e) ref return X[e, dz, dy, dx];
      inline proc  Y_(dx, dy, dz, e) ref return Y[e, dz, dy, dx];
      inline proc  D_(qx, qy, qz, e) ref return D[e, qz, qy, qx];

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for irep in 0..<run_reps {
            for e in 0..<NE {

              // MASS3DPA_0_CPU
              const MQ1 = MPA_Q1D;
              const MD1 = MPA_D1D;
              const MDQ = if MQ1 > MD1 then MQ1 else MD1;

              var sDQ: [0..<MDQ, 0..<MDQ] real;
              ref  Bsmem = sDQ;  // (MQ1, MD1));
              ref Btsmem = sDQ;  // (MD1, MQ1));

              var sm0: [0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              var sm1: [0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              ref Xsmem = sm0;  // (MD1, MD1, MD1)
              ref   DDQ = sm1;  // (MD1, MD1, MQ1)
              ref   DQQ = sm0;  // (MD1, MQ1, MQ1)
              ref   QQQ = sm1;  // (MQ1, MQ1, MQ1) 3?
              ref   QQD = sm0;  // (MQ1, MQ1, MD1)
              ref   QDD = sm1;  // (MQ1, MD1, MD1) 2?

              for dy in 0..<MPA_D1D {
                for dx in 0..<MPA_D1D{
                  // MASS3DPA_1
                  for dz in 0..<MPA_D1D do
                    Xsmem[dz, dy, dx] = X_(dx, dy, dz, e);
                }
                for dx in 0..<MPA_Q1D {
                  // MASS3DPA_2
                  Bsmem[dx, dy] = B_(dx, dy);
                }
              }

              for dy in 0..<MPA_D1D {
                for qx in 0..<MPA_Q1D {
                  // MASS3DPA_3
                  var u: [0..<MPA_D1D] real;
                  for dz in 0..<MPA_D1D do
                    u[dz] = 0;
                  for dx in 0..<MPA_D1D do
                    for dz in 0..<MPA_D1D do
                      u[dz] += Xsmem[dz, dy, dx] * Bsmem[qx, dx];
                  for dz in 0..<MPA_D1D do
                    DDQ[dz, dy, qx] = u[dz];
                }
              }

              for qy in 0..<MPA_Q1D {
                for qx in 0..<MPA_Q1D {
                  // MASS3DPA_4
                  var u: [0..<MPA_D1D] real;
                  for dz in 0..<MPA_D1D do
                    u[dz] = 0;
                  for dy in 0..<MPA_D1D do
                    for dz in 0..<MPA_D1D do
                      u[dz] += DDQ[dz, dy, qx] * Bsmem[qy, dy];
                  for dz in 0..<MPA_D1D do
                    DQQ[dz, qy, qx] = u[dz];
                }
              }

              for qy in 0..<MPA_Q1D {
                for qx in 0..<MPA_Q1D {
                  // MASS3DPA_5
                  var u: [0..<MPA_Q1D] real;
                  for dz in 0..<MPA_Q1D do
                    u[dz] = 0;
                  for dz in 0..<MPA_D1D do
                    for qz in 0..<MPA_Q1D do
                      u[qz] += DQQ[dz, qy, qx] * Bsmem[qz, dz];
                  for qz in 0..<MPA_Q1D do
                    QQQ[qz, qy, qx] = u[qz] * D_(qx, qy, qz, e);
                }
              }

              for d in 0..<MPA_D1D {
                for q in 0..<MPA_Q1D {
                  // MASS3DPA_6
                  Btsmem[d, q] = Bt_(q, d);
                }
              }

              for qy in 0..<MPA_Q1D {
                for dx in 0..<MPA_D1D {
                  // MASS3DPA_7
                  var u: [0..<MPA_Q1D] real;
                  for dz in 0..<MPA_Q1D do
                    u[dz] = 0;
                  for qx in 0..<MPA_Q1D do
                    for qz in 0..<MPA_Q1D do
                      u[qz] += QQQ[qz, qy, qx] * Btsmem[dx, qx];
                  for qz in 0..<MPA_Q1D do
                    QQD[qz, qy, dx] = u[qz];
                }
              }

              for dy in 0..<MPA_D1D {
                for dx in 0..<MPA_D1D {
                  // MASS3DPA_8
                  var u: [0..<MPA_Q1D] real;
                  for dz in 0..<MPA_Q1D do
                    u[dz] = 0;
                  for qy in 0..<MPA_Q1D do
                    for qz in 0..<MPA_Q1D do
                      u[qz] += QQD[qz, qy, dx] * Btsmem[dy, qy];
                  for qz in 0..<MPA_Q1D do
                    QDD[qz, dy, dx] = u[qz];
                }
              }

              for dy in 0..<MPA_D1D {
                for dx in 0..<MPA_D1D {
                  // MASS3DPA_9
                  var u: [0..<MPA_D1D] real;
                  for dz in 0..<MPA_D1D do
                    u[dz] = 0;
                  for qz in 0..<MPA_Q1D do
                    for dz in 0..<MPA_D1D do
                      u[dz] += QDD[qz, dy, dx] * Btsmem[dz, qz];
                  for dz in 0..<MPA_D1D do
                    Y_(dx, dy, dz, e) += u[dz];
                }
              }

            }  // element loop
          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for irep in 0..<run_reps {
            forall e in 0..<NE {

              // MASS3DPA_0_CPU
              const MQ1 = MPA_Q1D;
              const MD1 = MPA_D1D;
              const MDQ = if MQ1 > MD1 then MQ1 else MD1;

              var sDQ: [0..<MDQ, 0..<MDQ] real;
              ref  Bsmem = sDQ;  // (MQ1, MD1));
              ref Btsmem = sDQ;  // (MD1, MQ1));

              var sm0: [0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              var sm1: [0..<MDQ, 0..<MDQ, 0..<MDQ] real;
              ref Xsmem = sm0;  // (MD1, MD1, MD1)
              ref   DDQ = sm1;  // (MD1, MD1, MQ1)
              ref   DQQ = sm0;  // (MD1, MQ1, MQ1)
              ref   QQQ = sm1;  // (MQ1, MQ1, MQ1) 3?
              ref   QQD = sm0;  // (MQ1, MQ1, MD1)
              ref   QDD = sm1;  // (MQ1, MD1, MD1) 2?

              for dy in 0..<MPA_D1D {
                for dx in 0..<MPA_D1D{
                  // MASS3DPA_1
                  for dz in 0..<MPA_D1D do
                    Xsmem[dz, dy, dx] = X_(dx, dy, dz, e);
                }
                for dx in 0..<MPA_Q1D {
                  // MASS3DPA_2
                  Bsmem[dx, dy] = B_(dx, dy);
                }
              }

              for dy in 0..<MPA_D1D {
                for qx in 0..<MPA_Q1D {
                  // MASS3DPA_3
                  var u: [0..<MPA_D1D] real;
                  for dz in 0..<MPA_D1D do
                    u[dz] = 0;
                  for dx in 0..<MPA_D1D do
                    for dz in 0..<MPA_D1D do
                      u[dz] += Xsmem[dz, dy, dx] * Bsmem[qx, dx];
                  for dz in 0..<MPA_D1D do
                    DDQ[dz, dy, qx] = u[dz];
                }
              }

              for qy in 0..<MPA_Q1D {
                for qx in 0..<MPA_Q1D {
                  // MASS3DPA_4
                  var u: [0..<MPA_D1D] real;
                  for dz in 0..<MPA_D1D do
                    u[dz] = 0;
                  for dy in 0..<MPA_D1D do
                    for dz in 0..<MPA_D1D do
                      u[dz] += DDQ[dz, dy, qx] * Bsmem[qy, dy];
                  for dz in 0..<MPA_D1D do
                    DQQ[dz, qy, qx] = u[dz];
                }
              }

              for qy in 0..<MPA_Q1D {
                for qx in 0..<MPA_Q1D {
                  // MASS3DPA_5
                  var u: [0..<MPA_Q1D] real;
                  for dz in 0..<MPA_Q1D do
                    u[dz] = 0;
                  for dz in 0..<MPA_D1D do
                    for qz in 0..<MPA_Q1D do
                      u[qz] += DQQ[dz, qy, qx] * Bsmem[qz, dz];
                  for qz in 0..<MPA_Q1D do
                    QQQ[qz, qy, qx] = u[qz] * D_(qx, qy, qz, e);
                }
              }

              for d in 0..<MPA_D1D {
                for q in 0..<MPA_Q1D {
                  // MASS3DPA_6
                  Btsmem[d, q] = Bt_(q, d);
                }
              }

              for qy in 0..<MPA_Q1D {
                for dx in 0..<MPA_D1D {
                  // MASS3DPA_7
                  var u: [0..<MPA_Q1D] real;
                  for dz in 0..<MPA_Q1D do
                    u[dz] = 0;
                  for qx in 0..<MPA_Q1D do
                    for qz in 0..<MPA_Q1D do
                      u[qz] += QQQ[qz, qy, qx] * Btsmem[dx, qx];
                  for qz in 0..<MPA_Q1D do
                    QQD[qz, qy, dx] = u[qz];
                }
              }

              for dy in 0..<MPA_D1D {
                for dx in 0..<MPA_D1D {
                  // MASS3DPA_8
                  var u: [0..<MPA_Q1D] real;
                  for dz in 0..<MPA_Q1D do
                    u[dz] = 0;
                  for qy in 0..<MPA_Q1D do
                    for qz in 0..<MPA_Q1D do
                      u[qz] += QQD[qz, qy, dx] * Btsmem[dy, qy];
                  for qz in 0..<MPA_Q1D do
                    QDD[qz, dy, dx] = u[qz];
                }
              }

              for dy in 0..<MPA_D1D {
                for dx in 0..<MPA_D1D {
                  // MASS3DPA_9
                  var u: [0..<MPA_D1D] real;
                  for dz in 0..<MPA_D1D do
                    u[dz] = 0;
                  for qz in 0..<MPA_Q1D do
                    for dz in 0..<MPA_D1D do
                      u[dz] += QDD[qz, dy, dx] * Btsmem[dz, qz];
                  for dz in 0..<MPA_D1D do
                    Y_(dx, dy, dz, e) += u[dz];
                }
              }

            }  // element loop
          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(Y, MPA_D1D*MPA_D1D*MPA_D1D*m_NE);
    }
  }

  class PRESSURE: KernelBase {

    proc init() {
      super.init(KernelID.Apps_PRESSURE);

      setDefaultProblemSize(1000000);
      setDefaultReps(700);

      setActualProblemSize(getTargetProblemSize());

      setItsPerRep(2 * getActualProblemSize());
      setKernelsPerRep(2);
      setBytesPerRep((1*sizeof(Real_type) + 1*sizeof(Real_type)) * getActualProblemSize() +
                     (1*sizeof(Real_type) + 2*sizeof(Real_type)) * getActualProblemSize());
      setFLOPsPerRep((2 + 1) * getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var compression = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var bvc   = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var p_new = allocAndInitDataConst(Real_type, getActualProblemSize(), 0.0, vid);
      var e_old = allocAndInitData(Real_type, getActualProblemSize(), vid);
      var vnewc = allocAndInitData(Real_type, getActualProblemSize(), vid);

      var cls = initData(Real_type, vid);
      var p_cut = initData(Real_type, vid);
      var pmin = initData(Real_type, vid);
      var eosvmax = initData(Real_type, vid);

      const run_reps = getRunReps();
      const ibegin = 0;
      const iend = getActualProblemSize();

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            for i in ibegin..<iend do
              bvc[i] = cls * (compression[i] + 1.0);

            for i in ibegin..<iend {
              p_new[i] = bvc[i] * e_old[i];
              if abs(p_new[i]) <   p_cut then p_new[i] = 0.0;
              if     vnewc[i] >= eosvmax then p_new[i] = 0.0;
              if     p_new[i]  <    pmin then p_new[i] = pmin;
            }

          }

          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();

          for irep in 0..<run_reps {

            forall i in ibegin..<iend do
              bvc[i] = cls * (compression[i] + 1.0);

            forall i in ibegin..<iend {
              p_new[i] = bvc[i] * e_old[i];
              if abs(p_new[i]) <   p_cut then p_new[i] = 0.0;
              if     vnewc[i] >= eosvmax then p_new[i] = 0.0;
              if     p_new[i]  <    pmin then p_new[i] = pmin;
            }

          }

          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(p_new, getActualProblemSize());
    }
  }

  class VOL3D: KernelBase {

    var m_domain;
    var m_array_length: Index_type;

    proc init() {
      super.init(KernelID.Apps_VOL3D);

      setDefaultProblemSize(100*100*100);  // See rzmax in ADomain struct
      setDefaultReps(100);

      var rzmax = cbrt(getTargetProblemSize()):Index_type+1;
      m_domain = new ADomain(rzmax, /* ndims = */ 3);

      m_array_length = m_domain.nnalls;

      setActualProblemSize(m_domain.lpz+1 - m_domain.fpz);

      setItsPerRep(m_domain.lpz+1 - m_domain.fpz);
      setKernelsPerRep(1);
      // touched data size, not actual number of stores and loads
      setBytesPerRep((1*sizeof(Real_type) + 0*sizeof(Real_type)) *  getItsPerRep() +
                     (0*sizeof(Real_type) + 3*sizeof(Real_type)) * (getItsPerRep() + 1+m_domain.jp+m_domain.kp));
      setFLOPsPerRep(72 * (m_domain.lpz+1 - m_domain.fpz));

      checksum_scale_factor = 0.001 *
        (getDefaultProblemSize():Checksum_type/getActualProblemSize());

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var x = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var y = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);
      var z = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);

      var dx: Real_type = 0.3;
      var dy: Real_type = 0.2;
      var dz: Real_type = 0.1;
      setMeshPositions_3d(x, dx, y, dy, z, dz, m_domain);

      var vol = allocAndInitDataConst(Real_type, m_array_length, 0.0, vid);

      const vnormq: Real_type = 0.083333333333333333;  /* vnormq = 1/12 */

      // NDPTRSET(m_domain.jp, m_domain.kp, x,x0,x1,x2,x3,x4,x5,x6,x7);
      ref x0 = reindex(x,                 0..);
      ref x1 = reindex(x0[1..],           0..);
      ref x2 = reindex(x0[m_domain.jp..], 0..);
      ref x3 = reindex(x1[m_domain.jp..], 0..);
      ref x4 = reindex(x0[m_domain.kp..], 0..);
      ref x5 = reindex(x1[m_domain.kp..], 0..);
      ref x6 = reindex(x2[m_domain.kp..], 0..);
      ref x7 = reindex(x3[m_domain.kp..], 0..);

      // NDPTRSET(m_domain.jp, m_domain.kp, y,y0,y1,y2,y3,y4,y5,y6,y7);
      ref y0 = reindex(y,                 0..);
      ref y1 = reindex(y0[1..],           0..);
      ref y2 = reindex(y0[m_domain.jp..], 0..);
      ref y3 = reindex(y1[m_domain.jp..], 0..);
      ref y4 = reindex(y0[m_domain.kp..], 0..);
      ref y5 = reindex(y1[m_domain.kp..], 0..);
      ref y6 = reindex(y2[m_domain.kp..], 0..);
      ref y7 = reindex(y3[m_domain.kp..], 0..);

      // NDPTRSET(m_domain.jp, m_domain.kp, z,z0,z1,z2,z3,z4,z5,z6,z7);
      ref z0 = reindex(z,                 0..);
      ref z1 = reindex(z0[1..],           0..);
      ref z2 = reindex(z0[m_domain.jp..], 0..);
      ref z3 = reindex(z1[m_domain.jp..], 0..);
      ref z4 = reindex(z0[m_domain.kp..], 0..);
      ref z5 = reindex(z1[m_domain.kp..], 0..);
      ref z6 = reindex(z2[m_domain.kp..], 0..);
      ref z7 = reindex(z3[m_domain.kp..], 0..);

      const run_reps = getRunReps();
      const ibegin = m_domain.fpz;
      const iend = m_domain.lpz+1;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();
          for irep in 0..<run_reps {

            for i in ibegin..<iend {
              var x71: Real_type = x7[i] - x1[i];
              var x72: Real_type = x7[i] - x2[i];
              var x74: Real_type = x7[i] - x4[i];
              var x30: Real_type = x3[i] - x0[i];
              var x50: Real_type = x5[i] - x0[i];
              var x60: Real_type = x6[i] - x0[i];

              var y71: Real_type = y7[i] - y1[i];
              var y72: Real_type = y7[i] - y2[i];
              var y74: Real_type = y7[i] - y4[i];
              var y30: Real_type = y3[i] - y0[i];
              var y50: Real_type = y5[i] - y0[i];
              var y60: Real_type = y6[i] - y0[i];

              var z71: Real_type = z7[i] - z1[i];
              var z72: Real_type = z7[i] - z2[i];
              var z74: Real_type = z7[i] - z4[i];
              var z30: Real_type = z3[i] - z0[i];
              var z50: Real_type = z5[i] - z0[i];
              var z60: Real_type = z6[i] - z0[i];

              var xps: Real_type = x71 + x60;
              var yps: Real_type = y71 + y60;
              var zps: Real_type = z71 + z60;

              var cyz: Real_type = y72 * z30 - z72 * y30;
              var czx: Real_type = z72 * x30 - x72 * z30;
              var cxy: Real_type = x72 * y30 - y72 * x30;
              vol[i] = xps * cyz + yps * czx + zps * cxy;

              xps = x72 + x50;
              yps = y72 + y50;
              zps = z72 + z50;

              cyz = y74 * z60 - z74 * y60;
              czx = z74 * x60 - x74 * z60;
              cxy = x74 * y60 - y74 * x60;
              vol[i] += xps * cyz + yps * czx + zps * cxy;

              xps = x74 + x30;
              yps = y74 + y30;
              zps = z74 + z30;

              cyz = y74 * z60 - z74 * y60;
              czx = z74 * x60 - x74 * z60;
              cxy = x74 * y60 - y74 * x60;
              vol[i] += xps * cyz + yps * czx + zps * cxy;

              xps = x74 + x30;
              yps = y74 + y30;
              zps = z74 + z30;

              cyz = y71 * z50 - z71 * y50;
              czx = z71 * x50 - x71 * z50;
              cxy = x71 * y50 - y71 * x50;
              vol[i] += xps * cyz + yps * czx + zps * cxy;

              vol[i] *= vnormq;
            }

          }
          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();
          for irep in 0..<run_reps {

            forall i in ibegin..<iend {
              var x71: Real_type = x7[i] - x1[i];
              var x72: Real_type = x7[i] - x2[i];
              var x74: Real_type = x7[i] - x4[i];
              var x30: Real_type = x3[i] - x0[i];
              var x50: Real_type = x5[i] - x0[i];
              var x60: Real_type = x6[i] - x0[i];

              var y71: Real_type = y7[i] - y1[i];
              var y72: Real_type = y7[i] - y2[i];
              var y74: Real_type = y7[i] - y4[i];
              var y30: Real_type = y3[i] - y0[i];
              var y50: Real_type = y5[i] - y0[i];
              var y60: Real_type = y6[i] - y0[i];

              var z71: Real_type = z7[i] - z1[i];
              var z72: Real_type = z7[i] - z2[i];
              var z74: Real_type = z7[i] - z4[i];
              var z30: Real_type = z3[i] - z0[i];
              var z50: Real_type = z5[i] - z0[i];
              var z60: Real_type = z6[i] - z0[i];

              var xps: Real_type = x71 + x60;
              var yps: Real_type = y71 + y60;
              var zps: Real_type = z71 + z60;

              var cyz: Real_type = y72 * z30 - z72 * y30;
              var czx: Real_type = z72 * x30 - x72 * z30;
              var cxy: Real_type = x72 * y30 - y72 * x30;
              vol[i] = xps * cyz + yps * czx + zps * cxy;

              xps = x72 + x50;
              yps = y72 + y50;
              zps = z72 + z50;

              cyz = y74 * z60 - z74 * y60;
              czx = z74 * x60 - x74 * z60;
              cxy = x74 * y60 - y74 * x60;
              vol[i] += xps * cyz + yps * czx + zps * cxy;

              xps = x74 + x30;
              yps = y74 + y30;
              zps = z74 + z30;

              cyz = y74 * z60 - z74 * y60;
              czx = z74 * x60 - x74 * z60;
              cxy = x74 * y60 - y74 * x60;
              vol[i] += xps * cyz + yps * czx + zps * cxy;

              xps = x74 + x30;
              yps = y74 + y30;
              zps = z74 + z30;

              cyz = y71 * z50 - z71 * y50;
              czx = z71 * x50 - x71 * z50;
              cxy = x71 * y50 - y71 * x50;
              vol[i] += xps * cyz + yps * czx + zps * cxy;

              vol[i] *= vnormq;
            }

          }
          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(vol, m_array_length, checksum_scale_factor:Real_type);
    }
  }

  class COUPLE: KernelBase {

    var m_domain;

    var m_imin: Index_type;
    var m_imax: Index_type;
    var m_jmin: Index_type;
    var m_jmax: Index_type;
    var m_kmin: Index_type;
    var m_kmax: Index_type;

    inline proc zabs2(z) return z.re*z.re+z.im*z.im;

    proc init() {
      super.init(KernelID.Apps_COUPLE);

      setDefaultProblemSize(100*100*100);  // See rzmax in ADomain struct
      setDefaultReps(50);

      var rzmax: Index_type = cbrt(getTargetProblemSize()):Index_type+1;
      m_domain = new ADomain(rzmax, /* ndims = */ 3);

      m_imin = m_domain.imin;
      m_imax = m_domain.imax;
      m_jmin = m_domain.jmin;
      m_jmax = m_domain.jmax;
      m_kmin = m_domain.kmin;
      m_kmax = m_domain.kmax;

      setActualProblemSize(m_domain.n_real_zones);

      setItsPerRep(getActualProblemSize());
      setKernelsPerRep(1);
      setBytesPerRep((3*sizeof(Complex_type) + 5*sizeof(Complex_type)) * m_domain.n_real_zones);
      setFLOPsPerRep(0);

      setUsesFeature(FeatureID.Forall);

      setVariantDefined(VariantID.Base_Chpl);
      setVariantDefined(VariantID.Forall_Chpl);
    }

    override proc runVariant(vid:VariantID) {
      // setup
      var max_loop_index: Index_type = m_domain.lrn;

      var t0 = allocAndInitData(Complex_type, max_loop_index, vid);
      var t1 = allocAndInitData(Complex_type, max_loop_index, vid);
      var t2 = allocAndInitData(Complex_type, max_loop_index, vid);
      var denac = allocAndInitData(Complex_type, max_loop_index, vid);
      var denlw = allocAndInitData(Complex_type, max_loop_index, vid);

      var clight = 3.e+10:Real_type;
      var csound = 3.09e+7:Real_type;
      var omega0 = 0.9:Real_type;
      var omegar = 0.9:Real_type;
      var dt = 0.208:Real_type;
      var c10 = 0.25 * (clight/csound);
      var fratio = sqrt(omegar/omega0);
      var r_fratio = 1.0/fratio;
      var c20 = 0.25 * (clight/csound)*r_fratio;
      var ireal = (0.0 + 1.0i):Complex_type;

      const run_reps = getRunReps();

      const imin = m_imin;
      const imax = m_imax;
      const jmin = m_jmin;
      const jmax = m_jmax;
      const kmin = m_kmin;
      const kmax = m_kmax;

      // run
      select vid {

        when VariantID.Base_Chpl {
          startTimer();
          for irep in 0..<run_reps {

            for k in kmin..<kmax {
              for j in jmin..<jmax {

                var it0: Index_type    = (k*(jmax+1) + j)*(imax+1);
                var idenac: Index_type = (k*(jmax+2) + j)*(imax+2);

                for i in imin..<imax {

                  var c1: Complex_type = c10 * denac[idenac+i];
                  var c2: Complex_type = c20 * denlw[it0+i];

                  /* promote to doubles to avoid possible divide by zero */
                  var c1re: Real_type = c1.re; var c1im: Real_type = c1.im;
                  var c2re: Real_type = c2.re; var c2im: Real_type = c2.im;

                  /* lamda = sqrt(|c1|^2 + |c2|^2) uses doubles to avoid underflow. */
                  var zlam: Real_type = c1re*c1re + c1im*c1im +
                                        c2re*c2re + c2im*c2im + 1.0e-34;
                  zlam = sqrt(zlam);
                  var snlamt: Real_type = sin(zlam * dt * 0.5);
                  var cslamt: Real_type = cos(zlam * dt * 0.5);

                  var a0t: Complex_type = t0[it0+i];
                  var a1t: Complex_type = t1[it0+i];
                  var a2t: Complex_type = t2[it0+i] * fratio;

                  var r_zlam: Real_type = 1.0/zlam;
                  c1 *= r_zlam;
                  c2 *= r_zlam;
                  var zac1: Real_type = zabs2(c1);
                  var zac2: Real_type = zabs2(c2);

                  /* compute new A0 */
                  var z3: Complex_type = (c1*a1t + c2*a2t) * snlamt;
                  t0[it0+i] = a0t*cslamt - ireal*z3;

                  /* compute new A1  */
                  var r: Real_type = zac1*cslamt + zac2;
                  var z5: Complex_type = c2 * a2t;
                  var z4: Complex_type = conjg(c1) * z5 * (cslamt-1);
                  z3 = conjg(c1) * a0t * snlamt;
                  t1[it0+i] = a1t * r + z4 - ireal * z3;

                  /* compute new A2  */
                  r = zac1 + zac2 * cslamt;
                  z5 = c1 * a1t;
                  z4 = conjg(c2) * z5 * (cslamt-1);
                  z3 = conjg(c2) * a0t * snlamt;
                  t2[it0+i] = (a2t*r + z4 - ireal*z3) * r_fratio;

                } // i loop

              } // j loop
            }  // k loop

          }
          stopTimer();
        }

        when VariantID.Forall_Chpl {
          startTimer();
          for irep in 0..<run_reps {

            forall k in kmin..<kmax {
              for j in jmin..<jmax {

                var it0: Index_type    = (k*(jmax+1) + j)*(imax+1);
                var idenac: Index_type = (k*(jmax+2) + j)*(imax+2);

                for i in imin..<imax {

                  var c1: Complex_type = c10 * denac[idenac+i];
                  var c2: Complex_type = c20 * denlw[it0+i];

                  /* promote to doubles to avoid possible divide by zero */
                  var c1re: Real_type = c1.re; var c1im: Real_type = c1.im;
                  var c2re: Real_type = c2.re; var c2im: Real_type = c2.im;

                  /* lamda = sqrt(|c1|^2 + |c2|^2) uses doubles to avoid underflow. */
                  var zlam: Real_type = c1re*c1re + c1im*c1im +
                                        c2re*c2re + c2im*c2im + 1.0e-34;
                  zlam = sqrt(zlam);
                  var snlamt: Real_type = sin(zlam * dt * 0.5);
                  var cslamt: Real_type = cos(zlam * dt * 0.5);

                  var a0t: Complex_type = t0[it0+i];
                  var a1t: Complex_type = t1[it0+i];
                  var a2t: Complex_type = t2[it0+i] * fratio;

                  var r_zlam: Real_type = 1.0/zlam;
                  c1 *= r_zlam;
                  c2 *= r_zlam;
                  var zac1: Real_type = zabs2(c1);
                  var zac2: Real_type = zabs2(c2);

                  /* compute new A0 */
                  var z3: Complex_type = (c1*a1t + c2*a2t) * snlamt;
                  t0[it0+i] = a0t*cslamt - ireal*z3;

                  /* compute new A1  */
                  var r: Real_type = zac1*cslamt + zac2;
                  var z5: Complex_type = c2 * a2t;
                  var z4: Complex_type = conjg(c1) * z5 * (cslamt-1);
                  z3 = conjg(c1) * a0t * snlamt;
                  t1[it0+i] = a1t * r + z4 - ireal * z3;

                  /* compute new A2  */
                  r = zac1 + zac2 * cslamt;
                  z5 = c1 * a1t;
                  z4 = conjg(c2) * z5 * (cslamt-1);
                  z3 = conjg(c2) * a0t * snlamt;
                  t2[it0+i] = (a2t*r + z4 - ireal*z3) * r_fratio;

                } // i loop

              } // j loop
            }  // k loop

          }
          stopTimer();
        }

        otherwise halt();

      }

      // update checksum
      checksum[vid] += calcChecksum(t0, max_loop_index);
      checksum[vid] += calcChecksum(t1, max_loop_index);
      checksum[vid] += calcChecksum(t2, max_loop_index);
    }
  }

}
