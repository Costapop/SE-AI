/* Эллипсометрия тонкой плёнки: физика, синтез, фит и ИИ-модель на остатках.
   Конвенция: N = n - ik, фазовый множитель exp(-2iβ), ρ = r_p/r_s = tanΨ·exp(iΔ) (Аззам–Башара). */
(function (root) {
  'use strict';
  const HC = 1239.841984;
  const PI = Math.PI;

  // ---------------------------------------------------------------- комплексные числа [re, im]
  const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
  const cdiv = (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; };
  const cadd = (a, b) => [a[0] + b[0], a[1] + b[1]];
  const csub = (a, b) => [a[0] - b[0], a[1] - b[1]];
  const cabs2 = (a) => a[0] * a[0] + a[1] * a[1];
  function csqrt(a) { const r = Math.hypot(a[0], a[1]); let re = Math.sqrt(Math.max(0, (r + a[0]) / 2)); let im = Math.sqrt(Math.max(0, (r - a[0]) / 2)); if (a[1] < 0) im = -im; return [re, im]; }
  function cexp(a) { const e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; }

  // ---------------------------------------------------------------- дисперсия: Тауц–Лоренц + УФ-полюс (ε∞ = 1)
  function tlEps(E, A, E0, C, Eg) {
    const E2 = E * E, E02 = E0 * E0, C2 = C * C, Eg2 = Eg * Eg;
    const eps2 = E > Eg ? A * E0 * C * (E - Eg) * (E - Eg) / (((E2 - E02) * (E2 - E02) + C2 * E2) * E) : 0;
    const aln = (Eg2 - E02) * E2 + Eg2 * C2 - E02 * (E02 + 3 * Eg2);
    const aatan = (E2 - E02) * (E02 + Eg2) + Eg2 * C2;
    const alpha = Math.sqrt(4 * E02 - C2), gamma2 = E02 - C2 / 2;
    const zeta4 = (E2 - gamma2) * (E2 - gamma2) + alpha * alpha * C2 / 4;
    const eps1 = A * C * aln / (2 * PI * zeta4 * alpha * E0) * Math.log((E02 + Eg2 + alpha * Eg) / (E02 + Eg2 - alpha * Eg))
      - A * aatan / (PI * zeta4 * E0) * (PI - Math.atan((2 * Eg + alpha) / C) + Math.atan((-2 * Eg + alpha) / C))
      + 2 * A * E0 * Eg * (E2 - gamma2) / (PI * zeta4 * alpha) * (PI + 2 * Math.atan(2 * (gamma2 - Eg2) / (alpha * C)))
      - A * E0 * C * (E2 + Eg2) / (PI * zeta4 * E) * Math.log(Math.abs(E - Eg) / (E + Eg))
      + 2 * A * E0 * C * Eg / (PI * zeta4) * Math.log(Math.abs(E - Eg) * (E + Eg) / Math.sqrt((E02 - Eg2) * (E02 - Eg2) + Eg2 * C2));
    return [eps1, eps2];
  }
  /** Комплексный показатель плёнки в конвенции n - ik: массив [n, -k] на каждой λ. */
  function filmN(lam, m) {
    const out = new Array(lam.length);
    for (let i = 0; i < lam.length; i++) {
      const E = HC / lam[i];
      const [e1, e2] = tlEps(E, m.A, m.E0, m.C, m.Eg);
      const eps1 = 1 + e1 + m.Auv / (m.Euv * m.Euv - E * E);
      const N = csqrt([eps1, e2]);
      out[i] = [N[0], -Math.abs(N[1])];
    }
    return out;
  }
  function filmNK(lam, m) { const N = filmN(lam, m); return { n: N.map(v => v[0]), k: N.map(v => -v[1]) }; }

  // ---------------------------------------------------------------- подложки (полубесконечные)
  function nBK7(l) { const x = (l / 1000) ** 2; return Math.sqrt(1 + 1.03961212 * x / (x - 0.00600069867) + 0.231792344 * x / (x - 0.0200179144) + 1.01046945 * x / (x - 103.560653)); }
  function nSilica(l) { const x = (l / 1000) ** 2; return Math.sqrt(1 + 0.6961663 * x / (x - 0.0684043 ** 2) + 0.4079426 * x / (x - 0.1162414 ** 2) + 0.8974794 * x / (x - 9.896161 ** 2)); }
  const SI_TAB = [[400, 5.570, 0.387], [420, 5.010, 0.262], [440, 4.780, 0.183], [460, 4.610, 0.130], [480, 4.460, 0.098], [500, 4.298, 0.073],
    [520, 4.200, 0.058], [540, 4.120, 0.046], [560, 4.055, 0.038], [580, 3.995, 0.030], [600, 3.939, 0.024], [620, 3.910, 0.020],
    [640, 3.870, 0.017], [660, 3.840, 0.015], [680, 3.810, 0.013], [700, 3.775, 0.011], [720, 3.755, 0.010], [740, 3.735, 0.0086],
    [760, 3.720, 0.0076], [780, 3.705, 0.0068], [800, 3.692, 0.0060], [850, 3.664, 0.0045], [900, 3.641, 0.0033], [1000, 3.606, 0.0015]];
  function siNK(l) {
    const t = SI_TAB; if (l <= t[0][0]) return [t[0][1], t[0][2]]; if (l >= t[t.length - 1][0]) { const e = t[t.length - 1]; return [e[1], e[2]]; }
    let j = 1; while (t[j][0] < l) j++;
    const f = (l - t[j - 1][0]) / (t[j][0] - t[j - 1][0]);
    return [t[j - 1][1] + f * (t[j][1] - t[j - 1][1]), t[j - 1][2] + f * (t[j][2] - t[j - 1][2])];
  }
  const SUBSTRATES = {
    bk7: { name: 'Стекло BK7 (Зельмейер Schott)', fn: l => [nBK7(l), 0] },
    silica: { name: 'Плавленый кварц (Malitson 1965)', fn: l => [nSilica(l), 0] },
    si: { name: 'Кремний c-Si (табл. ≈ Aspnes–Studna 1983)', fn: l => siNK(l) },
    glass152: { name: 'Стекло n = 1.52 (без дисперсии)', fn: l => [1.52, 0] },
  };
  function substrateN(key, lam) { const f = SUBSTRATES[key].fn; return Array.from(lam).map(l => { const nk = f(l); return [nk[0], -nk[1]]; }); }

  // ---------------------------------------------------------------- справочник материалов (Тауц–Лоренц + УФ-полюс; типовые значения для аморфных плёнок)
  const MATERIALS = {
    tio2: { name: 'TiO₂ (аморфный)', A: 255.83, E0: 4.00, C: 1.77, Eg: 3.40, Auv: 137.65, Euv: 11.0 },
    ta2o5: { name: 'Ta₂O₅ (аморфный)', A: 321.34, E0: 5.30, C: 2.60, Eg: 4.20, Auv: 11.33, Euv: 12.0 },
    nb2o5: { name: 'Nb₂O₅ (аморфный)', A: 307.24, E0: 4.70, C: 2.30, Eg: 3.75, Auv: 69.07, Euv: 11.0 },
    hfo2: { name: 'HfO₂ (аморфный)', A: 396.97, E0: 6.20, C: 2.80, Eg: 5.20, Auv: 16.82, Euv: 13.0 },
    si3n4: { name: 'Si₃N₄ (LPCVD)', A: 156.51, E0: 7.20, C: 3.60, Eg: 4.60, Auv: 103.31, Euv: 13.0 },
    sio2: { name: 'SiO₂ (термический)', A: 295.50, E0: 10.35, C: 3.00, Eg: 8.70, Auv: 10.76, Euv: 16.0 },
    al2o3: { name: 'Al₂O₃ (аморфный)', A: 182.48, E0: 9.50, C: 3.50, Eg: 6.60, Auv: 11.45, Euv: 14.0 },
    zro2: { name: 'ZrO₂ (аморфный)', A: 449.00, E0: 5.95, C: 3.00, Eg: 5.00, Auv: 28.74, Euv: 13.0 },
  };

  // ---------------------------------------------------------------- слои
  function bruggeman(ea, eb, fa) {
    const fb = 1 - fa;
    const b = cadd([(3 * fa - 1) * ea[0], (3 * fa - 1) * ea[1]], [(3 * fb - 1) * eb[0], (3 * fb - 1) * eb[1]]);
    const disc = cadd(cmul(b, b), [8 * (ea[0] * eb[0] - ea[1] * eb[1]), 8 * (ea[0] * eb[1] + ea[1] * eb[0])]);
    const s = csqrt(disc);
    const r1 = [(b[0] + s[0]) / 4, (b[1] + s[1]) / 4], r2 = [(b[0] - s[0]) / 4, (b[1] - s[1]) / 4];
    return (r1[0] > 0) ? r1 : r2;
  }
  /** Слои сверху вниз: [{N: [[n,-k],...], d}]. delta > 0 — n растёт к поверхности. */
  function buildLayers(lam, d, m, delta, dRough, M) {
    const Nf = filmN(lam, m);
    const film = [];
    if (!delta) film.push({ N: Nf, d: d });
    else for (let j = M - 1; j >= 0; j--) { const z = (j + 0.5) / M; const s = 1 + delta * (z - 0.5); film.push({ N: Nf.map(v => [v[0] * s, v[1] * s]), d: d / M }); }
    const layers = [];
    if (dRough > 0) {
      const top = film[0].N;
      const Ne = top.map(v => { const eps = cmul(v, v); const e = bruggeman(eps, [1, 0], 0.5); const N = csqrt(e); return [Math.abs(N[0]), -Math.abs(N[1])]; });
      layers.push({ N: Ne, d: dRough });
    }
    return layers.concat(film);
  }

  // ---------------------------------------------------------------- отражение от стопки
  function cosT(N, s02) { const q = cdiv([s02, 0], cmul(N, N)); return csqrt([1 - q[0], -q[1]]); }
  function fresnelPS(Ni, ci, Nj, cj) {
    const a = cmul(Nj, ci), b = cmul(Ni, cj);
    const rp = cdiv(csub(a, b), cadd(a, b));
    const c = cmul(Ni, ci), e = cmul(Nj, cj);
    const rs = cdiv(csub(c, e), cadd(c, e));
    return [rp, rs];
  }
  /** r_p, r_s для каждой λ. layers сверху вниз; Nsub — массив [n,-k] подложки; phi — угол, град. */
  function rhoStack(lam, layers, Nsub, phi) {
    const s0 = Math.sin(phi * PI / 180), s02 = s0 * s0, c0 = [Math.cos(phi * PI / 180), 0], N0 = [1, 0];
    const L = layers.length, RP = new Array(lam.length), RS = new Array(lam.length);
    const cache = layers.map(() => null);
    for (let i = 0; i < lam.length; i++) {
      const Nb = Nsub[i], cb = cosT(Nb, s02);
      let Na = L ? layers[L - 1].N[i] : N0, ca = L ? cosT(Na, s02) : c0;
      let [rp, rs] = fresnelPS(Na, ca, Nb, cb);
      for (let j = L - 1; j >= 0; j--) {
        const Nj = layers[j].N[i], cj = ca;
        const bt = cmul([2 * PI * layers[j].d / lam[i], 0], cmul(Nj, cj));
        const X = cexp([2 * bt[1], -2 * bt[0]]);                        // exp(-2iβ)
        const Nu = j > 0 ? layers[j - 1].N[i] : N0, cu = j > 0 ? cosT(Nu, s02) : c0;
        const [rpu, rsu] = fresnelPS(Nu, cu, Nj, cj);
        rp = cdiv(cadd(rpu, cmul(rp, X)), cadd([1, 0], cmul(cmul(rpu, rp), X)));
        rs = cdiv(cadd(rsu, cmul(rs, X)), cadd([1, 0], cmul(cmul(rsu, rs), X)));
        ca = cu;
      }
      RP[i] = rp; RS[i] = rs;
    }
    return [RP, RS];
  }
  function psiDeltaFromR(RP, RS) {
    const n = RP.length, psi = new Float64Array(n), del = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const rp = RP[i], rs = RS[i];
      psi[i] = Math.atan(Math.sqrt(cabs2(rp) / cabs2(rs))) * 180 / PI;
      const re = rp[0] * rs[0] + rp[1] * rs[1], im = rp[1] * rs[0] - rp[0] * rs[1];
      let d = Math.atan2(im, re) * 180 / PI; if (d < 0) d += 360; del[i] = d;
    }
    return { psi, del };
  }
  function ncsFromR(RP, RS) {
    const n = RP.length, N = new Float64Array(n), C = new Float64Array(n), S = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const rp = RP[i], rs = RS[i], pp = cabs2(rp), ss = cabs2(rs), t = pp + ss;
      N[i] = (ss - pp) / t; C[i] = 2 * (rp[0] * rs[0] + rp[1] * rs[1]) / t; S[i] = 2 * (rp[1] * rs[0] - rp[0] * rs[1]) / t;
    }
    return { N, C, S };
  }
  /** Непрерывный линейный профиль n(z) = n·[1 + δ(z/d − ½)]: уравнение Риккати для локальной функции отражения
   *  (Furman & Tikhonravov, 1.1.19, 1.1.20), RK4 с шагом h (нм); EMA-слой сверху добавляется по Эйри.
   *  Возвращает [RP, RS] в конвенции ядра (r_p как в эллипсометрии, т.е. −r_p книги). */
  function rhoGraded(lam, Nsub, phi, d, m, delta, dRough, h) {
    h = h || 1.0;
    const s0 = Math.sin(phi * PI / 180), s02 = s0 * s0, c0 = [Math.cos(phi * PI / 180), 0], N0 = [1, 0];
    const Nf = filmN(lam, m);
    const RP = new Array(lam.length), RS = new Array(lam.length);
    const nSteps = Math.max(4, Math.ceil(d / h)), dz = d / nSteps;
    for (let i = 0; i < lam.length; i++) {
      const k = 2 * PI / lam[i], Nb = Nsub[i], cb = cosT(Nb, s02), Nc = Nf[i];
      const epsAt = (z) => { const sc = 1 + delta * (z / d - 0.5); const N = [Nc[0] * sc, Nc[1] * sc]; return cmul(N, N); };
      const out = {};
      for (const pol of ['s', 'p']) {
        const qa = pol === 's' ? c0[0] : 1 / c0[0];                             // внешняя среда — воздух
        const qs = pol === 's' ? cmul(Nb, cb) : cdiv(Nb, cb);
        const pref = [0, k / (2 * qa)];
        const f = (z, r) => {
          const eps = epsAt(z), om = csub([1, 0], r), op = cadd([1, 0], r), om2 = cmul(om, om), op2 = cmul(op, op);
          let term;
          if (pol === 's') term = csub([qa * qa * om2[0], qa * qa * om2[1]], cmul(csub(eps, [s02, 0]), op2));
          else { const g = cmul(csub([1, 0], cdiv([s02, 0], eps)), om2); term = csub([qa * qa * g[0], qa * qa * g[1]], cmul(eps, op2)); }
          return cmul(pref, term);
        };
        let r = cdiv(csub([qa, 0], qs), cadd([qa, 0], qs));                    // голая подложка (1.2.17)
        for (let j = 0; j < nSteps; j++) {
          const z = j * dz;
          const k1 = f(z, r), k2 = f(z + dz / 2, cadd(r, [k1[0] * dz / 2, k1[1] * dz / 2])), k3 = f(z + dz / 2, cadd(r, [k2[0] * dz / 2, k2[1] * dz / 2])), k4 = f(z + dz, cadd(r, [k3[0] * dz, k3[1] * dz]));
          r = [r[0] + dz / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), r[1] + dz / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])];
        }
        if (dRough > 0) {                                                         // EMA-слой поверх: r относительно среды ЭС, затем шаг Эйри
          const sc = 1 + delta * 0.5, Ntop = [Nc[0] * sc, Nc[1] * sc];
          const e = bruggeman(cmul(Ntop, Ntop), [1, 0], 0.5); let Ne = csqrt(e); Ne = [Math.abs(Ne[0]), -Math.abs(Ne[1])];
          const ce = cosT(Ne, s02), qe = pol === 's' ? cmul(Ne, ce) : cdiv(Ne, ce);
          const w = cmul([qa, 0], cdiv(csub([1, 0], r), cadd([1, 0], r)));       // w = v/u на верхней границе плёнки
          const re = cdiv(csub(qe, w), cadd(qe, w));                              // отражение плёнки, видимое из среды ЭС
          const rae = cdiv(csub([qa, 0], qe), cadd([qa, 0], qe));                // Френель воздух → ЭС (конвенция книги)
          const bt = cmul([2 * PI * dRough / lam[i], 0], cmul(Ne, ce)); const X = cexp([2 * bt[1], -2 * bt[0]]);
          r = cdiv(cadd(rae, cmul(re, X)), cadd([1, 0], cmul(cmul(rae, re), X)));
        }
        out[pol] = r;
      }
      RS[i] = out.s; RP[i] = [-out.p[0], -out.p[1]];                             // знак p: конвенция эллипсометрии
    }
    return [RP, RS];
  }
  /** Ψ, Δ модели на сетке lam (без полосы). */
  function modelPsiDelta(lam, Nsub, phi, d, m, delta, dRough, M) {
    const [RP, RS] = rhoStack(lam, buildLayers(lam, d, m, delta || 0, dRough || 0, M || 20), Nsub, phi);
    return psiDeltaFromR(RP, RS);
  }
  /** Чистый спектр «измерения»: полоса (свёртка N, C, S на мелкой сетке), сдвиг шкалы, фактический угол.
   *  opt.continuous — непрерывный профиль градиента (Риккати, шаг opt.h, по умолчанию 0.5 нм) вместо лестницы из 40 подслоёв. */
  function synthesizeClean(lam, subKey, phiTrue, d, m, delta, dRough, opt) {
    opt = opt || {};
    const bw = opt.bw || 0, off = opt.lamOffset || 0, cont = !!(opt.continuous && delta), h = opt.h || 0.5, step = cont ? 1.0 : 0.5;
    const rps = (lf, Nsub) => cont ? rhoGraded(lf, Nsub, phiTrue, d, m, delta, dRough, h) : rhoStack(lf, buildLayers(lf, d, m, delta, dRough, 40), Nsub, phiTrue);
    let psi, del;
    if (bw > 0) {
      const lo = lam[0] - 4 * bw - 2, hi = lam[lam.length - 1] + 4 * bw + 2;
      const nf = Math.floor((hi - lo) / step) + 1; const lf = new Array(nf); for (let i = 0; i < nf; i++) lf[i] = lo + i * step + off;
      const [RP, RS] = rps(lf, substrateN(subKey, lf));
      const { N, C, S } = ncsFromR(RP, RS);
      const sg = bw / 2.3548, half = Math.max(1, Math.ceil(3.5 * sg / step)); const ker = []; let ks = 0;
      for (let j = -half; j <= half; j++) { const w = Math.exp(-0.5 * (j * step / sg) ** 2); ker.push(w); ks += w; }
      const conv = (y) => { const out = new Float64Array(nf); for (let i = 0; i < nf; i++) { let s = 0; for (let j = -half; j <= half; j++) { const ii = Math.min(nf - 1, Math.max(0, i + j)); s += y[ii] * ker[j + half]; } out[i] = s / ks; } return out; };
      const Nc = conv(N), Cc = conv(C), Sc = conv(S);
      psi = new Float64Array(lam.length); del = new Float64Array(lam.length);
      for (let i = 0; i < lam.length; i++) {
        const x = (lam[i] + off - (lo + off)) / step; const j = Math.min(nf - 2, Math.max(0, Math.floor(x))); const f = x - j;
        const Ni = Nc[j] + f * (Nc[j + 1] - Nc[j]), Ci = Cc[j] + f * (Cc[j + 1] - Cc[j]), Si = Sc[j] + f * (Sc[j + 1] - Sc[j]);
        psi[i] = 0.5 * Math.acos(Math.max(-1, Math.min(1, Ni))) * 180 / PI;
        let dd = Math.atan2(Si, Ci) * 180 / PI; if (dd < 0) dd += 360; del[i] = dd;
      }
    } else {
      const la = lam.map(l => l + off);
      const r = psiDeltaFromR(...rps(la, substrateN(subKey, la))); psi = r.psi; del = r.del;
    }
    return { psi, del };
  }
  /** Дрейф калибровки и шум поверх чистого спектра (возвращает новые массивы). */
  function applyNoise(clean, lam, opt, rng) {
    const n = lam.length, t = (l) => (l - (lam[0] + lam[n - 1]) / 2) / ((lam[n - 1] - lam[0]) / 2);
    const dp = opt.driftPsi || [0, 0], dd = opt.driftDel || [0, 0];
    const psi = new Float64Array(n), del = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      psi[i] = clean.psi[i] + dp[0] + dp[1] * t(lam[i]) + (opt.sigPsi ? opt.sigPsi[i] * rng.normal() : 0);
      let v = clean.del[i] + dd[0] + dd[1] * t(lam[i]) + (opt.sigDel ? opt.sigDel[i] * rng.normal() : 0);
      v = v % 360; if (v < 0) v += 360; del[i] = v;
    }
    return { psi, del };
  }
  /** Синтез «измерения» = чистый спектр + дрейф + шум. */
  function synthesize(lam, subKey, phiTrue, d, m, delta, dRough, opt, rng) {
    return applyNoise(synthesizeClean(lam, subKey, phiTrue, d, m, delta, dRough, opt), lam, opt || {}, rng);
  }

  // ---------------------------------------------------------------- генератор случайных чисел (детерминированный)
  function makeRng(seed) {
    let a = seed >>> 0;
    const uniform = () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
    let spare = null;
    const normal = () => { if (spare !== null) { const s = spare; spare = null; return s; } let u, v, s; do { u = 2 * uniform() - 1; v = 2 * uniform() - 1; s = u * u + v * v; } while (s >= 1 || s === 0); const m = Math.sqrt(-2 * Math.log(s) / s); spare = v * m; return u * m; };
    return { uniform, normal, range: (lo, hi) => lo + (hi - lo) * uniform(), int: (n) => Math.floor(uniform() * n) };
  }

  // ---------------------------------------------------------------- фит (Левенберг–Марквардт, численный якобиан, границы)
  function solveLinear(A, b) {
    const n = b.length, M = A.map((r, i) => r.concat([b[i]]));
    for (let c = 0; c < n; c++) {
      let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      [M[c], M[p]] = [M[p], M[c]];
      const piv = M[c][c]; if (Math.abs(piv) < 1e-300) return null;
      for (let r = 0; r < n; r++) { if (r === c) continue; const f = M[r][c] / piv; if (f) for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
    }
    return M.map((r, i) => r[n] / r[i]);
  }
  function invert(A) {
    const n = A.length, M = A.map((r, i) => r.concat(Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))));
    for (let c = 0; c < n; c++) {
      let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      [M[c], M[p]] = [M[p], M[c]];
      const piv = M[c][c]; if (Math.abs(piv) < 1e-300) return null;
      for (let k = 0; k < 2 * n; k++) M[c][k] /= piv;
      for (let r = 0; r < n; r++) { if (r === c) continue; const f = M[r][c]; if (f) for (let k = 0; k < 2 * n; k++) M[r][k] -= f * M[c][k]; }
    }
    return M.map(r => r.slice(n));
  }
  function lmFit(residFn, x0, lo, hi, scale, maxIter) {
    maxIter = maxIter || 60;
    const p = x0.length, clamp = (x) => x.map((v, k) => Math.min(hi[k], Math.max(lo[k], v)));
    let x = clamp(x0.slice()), r = residFn(x), cost = r.reduce((s, v) => s + v * v, 0), lam = 1e-3, nfev = 1;
    let J = null;
    for (let it = 0; it < maxIter; it++) {
      const m = r.length; J = Array.from({ length: p }, () => new Float64Array(m));
      for (let k = 0; k < p; k++) {
        const h = 1e-4 * scale[k]; const xp = x.slice(); xp[k] = Math.min(hi[k], x[k] + h); const hh = xp[k] - x[k] || h;
        const rp = residFn(xp); nfev++; for (let i = 0; i < m; i++) J[k][i] = (rp[i] - r[i]) / hh;
      }
      const JTJ = Array.from({ length: p }, () => new Array(p).fill(0)), JTr = new Array(p).fill(0);
      for (let a = 0; a < p; a++) { for (let b = a; b < p; b++) { let s = 0; for (let i = 0; i < m; i++) s += J[a][i] * J[b][i]; JTJ[a][b] = s; JTJ[b][a] = s; } let s = 0; for (let i = 0; i < m; i++) s += J[a][i] * r[i]; JTr[a] = s; }
      let improved = false;
      for (let tries = 0; tries < 12; tries++) {
        const Aug = JTJ.map((row, a) => row.map((v, b) => (a === b ? v * (1 + lam) + 1e-12 : v)));
        const dx = solveLinear(Aug, JTr.map(v => -v)); if (!dx) { lam *= 10; continue; }
        const xn = clamp(x.map((v, k) => v + dx[k])); const rn = residFn(xn); nfev++;
        const cn = rn.reduce((s, v) => s + v * v, 0);
        if (cn < cost) { const rel = (cost - cn) / Math.max(cost, 1e-300); x = xn; r = rn; cost = cn; lam = Math.max(lam * 0.3, 1e-9); improved = true; if (rel < 1e-9) it = maxIter; break; }
        lam *= 5;
      }
      if (!improved) break;
    }
    const m = r.length, nu = Math.max(1, m - p), chi2 = cost / nu;
    // ковариация по последнему якобиану
    let cov = null;
    if (J) { const JTJ = Array.from({ length: p }, () => new Array(p).fill(0)); for (let a = 0; a < p; a++) for (let b = 0; b < p; b++) { let s = 0; for (let i = 0; i < m; i++) s += J[a][i] * J[b][i]; JTJ[a][b] = s; } const inv = invert(JTJ); if (inv) cov = inv.map(row => row.map(v => v * chi2)); }
    return { x, resid: r, chi2, cov, nfev, err: cov ? cov.map((row, k) => Math.sqrt(Math.max(0, row[k]))) : x.map(() => NaN) };
  }

  // ---------------------------------------------------------------- модели фита
  const wrap = (v) => { let w = ((v + 180) % 360 + 360) % 360 - 180; return w; };
  const PARAM_NAMES = ['d', 'A', 'Auv', 'Eg', 'dRough', 'delta'];
  function makeResidFn(data, cfg, free, fixedVals, M) {
    const Nsub = cfg.Nsub, lam = data.lam, n = lam.length;
    return (x) => {
      const P = Object.assign({}, fixedVals); free.forEach((k, i) => P[k] = x[i]);
      const m = Object.assign({}, cfg.mat, { A: P.A, Auv: P.Auv, Eg: P.Eg });
      const r = modelPsiDelta(lam, Nsub, cfg.phi, P.d, m, P.delta || 0, P.dRough || 0, M || 20);
      const out = new Float64Array(2 * n);
      for (let i = 0; i < n; i++) { out[i] = (data.psi[i] - r.psi[i]) / data.sigPsi[i]; out[n + i] = wrap(data.del[i] - r.del[i]) / data.sigDel[i]; }
      return out;
    };
  }
  function boundsFor(k, cfg, prior) {
    const m = cfg.mat;
    switch (k) {
      case 'd': return [prior.dMin, prior.dMax, 5];
      case 'A': return [0.2 * m.A, 4 * m.A, 0.05 * m.A];
      case 'Auv': return [0, 3 * Math.abs(m.Auv) + 100, 0.05 * Math.abs(m.Auv) + 5];
      case 'Eg': return [Math.max(1.5, m.Eg - 1.0), Math.min(m.E0 - 0.15, m.Eg + 1.0), 0.05];
      case 'dRough': return [0, 12, 0.5];
      case 'delta': return [-0.20, 0.20, 0.01];
    }
  }
  /** Опорный фит однородной плёнки: грубый перебор d + ЛМ по (d, A, Auv, Eg). */
  function referenceFit(data, cfg, prior) {
    const fixed = { dRough: 0, delta: 0 };
    const base = { d: 0, A: cfg.mat.A, Auv: cfg.mat.Auv, Eg: cfg.mat.Eg };
    // грубая кривая χ²(d) при справочной дисперсии → несколько лучших локальных минимумов → ЛМ из каждого
    const rf = makeResidFn(data, cfg, ['d'], Object.assign({}, fixed, base), 1);
    const step = Math.max(0.5, (prior.dMax - prior.dMin) / 120), grid = [], cost = [];
    for (let dg = prior.dMin; dg <= prior.dMax + 1e-9; dg += step) { const r = rf([dg]); let c = 0; for (let i = 0; i < r.length; i++) c += r[i] * r[i]; grid.push(dg); cost.push(c); }
    const cands = [];
    for (let i = 0; i < grid.length; i++) if ((i === 0 || cost[i] <= cost[i - 1]) && (i === grid.length - 1 || cost[i] <= cost[i + 1])) cands.push(i);
    cands.sort((a, b) => cost[a] - cost[b]);
    const starts = cands.slice(0, 3).map(i => grid[i]);
    const free = ['d', 'A', 'Auv', 'Eg'];
    const b = free.map(k => boundsFor(k, cfg, prior));
    const margin = 0.3 * (prior.dMax - prior.dMin); b[0] = [prior.dMin - margin, prior.dMax + margin, b[0][2]];   // ЛМ может выйти за априорный разброс
    let res = null;
    for (const d0 of starts) {
      const r = lmFit(makeResidFn(data, cfg, free, fixed, 1), [d0, base.A, base.Auv, base.Eg], b.map(v => v[0]), b.map(v => v[1]), b.map(v => v[2]), 80);
      if (!res || r.chi2 < res.chi2) res = r;
    }
    const P = {}; free.forEach((k, i) => P[k] = res.x[i]); P.dRough = 0; P.delta = 0;
    const E = {}; free.forEach((k, i) => E[k] = res.err[i]);
    return { P, E, resid: res.resid, chi2: res.chi2, free };
  }
  /** Фит в выбранной модели дефектов. */
  function finalFit(data, cfg, prior, flags, p0) {
    const free = ['d', 'A', 'Auv', 'Eg']; if (flags.rough) free.push('dRough'); if (flags.grad) free.push('delta');
    const fixed = { dRough: 0, delta: 0 };
    const b = free.map(k => boundsFor(k, cfg, prior));
    const margin = 0.3 * (prior.dMax - prior.dMin); b[0] = [prior.dMin - margin, prior.dMax + margin, b[0][2]];
    const x0 = free.map(k => Math.min(b[free.indexOf(k)][1], Math.max(b[free.indexOf(k)][0], p0[k])));
    let res = lmFit(makeResidFn(data, cfg, free, fixed, 40), x0, b.map(v => v[0]), b.map(v => v[1]), b.map(v => v[2]), 80);
    const iD = free.indexOf('delta');
    if (iD >= 0 && Math.abs(res.x[iD]) > 0.04) res = lmFit(makeResidFn(data, cfg, free, fixed, 80), res.x, b.map(v => v[0]), b.map(v => v[1]), b.map(v => v[2]), 40);   // |δ| > 4 %: 80 подслоёв (ошибка лестницы ≤ 0.4σ)
    const P = { dRough: 0, delta: 0 }; free.forEach((k, i) => P[k] = res.x[i]);
    const E = {}; free.forEach((k, i) => E[k] = res.err[i]);
    return { P, E, resid: res.resid, chi2: res.chi2, free };
  }

  // ---------------------------------------------------------------- якоря (полуволновые/четвертьволновые точки) и каналы остатка
  function unwrapDeg(a) { const out = new Float64Array(a.length); out[0] = a[0]; for (let i = 1; i < a.length; i++) { let d = a[i] - a[i - 1]; d = ((d + 180) % 360 + 360) % 360 - 180; out[i] = out[i - 1] + d; } return out; }
  function interp(x, xs, ys) { if (x <= xs[0]) return ys[0]; if (x >= xs[xs.length - 1]) return ys[ys.length - 1]; let j = 1; while (xs[j] < x) j++; const f = (x - xs[j - 1]) / (xs[j] - xs[j - 1]); return ys[j - 1] + f * (ys[j] - ys[j - 1]); }
  function anchors(lam, del, delSub) {
    const diff = new Float64Array(lam.length); for (let i = 0; i < lam.length; i++) diff[i] = wrap(del[i] - delSub[i]);
    let ph = unwrapDeg(diff); if (ph[ph.length - 1] < ph[0]) ph = ph.map(v => -v);
    const HW = [], QW = [];
    for (let t = Math.ceil(ph[0] / 360 + 1e-9) * 360; t < ph[ph.length - 1]; t += 360) { HW.push(interp(t, ph, lam)); if (t + 180 < ph[ph.length - 1]) QW.push(interp(t + 180, ph, lam)); }
    return { HW, QW };
  }
  function channels(lam, resid, HW, QW) {
    const n = lam.length, rp = Array.from(resid.slice(0, n)), rd = Array.from(resid.slice(n));
    const at = (xs, ys) => xs.map(x => interp(x, lam, ys));
    const psiHW = at(HW, rp), delHW = at(HW, rd), delQW = at(QW, rd);
    const mean = a => a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN, mabs = a => mean(a.map(Math.abs));
    const rms = a => Math.sqrt(a.reduce((s, v) => s + v * v, 0) / a.length);
    return { psiHW: mean(psiHW), psiHWsd: Math.sqrt(Math.max(0, mean(psiHW.map(v => v * v)) - mean(psiHW) ** 2)), delHW: mabs(delHW), delQW: mabs(delQW), rmsPsi: rms(rp), rmsDel: rms(rd), nHW: HW.length };
  }

  // ---------------------------------------------------------------- ИИ: мультиномиальная логистическая регрессия (softmax, Adam) и гребневая регрессия
  function standardizer(X) {
    const n = X.length, D = X[0].length, mu = new Float64Array(D), sd = new Float64Array(D);
    for (let j = 0; j < D; j++) { let s = 0; for (let i = 0; i < n; i++) s += X[i][j]; mu[j] = s / n; }
    for (let j = 0; j < D; j++) { let s = 0; for (let i = 0; i < n; i++) { const v = X[i][j] - mu[j]; s += v * v; } sd[j] = Math.sqrt(s / n) || 1; }
    return { mu, sd, apply: (x) => { const z = new Float64Array(D); for (let j = 0; j < D; j++) z[j] = (x[j] - mu[j]) / sd[j]; return z; } };
  }
  function trainSoftmax(X, y, K, opt) {
    opt = opt || {}; const l2 = opt.l2 ?? 2e-3, epochs = opt.epochs ?? 250, lr = opt.lr ?? 0.05;
    const st = standardizer(X), Z = X.map(st.apply), n = Z.length, D = Z[0].length;
    const W = Array.from({ length: K }, () => new Float64Array(D)), b = new Float64Array(K);
    const mW = Array.from({ length: K }, () => new Float64Array(D)), vW = Array.from({ length: K }, () => new Float64Array(D)), mb = new Float64Array(K), vb = new Float64Array(K);
    const b1 = 0.9, b2 = 0.999, eps = 1e-8;
    const P = Array.from({ length: n }, () => new Float64Array(K));
    for (let ep = 1; ep <= epochs; ep++) {
      // прямой проход
      for (let i = 0; i < n; i++) { const z = Z[i]; let mx = -Infinity; for (let k = 0; k < K; k++) { let s = b[k]; const w = W[k]; for (let j = 0; j < D; j++) s += w[j] * z[j]; P[i][k] = s; if (s > mx) mx = s; } let sum = 0; for (let k = 0; k < K; k++) { P[i][k] = Math.exp(P[i][k] - mx); sum += P[i][k]; } for (let k = 0; k < K; k++) P[i][k] /= sum; }
      // градиенты
      const gW = Array.from({ length: K }, () => new Float64Array(D)), gb = new Float64Array(K);
      for (let i = 0; i < n; i++) { const z = Z[i]; for (let k = 0; k < K; k++) { const g = (P[i][k] - (y[i] === k ? 1 : 0)) / n; gb[k] += g; const gw = gW[k]; for (let j = 0; j < D; j++) gw[j] += g * z[j]; } }
      for (let k = 0; k < K; k++) { const gw = gW[k], w = W[k]; for (let j = 0; j < D; j++) { const g = gw[j] + l2 * w[j]; mW[k][j] = b1 * mW[k][j] + (1 - b1) * g; vW[k][j] = b2 * vW[k][j] + (1 - b2) * g * g; w[j] -= lr * (mW[k][j] / (1 - b1 ** ep)) / (Math.sqrt(vW[k][j] / (1 - b2 ** ep)) + eps); }
        mb[k] = b1 * mb[k] + (1 - b1) * gb[k]; vb[k] = b2 * vb[k] + (1 - b2) * gb[k] * gb[k]; b[k] -= lr * (mb[k] / (1 - b1 ** ep)) / (Math.sqrt(vb[k] / (1 - b2 ** ep)) + eps); }
    }
    const predictProba = (x) => { const z = st.apply(x); const s = new Float64Array(K); let mx = -Infinity; for (let k = 0; k < K; k++) { let v = b[k]; for (let j = 0; j < D; j++) v += W[k][j] * z[j]; s[k] = v; if (v > mx) mx = v; } let sum = 0; for (let k = 0; k < K; k++) { s[k] = Math.exp(s[k] - mx); sum += s[k]; } return Array.from(s, v => v / sum); };
    return { predictProba, W, b, st };
  }
  function ridge(X, y, alpha) {
    const st = standardizer(X), Z = X.map(st.apply), n = Z.length, D = Z[0].length;
    const ym = y.reduce((s, v) => s + v, 0) / n;
    const G = Array.from({ length: D }, () => new Array(D).fill(0)), g = new Array(D).fill(0);
    for (let i = 0; i < n; i++) { const z = Z[i], yc = y[i] - ym; for (let a = 0; a < D; a++) { const za = z[a]; if (!za) continue; g[a] += za * yc; const row = G[a]; for (let b = a; b < D; b++) row[b] += za * z[b]; } }
    for (let a = 0; a < D; a++) { for (let b = 0; b < a; b++) G[a][b] = G[b][a]; G[a][a] += alpha; }
    const w = solveLinear(G, g) || new Array(D).fill(0);
    return { predict: (x) => { const z = st.apply(x); let s = ym; for (let j = 0; j < D; j++) s += w[j] * z[j]; return s; } };
  }

  // ---------------------------------------------------------------- классы дефектов
  const CLASSES = [{ g: 0, r: 0, name: 'однородная плёнка' }, { g: 1, r: 0, name: 'градиент, n растёт к поверхности' }, { g: -1, r: 0, name: 'градиент, n падает к поверхности' },
    { g: 0, r: 1, name: 'шероховатость' }, { g: 1, r: 1, name: 'градиент (n растёт) + шероховатость' }, { g: -1, r: 1, name: 'градиент (n падает) + шероховатость' }];

  /** Вектор признаков: остаток опорного фита + 6 чисел каналов якорей (Ψ в HW, разброс, |Δ| в HW и QW, RMS Ψ/Δ). */
  function featureVector(dd, cfg, resid) {
    const subDel = cfg.subDel || (cfg.subDel = modelPsiDelta(dd.lam, cfg.Nsub, cfg.phi, 0, cfg.mat, 0, 0, 1).del);
    const { HW, QW } = anchors(dd.lam, dd.del, subDel);
    let ch = { psiHW: 0, psiHWsd: 0, delHW: 0, delQW: 0, rmsPsi: 0, rmsDel: 0 };
    if (HW.length >= 2 && QW.length >= 1) ch = channels(dd.lam, resid, HW, QW);
    const f = (v) => (isFinite(v) ? v : 0);
    return Array.from(resid).concat([f(ch.psiHW), f(ch.psiHWsd), f(ch.delHW), f(ch.delQW), f(ch.rmsPsi), f(ch.rmsDel)]);
  }
  /** Один обучающий пример: истина со случайными дефектами и мешающими факторами → опорный фит → остаток. */
  function trainingExample(data, cfg, prior, nuis, rng, cls) {
    const c = CLASSES[cls];
    const dm = nuis.deltaMax, small = Math.min(dm, 0.02);
    const mag = (dm > small + 1e-9 && rng.uniform() < 0.5) ? rng.range(small, dm) : rng.range(0.0002, small);
    const delta = c.g * mag, dRough = c.r ? rng.range(0.02, nuis.roughMax) : 0;
    const d = rng.range(prior.dMin, prior.dMax);
    const m = Object.assign({}, cfg.mat, { A: cfg.mat.A * (1 + rng.range(-0.05, 0.05)), Auv: cfg.mat.Auv * (1 + rng.range(-0.05, 0.05)) + rng.range(-2, 2), E0: cfg.mat.E0 + rng.range(-0.05, 0.05), C: cfg.mat.C + rng.range(-0.06, 0.06), Eg: cfg.mat.Eg + rng.range(-0.05, 0.05) });
    if (m.Auv < 0) m.Auv = 0;
    const opt = { bw: Math.max(0, nuis.bw + rng.range(-0.3, 0.3)), lamOffset: rng.range(-0.1, 0.1), driftPsi: [rng.range(-0.01, 0.01), rng.range(-0.015, 0.015)], driftDel: [rng.range(-0.03, 0.03), rng.range(-0.04, 0.04)], sigPsi: data.sigPsi, sigDel: data.sigDel };
    const s = synthesize(data.lam, cfg.subKey, cfg.phi + rng.range(-nuis.phiErr, nuis.phiErr), d, m, delta, dRough, opt, rng);
    const dd = { lam: data.lam, psi: s.psi, del: s.del, sigPsi: data.sigPsi, sigDel: data.sigDel };
    const fit = referenceFit(dd, cfg, prior);
    return { x: featureVector(dd, cfg, fit.resid), y: cls, delta, dRough, d, dFit: fit.P.d, chi2: fit.chi2 };
  }

  const api = { HC, MATERIALS, SUBSTRATES, CLASSES, filmN, filmNK, substrateN, buildLayers, rhoStack, psiDeltaFromR, modelPsiDelta, rhoGraded, synthesizeClean, applyNoise, synthesize, makeRng, lmFit, referenceFit, finalFit, anchors, channels, featureVector, trainSoftmax, ridge, trainingExample, wrap, interp, tlEps };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.EllipCore = api;
})(typeof window !== 'undefined' ? window : globalThis);
