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
  /** Хвост поглощения ниже края (дефект «поглощение»): второй осциллятор Тауца–Лоренца с теми же E₀ и C, что у материала,
   *  порогом E_t = TAIL_ET и силой m.Ak (эВ). ε₂ хвоста ∝ (E − E_t)²/E · L(E) растёт к синему краю; ε₁ — то же КК-выражение
   *  Джеллисона–Модина, поэтому хвост причинен. При m.Ak = 0 (или undefined) дисперсия совпадает с чистым материалом. */
  const TAIL_ET = 1.2;
  /** Комплексный показатель плёнки в конвенции n - ik: массив [n, -k] на каждой λ. */
  function filmN(lam, m) {
    const out = new Array(lam.length), Ak = m.Ak || 0;
    for (let i = 0; i < lam.length; i++) {
      const E = HC / lam[i];
      const [e1, e2] = tlEps(E, m.A, m.E0, m.C, m.Eg);
      let eps1 = 1 + e1 + m.Auv / (m.Euv * m.Euv - E * E), eps2 = e2;
      if (Ak) { const [t1, t2] = tlEps(E, Ak, m.E0, m.C, TAIL_ET); eps1 += t1; eps2 += t2; }
      const N = csqrt([eps1, eps2]);
      out[i] = [N[0], -Math.abs(N[1])];
    }
    return out;
  }
  function filmNK(lam, m) { const N = filmN(lam, m); return { n: N.map(v => v[0]), k: N.map(v => -v[1]) }; }
  /** k хвоста поглощения при 400 нм (k плёнки с хвостом минус k чистого материала) — характерный размер дефекта «поглощение». */
  function tailK400(m, Ak) {
    const k1 = filmNK([400], Object.assign({}, m, { Ak: Ak })).k[0], k0 = filmNK([400], Object.assign({}, m, { Ak: 0 })).k[0];
    return k1 - k0;
  }
  /** Обратное к tailK400: сила хвоста A_k (эВ), дающая заданный k хвоста при 400 нм (монотонная зависимость; Ньютон). */
  function tailAk(m, k400) {
    if (!(k400 > 0)) return 0;
    let Ak = k400 / Math.max(1e-12, tailK400(m, 1e-3) / 1e-3);            // линейное приближение (k ≈ ε₂/2n)
    for (let it = 0; it < 12; it++) {
      const f = tailK400(m, Ak) - k400; if (Math.abs(f) < 1e-9 * Math.max(k400, 1e-6)) break;
      const h = Math.max(1e-6, 1e-3 * Ak), df = (tailK400(m, Ak + h) - tailK400(m, Ak)) / h;
      const step = f / df; Ak = Math.max(0.1 * Ak, Ak - step);
    }
    return Ak;
  }

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
    // материалы с краем поглощения в измеряемом диапазоне: параметры подобраны по табличным n, k из литературы (см. README, разд. 3.3)
    sinx: { name: 'SiNₓ:H (PECVD, обогащённый Si)', A: 52.04, E0: 7.38, C: 3.26, Eg: 2.17, Auv: 79.35, Euv: 13.0 },   // Vogt 2015 (n = 2.13 при 633 нм): k(400) ≈ 0.036, k(500) ≈ 0.005
    ceo2: { name: 'CeO₂ (поликристаллический)', A: 65.58, E0: 3.92, C: 0.89, Eg: 2.81, Auv: 414.58, Euv: 13.0 },    // ALD, arXiv:1705.04071: k(400) ≈ 0.03, k ≥ 440 нм ≈ 0
    wo3: { name: 'WO₃ (тонкая плёнка)', A: 57.22, E0: 4.41, C: 1.20, Eg: 3.26, Auv: 311.01, Euv: 13.0 },              // Kulikova 2020, Opt. Express 28, 32049: k(350) ≈ 0.023, k ≥ 380 нм ≈ 0
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
    // рекурсия Эйри снизу вверх без промежуточных аллокаций: N = n − ik, X = exp(−2iβ), r_p — эллипсометрический знак
    const s0 = Math.sin(phi * PI / 180), s02 = s0 * s0, c0 = Math.cos(phi * PI / 180);
    const L = layers.length, n = lam.length, RP = new Array(n), RS = new Array(n);
    const cr_ = new Float64Array(L + 1), ci_ = new Float64Array(L + 1), nr_ = new Float64Array(L + 1), ni_ = new Float64Array(L + 1);
    for (let i = 0; i < n; i++) {
      const k = 2 * PI / lam[i];
      // косинусы: индекс L — подложка, 0..L−1 — слои сверху вниз; внешняя среда: (1, cos φ0)
      for (let j = 0; j <= L; j++) {
        const Nv = j < L ? layers[j].N[i] : Nsub[i], Nr = Nv[0], Ni = Nv[1];
        const N2r = Nr * Nr - Ni * Ni, N2i = 2 * Nr * Ni, den = N2r * N2r + N2i * N2i;
        const ar = 1 - s02 * N2r / den, ai = s02 * N2i / den, mod = Math.hypot(ar, ai);
        let cre = Math.sqrt(Math.max(0, (mod + ar) / 2)), cim = Math.sqrt(Math.max(0, (mod - ar) / 2)); if (ai < 0) cim = -cim;
        cr_[j] = cre; ci_[j] = cim; nr_[j] = Nr; ni_[j] = Ni;
      }
      // старт: граница нижний слой (или воздух) | подложка
      let ur, ui, ucr, uci;
      if (L) { ur = nr_[L - 1]; ui = ni_[L - 1]; ucr = cr_[L - 1]; uci = ci_[L - 1]; } else { ur = 1; ui = 0; ucr = c0; uci = 0; }
      const br = nr_[L], bi = ni_[L], bcr = cr_[L], bci = ci_[L];
      // r_p = (Nb cu − Nu cb)/(Nb cu + Nu cb), r_s = (Nu cu − Nb cb)/(Nu cu + Nb cb)
      let ar = br * ucr - bi * uci, ai = br * uci + bi * ucr, cr = ur * bcr - ui * bci, ci = ur * bci + ui * bcr;
      let dr = ar + cr, di = ai + ci, dm = dr * dr + di * di, nr = ar - cr, ni = ai - ci;
      let rpr = (nr * dr + ni * di) / dm, rpi = (ni * dr - nr * di) / dm;
      let er = ur * ucr - ui * uci, ei = ur * uci + ui * ucr, fr = br * bcr - bi * bci, fi = br * bci + bi * bcr;
      dr = er + fr; di = ei + fi; dm = dr * dr + di * di; nr = er - fr; ni = ei - fi;
      let rsr = (nr * dr + ni * di) / dm, rsi = (ni * dr - nr * di) / dm;
      for (let j = L - 1; j >= 0; j--) {
        const Nr = nr_[j], Ni = ni_[j], cjr = cr_[j], cji = ci_[j], d = layers[j].d;
        // X = exp(−2iβ), β = k d N cos γ → −2iβ = 2 k d Im(Ncosγ) − 2i k d Re(Ncosγ)
        const ncr = Nr * cjr - Ni * cji, nci = Nr * cji + Ni * cjr, ex = Math.exp(2 * k * d * nci), ang = -2 * k * d * ncr;
        const Xr = ex * Math.cos(ang), Xi = ex * Math.sin(ang);
        let Ur, Ui, Ucr, Uci; if (j > 0) { Ur = nr_[j - 1]; Ui = ni_[j - 1]; Ucr = cr_[j - 1]; Uci = ci_[j - 1]; } else { Ur = 1; Ui = 0; Ucr = c0; Uci = 0; }
        // Френель верхняя среда U → слой j
        ar = Nr * Ucr - Ni * Uci; ai = Nr * Uci + Ni * Ucr; cr = Ur * cjr - Ui * cji; ci = Ur * cji + Ui * cjr;
        dr = ar + cr; di = ai + ci; dm = dr * dr + di * di; nr = ar - cr; ni = ai - ci;
        const rpur = (nr * dr + ni * di) / dm, rpui = (ni * dr - nr * di) / dm;
        er = Ur * Ucr - Ui * Uci; ei = Ur * Uci + Ui * Ucr; fr = Nr * cjr - Ni * cji; fi = Nr * cji + Ni * cjr;
        dr = er + fr; di = ei + fi; dm = dr * dr + di * di; nr = er - fr; ni = ei - fi;
        const rsur = (nr * dr + ni * di) / dm, rsui = (ni * dr - nr * di) / dm;
        // r ← (r_u + r X)/(1 + r_u r X)
        let rXr = rpr * Xr - rpi * Xi, rXi = rpr * Xi + rpi * Xr;
        nr = rpur + rXr; ni = rpui + rXi; dr = 1 + (rpur * rXr - rpui * rXi); di = rpur * rXi + rpui * rXr; dm = dr * dr + di * di;
        rpr = (nr * dr + ni * di) / dm; rpi = (ni * dr - nr * di) / dm;
        rXr = rsr * Xr - rsi * Xi; rXi = rsr * Xi + rsi * Xr;
        nr = rsur + rXr; ni = rsui + rXi; dr = 1 + (rsur * rXr - rsui * rXi); di = rsur * rXi + rsui * rXr; dm = dr * dr + di * di;
        rsr = (nr * dr + ni * di) / dm; rsi = (ni * dr - nr * di) / dm;
      }
      RP[i] = [rpr, rpi]; RS[i] = [rsr, rsi];
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
  const PARAM_NAMES = ['d', 'A', 'Auv', 'Eg', 'dRough', 'delta', 'Ak'];
  const K400_MAX = 0.05;                                                  // верхняя граница k хвоста при 400 нм в окончательном фите
  function makeResidFn(data, cfg, free, fixedVals, M) {
    const Nsub = cfg.Nsub, lam = data.lam, n = lam.length;
    return (x) => {
      const P = Object.assign({}, fixedVals); free.forEach((k, i) => P[k] = x[i]);
      const m = Object.assign({}, cfg.mat, { A: P.A, Auv: P.Auv, Eg: P.Eg, Ak: P.Ak || 0 });
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
      case 'Eg': return [Math.max(0.8, m.Eg - 1.0), Math.min(m.E0 - 0.15, m.Eg + 1.0), 0.05];
      case 'dRough': return [0, 12, 0.5];
      case 'delta': return [-0.20, 0.20, 0.01];
      case 'Ak': { const s = tailAk(m, 1e-3); return [0, tailAk(m, K400_MAX), s]; }      // k хвоста при 400 нм от 0 до K400_MAX
    }
  }
  /** Опорный фит однородной плёнки: грубый перебор d при справочной дисперсии → три лучших локальных минимума χ²(d) →
   *  ЛМ (аналитический якобиан, масштабированные параметры) по (d, A, Auv, Eg) из каждого; лучший по χ². */
  function referenceFit(data, cfg, prior) {
    const fixed = { dRough: 0, delta: 0, Ak: 0 };
    const base = { A: cfg.mat.A, Auv: cfg.mat.Auv, Eg: cfg.mat.Eg };
    const scan = makeResidJacFn(data, cfg, ['d'], Object.assign({}, fixed, base), 1);
    const step = Math.max(0.5, (prior.dMax - prior.dMin) / 120), grid = [], cost = [];
    for (let dg = prior.dMin; dg <= prior.dMax + 1e-9; dg += step) { const r = scan([dg], false).r; let c = 0; for (let i = 0; i < r.length; i++) c += r[i] * r[i]; grid.push(dg); cost.push(c); }
    const cands = [];
    for (let i = 0; i < grid.length; i++) if ((i === 0 || cost[i] <= cost[i - 1]) && (i === grid.length - 1 || cost[i] <= cost[i + 1])) cands.push(i);
    cands.sort((a, b) => cost[a] - cost[b]);
    const starts = cands.slice(0, 3).map(i => grid[i]), startCosts = cands.slice(0, 3).map(i => cost[i]);
    const free = ['d', 'A', 'Auv', 'Eg'];
    const b = free.map(k => boundsFor(k, cfg, prior));
    const margin = 0.3 * (prior.dMax - prior.dMin); b[0] = [prior.dMin - margin, prior.dMax + margin, b[0][2]];
    const fn = makeResidJacFn(data, cfg, free, fixed, 1);
    let res = null, nfev = 0;
    for (const d0 of starts) {
      const r = lmFitJac(fn, [d0, base.A, base.Auv, base.Eg], b.map(v => v[0]), b.map(v => v[1]), 80); nfev += r.nfev;
      if (!res || r.chi2 < res.chi2) res = r;
    }
    const P = {}; free.forEach((k, i) => P[k] = res.x[i]); P.dRough = 0; P.delta = 0; P.Ak = 0; P.k400 = 0;
    const E = {}; free.forEach((k, i) => E[k] = res.err[i]);
    return { P, E, resid: res.resid, chi2: res.chi2, free, starts, startCosts, nfev };
  }
  /** Фит в выбранной модели дефектов flags = {grad, rough, abs}: несколько стартов (толщины из перебора опорного фита × знак градиента),
   *  лучший по χ² дотягивается до сходимости; при |δ| > 4 % — повтор с 80 подслоями. Поглощение — сила хвоста A_k (старт из p0.Ak). */
  function finalFit(data, cfg, prior, flags, p0) {
    const free = ['d', 'A', 'Auv', 'Eg']; if (flags.rough) free.push('dRough'); if (flags.grad) free.push('delta'); if (flags.abs) free.push('Ak');
    const fixed = { dRough: 0, delta: 0, Ak: 0 };
    const b = free.map(k => boundsFor(k, cfg, prior));
    const margin = 0.3 * (prior.dMax - prior.dMin); b[0] = [prior.dMin - margin, prior.dMax + margin, b[0][2]];
    const lo = b.map(v => v[0]), hi = b.map(v => v[1]);
    const clampP = (P) => free.map((k, i) => Math.min(hi[i], Math.max(lo[i], P[k] === undefined ? 0 : P[k])));
    // альтернативные старты: толщины из перебора опорного фита, если их грубый χ² не хуже лучшего более чем в 5 раз; знак δ — если |δ| ≤ 1 %
    const sc = p0.startCosts || [], best0 = sc.length ? Math.min(...sc) : 0;
    const alt = (p0.starts || []).filter((v, i) => Math.abs(v - p0.d) > 1.0 && (!sc.length || sc[i] <= 5 * best0));
    const dStarts = [p0.d].concat(alt).slice(0, 3);
    const d0 = flags.grad ? (p0.delta || 0) : 0;
    const deltaStarts = flags.grad ? (Math.abs(d0) <= 0.01 ? (Math.abs(d0) > 1e-4 ? [d0, -d0] : [0.004, -0.004]) : [d0]) : [0];
    // разведка бассейнов на грубой модели (10 подслоёв), затем лучший старт дотягивается на 40 подслоях
    const fn10 = makeResidJacFn(data, cfg, free, fixed, flags.grad ? 10 : 40), fn40 = makeResidJacFn(data, cfg, free, fixed, 40);
    let best = null, nfev = 0, nStarts = 0;
    let base = Object.assign({}, p0, { dRough: Math.max(0.1, p0.dRough || 0.1), Ak: 0 });
    if (flags.abs) {
      // Поглощение: опорный фит без хвоста искажает дисперсию (A, E_g растут, компенсируя затухание полос), и старт из него
      // ведёт в ложный минимум. Раунд А: старты по силе хвоста (оценка ИИ и k₄₀₀ = 10⁻⁴, 10⁻³, 10⁻²) × дисперсия
      // (опорный фит / справочник) при основной толщине; лучший по χ² задаёт дисперсию и хвост для раунда Б.
      const kEst = p0.Ak > 0 ? tailK400(cfg.mat, p0.Ak) : 0, ks = [];
      for (const k of [kEst, 1e-4, 1e-3, 1e-2]) if (k >= 3e-5 && !ks.some(v => Math.abs(Math.log(v / k)) < 0.5)) ks.push(k);
      const disp = [{ A: p0.A, Auv: p0.Auv, Eg: p0.Eg }, { A: cfg.mat.A, Auv: cfg.mat.Auv, Eg: cfg.mat.Eg }];
      let bestA = null;
      for (const dp of disp) for (const k of ks) {
        const x0 = clampP(Object.assign({}, base, dp, { delta: deltaStarts[0], Ak: tailAk(cfg.mat, k) }));
        const r = lmFitJac(fn10, x0, lo, hi, 20); nfev += r.nfev; nStarts++;
        if (!bestA || r.chi2 < bestA.chi2) bestA = r;
      }
      const PA = {}; free.forEach((k, i) => PA[k] = bestA.x[i]);
      base = Object.assign({}, base, { A: PA.A, Auv: PA.Auv, Eg: PA.Eg, Ak: PA.Ak, dRough: Math.max(0.1, PA.dRough || 0.1) });
      if (dStarts.length === 1 && deltaStarts.length === 1) best = bestA;
    }
    const cands = [];
    for (const dS of dStarts) for (const dl of deltaStarts) cands.push(clampP(Object.assign({}, base, { d: dS, delta: dl })));
    if (cands.length > 1) cands.forEach((x0) => { const r = lmFitJac(fn10, x0, lo, hi, 25); nfev += r.nfev; nStarts++; if (!best || r.chi2 < best.chi2) best = r; });
    let res = lmFitJac(fn40, best ? best.x : cands[0], lo, hi, 80); nfev += res.nfev; nStarts = Math.max(1, nStarts);
    const iD = free.indexOf('delta');
    if (iD >= 0 && Math.abs(res.x[iD]) > 0.04) { const r80 = lmFitJac(makeResidJacFn(data, cfg, free, fixed, 80), res.x, lo, hi, 40); nfev += r80.nfev; res = r80; }   // 80 подслоёв
    const P = { dRough: 0, delta: 0, Ak: 0 }; free.forEach((k, i) => P[k] = res.x[i]);
    const E = {}; free.forEach((k, i) => E[k] = res.err[i]);
    if (flags.abs) { const h = Math.max(1e-6, 1e-3 * P.Ak), dk = (tailK400(cfg.mat, P.Ak + h) - tailK400(cfg.mat, P.Ak)) / h; P.k400 = tailK400(cfg.mat, P.Ak); E.k400 = E.Ak * dk; } else { P.k400 = 0; }
    return { P, E, resid: res.resid, chi2: res.chi2, free, nfev, starts: nStarts };
  }

  // ---------------------------------------------------------------- аналитический якобиан (Фурман–Тихонравов, разд. 1.4.1: формулы 1.4.9–1.4.16) и ЛМ с масштабированием
  function ccos(a) { return [Math.cos(a[0]) * Math.cosh(a[1]), -Math.sin(a[0]) * Math.sinh(a[1])]; }
  function csin(a) { return [Math.sin(a[0]) * Math.cosh(a[1]), Math.cos(a[0]) * Math.sinh(a[1])]; }
  function mmul(A, B) {
    return [[cadd(cmul(A[0][0], B[0][0]), cmul(A[0][1], B[1][0])), cadd(cmul(A[0][0], B[0][1]), cmul(A[0][1], B[1][1]))],
            [cadd(cmul(A[1][0], B[0][0]), cmul(A[1][1], B[1][0])), cadd(cmul(A[1][0], B[0][1]), cmul(A[1][1], B[1][1]))]];
  }
  const trace2 = (G, X) => cadd(cadd(cmul(G[0][0], X[0][0]), cmul(G[0][1], X[1][0])), cadd(cmul(G[1][0], X[0][1]), cmul(G[1][1], X[1][1])));
  const I2 = [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];

  /** Показатель плёнки N = n − ik и его производные по A, Auv, Eg, Ak (комплексные): N = sqrt(conj ε), dN/dp = conj(∂ε/∂p)/(2N).
   *  ∂ε/∂A, ∂ε/∂Auv и ∂ε/∂Ak аналитические (ε_TL линейна по силе осциллятора), ∂ε/∂Eg — центральная разность по скалярной функции. */
  function filmNJac(lam, m) {
    const n = lam.length, N = new Array(n), dA = new Array(n), dAuv = new Array(n), dEg = new Array(n), dAk = new Array(n), h = 1e-5, Ak = m.Ak || 0;
    for (let i = 0; i < n; i++) {
      const E = HC / lam[i], [e1, e2] = tlEps(E, m.A, m.E0, m.C, m.Eg), uvd = m.Euv * m.Euv - E * E;
      const [t1, t2] = tlEps(E, 1, m.E0, m.C, TAIL_ET);                                   // хвост единичной силы: ∂ε/∂Ak
      let Nv = csqrt([1 + e1 + m.Auv / uvd + (Ak ? Ak * t1 : 0), -(e2 + (Ak ? Ak * t2 : 0))]); Nv = [Math.abs(Nv[0]), -Math.abs(Nv[1])]; N[i] = Nv;
      const g = cdiv([0.5, 0], Nv);
      dA[i] = cmul(g, [e1 / m.A, -e2 / m.A]);
      dAuv[i] = cmul(g, [1 / uvd, 0]);
      dAk[i] = cmul(g, [t1, -t2]);
      const [p1, p2] = tlEps(E, m.A, m.E0, m.C, m.Eg + h), [q1, q2] = tlEps(E, m.A, m.E0, m.C, m.Eg - h);
      dEg[i] = cmul(g, [(p1 - q1) / (2 * h), -(p2 - q2) / (2 * h)]);
    }
    return { N, dN: { A: dA, Auv: dAuv, Eg: dEg, Ak: dAk } };
  }

  /** Слои сверху вниз с чувствительностями: {N, d, dN: {param: [complex]}, dNs: {param: scalar}, dd: {param: number}};
   *  чувствительность по параметру = dNs[param] · dN[param][i]. Геометрия совпадает с buildLayers. */
  function buildLayersJac(lam, P, mat, M) {
    const fj = filmNJac(lam, Object.assign({}, mat, { A: P.A, Auv: P.Auv, Eg: P.Eg, Ak: P.Ak || 0 }));
    const Nf = fj.N, delta = P.delta || 0, d = P.d, dRough = P.dRough || 0, n = lam.length, sub = delta ? M : 1;
    const film = [];
    for (let j = sub - 1; j >= 0; j--) {
      const z = sub === 1 ? 0.5 : (j + 0.5) / sub, sc = 1 + delta * (z - 0.5);
      film.push({ N: sc === 1 ? Nf : Nf.map(v => [v[0] * sc, v[1] * sc]), d: d / sub, dN: { A: fj.dN.A, Auv: fj.dN.Auv, Eg: fj.dN.Eg, Ak: fj.dN.Ak, delta: Nf }, dNs: { A: sc, Auv: sc, Eg: sc, Ak: sc, delta: z - 0.5 }, dd: { d: 1 / sub } });
    }
    const layers = [];
    if (dRough > 0) {
      const top = film[0], N = new Array(n), dN = { A: new Array(n), Auv: new Array(n), Eg: new Array(n), Ak: new Array(n), delta: new Array(n) };
      for (let i = 0; i < n; i++) {
        const Nt = top.N[i], eps = cmul(Nt, Nt);
        const b = [0.5 * (eps[0] + 1), 0.5 * eps[1]], sq = csqrt(cadd(cmul(b, b), [8 * eps[0], 8 * eps[1]]));
        let sgn = 1, e = [(b[0] + sq[0]) / 4, (b[1] + sq[1]) / 4]; if (e[0] <= 0) { sgn = -1; e = [(b[0] - sq[0]) / 4, (b[1] - sq[1]) / 4]; }
        let Ne = csqrt(e); Ne = [Math.abs(Ne[0]), -Math.abs(Ne[1])]; N[i] = Ne;
        const dEdEt = cadd([0.125, 0], cdiv([sgn * (b[0] + 8) / 8, sgn * b[1] / 8], sq));           // dε_ema/dε_top
        const chain = cmul(cmul(dEdEt, cdiv([0.5, 0], Ne)), [2 * Nt[0], 2 * Nt[1]]);               // dN_e/dN_top
        for (const key of ['A', 'Auv', 'Eg', 'Ak', 'delta']) { const v = top.dN[key][i], sc = top.dNs[key]; dN[key][i] = cmul(chain, [v[0] * sc, v[1] * sc]); }
      }
      layers.push({ N, d: dRough, dN, dNs: { A: 1, Auv: 1, Eg: 1, Ak: 1, delta: 1 }, dd: { dRough: 1 } });
    }
    return layers.concat(film);
  }

  /** r_p, r_s и ∂r/∂x по списку параметров free: матричный метод Абелеса; D_j (1.4.9), N_j (∂M_j/∂n_j), Ψ_r (1.4.10),
   *  ∂r/∂x_j = Tr(Ψ_r^{j−1} M^{j+1} X_j) (1.4.12, 1.4.16) и цепное правило к глобальным параметрам. r_p — в эллипсометрической конвенции.
   *  Реализация без промежуточных аллокаций: комплексные 2×2 матрицы — блоки по 8 чисел в Float64Array. */
  function rhoStackJac(lam, layers, Nsub, phi, free) {
    const s0 = Math.sin(phi * PI / 180), s02 = s0 * s0, c0 = Math.cos(phi * PI / 180);
    const nL = layers.length, nP = free.length, n = lam.length;
    const RP = new Array(n), RS = new Array(n), dRP = free.map(() => new Array(n)), dRS = free.map(() => new Array(n));
    const Mb = new Float64Array(8 * Math.max(1, nL)), Db = new Float64Array(8 * Math.max(1, nL)), Nb_ = new Float64Array(8 * Math.max(1, nL));
    const Pb = new Float64Array(8 * (nL + 1)), Bb = new Float64Array(8 * Math.max(1, nL)), T = new Float64Array(8), G = new Float64Array(8), PsiR = new Float64Array(8);
    const ddKey = layers.map(L => free.map(k => (L.dd && L.dd[k]) ? L.dd[k] : 0)), dNKey = layers.map(L => free.map(k => (L.dN && L.dN[k]) ? L.dN[k] : null)), dNs = layers.map(L => free.map(k => (L.dNs && L.dNs[k] !== undefined) ? L.dNs[k] : 1));
    // C = A × B для блоков (offsets a, b, c)
    const mm = (A, a, B, b, Cm, c) => {
      const a11r = A[a], a11i = A[a + 1], a12r = A[a + 2], a12i = A[a + 3], a21r = A[a + 4], a21i = A[a + 5], a22r = A[a + 6], a22i = A[a + 7];
      const b11r = B[b], b11i = B[b + 1], b12r = B[b + 2], b12i = B[b + 3], b21r = B[b + 4], b21i = B[b + 5], b22r = B[b + 6], b22i = B[b + 7];
      Cm[c] = a11r * b11r - a11i * b11i + a12r * b21r - a12i * b21i; Cm[c + 1] = a11r * b11i + a11i * b11r + a12r * b21i + a12i * b21r;
      Cm[c + 2] = a11r * b12r - a11i * b12i + a12r * b22r - a12i * b22i; Cm[c + 3] = a11r * b12i + a11i * b12r + a12r * b22i + a12i * b22r;
      Cm[c + 4] = a21r * b11r - a21i * b11i + a22r * b21r - a22i * b21i; Cm[c + 5] = a21r * b11i + a21i * b11r + a22r * b21i + a22i * b21r;
      Cm[c + 6] = a21r * b12r - a21i * b12i + a22r * b22r - a22i * b22i; Cm[c + 7] = a21r * b12i + a21i * b12r + a22r * b22i + a22i * b22r;
    };
    let trR = 0, trI = 0;
    const tr = (A, a, X, x) => {   // Tr(A X) = a11 x11 + a12 x21 + a21 x12 + a22 x22
      trR = A[a] * X[x] - A[a + 1] * X[x + 1] + A[a + 2] * X[x + 4] - A[a + 3] * X[x + 5] + A[a + 4] * X[x + 2] - A[a + 5] * X[x + 3] + A[a + 6] * X[x + 6] - A[a + 7] * X[x + 7];
      trI = A[a] * X[x + 1] + A[a + 1] * X[x] + A[a + 2] * X[x + 5] + A[a + 3] * X[x + 4] + A[a + 4] * X[x + 3] + A[a + 5] * X[x + 2] + A[a + 6] * X[x + 7] + A[a + 7] * X[x + 6];
    };
    const drR = new Float64Array(nP), drI = new Float64Array(nP);
    for (let i = 0; i < n; i++) {
      const k = 2 * PI / lam[i], Nbr = Nsub[i][0], Nbi = Nsub[i][1];
      // cos γ подложки: sqrt(1 − s0²/Nb²)
      let N2r = Nbr * Nbr - Nbi * Nbi, N2i = 2 * Nbr * Nbi, den = N2r * N2r + N2i * N2i, qr = s02 * N2r / den, qi = -s02 * N2i / den, ar = 1 - qr, ai = -qi, mod = Math.hypot(ar, ai);
      let cbr = Math.sqrt(Math.max(0, (mod + ar) / 2)), cbi = Math.sqrt(Math.max(0, (mod - ar) / 2)); if (ai < 0) cbi = -cbi;
      for (let pol = 0; pol < 2; pol++) {
        const isS = pol === 0, qa = isS ? c0 : 1 / c0;
        let qsr, qsi; if (isS) { qsr = Nbr * cbr - Nbi * cbi; qsi = Nbr * cbi + Nbi * cbr; } else { const dd = cbr * cbr + cbi * cbi; qsr = (Nbr * cbr + Nbi * cbi) / dd; qsi = (Nbi * cbr - Nbr * cbi) / dd; }
        for (let a = 0; a < nL; a++) {
          const Nr = layers[a].N[i][0], Ni = layers[a].N[i][1], dj = layers[a].d;
          N2r = Nr * Nr - Ni * Ni; N2i = 2 * Nr * Ni; den = N2r * N2r + N2i * N2i; qr = s02 * N2r / den; qi = -s02 * N2i / den; ar = 1 - qr; ai = -qi; mod = Math.hypot(ar, ai);
          let cr = Math.sqrt(Math.max(0, (mod + ar) / 2)), ci = Math.sqrt(Math.max(0, (mod - ar) / 2)); if (ai < 0) ci = -ci;    // cos γ_j
          const ncr = Nr * cr - Ni * ci, nci = Nr * ci + Ni * cr;                                        // N cos γ
          let qjr, qji; if (isS) { qjr = ncr; qji = nci; } else { const dd = cr * cr + ci * ci; qjr = (Nr * cr + Ni * ci) / dd; qji = (Ni * cr - Nr * ci) / dd; }
          const phr = k * dj * ncr, phi_ = k * dj * nci;
          const cpr = Math.cos(phr) * Math.cosh(phi_), cpi = -Math.sin(phr) * Math.sinh(phi_), spr = Math.sin(phr) * Math.cosh(phi_), spi = Math.cos(phr) * Math.sinh(phi_);
          const qq = qjr * qjr + qji * qji, ioqr = qji / qq, ioqi = qjr / qq;                             // i/q = (i q*)/|q|² → (qji + i qjr)/|q|²
          const iqr = -qji, iqi = qjr;                                                                    // i q
          const o = 8 * a;
          // M_j = [[cos φ, (i/q) sin φ], [i q sin φ, cos φ]]
          Mb[o] = cpr; Mb[o + 1] = cpi; Mb[o + 2] = ioqr * spr - ioqi * spi; Mb[o + 3] = ioqr * spi + ioqi * spr;
          Mb[o + 4] = iqr * spr - iqi * spi; Mb[o + 5] = iqr * spi + iqi * spr; Mb[o + 6] = cpr; Mb[o + 7] = cpi;
          // D_j = k N cos γ · [[−sin φ, (i/q) cos φ], [i q cos φ, −sin φ]]
          const dpr = k * ncr, dpi = k * nci;
          const x12r = ioqr * cpr - ioqi * cpi, x12i = ioqr * cpi + ioqi * cpr, x21r = iqr * cpr - iqi * cpi, x21i = iqr * cpi + iqi * cpr;
          Db[o] = -(dpr * spr - dpi * spi); Db[o + 1] = -(dpr * spi + dpi * spr); Db[o + 2] = dpr * x12r - dpi * x12i; Db[o + 3] = dpr * x12i + dpi * x12r;
          Db[o + 4] = dpr * x21r - dpi * x21i; Db[o + 5] = dpr * x21i + dpi * x21r; Db[o + 6] = Db[o]; Db[o + 7] = Db[o + 1];
          // N_j = f1 D_j + c [[0, −(i/q²) sin φ], [i sin φ, 0]],  f1 = d/(N cos²γ),  c = 1/cosγ (s) или (1 − tan²γ)/cosγ = (2cos²γ − 1)/cos³γ (p)
          const c2r = cr * cr - ci * ci, c2i = 2 * cr * ci;                                               // cos²γ
          const nc2r = Nr * c2r - Ni * c2i, nc2i = Nr * c2i + Ni * c2r, nc2m = nc2r * nc2r + nc2i * nc2i;
          const f1r = dj * nc2r / nc2m, f1i = -dj * nc2i / nc2m;
          let cfr, cfi;
          if (isS) { const cm = cr * cr + ci * ci; cfr = cr / cm; cfi = -ci / cm; }
          else { const numr = 2 * c2r - 1, numi = 2 * c2i, c3r = c2r * cr - c2i * ci, c3i = c2r * ci + c2i * cr, c3m = c3r * c3r + c3i * c3i; cfr = (numr * c3r + numi * c3i) / c3m; cfi = (numi * c3r - numr * c3i) / c3m; }
          // −(i/q²) sin φ = −i sin φ / q²;  q² = qjr² − qji² + 2i qjr qji
          const q2r = qjr * qjr - qji * qji, q2i = 2 * qjr * qji, q2m = q2r * q2r + q2i * q2i;
          const isr = -spi, isi = spr;                                                                     // i sin φ
          const t12r = -(isr * q2r + isi * q2i) / q2m, t12i = -(isi * q2r - isr * q2i) / q2m;             // −(i sin φ)/q²
          const T12r = cfr * t12r - cfi * t12i, T12i = cfr * t12i + cfi * t12r, T21r = cfr * isr - cfi * isi, T21i = cfr * isi + cfi * isr;
          for (let e = 0; e < 8; e += 2) { Nb_[o + e] = f1r * Db[o + e] - f1i * Db[o + e + 1]; Nb_[o + e + 1] = f1r * Db[o + e + 1] + f1i * Db[o + e]; }
          Nb_[o + 2] += T12r; Nb_[o + 3] += T12i; Nb_[o + 4] += T21r; Nb_[o + 5] += T21i;
        }
        // префиксы P_a = M_0 ⋯ M_{a−1}
        Pb[0] = 1; Pb[1] = 0; Pb[2] = 0; Pb[3] = 0; Pb[4] = 0; Pb[5] = 0; Pb[6] = 1; Pb[7] = 0;
        for (let a = 0; a < nL; a++) mm(Pb, 8 * a, Mb, 8 * a, Pb, 8 * (a + 1));
        const mo = 8 * nL, m11r = Pb[mo], m11i = Pb[mo + 1], m12r = Pb[mo + 2], m12i = Pb[mo + 3], m21r = Pb[mo + 4], m21i = Pb[mo + 5], m22r = Pb[mo + 6], m22i = Pb[mo + 7];
        const qaqsr = qa * qsr, qaqsi = qa * qsi;
        const dnr = qa * m11r + (qsr * m22r - qsi * m22i) + (qaqsr * m12r - qaqsi * m12i) + m21r, dni = qa * m11i + (qsr * m22i + qsi * m22r) + (qaqsr * m12i + qaqsi * m12r) + m21i;
        const nur = qa * m11r + (qaqsr * m12r - qaqsi * m12i) - (qsr * m22r - qsi * m22i) - m21r, nui = qa * m11i + (qaqsr * m12i + qaqsi * m12r) - (qsr * m22i + qsi * m22r) - m21i;
        const dm = dnr * dnr + dni * dni, rr = (nur * dnr + nui * dni) / dm, ri = (nui * dnr - nur * dni) / dm;
        // f = t/(2 q_a) = 1/den
        const fr = dnr / dm, fi = -dni / dm, omr_r = 1 - rr, omr_i = -ri, mop_r = -1 - rr, mop_i = -ri;
        // Ψ_r = f [[q_a(1−r), −(1+r)], [q_a q_s (1−r), −q_s (1+r)]]
        let tr_ = qa * omr_r, ti_ = qa * omr_i; PsiR[0] = fr * tr_ - fi * ti_; PsiR[1] = fr * ti_ + fi * tr_;
        PsiR[2] = fr * mop_r - fi * mop_i; PsiR[3] = fr * mop_i + fi * mop_r;
        tr_ = qaqsr * omr_r - qaqsi * omr_i; ti_ = qaqsr * omr_i + qaqsi * omr_r; PsiR[4] = fr * tr_ - fi * ti_; PsiR[5] = fr * ti_ + fi * tr_;
        tr_ = qsr * mop_r - qsi * mop_i; ti_ = qsr * mop_i + qsi * mop_r; PsiR[6] = fr * tr_ - fi * ti_; PsiR[7] = fr * ti_ + fi * tr_;
        // суффиксы B_a = M_{a+1} ⋯ M_{nL−1}
        if (nL) { const bo = 8 * (nL - 1); Bb[bo] = 1; Bb[bo + 1] = 0; Bb[bo + 2] = 0; Bb[bo + 3] = 0; Bb[bo + 4] = 0; Bb[bo + 5] = 0; Bb[bo + 6] = 1; Bb[bo + 7] = 0; for (let a = nL - 2; a >= 0; a--) mm(Mb, 8 * (a + 1), Bb, 8 * (a + 1), Bb, 8 * a); }
        drR.fill(0); drI.fill(0);
        for (let a = 0; a < nL; a++) {
          mm(PsiR, 0, Pb, 8 * a, T, 0); mm(Bb, 8 * a, T, 0, G, 0);
          tr(G, 0, Db, 8 * a); const ddr = trR, ddi = trI;
          tr(G, 0, Nb_, 8 * a); const dnr_ = trR, dni_ = trI;
          const ddk = ddKey[a], dNk = dNKey[a], dNsA = dNs[a];
          for (let p = 0; p < nP; p++) {
            if (ddk[p]) { drR[p] += ddr * ddk[p]; drI[p] += ddi * ddk[p]; }
            if (dNk[p]) { const v = dNk[p][i], sc = dNsA[p]; drR[p] += (dnr_ * v[0] - dni_ * v[1]) * sc; drI[p] += (dnr_ * v[1] + dni_ * v[0]) * sc; }
          }
        }
        if (isS) { RS[i] = [rr, ri]; for (let p = 0; p < nP; p++) dRS[p][i] = [drR[p], drI[p]]; }
        else { RP[i] = [-rr, -ri]; for (let p = 0; p < nP; p++) dRP[p][i] = [-drR[p], -drI[p]]; }
      }
    }
    return { RP, RS, dRP, dRS };
  }

  /** Ψ, Δ (град) и их производные: dΨ = ½ sin2Ψ · Re(ρ'/ρ), dΔ = Im(ρ'/ρ), ρ = r_p/r_s. */
  function psiDeltaJac(RP, RS, dRP, dRS) {
    const n = RP.length, nP = dRP.length, deg = 180 / PI;
    const psi = new Float64Array(n), del = new Float64Array(n), dpsi = dRP.map(() => new Float64Array(n)), ddel = dRP.map(() => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      const rp = RP[i], rs = RS[i], tanPsi = Math.sqrt((rp[0] * rp[0] + rp[1] * rp[1]) / (rs[0] * rs[0] + rs[1] * rs[1]));
      psi[i] = Math.atan(tanPsi) * deg;
      let d = Math.atan2(rp[1] * rs[0] - rp[0] * rs[1], rp[0] * rs[0] + rp[1] * rs[1]) * deg; if (d < 0) d += 360; del[i] = d;
      const s2 = 2 * tanPsi / (1 + tanPsi * tanPsi);
      for (let p = 0; p < nP; p++) { const g = csub(cdiv(dRP[p][i], rp), cdiv(dRS[p][i], rs)); dpsi[p][i] = 0.5 * s2 * g[0] * deg; ddel[p][i] = g[1] * deg; }
    }
    return { psi, del, dpsi, ddel };
  }

  /** Невязка и аналитический якобиан для ЛМ: fn(x, needJ) → {r, J?, psi, del}; при needJ = false — только невязка (быстрый путь). */
  function makeResidJacFn(data, cfg, free, fixedVals, M) {
    const lam = data.lam, n = lam.length, Mm = M || 40;
    return (x, needJ) => {
      const P = Object.assign({}, fixedVals); free.forEach((k, i) => P[k] = x[i]);
      const r = new Float64Array(2 * n);
      if (needJ === false) {
        const m = Object.assign({}, cfg.mat, { A: P.A, Auv: P.Auv, Eg: P.Eg, Ak: P.Ak || 0 });
        const [RP, RS] = rhoStack(lam, buildLayers(lam, P.d, m, P.delta || 0, P.dRough || 0, Mm), cfg.Nsub, cfg.phi);
        const { psi, del } = psiDeltaFromR(RP, RS);
        for (let i = 0; i < n; i++) { r[i] = (data.psi[i] - psi[i]) / data.sigPsi[i]; r[n + i] = wrap(data.del[i] - del[i]) / data.sigDel[i]; }
        return { r, psi, del };
      }
      const layers = buildLayersJac(lam, P, cfg.mat, Mm);
      const { RP, RS, dRP, dRS } = rhoStackJac(lam, layers, cfg.Nsub, cfg.phi, free);
      const { psi, del, dpsi, ddel } = psiDeltaJac(RP, RS, dRP, dRS);
      const J = free.map(() => new Float64Array(2 * n));
      for (let i = 0; i < n; i++) {
        r[i] = (data.psi[i] - psi[i]) / data.sigPsi[i]; r[n + i] = wrap(data.del[i] - del[i]) / data.sigDel[i];
        for (let p = 0; p < free.length; p++) { J[p][i] = -dpsi[p][i] / data.sigPsi[i]; J[p][n + i] = -ddel[p][i] / data.sigDel[i]; }
      }
      return { r, J, psi, del };
    };
  }

  /** ЛМ с аналитическим якобианом. Параметры масштабированы преобразованием x = lo + (hi − lo)(1 + sin θ)/2:
   *  θ безразмерны и O(1), границы соблюдаются гладко (без отсечения). Демпфирование по Нильсену (gain ratio),
   *  диагональное масштабирование Марквардта с нижним порогом; пробные шаги — без якобиана. Ковариация — по J_x в решении. */
  function lmFitJac(residJac, x0, lo, hi, maxIter) {
    maxIter = maxIter || 60;
    const p = x0.length, W = lo.map((v, k) => hi[k] - v);
    const toX = (th) => th.map((v, k) => lo[k] + W[k] * (1 + Math.sin(v)) / 2);
    const toTh = (x) => x.map((v, k) => Math.asin(Math.max(-0.999999, Math.min(0.999999, 2 * (v - lo[k]) / W[k] - 1))));
    let th = toTh(x0), x = toX(th), ev = residJac(x, true), r = ev.r, Jx = ev.J;
    const m = r.length;
    let F = 0; for (let i = 0; i < m; i++) F += r[i] * r[i]; F *= 0.5;
    let lam = -1, nu = 2, nfev = 1, nit = 0;
    const JTJ = Array.from({ length: p }, () => new Array(p).fill(0)), g = new Array(p).fill(0), D = new Array(p).fill(0);
    for (let it = 0; it < maxIter; it++) {
      nit++;
      const sc = th.map((v, k) => W[k] / 2 * Math.cos(v));                     // dx/dθ
      let maxD = 0;
      for (let a = 0; a < p; a++) {
        for (let b = a; b < p; b++) { let s = 0; const Ja = Jx[a], Jb = Jx[b]; for (let i = 0; i < m; i++) s += Ja[i] * Jb[i]; s *= sc[a] * sc[b]; JTJ[a][b] = s; JTJ[b][a] = s; }
        let s = 0; const Ja = Jx[a]; for (let i = 0; i < m; i++) s += Ja[i] * r[i]; g[a] = s * sc[a];
        maxD = Math.max(maxD, JTJ[a][a]);
      }
      for (let a = 0; a < p; a++) D[a] = Math.max(JTJ[a][a], 1e-8 * maxD, 1e-300);
      if (lam < 0) lam = 1e-3;
      let accepted = false;
      for (let tries = 0; tries < 15; tries++) {
        const Aug = JTJ.map((row, a) => row.map((v, b) => (a === b ? v + lam * D[a] : v)));
        const dth = solveLinear(Aug, g.map(v => -v));
        if (!dth) { lam *= nu; nu *= 2; continue; }
        // геодезическое ускорение (Transtrum–Sethna): вторая производная невязки вдоль шага по разностям, поправка ½a
        let vnorm = 0; for (let k = 0; k < p; k++) vnorm += dth[k] * dth[k]; vnorm = Math.sqrt(vnorm);
        if (vnorm > 1e-12) {
          const hh = Math.min(0.1, 0.1 / vnorm);
          const rP = residJac(toX(th.map((v, k) => v + hh * dth[k])), false).r, rM = residJac(toX(th.map((v, k) => v - hh * dth[k])), false).r; nfev += 2;
          const rhs = new Array(p).fill(0);
          for (let a = 0; a < p; a++) { let sacc = 0; const Ja = Jx[a]; for (let i = 0; i < m; i++) sacc += Ja[i] * (rP[i] - 2 * r[i] + rM[i]); rhs[a] = -sacc * sc[a] / (hh * hh); }
          const acc = solveLinear(Aug, rhs);
          if (acc) { let anorm = 0; for (let k = 0; k < p; k++) anorm += acc[k] * acc[k]; anorm = Math.sqrt(anorm); if (anorm < 0.75 * vnorm) for (let k = 0; k < p; k++) dth[k] += 0.5 * acc[k]; }
        }
        let maxStep = 0; for (let k = 0; k < p; k++) { dth[k] = Math.max(-1.0, Math.min(1.0, dth[k])); maxStep = Math.max(maxStep, Math.abs(dth[k])); }
        const thn = th.map((v, k) => v + dth[k]), xn = toX(thn), evn = residJac(xn, false); nfev++;
        let Fn = 0; for (let i = 0; i < m; i++) Fn += evn.r[i] * evn.r[i]; Fn *= 0.5;
        let pred = 0; for (let k = 0; k < p; k++) pred += dth[k] * (lam * D[k] * dth[k] - g[k]); pred *= 0.5;
        const rho = pred > 0 ? (F - Fn) / pred : (Fn < F ? 1 : -1);
        if (Fn < F) {
          const rel = (F - Fn) / Math.max(F, 1e-300);
          th = thn; x = xn; r = evn.r; F = Fn; accepted = true;
          lam *= Math.max(1 / 3, 1 - Math.pow(2 * rho - 1, 3)); nu = 2;
          if (rel < 1e-8 || maxStep < 1e-7) { it = maxIter; } else { Jx = residJac(x, true).J; }
          break;
        }
        lam *= nu; nu *= 2;
      }
      if (!accepted) break;
    }
    Jx = residJac(x, true).J;                                                  // якобиан в решении — для ковариации
    const nuDof = Math.max(1, m - p), chi2 = 2 * F / nuDof;
    const JTJx = Array.from({ length: p }, () => new Array(p).fill(0));
    for (let a = 0; a < p; a++) for (let b = a; b < p; b++) { let s = 0; for (let i = 0; i < m; i++) s += Jx[a][i] * Jx[b][i]; JTJx[a][b] = s; JTJx[b][a] = s; }
    const inv = invert(JTJx); const cov = inv ? inv.map(row => row.map(v => v * chi2)) : null;
    return { x, resid: r, chi2, cov, nfev, nit, err: cov ? cov.map((row, k) => Math.sqrt(Math.max(0, row[k]))) : x.map(() => NaN) };
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
  /** Каналы остатка по якорям. В полуволновых точках однородная прозрачная плёнка «отсутствует» и Ψ = Ψ подложки при любых n и d,
   *  поэтому Ψ-остаток в HW-якорях не зависит от того, как опорный фит подобрал дисперсию, и читает дефекты напрямую:
   *  градиент сдвигает Ψ в HW-якорях примерно одинаково по всему диапазону (канал psiHW — среднее), поглощение — гасит контраст
   *  полос тем сильнее, чем ближе к синему краю (∝ k(λ)·d/λ), поэтому его канал — перепад Ψ-остатка между синей и красной половинами
   *  якорей: psiHWtrend = ⟨r_Ψ(HW)⟩_синие − ⟨r_Ψ(HW)⟩_красные; в четвертьволновых точках контраст гасится с противоположной стороны
   *  (psiQWtrend — тот же перепад по QW-якорям). Половины — по середине диапазона λ; если одна пуста — по номеру якоря. */
  function channels(lam, resid, HW, QW) {
    const n = lam.length, rp = Array.from(resid.slice(0, n)), rd = Array.from(resid.slice(n));
    const at = (xs, ys) => xs.map(x => interp(x, lam, ys));
    const psiHW = at(HW, rp), delHW = at(HW, rd), delQW = at(QW, rd), psiQW = at(QW, rp);
    const mean = a => a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN, mabs = a => mean(a.map(Math.abs));
    const rms = a => Math.sqrt(a.reduce((s, v) => s + v * v, 0) / a.length);
    const mid = (lam[0] + lam[n - 1]) / 2;
    const trend = (xs, ys) => {
      if (xs.length < 2) return 0;
      let blue = [], red = []; xs.forEach((x, j) => (x < mid ? blue : red).push(ys[j]));
      if (!blue.length || !red.length) { const h = Math.floor(xs.length / 2); blue = ys.slice(0, h); red = ys.slice(h); }
      return mean(blue) - mean(red);
    };
    return { psiHW: mean(psiHW), psiHWsd: Math.sqrt(Math.max(0, mean(psiHW.map(v => v * v)) - mean(psiHW) ** 2)), delHW: mabs(delHW), delQW: mabs(delQW), rmsPsi: rms(rp), rmsDel: rms(rd),
      psiHWtrend: trend(HW, psiHW), psiQWtrend: trend(QW, psiQW), nHW: HW.length };
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

  // ---------------------------------------------------------------- классы дефектов: градиент g ∈ {0, +1, −1} × шероховатость r ∈ {0, 1} × поглощение a ∈ {0, 1}
  const CLASSES = (() => {
    const base = [{ g: 0, r: 0, name: 'однородная плёнка' }, { g: 1, r: 0, name: 'градиент, n растёт к поверхности' }, { g: -1, r: 0, name: 'градиент, n падает к поверхности' },
      { g: 0, r: 1, name: 'шероховатость' }, { g: 1, r: 1, name: 'градиент (n растёт) + шероховатость' }, { g: -1, r: 1, name: 'градиент (n падает) + шероховатость' }];
    const abs = base.map(c => ({ g: c.g, r: c.r, a: 1, name: c.g === 0 && c.r === 0 ? 'поглощение' : c.name.replace('градиент, n растёт к поверхности', 'градиент (n растёт)').replace('градиент, n падает к поверхности', 'градиент (n падает)') + ' + поглощение' }));
    return base.map(c => Object.assign({ a: 0 }, c)).concat(abs);
  })();
  /** Индекс класса по набору дефектов. */
  const classIndex = (g, r, a) => CLASSES.findIndex(c => c.g === g && c.r === r && c.a === a);
  const N_CHANNELS = 8;

  /** Вектор признаков: остаток опорного фита + 8 чисел каналов якорей (Ψ в HW, разброс, |Δ| в HW и QW, RMS Ψ/Δ,
   *  перепад Ψ-остатка синие − красные якоря HW и QW — затухание контраста полос к синему краю). */
  function featureVector(dd, cfg, resid) {
    const subDel = cfg.subDel || (cfg.subDel = modelPsiDelta(dd.lam, cfg.Nsub, cfg.phi, 0, cfg.mat, 0, 0, 1).del);
    const { HW, QW } = anchors(dd.lam, dd.del, subDel);
    let ch = { psiHW: 0, psiHWsd: 0, delHW: 0, delQW: 0, rmsPsi: 0, rmsDel: 0, psiHWtrend: 0, psiQWtrend: 0 };
    if (HW.length >= 2 && QW.length >= 1) ch = channels(dd.lam, resid, HW, QW);
    const f = (v) => (isFinite(v) ? v : 0);
    return Array.from(resid).concat([f(ch.psiHW), f(ch.psiHWsd), f(ch.delHW), f(ch.delQW), f(ch.rmsPsi), f(ch.rmsDel), f(ch.psiHWtrend), f(ch.psiQWtrend)]);
  }
  /** Случайный размер дефекта в двух шкалах: с вероятностью ½ в тонкой [lo, small], иначе в [small, max] (если max > small). */
  function twoScale(rng, lo, small, max) { small = Math.min(small, max); return (max > small + 1e-12 && rng.uniform() < 0.5) ? rng.range(small, max) : rng.range(lo, small); }
  /** Один обучающий пример: истина со случайными дефектами и мешающими факторами → опорный фит → остаток.
   *  nuis = {deltaMax, roughMax, kMax, bw, phiErr}; kMax — верхняя граница k хвоста поглощения при 400 нм. */
  function trainingExample(data, cfg, prior, nuis, rng, cls) {
    const c = CLASSES[cls];
    const delta = c.g * twoScale(rng, 0.0002, 0.02, nuis.deltaMax), dRough = c.r ? rng.range(0.02, nuis.roughMax) : 0;
    const k400 = c.a ? twoScale(rng, 2e-5, 1e-3, nuis.kMax || 5e-3) : 0;
    const d = rng.range(prior.dMin, prior.dMax);
    const m = Object.assign({}, cfg.mat, { A: cfg.mat.A * (1 + rng.range(-0.05, 0.05)), Auv: cfg.mat.Auv * (1 + rng.range(-0.05, 0.05)) + rng.range(-2, 2), E0: cfg.mat.E0 + rng.range(-0.05, 0.05), C: cfg.mat.C + rng.range(-0.06, 0.06), Eg: cfg.mat.Eg + rng.range(-0.05, 0.05) });
    if (m.Auv < 0) m.Auv = 0;
    if (k400 > 0) m.Ak = tailAk(m, k400);
    const opt = { bw: Math.max(0, nuis.bw + rng.range(-0.3, 0.3)), lamOffset: rng.range(-0.1, 0.1), driftPsi: [rng.range(-0.01, 0.01), rng.range(-0.015, 0.015)], driftDel: [rng.range(-0.03, 0.03), rng.range(-0.04, 0.04)], sigPsi: data.sigPsi, sigDel: data.sigDel };
    const s = synthesize(data.lam, cfg.subKey, cfg.phi + rng.range(-nuis.phiErr, nuis.phiErr), d, m, delta, dRough, opt, rng);
    const dd = { lam: data.lam, psi: s.psi, del: s.del, sigPsi: data.sigPsi, sigDel: data.sigDel };
    const fit = referenceFit(dd, cfg, prior);
    return { x: featureVector(dd, cfg, fit.resid), y: cls, delta, dRough, k400, d, dFit: fit.P.d, chi2: fit.chi2 };
  }

  const api = { HC, TAIL_ET, K400_MAX, N_CHANNELS, MATERIALS, SUBSTRATES, CLASSES, classIndex, filmN, filmNK, tailK400, tailAk, substrateN, buildLayers, rhoStack, psiDeltaFromR, modelPsiDelta, rhoGraded, synthesizeClean, applyNoise, synthesize, makeRng, lmFit, lmFitJac, makeResidJacFn, rhoStackJac, buildLayersJac, filmNJac, psiDeltaJac, referenceFit, finalFit, anchors, channels, featureVector, trainSoftmax, ridge, trainingExample, wrap, interp, tlEps };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.EllipCore = api;
})(typeof window !== 'undefined' ? window : globalThis);
