/* Проверка физики ядра по книге Sh. A. Furman, A. V. Tikhonravov, «Basics of Optics of Multilayer Systems» (1992):
   1) матричный метод Абелеса (формулы 1.2.9, 1.2.13, 1.2.14) — независимая реализация;
   2) воспроизведение таблицы 1.2 (17-слойный узкополосный фильтр S (HL)^4 2H (LH)^4, λ0 = 500 нм);
   3) сравнение рекурсии Эйри ядра (rhoStack) с матричным методом при наклонном падении, поглощающих слоях и подложке;
   4) уравнение Риккати для локальной функции отражения (1.1.19, 1.1.20): непрерывный градиент n(z) против лестницы подслоёв;
   5) закон сохранения энергии R + T = 1 для непоглощающей системы.
   Запуск: node test/verify_book.js */
const E = require("../src/ellip_core.js");
const PI = Math.PI;
// ---- комплексная арифметика (конвенция книги: exp(+iωt), N = n − ik)
const C = {
  add: (a, b) => [a[0] + b[0], a[1] + b[1]], sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
  mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
  div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
  scale: (a, s) => [a[0] * s, a[1] * s], abs2: (a) => a[0] * a[0] + a[1] * a[1], abs: (a) => Math.hypot(a[0], a[1]),
  sqrt: (a) => { const r = Math.hypot(a[0], a[1]); let re = Math.sqrt(Math.max(0, (r + a[0]) / 2)), im = Math.sqrt(Math.max(0, (r - a[0]) / 2)); if (a[1] < 0) im = -im; return [re, im]; },
  cos: (a) => [Math.cos(a[0]) * Math.cosh(a[1]), -Math.sin(a[0]) * Math.sinh(a[1])],
  sin: (a) => [Math.sin(a[0]) * Math.cosh(a[1]), Math.cos(a[0]) * Math.sinh(a[1])],
  I: [0, 1], ONE: [1, 0],
};
const cosGamma = (N, alpha2) => C.sqrt(C.sub(C.ONE, C.div([alpha2, 0], C.mul(N, N))));   // cos γ_j = sqrt(1 − α²/n_j²), главная ветвь

/** Матричный метод книги. layers: [{N: [n,−k], d}] сверху вниз (как в ядре); возвращает {rs, rp, ts, tp, Rs, Rp, Ts, Tp} на одной λ. */
function abeles(lam, layers, Nsub, phiDeg, Nout) {
  Nout = Nout || [1, 0];
  const alpha = Nout[0] * Math.sin(phiDeg * PI / 180), alpha2 = alpha * alpha, k = 2 * PI / lam;
  const out = {};
  for (const pol of ["s", "p"]) {
    const q = (N) => { const c = cosGamma(N, alpha2); return pol === "s" ? C.mul(N, c) : C.div(N, c); };
    // M = M_m ... M_1, сомножители от внешней среды к подложке; слои у нас перечислены сверху вниз
    let M = [[C.ONE, [0, 0]], [[0, 0], C.ONE]];
    for (const L of layers) {
      const c = cosGamma(L.N, alpha2), phi = C.scale(C.mul(L.N, c), k * L.d), qj = q(L.N);
      const cp = C.cos(phi), sp = C.sin(phi);
      const Mj = [[cp, C.div(C.mul(C.I, sp), qj)], [C.mul(C.mul(C.I, qj), sp), cp]];
      M = [[C.add(C.mul(M[0][0], Mj[0][0]), C.mul(M[0][1], Mj[1][0])), C.add(C.mul(M[0][0], Mj[0][1]), C.mul(M[0][1], Mj[1][1]))],
           [C.add(C.mul(M[1][0], Mj[0][0]), C.mul(M[1][1], Mj[1][0])), C.add(C.mul(M[1][0], Mj[0][1]), C.mul(M[1][1], Mj[1][1]))]];
    }
    const qa = q(Nout), qs = q(Nsub);
    const [m11, m12, m21, m22] = [M[0][0], M[0][1], M[1][0], M[1][1]];
    const den = C.add(C.add(C.mul(qa, m11), C.mul(qs, m22)), C.add(C.mul(C.mul(qa, qs), m12), m21));
    const num = C.sub(C.add(C.mul(qa, m11), C.mul(C.mul(qa, qs), m12)), C.add(C.mul(qs, m22), m21));
    const r = C.div(num, den), t = C.div(C.scale(qa, 2), den);
    out["r" + pol] = r; out["t" + pol] = t; out["R" + pol] = C.abs2(r); out["T" + pol] = qs[0] / qa[0] * C.abs2(t);
  }
  return out;
}

/** Уравнение Риккати (1.1.19)/(1.1.20) для локальной функции отражения; профиль eps(z), z от подложки (0) к поверхности (za). RK4. */
function riccati(lam, epsOf, za, Nsub, phiDeg, h) {
  const alpha = Math.sin(phiDeg * PI / 180), alpha2 = alpha * alpha, k = 2 * PI / lam;
  const out = {};
  for (const pol of ["s", "p"]) {
    const ca = Math.sqrt(1 - alpha2), qa = pol === "s" ? [ca, 0] : [1 / ca, 0];
    const cs = cosGamma(Nsub, alpha2), qs = pol === "s" ? C.mul(Nsub, cs) : C.div(Nsub, cs);
    const f = (z, r) => {        // dr/dz = (ik/2qa) [ qa² (1−r)² (·) − (·)(1+r)² ]
      const eps = epsOf(z), m = C.mul(C.sub(C.ONE, r), C.sub(C.ONE, r)), p = C.mul(C.add(C.ONE, r), C.add(C.ONE, r));
      let term;
      if (pol === "s") term = C.sub(C.scale(m, qa[0] * qa[0]), C.mul(C.sub(eps, [alpha2, 0]), p));
      else term = C.sub(C.scale(C.mul(C.sub(C.ONE, C.div([alpha2, 0], eps)), m), qa[0] * qa[0]), C.mul(eps, p));
      return C.mul([0, k / (2 * qa[0])], term);
    };
    let r = C.div(C.sub(qa, qs), C.add(qa, qs));
    const n = Math.ceil(za / h), dz = za / n;
    for (let i = 0; i < n; i++) {
      const z = i * dz;
      const k1 = f(z, r), k2 = f(z + dz / 2, C.add(r, C.scale(k1, dz / 2))), k3 = f(z + dz / 2, C.add(r, C.scale(k2, dz / 2))), k4 = f(z + dz, C.add(r, C.scale(k3, dz)));
      r = C.add(r, C.scale(C.add(C.add(k1, C.scale(k2, 2)), C.add(C.scale(k3, 2), k4)), dz / 6));
    }
    out["r" + pol] = r;
  }
  return out;
}
const psiDelta = (rp, rs) => { const psi = Math.atan(Math.sqrt(C.abs2(rp) / C.abs2(rs))) * 180 / PI; const pr = C.mul(rp, [rs[0], -rs[1]]); let d = Math.atan2(pr[1], pr[0]) * 180 / PI; if (d < 0) d += 360; return [psi, d]; };
const wrap = (v) => ((v + 180) % 360 + 360) % 360 - 180;
let failures = 0; const check = (ok, msg) => { console.log((ok ? "  ok   " : "  FAIL ") + msg); if (!ok) failures++; };

// ================================================================ 1–2. Таблица 1.2 книги
console.log("1. Таблица 1.2 (17 слоёв S (HL)^4 xH (LH)^4, n_s = 1.52, n_H = 2.3 − 0.0002i, n_L = 1.35 − 0.0002i, λ0 = 500 нм, нормальное падение)");
const nH = [2.3, -0.0002], nL = [1.35, -0.0002], ns = [1.52, 0];
const dH = 500 / (4 * 2.3), dL = 500 / (4 * 1.35);
const bottomUp = [];                                    // от подложки к воздуху
for (let i = 0; i < 4; i++) bottomUp.push({ N: nH, d: dH }, { N: nL, d: dL });
bottomUp.push({ N: nH, d: 4 * dH });                    // спейсер с оптической толщиной λ0 (см. примечание в отчёте)
for (let i = 0; i < 4; i++) bottomUp.push({ N: nL, d: dL }, { N: nH, d: dH });
const topDown = bottomUp.slice().reverse();
const table = [[480, 0.146, 99.721, 0.133], [485, 0.238, 99.616, 0.146], [492, 0.780, 98.982, 0.238], [497, 5.196, 93.786, 1.018], [499, 31.286, 63.049, 5.665],
  [500, 83.401, 1.597, 15.002], [501, 31.429, 62.827, 5.744], [503, 5.312, 93.624, 1.064], [508, 0.829, 98.917, 0.254], [515, 0.264, 99.585, 0.151], [520, 0.165, 99.702, 0.133]];
let maxDev = 0, maxDevT = 0, maxDevApp = 0;
console.log("   λ, нм | T книга / расчёт | R книга / расчёт (свет из подложки) | A книга / расчёт | R прямого хода: матрица / ядро");
for (const [lam, Tb, Rb, Ab] of table) {
  const rev = abeles(lam, bottomUp, [1, 0], 0, ns);      // обратный ход: внешняя среда 1.52, выход в воздух
  const fwd = abeles(lam, topDown, ns, 0);               // прямой ход: из воздуха
  const T = 100 * rev.Ts, R = 100 * rev.Rs, A = 100 - T - R;
  const [RP, RS] = E.rhoStack([lam], topDown.map(L => ({ N: [L.N], d: L.d })), [ns], 1e-6);
  const Rapp = 100 * C.abs2(RS[0]);
  maxDev = Math.max(maxDev, Math.abs(T - Tb), Math.abs(R - Rb), Math.abs(A - Ab)); maxDevT = Math.max(maxDevT, Math.abs(100 * fwd.Ts - Tb)); maxDevApp = Math.max(maxDevApp, Math.abs(Rapp - 100 * fwd.Rs));
  console.log(`   ${lam}   | ${Tb.toFixed(3)} / ${T.toFixed(3)} | ${Rb.toFixed(3)} / ${R.toFixed(3)} | ${Ab.toFixed(3)} / ${A.toFixed(3)} | ${(100 * fwd.Rs).toFixed(3)} / ${Rapp.toFixed(3)}`);
}
check(maxDev < 0.002, `матричный метод воспроизводит таблицу 1.2 (T, R, A) с точностью ${maxDev.toFixed(4)} % — при свете со стороны подложки`);
check(maxDevT < 0.002, `T не зависит от направления хода (разд. 1.3.2): прямой ход даёт те же T с точностью ${maxDevT.toFixed(4)} %`);
check(maxDevApp < 1e-9, `рекурсия Эйри ядра совпадает с матричным методом по R прямого хода: макс. отклонение ${maxDevApp.toExponential(1)} %`);

// ================================================================ 3. Ядро против матричного метода при наклонном падении
console.log("\n2. rhoStack против матричного метода: TiO2 1000 нм с градиентом +5 % (40 подслоёв) и EMA-слоем 2 нм на c-Si, 65°, 400–800 нм");
{
  const lam = []; for (let l = 400; l <= 800; l += 4) lam.push(l);
  const m = E.MATERIALS.tio2, Nsub = E.substrateN("si", lam);
  const layers = E.buildLayers(lam, 1000, m, 0.05, 2.0, 40);
  const [RP, RS] = E.rhoStack(lam, layers, Nsub, 65);
  let dRs = 0, dRpSign = 0, dPsi = 0, dDel = 0;
  for (let i = 0; i < lam.length; i++) {
    const a = abeles(lam[i], layers.map(L => ({ N: L.N[i], d: L.d })), Nsub[i], 65);
    dRs = Math.max(dRs, C.abs(C.sub(a.rs, RS[i])));
    dRpSign = Math.max(dRpSign, C.abs(C.add(a.rp, RP[i])));           // r_p ядра = −r_p книги (определение через H, а не через E_танг)
    const [ps1, d1] = psiDelta(RP[i], RS[i]); const [ps2, d2] = psiDelta(C.scale(a.rp, -1), a.rs);
    dPsi = Math.max(dPsi, Math.abs(ps1 - ps2)); dDel = Math.max(dDel, Math.abs(wrap(d1 - d2)));
  }
  check(dRs < 1e-10, `r_s совпадает: макс. |Δr_s| = ${dRs.toExponential(1)}`);
  check(dRpSign < 1e-10, `r_p ядра равен −r_p книги (эллипсометрическая конвенция знака p): макс. невязка ${dRpSign.toExponential(1)}`);
  check(dPsi < 1e-8 && dDel < 1e-8, `Ψ, Δ совпадают после учёта знака: max|ΔΨ| = ${dPsi.toExponential(1)}°, max|ΔΔ| = ${dDel.toExponential(1)}°`);
  // Δ книги отличается на 180°
  const [, dCore] = psiDelta(RP[10], RS[10]); const a = abeles(lam[10], layers.map(L => ({ N: L.N[10], d: L.d })), Nsub[10], 65); const [, dBook] = psiDelta(a.rp, a.rs);
  console.log(`   при λ = ${lam[10]} нм: Δ ядра ${dCore.toFixed(2)}°, Δ в конвенции книги ${dBook.toFixed(2)}° (разница 180°)`);
}

// ================================================================ 4. Риккати: непрерывный градиент против лестницы
console.log("\n3. Непрерывный линейный профиль n(z) (уравнение Риккати, RK4, шаг 0.25 нм) против лестницы подслоёв; TiO2 1000 нм на BK7, 65°");
{
  const lam = []; for (let l = 400; l <= 800; l += 8) lam.push(l);
  const m = E.MATERIALS.tio2, Nsub = E.substrateN("bk7", lam), Nf = E.filmN(lam, m);
  const sigP = lam.map(l => 0.015 + 0.02 * Math.max(0, (650 - l) / 250) ** 2), sigD = sigP.map(v => 2.2 * v);
  for (const delta of [0.01, 0.05, 0.12]) {
    const refPD = [];
    const t0 = Date.now();
    for (let i = 0; i < lam.length; i++) {
      const N0 = Nf[i], d = 1000;
      const epsOf = (z) => { const s = 1 + delta * (z / d - 0.5); const N = C.scale(N0, s); return C.mul(N, N); };
      const r = riccati(lam[i], epsOf, d, Nsub[i], 65, 0.25);
      refPD.push(psiDelta(C.scale(r.rp, -1), r.rs));
    }
    const ms = Date.now() - t0;
    const row = [];
    for (const M of [10, 20, 40, 80, 160]) {
      const pd = E.modelPsiDelta(lam, Nsub, 65, 1000, m, delta, 0, M);
      let sp = 0, sd = 0; for (let i = 0; i < lam.length; i++) { sp += ((pd.psi[i] - refPD[i][0]) / sigP[i]) ** 2; sd += (wrap(pd.del[i] - refPD[i][1]) / sigD[i]) ** 2; }
      row.push(`M=${M}: Ψ ${Math.sqrt(sp / lam.length).toFixed(2)}σ, Δ ${Math.sqrt(sd / lam.length).toFixed(2)}σ`);
    }
    console.log(`   δ = ${(100 * delta).toFixed(0)} % (Риккати ${ms} мс на ${lam.length} λ): RMS отклонения лестницы от непрерывного профиля — ${row.join("; ")}`);
  }
  // сходимость Риккати по шагу
  const N0 = Nf[20], epsOf = (z) => { const s = 1 + 0.05 * (z / 1000 - 0.5); const N = C.scale(N0, s); return C.mul(N, N); };
  const r1 = riccati(lam[20], epsOf, 1000, Nsub[20], 65, 0.5), r2 = riccati(lam[20], epsOf, 1000, Nsub[20], 65, 0.125);
  check(C.abs(C.sub(r1.rp, r2.rp)) < 1e-6 && C.abs(C.sub(r1.rs, r2.rs)) < 1e-6, `сходимость RK4 по шагу: |r(0.5 нм) − r(0.125 нм)| = ${Math.max(C.abs(C.sub(r1.rp, r2.rp)), C.abs(C.sub(r1.rs, r2.rs))).toExponential(1)}`);
  // rhoGraded ядра (h = 0.5) против независимого Риккати теста (h = 0.25), с EMA-слоем
  {
    const lam2 = [420, 550, 700], Ns2 = E.substrateN("bk7", lam2), Nf2 = E.filmN(lam2, m);
    const [RPg, RSg] = E.rhoGraded(lam2, Ns2, 65, 1000, m, 0.07, 0, 0.5);
    let dm = 0;
    for (let i = 0; i < lam2.length; i++) {
      const epsOf = (z) => { const sc = 1 + 0.07 * (z / 1000 - 0.5); const N = C.scale(Nf2[i], sc); return C.mul(N, N); };
      const r = riccati(lam2[i], epsOf, 1000, Ns2[i], 65, 0.25);
      dm = Math.max(dm, C.abs(C.add(RPg[i], r.rp)), C.abs(C.sub(RSg[i], r.rs)));
    }
    check(dm < 1e-4, `rhoGraded ядра совпадает с независимой реализацией Риккати: max|Δr| = ${dm.toExponential(1)}`);
  }
  // однородная плёнка: Риккати должен совпасть с Эйри точно
  const eps0 = C.mul(N0, N0); const r0 = riccati(lam[20], () => eps0, 1000, Nsub[20], 65, 0.25);
  const [RP, RS] = E.rhoStack([lam[20]], [{ N: [N0], d: 1000 }], [Nsub[20]], 65);
  check(C.abs(C.add(r0.rp, RP[0])) < 1e-8 && C.abs(C.sub(r0.rs, RS[0])) < 1e-8, `однородная плёнка: Риккати = Эйри с точностью ${Math.max(C.abs(C.add(r0.rp, RP[0])), C.abs(C.sub(r0.rs, RS[0]))).toExponential(1)}`);
}

// ================================================================ 5. Сохранение энергии и предельные случаи
console.log("\n4. Сохранение энергии, френелевские пределы, ветвь корня");
{
  const lam = 550, m = E.MATERIALS.ta2o5, N0 = E.filmN([lam], m)[0];
  const a = abeles(lam, [{ N: [N0[0], 0], d: 700 }], [1.52, 0], 60);
  check(Math.abs(a.Rs + a.Ts - 1) < 1e-12 && Math.abs(a.Rp + a.Tp - 1) < 1e-12, `непоглощающая плёнка на стекле, 60°: R_s + T_s − 1 = ${(a.Rs + a.Ts - 1).toExponential(1)}, R_p + T_p − 1 = ${(a.Rp + a.Tp - 1).toExponential(1)}`);
  const Nuv = E.filmN([340], E.MATERIALS.tio2)[0];
  const b = abeles(340, [{ N: Nuv, d: 300 }], E.substrateN("bk7", [340])[0], 60);
  check(Nuv[1] < 0 && b.Rs + b.Ts < 1 && b.Rp + b.Tp < 1, `TiO2 при 340 нм: N = ${Nuv[0].toFixed(3)} − ${(-Nuv[1]).toFixed(3)}i (k > 0 в конвенции n − ik); плёнка 300 нм на BK7, 60°: A_s = ${(1 - b.Rs - b.Ts).toFixed(4)}, A_p = ${(1 - b.Rp - b.Tp).toFixed(4)} — положительны (знак мнимой части верен)`);
  // EMA Бруггемана: корень удовлетворяет уравнению и физичен
  const eps = C.mul(N0, N0), b0 = C.scale(C.add(C.scale(eps, 0.5), [0.5, 0]), 1); // (3f−1)ε_a + (3f−1)ε_b при f = 0.5
  const Nema = E.buildLayers([lam], 100, m, 0, 5, 1)[0].N[0], eEma = C.mul(Nema, Nema);
  const resid = C.add(C.scale(C.div(C.sub(eps, eEma), C.add(eps, C.scale(eEma, 2))), 0.5), C.scale(C.div(C.sub(C.ONE, eEma), C.add(C.ONE, C.scale(eEma, 2))), 0.5));
  check(C.abs(resid) < 1e-12 && eEma[0] > 1 && eEma[0] < eps[0], `слой Бруггемана 50/50: невязка уравнения ЭС ${C.abs(resid).toExponential(1)}, ε_ema = ${eEma[0].toFixed(4)} между 1 и ε_плёнки = ${eps[0].toFixed(4)} (n_ema = ${Nema[0].toFixed(4)})`);
  // голая подложка: Френель; угол Брюстера BK7 при 550 нм
  const nb = E.substrateN("bk7", [550])[0]; const brew = Math.atan(nb[0]) * 180 / PI;
  const [RPb, RSb] = E.rhoStack([550], [], [nb], brew);
  check(C.abs(RPb[0]) < 1e-6, `угол Брюстера BK7 (${brew.toFixed(2)}°): |r_p| = ${C.abs(RPb[0]).toExponential(1)}`);
  const [RP1, RS1] = E.rhoStack([550], [], [nb], 45), [RP2, RS2] = E.rhoStack([550], [], [nb], 70);
  const d45 = psiDelta(RP1[0], RS1[0])[1], d70 = psiDelta(RP2[0], RS2[0])[1];
  check(Math.abs(d45 - 180) < 1e-9 && Math.abs(d70) < 1e-9, `Δ голого диэлектрика: ${d45.toFixed(1)}° ниже Брюстера, ${d70.toFixed(1)}° выше (эллипсометрическая конвенция)`);
  const [RPn, RSn] = E.rhoStack([550], [], [nb], 1e-9);
  check(C.abs(C.add(RPn[0], RSn[0])) < 1e-9, `нормальное падение: r_p = −r_s (ρ = −1, Δ = 180°), как и должно быть при определении r_p через H`);
  // ветвь корня: затухание в поглощающей среде и в плёнке
  const Nsi = E.substrateN("si", [400])[0]; const c = cosGamma(Nsi, Math.sin(65 * PI / 180) ** 2); const qz = C.mul(Nsi, c);
  check(qz[1] < 0, `Im(N cos γ) = ${qz[1].toFixed(4)} < 0 для c-Si при 400 нм: прошедшая волна затухает (требование книги для exp(+iωt))`);
}

console.log(failures ? `\nОШИБОК: ${failures}` : "\nВсе проверки по книге пройдены");
process.exit(failures ? 1 : 0);
