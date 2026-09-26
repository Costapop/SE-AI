const E = require("../src/ellip_core.js"); const fs = require("fs");
const ex = JSON.parse(fs.readFileSync(require("path").join(__dirname, "../examples/examples.json")));
const NTR = +(process.argv[2] || 500);
for (const e of ex) {
  const lam = e.lam, data = { lam, psi: e.psi, del: e.del, sigPsi: e.sigPsi, sigDel: e.sigDel };
  const cfg = { mat: E.MATERIALS[e.mat], subKey: e.sub, Nsub: E.substrateN(e.sub, lam), phi: e.phi };
  const prior = { dMin: e.dNom * (1 - e.dRange / 100), dMax: e.dNom * (1 + e.dRange / 100) };
  let t = Date.now();
  const rf = E.referenceFit(data, cfg, prior);
  const tRef = Date.now() - t;
  const subDel = E.modelPsiDelta(lam, cfg.Nsub, cfg.phi, 0, cfg.mat, 0, 0, 1).del;
  const { HW, QW } = E.anchors(lam, data.del, subDel); const ch = E.channels(lam, rf.resid, HW, QW);
  const rng = E.makeRng(77); const X = [], y = [], sz = [];
  t = Date.now();
  for (let i = 0; i < NTR; i++) { const s = E.trainingExample(data, cfg, prior, { deltaMax: 0.05, roughMax: 3, bw: e.bw, phiErr: 0.03 }, rng, i % 6); X.push(s.x); y.push(s.y); sz.push([s.delta, s.dRough]); }
  const tGen = Date.now() - t;
  const idx = X.map((_, i) => i); for (let i = idx.length - 1; i > 0; i--) { const j = rng.int(i + 1); [idx[i], idx[j]] = [idx[j], idx[i]]; }
  const nte = Math.floor(NTR / 4), te = idx.slice(0, nte), tr = idx.slice(nte);
  const model = E.trainSoftmax(tr.map(i => X[i]), tr.map(i => y[i]), 6, { epochs: 250 });
  let acc = 0; for (const i of te) { const p = model.predictProba(X[i]); if (p.indexOf(Math.max(...p)) === y[i]) acc++; }
  const full = E.trainSoftmax(X, y, 6, { epochs: 250 });
  const rd = E.ridge(X, sz.map(v => 100 * v[0]), 30), rr = E.ridge(X, sz.map(v => v[1]), 30);
  const xr = E.featureVector(data, cfg, rf.resid); const probs = full.predictProba(xr); const k = probs.indexOf(Math.max(...probs)); const cls = E.CLASSES[k];
  const p0 = Object.assign({}, rf.P, { dRough: Math.max(0.1, rr.predict(xr)), delta: rd.predict(xr) / 100 });
  t = Date.now();
  let ff = E.finalFit(data, cfg, prior, { grad: cls.g !== 0, rough: cls.r === 1 }, p0);
  let flags = { grad: cls.g !== 0, rough: cls.r === 1 }, veto = "";
  for (let pass = 0; pass < 2; pass++) {
    let changed = false;
    if (flags.rough && ff.P.dRough < 3 * (ff.E.dRough || 0.01)) { flags.rough = false; changed = true; veto += " EMA незначим →"; }
    if (flags.grad && Math.abs(ff.P.delta) < 3 * (ff.E.delta || 0.0001)) { flags.grad = false; changed = true; veto += " δ незначим →"; }
    if (!changed) break;
    ff = E.finalFit(data, cfg, prior, flags, p0);
  }
  const tFin = Date.now() - t;
  const T = e.truth;
  const trueCls = E.CLASSES.findIndex(c => c.g === Math.sign(T.delta) && c.r === (T.dRough > 0 ? 1 : 0));
  console.log(`\n${e.title}: опорный фит ${tRef} мс, χ²/ν = ${rf.chi2.toFixed(1)}, d = ${rf.P.d.toFixed(1)}; каналы Ψ_HW ${ch.psiHW.toFixed(1)}σ, |Δ|_HW ${ch.delHW.toFixed(1)}σ, |Δ|_QW ${ch.delQW.toFixed(1)}σ (якорей ${ch.nHW})`);
  console.log(`   обучение: ${NTR} примеров за ${(tGen / 1000).toFixed(0)} с, точность на отложенных ${(100 * acc / nte).toFixed(1)} %`);
  console.log(`   вероятности: ${probs.map((p, i) => E.CLASSES[i].name + " " + (100 * p).toFixed(1) + "%").join("; ")}`);
  console.log(`   выбрано: «${cls.name}» ${k === trueCls ? "✓ верно" : "✗ (истина: " + E.CLASSES[trueCls].name + ")"}; ridge: δ ${rd.predict(xr).toFixed(2)} %, EMA ${rr.predict(xr).toFixed(2)} нм`);
  console.log(`   финальный фит ${tFin} мс${veto ? " (вето:" + veto + " упрощено)" : ""}: χ²/ν = ${ff.chi2.toFixed(2)}; d = ${ff.P.d.toFixed(2)} ± ${(ff.E.d || 0).toFixed(2)} (истина ${T.d}); δ = ${(100 * ff.P.delta).toFixed(3)} ± ${(100 * (ff.E.delta || 0)).toFixed(3)} % (истина ${(100 * T.delta).toFixed(2)}); EMA = ${ff.P.dRough.toFixed(3)} ± ${(ff.E.dRough || 0).toFixed(3)} (истина ${T.dRough})`);
  const nk = E.filmNK([550], Object.assign({}, cfg.mat, { A: ff.P.A, Auv: ff.P.Auv, Eg: ff.P.Eg }));
  console.log(`   n(550) = ${nk.n[0].toFixed(4)} (истина ${T.n["550"].toFixed(4)})`);
}
