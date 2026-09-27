/* Дымовой тест интерфейса в jsdom: встроенный образец, симуляция с большим градиентом, симуляция с поглощением,
   образец с задней стороной подложки, очистка панелей. Запуск: npm install && node test/ui_smoke.js   (≈12–15 мин на одном ядре; jsdom втрое медленнее браузера) */
const { JSDOM } = require("jsdom"); const fs = require("fs"), path = require("path");
let html = fs.readFileSync(path.join(__dirname, "../index.html"), "utf8");
html = html.replace(/<script src="https:\/\/cdnjs[^"]*chart[^"]*"><\/script>/, '<script>window.Chart = class { constructor(){} destroy(){} resize(){} };</script>')
           .replace(/<script src="https:\/\/cdnjs[^"]*xlsx[^"]*"><\/script>/, '<script>window.XLSX = {};</script>');
const dom = new JSDOM(html, { runScripts: "dangerously", pretendToBeVisual: true }); const w = dom.window; w.HTMLCanvasElement.prototype.getContext = () => ({});
const $ = (id) => w.document.getElementById(id); const txt = (el) => el.textContent.replace(/\s+/g, " ").trim();
let failures = 0; const check = (ok, msg) => { console.log((ok ? "  ok   " : "  FAIL ") + msg); if (!ok) failures++; };
w.addEventListener("error", (e) => { console.log("window error:", e.message); failures++; });
async function runAll(nTrain) {
  $("nTrain").value = String(nTrain); $("nTrain").dispatchEvent(new w.Event("input"));
  const t = Date.now(); $("btnAll").click();
  while ($("globalStatus").textContent !== "Готово." && Date.now() - t < 600000) await new Promise(r => setTimeout(r, 200));
  const t2 = Date.now(); while (/Анализ причин отклонения…/.test($("causeOut").textContent) && Date.now() - t2 < 300000) await new Promise(r => setTimeout(r, 200));
  return (Date.now() - t) / 1000;
}
const truthRows = () => Object.fromEntries([...$("truthOut").querySelectorAll("tr")].slice(1).map(r => [...r.children].map(c => c.textContent.trim())).map(r => [r[0], r]));
(async () => {
  await new Promise(r => setTimeout(r, 50));
  console.log("1. Встроенный образец:", $("exSel").selectedOptions[0].textContent);
  let sec = await runAll(200); const tr = truthRows();
  check(/совпадает/.test(tr["Набор дефектов"][2]) && !/не совпадает/.test(tr["Набор дефектов"][2]), `набор дефектов совпал (${sec.toFixed(0)} с)`);
  check(/Анализ причин/.test(txt($("causeOut"))), "разбор причин отображён");
  console.log("2. Симуляция: TiO2/BK7, δ = +8 %, слой 1 нм");
  $("simD").value = "1017.3"; $("simDelta").value = "8"; $("simRough").value = "1.0"; $("simDphi").value = "0.02"; $("simSeed").value = "9"; $("deltaMax").value = "10";
  $("btnSim").click(); { const t0 = Date.now(); while (!/Сгенерировано|неверно|Нужно|ограничен/.test($("simStatus").textContent) && Date.now() - t0 < 60000) await new Promise(r => setTimeout(r, 100)); }
  check(/Сгенерировано/.test($("simStatus").textContent), "спектр сгенерирован (непрерывный профиль, Риккати)");
  check(/Сначала диагноз/.test(txt($("finOut"))) && $("causeOut").innerHTML === "", "панели шагов 2–4 очищены после генерации");
  sec = await runAll(200); const tr2 = truthRows();
  check(!/не совпадает/.test(tr2["Набор дефектов"][2]), `набор дефектов совпал (${sec.toFixed(0)} с)`);
  const dFound = parseFloat((tr2["δ, %"][2].match(/найдено (-?[\d.]+)/) || [])[1]);
  check(Math.abs(dFound - 8) < 0.3, `градиент найден: ${dFound} % (истина 8 %)`);
  console.log("3. Симуляция: HfO2/Si 65°, поглощение k₄₀₀ = 1.5·10⁻³ + слой 1 нм, без градиента");
  $("matSel").value = "hfo2"; $("matSel").dispatchEvent(new w.Event("change")); $("subSel").value = "si"; $("phi").value = "65"; $("dNom").value = "950"; $("dNom").dispatchEvent(new w.Event("change"));
  $("simD").value = "933"; $("simDelta").value = "0"; $("simRough").value = "1.0"; $("simK").value = "1.5"; $("simK").dispatchEvent(new w.Event("input")); $("simSeed").value = "11";
  check(/A_k/.test($("simKInfo").textContent), "подсказка к полю поглощения показывает k хвоста и A_k: " + $("simKInfo").textContent);
  $("btnSim").click(); { const t0 = Date.now(); while (!/Сгенерировано|неверно|Нужно|ограничен/.test($("simStatus").textContent) && Date.now() - t0 < 60000) await new Promise(r => setTimeout(r, 100)); }
  check(/Хвост поглощения/.test($("simStatus").textContent), "спектр с поглощением сгенерирован: " + $("simStatus").textContent.replace(/^.*Хвост/, "Хвост"));
  sec = await runAll(500); const tr3 = truthRows();   // с 12 классами при 300 примерах классификатор ещё добавляет ложный градиент
  check(!/не совпадает/.test(tr3["Набор дефектов"][2]), `набор дефектов совпал: «${tr3["Набор дефектов"][1]}» (${sec.toFixed(0)} с)`);
  const kFound = parseFloat((tr3["Поглощение: k хвоста при 400 нм"][2].match(/найдено ([\d.]+)/) || [])[1]);
  check(Math.abs(kFound - 1.5) < 0.4, `поглощение найдено: k₄₀₀ = ${kFound}·10⁻³ (истина 1.5·10⁻³)`);
  check(/k₄₀₀/.test(txt($("causeOut"))) && /P\(поглощение\)/.test(txt($("diagOut"))), "поглощение показано в диагнозе и в анализе причин");
  console.log("4. Встроенный образец 09: задняя сторона подложки (задано 0.95 ± 0.02, истина 0.96)");
  $("exSel").value = "8"; $("exSel").dispatchEvent(new w.Event("change")); await new Promise(r => setTimeout(r, 100));
  check($("ckBack").checked && $("backF").value === "0.95" && !$("backF").disabled && /задняя сторона подложки: доля 0.95/.test($("dataInfo").textContent), "флажок задней стороны и доля подставлены из образца: " + $("dataInfo").textContent);
  sec = await runAll(400); const tr4 = truthRows();
  check(!/не совпадает/.test(tr4["Набор дефектов"][2]), `набор дефектов совпал: «${tr4["Набор дефектов"][1]}» (${sec.toFixed(0)} с)`);
  check(/0\.960/.test(tr4["Задняя сторона подложки, доля"][1]) && /в модели 0\.950/.test(tr4["Задняя сторона подложки, доля"][2]), "таблица истины показывает истинную и заданную долю: " + tr4["Задняя сторона подложки, доля"].join(" | "));
  check(/доля задней стороны как у образца/.test(txt($("causeOut"))), "в бюджете отклонения есть строка про долю задней стороны");
  const dErr = parseFloat((tr4["d, нм"][2].match(/ошибка (-?[\d.]+)/) || [])[1]);
  check(Math.abs(dErr) < 1.5, `толщина найдена с ошибкой ${dErr} нм`);
  console.log("5. Снятие флажка задней стороны сбрасывает шаги, установка обратно — тоже");
  $("ckBack").click(); await new Promise(r => setTimeout(r, 50));
  check(/Нажмите «Опорный фит»/.test(txt($("refOut"))) && $("backF").disabled && /задней стороны нет/.test($("dataInfo").textContent), "панели сброшены, поля задней стороны отключены");
  $("ckBack").click(); await new Promise(r => setTimeout(r, 50));
  check(!$("backF").disabled && /доля 0.95/.test($("dataInfo").textContent), "поля снова активны");
  console.log("6. Повторный шаг 1 очищает шаги 2–4");
  $("btn1").click(); await new Promise(r => setTimeout(r, 100));
  check(/Сначала опорный фит/.test(txt($("trainOut"))) && /Сначала обучение/.test(txt($("diagOut"))) && /Сначала диагноз/.test(txt($("finOut"))) && $("causeOut").innerHTML === "" && $("btn3").disabled && $("btn4").disabled, "панели очищены, кнопки 3–4 отключены");
  console.log(failures ? `\nОШИБОК: ${failures}` : "\nВсе проверки пройдены");
  process.exit(failures ? 1 : 0);
})();
