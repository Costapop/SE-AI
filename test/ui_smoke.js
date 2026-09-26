/* Дымовой тест интерфейса в jsdom: встроенный образец, симуляция с большим градиентом, очистка панелей.
   Запуск: npm install && node test/ui_smoke.js   (≈3–4 мин на одном ядре) */
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
  $("btnSim").click(); check(/Сгенерировано/.test($("simStatus").textContent), "спектр сгенерирован");
  check(/Сначала диагноз/.test(txt($("finOut"))) && $("causeOut").innerHTML === "", "панели шагов 2–4 очищены после генерации");
  sec = await runAll(200); const tr2 = truthRows();
  check(!/не совпадает/.test(tr2["Набор дефектов"][2]), `набор дефектов совпал (${sec.toFixed(0)} с)`);
  const dFound = parseFloat((tr2["δ, %"][2].match(/найдено (-?[\d.]+)/) || [])[1]);
  check(Math.abs(dFound - 8) < 0.3, `градиент найден: ${dFound} % (истина 8 %)`);
  console.log("3. Повторный шаг 1 очищает шаги 2–4");
  $("btn1").click(); await new Promise(r => setTimeout(r, 100));
  check(/Сначала опорный фит/.test(txt($("trainOut"))) && /Сначала обучение/.test(txt($("diagOut"))) && /Сначала диагноз/.test(txt($("finOut"))) && $("causeOut").innerHTML === "" && $("btn3").disabled && $("btn4").disabled, "панели очищены, кнопки 3–4 отключены");
  console.log(failures ? `\nОШИБОК: ${failures}` : "\nВсе проверки пройдены");
  process.exit(failures ? 1 : 0);
})();
