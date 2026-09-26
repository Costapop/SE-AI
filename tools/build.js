#!/usr/bin/env node
/* Сборка одностраничного приложения: src/app_template.html + src/ellip_core.js + examples/examples.json → index.html */
const fs = require("fs"), path = require("path");
const root = path.join(__dirname, "..");
const tpl = fs.readFileSync(path.join(root, "src/app_template.html"), "utf8");
const core = fs.readFileSync(path.join(root, "src/ellip_core.js"), "utf8");
const examples = fs.readFileSync(path.join(root, "examples/examples.json"), "utf8").trim();
if (!tpl.includes("/*__CORE__*/") || !tpl.includes("/*__EXAMPLES__*/[]")) throw new Error("в шаблоне нет плейсхолдеров");
const html = tpl.replace("/*__CORE__*/\n", core + "\n").replace("/*__EXAMPLES__*/[]", examples);
fs.writeFileSync(path.join(root, "index.html"), html);
console.log(`index.html: ${Buffer.byteLength(html)} байт`);
