// Run: node calculations/code/test_latex_cell_fix.cjs
const fs = require("node:fs");
const vm = require("node:vm");
const assert = require("node:assert/strict");
const path = require("node:path");
const root = path.join(process.env.USERPROFILE, ".cache/opencode/npm/opencode-latex@0.2.0/1793138687352/node_modules");
const ts = require(path.join(root, "typescript"));
// Test the installed pure formatter without loading Bun-only native rendering.
const bundle = fs.readFileSync(path.join(root, "opentui-math/dist/chunk-nf3q6tfq.js"), "utf8");
assert(bundle.includes("// src/renderable.ts"));
const math = {};
vm.runInNewContext(bundle.split("// src/renderable.ts")[0].replace("var __require = import.meta.require;", ""), math);
const fraction = math.renderLatexToString("\\frac{a+b}{2}", { displayMode: true });
assert(fraction.includes("\n") && fraction.includes("2") && !fraction.includes("\\frac"));
const source = fs.readFileSync(path.join(root, "opencode-latex/src/tui.ts"), "utf8").replace(/^import .*$/gm, "");
const out = ts.transpileModule(source, { compilerOptions: { module: ts.ModuleKind.CommonJS }, reportDiagnostics: true });
assert.equal(out.diagnostics.length, 0);
let callback, commands, alert, destroyed = false;
class TextRenderable {
  constructor(ctx, options) { this.options = options; }
  destroy() { destroyed = true; }
}
const sandbox = {
  exports: {}, Error, TextRenderable, renderLatexToString: math.renderLatexToString,
  Plugin: { define: x => x },
  GraphicalLatexRenderable: class { constructor() { throw Error("Legacy renderer must not run in cells mode"); } },
};
vm.runInNewContext(out.outputText, sandbox);
sandbox.exports.default.setup({
  options: {},
  storage: { store: () => [{ graphicsMode: "cells" }, () => {}] },
  keymap: { layer: f => { commands = f().commands; } },
  markdown: { registerCodeBlockRenderer: (language, f) => { callback = f; } },
  ui: { dialog: { alert: async x => { alert = x; } }, toast: { show: () => {} } },
});
const fallback = { defaultRender: () => "fallback" };
const rendered = callback({ text: "\\frac{a+b}{2}" }, fallback);
assert.equal(rendered.options.content, fraction);
assert.equal(rendered.options.wrapMode, "none");
assert.equal(callback({ text: "\\frac{" }, fallback), "fallback");
commands.find(x => x.slash.name === "latex-diagnostic").run().then(() => {
  assert(destroyed);
  assert(alert.message.includes("constructor succeeded"));
  console.log("PASS: real formula formatting, cell dispatch, invalid-input fallback, diagnostic cleanup; host TextRenderable mocked.");
  console.log(fraction);
}).catch(error => { console.error(error); process.exitCode = 1; });
