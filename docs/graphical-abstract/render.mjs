// Render graphical-abstract.html to a 2400 x 1200 PNG with the pre-installed Chromium.
// Usage: NODE_PATH=$(npm root -g) node render.mjs   (needs the `playwright` package)
import { createRequire } from "node:module";
import { fileURLToPath, pathToFileURL } from "node:url";
import path from "node:path";

const here = path.dirname(fileURLToPath(import.meta.url));
// ESM ignores NODE_PATH, so resolve playwright through require (local install or NODE_PATH).
const { chromium } = createRequire(import.meta.url)("playwright");
const browser = await chromium.launch();
const page = await browser.newPage({ viewport: { width: 1600, height: 800 }, deviceScaleFactor: 1.5 });
await page.goto(pathToFileURL(path.join(here, "graphical-abstract.html")).href, { waitUntil: "networkidle" });
await page.evaluate(() => document.fonts.ready);
await page.screenshot({ path: path.join(here, "graphical-abstract.png"), clip: { x: 0, y: 0, width: 1600, height: 800 } });
await browser.close();
console.log("wrote graphical-abstract.png");
