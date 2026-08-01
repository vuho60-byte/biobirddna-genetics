#!/usr/bin/env node
import { readFile, writeFile } from "node:fs/promises";
import { resolve } from "node:path";
import { ENGINE_VERSION, BioBirdDNAError, dispatch } from "../src/index.mjs";

function usage() { process.stderr.write("Usage: biobirddna-genetics <cross|infer|pigment|rank-mates|normalize-records|simulate-recovery> --input <file> [--output <file>] [--pretty]\n"); }
const [operation, ...args] = process.argv.slice(2);
const valueAfter = (flag) => { const index = args.indexOf(flag); return index >= 0 ? args[index + 1] : undefined; };
try {
  const input = valueAfter("--input"), output = valueAfter("--output"); if (!operation || !input) { usage(); process.exitCode = 2; } else {
    const request = JSON.parse(await readFile(resolve(input), "utf8")); const envelope = { ok: true, operation, engine_version: ENGINE_VERSION, result: dispatch(operation, request) }; const json = `${JSON.stringify(envelope, null, args.includes("--pretty") || output ? 2 : undefined)}\n`;
    if (output) await writeFile(resolve(output), json, "utf8"); else process.stdout.write(json);
  }
} catch (error) { const code = error instanceof BioBirdDNAError ? error.message : "input_or_runtime_error"; process.stderr.write(`${JSON.stringify({ ok: false, operation, error: code })}\n`); process.exitCode = 2; }
