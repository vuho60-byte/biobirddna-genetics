#!/usr/bin/env node
import { readFile, stat, writeFile } from "node:fs/promises";
import { resolve } from "node:path";
import { createFeedbackArtifact } from "../src/index.mjs";

const MAX_INPUT_BYTES = 1024 * 1024;
const usage = () => process.stderr.write("Usage: biobirddna-feedback --input <file> --expected-integrity <sha256> [--output <file>]\n");
function options(args) {
  const known = new Set(["--input", "--expected-integrity", "--output"]), result = {};
  for (let index = 0; index < args.length; index += 2) {
    const flag = args[index], value = args[index + 1];
    if (!known.has(flag) || !value || result[flag] !== undefined) throw new Error("feedback_cli_arguments_invalid");
    result[flag] = value;
  }
  return result;
}
try {
  const args = options(process.argv.slice(2));
  if (!args["--input"] || !args["--expected-integrity"]) { usage(); process.exitCode = 2; }
  else {
    const inputPath = resolve(args["--input"]), metadata = await stat(inputPath);
    if (!metadata.isFile()) throw new Error("feedback_input_not_regular_file");
    if (metadata.size > MAX_INPUT_BYTES) throw new Error("feedback_input_too_large");
    const bytes = await readFile(inputPath);
    let envelope;
    try { envelope = JSON.parse(bytes.toString("utf8")); } catch { throw new Error("feedback_input_invalid_json"); }
    if (!envelope || typeof envelope !== "object" || Array.isArray(envelope)) throw new Error("feedback_input_root_must_be_object");
    const artifact = createFeedbackArtifact(envelope, args["--expected-integrity"]), json = `${JSON.stringify(artifact, null, args["--output"] ? 2 : undefined)}\n`;
    if (args["--output"]) await writeFile(resolve(args["--output"]), json, "utf8"); else process.stdout.write(json);
  }
} catch (error) {
  process.stderr.write(`${JSON.stringify({ ok: false, error: error instanceof Error ? error.message : "feedback_runtime_error" })}\n`);
  process.exitCode = 2;
}
