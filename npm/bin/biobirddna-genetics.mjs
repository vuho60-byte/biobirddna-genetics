#!/usr/bin/env node
import { readFile, writeFile } from "node:fs/promises";
import { resolve } from "node:path";
import { ENGINE_VERSION, BioBirdDNAError, dispatch } from "../src/index.mjs";

function usage() { process.stderr.write("Usage: biobirddna-genetics <cross|infer|pigment|rank-mates|normalize-records|simulate-recovery> --input <file> [--output <file>] [--pretty]\n"); }
// Human-readable Vietnamese messages for the most common machine error slugs.
// Keys match the slug before the first ":". The machine `error` field stays
// unchanged; `message_vi` is additive. Mirrored verbatim from
// src/biobirddna_genetics/cli.py — keep both maps in sync.
const MESSAGES_VI = {
  input_file_missing: "Không tìm thấy file input — kiểm tra lại đường dẫn.",
  input_json_invalid: "File input không phải JSON hợp lệ — kiểm tra lại cú pháp.",
  input_json_non_finite: "Input chứa số không hữu hạn (NaN/Infinity) — chỉ dùng số thường.",
  input_root_must_be_object: "Nội dung JSON phải là một object ở cấp cao nhất.",
  operation_unknown: "Thao tác không được hỗ trợ — dùng: cross, infer, pigment, rank-mates, normalize-records, simulate-recovery.",
  male_and_female_genotype_required: "Thiếu kiểu gene bố hoặc mẹ — cần cả male_genotype và female_genotype.",
  offspring_count_invalid: "Số lượng chim con không hợp lệ — dùng số nguyên không âm.",
  offspring_sex_invalid: "Giới tính chim con không hợp lệ — dùng male, female hoặc unknown.",
  internal_error: "Lỗi nội bộ của công cụ, không phải do dữ liệu bạn nhập — báo lại cho người phát triển kèm traceback.",
};
const MESSAGE_VI_FALLBACK = "Dữ liệu đầu vào không hợp lệ — xem mã lỗi ở trường error và đối chiếu references/data-contracts.md.";
const messageVi = (code) => Object.prototype.hasOwnProperty.call(MESSAGES_VI, code.split(":", 1)[0]) ? MESSAGES_VI[code.split(":", 1)[0]] : MESSAGE_VI_FALLBACK;
const [operation, ...args] = process.argv.slice(2);
const valueAfter = (flag) => { const index = args.indexOf(flag); return index >= 0 ? args[index + 1] : undefined; };
const EXIT_INPUT_ERROR = 2, EXIT_INTERNAL_ERROR = 70;
try {
  const input = valueAfter("--input"), output = valueAfter("--output"); if (!operation || !input) { usage(); process.exitCode = EXIT_INPUT_ERROR; } else {
    const inputPath = resolve(input);
    let raw; try { raw = await readFile(inputPath, "utf8"); } catch (error) { if (error?.code === "ENOENT") throw new BioBirdDNAError(`input_file_missing:${inputPath}`); throw error; }
    let request; try { request = JSON.parse(raw); } catch { throw new BioBirdDNAError(`input_json_invalid:${inputPath}`); }
    const envelope = { ok: true, operation, engine_version: ENGINE_VERSION, result: dispatch(operation, request) }; const json = `${JSON.stringify(envelope, null, args.includes("--pretty") || output ? 2 : undefined)}\n`;
    if (output) await writeFile(resolve(output), json, "utf8"); else process.stdout.write(json);
  }
} catch (error) {
  // Same exit-code contract as the Python CLI (see references/data-contracts.md,
  // "CLI error envelope & exit codes"): a BioBirdDNAError slug means bad input
  // (exit 2); anything else is an engine bug, reported as internal_error:<Type>
  // with the stack preserved and exit 70, so wrappers never present engine
  // bugs as user input errors.
  const isInputError = error instanceof BioBirdDNAError;
  const code = isInputError ? error.message : `internal_error:${error?.name || error?.constructor?.name || "Error"}`;
  if (!isInputError) process.stderr.write(`${error?.stack ?? String(error)}\n`);
  process.stderr.write(`${JSON.stringify({ ok: false, operation: operation ?? null, error: code, message_vi: messageVi(code) })}\n`);
  process.exitCode = isInputError ? EXIT_INPUT_ERROR : EXIT_INTERNAL_ERROR;
}
