import assert from "node:assert/strict";
import { execFile } from "node:child_process";
import { access, mkdir, mkdtemp, readFile, rm, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

const execFileAsync = promisify(execFile);
const packageRoot = dirname(dirname(fileURLToPath(import.meta.url)));
const npmCli = join(dirname(process.execPath), "node_modules", "npm", "bin", "npm-cli.js");
const execNpm = (args, options) => execFileAsync(process.execPath, [npmCli, ...args], options);
const temp = await mkdtemp(join(tmpdir(), "biobirddna-pack-smoke-"));
try {
  const { stdout } = await execNpm(["pack", "--json", "--pack-destination", temp], { cwd: packageRoot });
  const [{ filename, files }] = JSON.parse(stdout);
  assert.ok(files.every((entry) => !entry.path.startsWith("test/")), "test files must not ship");
  assert.ok(files.some((entry) => entry.path === "bin/biobirddna-genetics.mjs"), "bin must ship");
  assert.ok(files.some((entry) => entry.path === "bin/biobirddna-feedback.mjs"), "feedback bin must ship");
  const consumer = join(temp, "consumer with spaces");
  await mkdir(consumer);
  await writeFile(join(consumer, "package.json"), JSON.stringify({ name: "consumer", private: true, type: "module" }), "utf8");
  await execNpm(["install", "--ignore-scripts", "--no-audit", "--no-fund", join(temp, filename)], { cwd: consumer });
  const input = join(consumer, "input.json"), output = join(consumer, "output.json");
  await writeFile(input, JSON.stringify({ model: { model_id: "M1_AR", mode: "autosomal_recessive" }, male_genotype: "N/m", female_genotype: "m/m" }));
  const installedBin = join(consumer, "node_modules", "@biobirddna", "genetics", "bin", "biobirddna-genetics.mjs");
  if (process.platform === "win32") {
    await access(join(consumer, "node_modules", ".bin", "biobirddna-genetics.cmd"));
    await execFileAsync(process.execPath, [installedBin, "cross", "--input", input, "--output", output], { cwd: consumer });
  } else {
    await execFileAsync(installedBin, ["cross", "--input", input, "--output", output], { cwd: consumer });
  }
  assert.equal(JSON.parse(await readFile(output, "utf8")).ok, true);
  const feedbackInput = join(consumer, "feedback input.json"), feedbackOutput = join(consumer, "feedback output.json"), feedbackBin = join(consumer, "node_modules", "@biobirddna", "genetics", "bin", "biobirddna-feedback.mjs"), integrity = "c".repeat(64);
  await writeFile(feedbackInput, JSON.stringify({ schema_version: "1.0.0", run_id: "CM-RUN-PACK", operation: "cross", engine_version: "0.2.0", input_sha256: "a".repeat(64), output_sha256: "b".repeat(64), package: { name: "@biobirddna/genetics", version: "0.2.0", integrity_sha256: integrity }, result: { evidence: { evidence_grade: "G0" }, limits: { decision_support_only: true, gene_status: "not_assessed" } }, human_review_status: "pending" }));
  await execFileAsync(process.execPath, [feedbackBin, "--input", feedbackInput, "--expected-integrity", integrity, "--output", feedbackOutput], { cwd: consumer });
  assert.equal(JSON.parse(await readFile(feedbackOutput, "utf8")).schema_version, "biobirddna_feedback_v1");
  process.stdout.write("pack-install smoke passed\n");
} finally {
  await rm(temp, { recursive: true, force: true });
}
