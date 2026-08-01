#!/usr/bin/env node
import assert from "node:assert/strict";
import { execFile } from "node:child_process";
import { access, cp, mkdtemp, mkdir, readFile, rm, stat, writeFile } from "node:fs/promises";
import { tmpdir } from "node:os";
import { basename, dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { promisify } from "node:util";

const execFileAsync = promisify(execFile);
const skillRoot = dirname(dirname(fileURLToPath(import.meta.url)));
const evidenceArgument = process.argv.indexOf("--evidence");
const evidencePath = evidenceArgument === -1 ? undefined : process.argv[evidenceArgument + 1];
const installerArgument = process.argv.indexOf("--installer");
const requestedInstaller = installerArgument === -1 ? undefined : process.argv[installerArgument + 1];
const tarballArgument = process.argv.indexOf("--tarball");
const retainedTarball = tarballArgument === -1 ? undefined : process.argv[tarballArgument + 1];
const PINNED_INTEGRITY = "sha512-KX/kxHHTEwbPiHxC86eTrm0Z/mgx+l1/JYgTSDC4UgEs0msLcUgRhpmNoiS1wd1oMHzxKD4kJXk+nSHqToJYmg==";

function requiredArgumentValue(index, name) {
  if (index !== -1 && !process.argv[index + 1]) {
    throw new Error(`missing_${name}_value`);
  }
}

function nativeInstallerRelativePath(platform = process.platform) {
  return platform === "win32" ? "scripts/install_skill.ps1" : "scripts/install_skill.sh";
}

async function resolveNpmCli(environment = process.env) {
  const candidates = [
    environment.npm_execpath,
    process.platform === "win32"
      ? join(dirname(process.execPath), "node_modules", "npm", "bin", "npm-cli.js")
      : join(dirname(dirname(process.execPath)), "lib", "node_modules", "npm", "bin", "npm-cli.js"),
  ].filter((candidate, index, values) => candidate && values.indexOf(candidate) === index);
  for (const candidate of candidates) {
    try {
      await access(candidate);
      return candidate;
    } catch (error) {
      if (error?.code !== "ENOENT" && error?.code !== "ENOTDIR") throw error;
    }
  }
  throw new Error("npm_cli_not_found");
}

async function run(executable, args, options = {}) {
  const completed = await execFileAsync(executable, args, {
    ...options,
    maxBuffer: 8 * 1024 * 1024,
    windowsHide: true,
  });
  return completed.stdout.trim();
}

async function runNpm(args, options) {
  return run(process.execPath, [await resolveNpmCli(), ...args], options);
}

async function exists(path) {
  try {
    await access(path);
    return true;
  } catch (error) {
    if (error?.code === "ENOENT") return false;
    throw error;
  }
}

function safeErrorCode(error) {
  return typeof error?.code === "string" && /^[A-Z0-9_]+$/i.test(error.code)
    ? error.code
    : "native_installer_smoke_failed";
}

function credentialFree(value) {
  return !/(?:-----BEGIN [A-Z ]*PRIVATE KEY-----|\bgh[pousr]_[A-Za-z0-9]{20,}\b|\beyJ[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}\.[A-Za-z0-9_-]{10,}\b|(?:token|secret|password)\s*[:=])/i.test(value);
}

async function writeEvidence(path, evidence) {
  if (!path) return;
  await mkdir(dirname(path), { recursive: true });
  await writeFile(path, `${JSON.stringify(evidence, null, 2)}\n`, "utf8");
}

async function main() {
  requiredArgumentValue(evidenceArgument, "evidence");
  requiredArgumentValue(installerArgument, "installer");
  requiredArgumentValue(tarballArgument, "tarball");
  const expectedInstaller = nativeInstallerRelativePath();
  if (requestedInstaller && requestedInstaller !== expectedInstaller) {
    throw new Error("installer_platform_mismatch");
  }
  const python = process.env.PYTHON ?? (process.platform === "win32" ? "python" : "python3");
  const workspace = await mkdtemp(join(tmpdir(), "biobirddna native Việt "));
  const cleanRoot = join(workspace, "clean-release");
  const childHome = join(workspace, "isolated-home");
  const childTemp = await mkdtemp(join(tmpdir(), "biobirddna-native-child-"));
  const childEnvironment = {
    PATH: process.env.PATH,
    SystemRoot: process.env.SystemRoot,
    SYSTEMROOT: process.env.SYSTEMROOT,
    ComSpec: process.env.ComSpec,
    PATHEXT: process.env.PATHEXT,
    HOME: childHome,
    USERPROFILE: childHome,
    TEMP: childTemp,
    TMP: childTemp,
    TMPDIR: childTemp,
    npm_config_cache: join(workspace, "npm-cache"),
    npm_config_userconfig: join(workspace, "npmrc"),
    PYTHONDONTWRITEBYTECODE: "1",
    PYTHONUTF8: "1",
  };
  const installer = join(cleanRoot, expectedInstaller);
  const evidence = {
    status: "failed",
    phase: "start",
    source: { clean_release_closure: false, vcs_state_copied: false, installed_git_absent: false, committed_posix_installer_mode: "100755" },
    environment: { isolated: true, inherited_user_state: false, temporary_home: true, temporary_npm_cache: true, temporary_npm_userconfig: true },
    installer: { platform: process.platform, path: expectedInstaller, executed: false },
    runtime: { runner: { os: process.platform, image: process.env.BIOBIRDDNA_RUNNER_IMAGE ?? "local-unverified" }, node: { version: process.version }, npm: { version: "unavailable" }, package: { name: "unavailable", version: "unavailable" } },
    test_counts: { passed: 0, failed: 0, skipped: 0 },
    gates: { zero_skipped_native_gates: false },
    skill_validation: { passed: false },
    npm: { install_clean: false, cli_probes: { genetics: false, feedback: false }, tarball: { retained: false } },
    paths: { space_and_non_ascii: false },
    executable_bits: { passed: process.platform === "win32", not_applicable: process.platform === "win32", committed_mode: "unverified" },
    cleanup: { passed: false },
    credential_free_log: false,
    log_bytes: 0,
  };

  try {
    evidence.phase = "clean-release-closure";
    await mkdir(childHome, { recursive: true });
    await mkdir(childEnvironment.npm_config_cache, { recursive: true });
    await cp(skillRoot, cleanRoot, { recursive: true, filter: (sourcePath) => ![".git", ".gitignore"].includes(basename(sourcePath)) });
    evidence.source.clean_release_closure = true;
    assert.equal(await exists(join(cleanRoot, ".git")), false, "clean release closure must exclude checkout VCS state");
    evidence.phase = "source-access";
    await access(installer);
    const gitEntry = await run("git", ["ls-tree", "HEAD", "scripts/install_skill.sh"], { cwd: skillRoot, env: childEnvironment });
    const committedMode = /^([0-7]{6})\s+blob\s+[0-9a-f]{40}\s+scripts\/install_skill\.sh$/m.exec(gitEntry)?.[1];
    assert.equal(committedMode, "100755", "POSIX installer must be executable in the committed release candidate");
    evidence.source.committed_posix_installer_mode = committedMode;
    evidence.executable_bits.committed_mode = committedMode;
    evidence.phase = "full-skill-validation";
    const validation = JSON.parse(await run(python, ["-X", "utf8", "-B", join(cleanRoot, "scripts", "validate_package.py")], { cwd: cleanRoot, env: childEnvironment }));
    assert.equal(validation.tests_passed, true, "Python release tests must pass");
    assert.equal(validation.examples_verified, true, "release examples must verify");
    assert.equal(validation.npm_tests_passed, true, "npm release tests must pass");
    assert.equal(validation.npm_pack_verified, true, "npm package smoke must pass");
    assert.equal(validation.package_verified, true, "release inventory must verify");
    const pythonPassed = Number(/Ran (\d+) tests/.exec(validation.test_stderr)?.[1] ?? 0);
    const npmPassed = Number(/ℹ pass (\d+)/.exec(validation.npm_test_stdout)?.[1] ?? 0);
    const npmSkipped = Number(/ℹ skipped (\d+)/.exec(validation.npm_test_stdout)?.[1] ?? 0);
    const pythonSkipped = (validation.test_stderr.match(/\.\.\. skipped(?:\s|$)/gm) ?? []).length;
    evidence.test_counts = { passed: pythonPassed + npmPassed, failed: 0, skipped: pythonSkipped + npmSkipped };
    assert.ok(evidence.test_counts.passed > 0, "release validation must report executed tests");
    assert.equal(evidence.test_counts.skipped, 0, "release validation must not skip native gates");
    evidence.gates.zero_skipped_native_gates = evidence.test_counts.skipped === 0;
    evidence.skill_validation.passed = true;

    evidence.phase = "native-installer";
    const codexHome = join(workspace, "CODEX HOME");
    const installRoot = join(codexHome, "skills");
    const installerEnvironment = { ...childEnvironment, CODEX_HOME: codexHome };
    if (process.platform === "win32") {
      await run(process.env.POWERSHELL ?? "powershell.exe", [
        "-NoLogo", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", installer,
        "-DestinationRoot", installRoot, "-PythonExecutable", python,
      ], { cwd: cleanRoot, env: installerEnvironment });
    } else {
      await run(installer, [installRoot], { cwd: cleanRoot, env: installerEnvironment });
    }
    evidence.installer.executed = true;
    evidence.phase = "installed-skill-access";
    const installedSkill = join(installRoot, "biobirddna-genetics");
    await access(installedSkill);
    evidence.source.installed_git_absent = !(await exists(join(installedSkill, ".git")));
    evidence.paths.space_and_non_ascii = workspace.includes("Việt") && await exists(installedSkill);

    if (process.platform !== "win32") {
      const sourceMode = (await stat(installer)).mode;
      const installedMode = (await stat(join(installedSkill, "scripts", "install_skill.sh"))).mode;
      evidence.executable_bits.passed = Boolean((sourceMode & 0o111) && (installedMode & 0o111));
      assert.equal(evidence.executable_bits.passed, true, "POSIX installer must remain executable after install");
    }

    evidence.phase = "npm-pack";
    const packageRoot = join(installedSkill, "npm");
    const packageMetadata = JSON.parse(await readFile(join(packageRoot, "package.json"), "utf8"));
    evidence.runtime.package = { name: packageMetadata.name, version: packageMetadata.version };
    assert.deepEqual(evidence.runtime.package, { name: "@biobirddna/genetics", version: "0.2.0" }, "installed package identity must match the reviewed candidate");
    evidence.runtime.npm.version = await runNpm(["--version"], { cwd: packageRoot, env: childEnvironment });
    const packed = JSON.parse(await runNpm(["pack", "--json", "--pack-destination", dirname(retainedTarball)], { cwd: packageRoot, env: childEnvironment }));
    const [{ filename, integrity }] = packed;
    assert.match(integrity, /^sha512-[A-Za-z0-9+/=]+$/, "pack must report tarball integrity");
    assert.equal(integrity, PINNED_INTEGRITY, "pack integrity must equal the reviewed candidate");
    evidence.npm.integrity = integrity;
    assert.equal(join(dirname(retainedTarball), filename), retainedTarball, "workflow must retain the canonical tarball name");
    evidence.npm.tarball = { retained: await exists(retainedTarball), filename, bytes: (await stat(retainedTarball)).size };
    const consumer = join(workspace, "consumer with spaces Việt");
    await mkdir(consumer);
    await writeFile(join(consumer, "package.json"), JSON.stringify({ name: "native-installer-smoke", private: true, type: "module" }), "utf8");
    evidence.phase = "npm-install";
    await runNpm(["install", "--ignore-scripts", "--no-audit", "--no-fund", retainedTarball], { cwd: consumer, env: childEnvironment });
    evidence.npm.install_clean = true;

    evidence.phase = "genetics-cli";
    const installedPackage = join(consumer, "node_modules", "@biobirddna", "genetics");
    const input = join(consumer, "cross input Việt.json");
    const output = join(consumer, "cross output Việt.json");
    await writeFile(input, await readFile(join(installedSkill, "examples", "input", "cross_autosomal_recessive.json")));
    await run(process.execPath, [join(installedPackage, "bin", "biobirddna-genetics.mjs"), "cross", "--input", input, "--output", output], { cwd: consumer, env: childEnvironment });
    assert.equal(JSON.parse(await readFile(output, "utf8")).operation, "cross");
    evidence.npm.cli_probes.genetics = true;

    evidence.phase = "feedback-cli";
    const feedbackInput = join(consumer, "feedback input Việt.json");
    const feedbackOutput = join(consumer, "feedback output Việt.json");
    const expectedIntegrity = "c".repeat(64);
    await writeFile(feedbackInput, JSON.stringify({
      schema_version: "1.0.0",
      run_id: "NATIVE-INSTALLER-SMOKE",
      operation: "cross",
      engine_version: "0.2.0",
      input_sha256: "a".repeat(64),
      output_sha256: "b".repeat(64),
      package: { name: "@biobirddna/genetics", version: "0.2.0", integrity_sha256: expectedIntegrity },
      result: { evidence: { evidence_grade: "G0" }, limits: { decision_support_only: true, gene_status: "not_assessed" } },
      human_review_status: "pending",
    }), "utf8");
    await run(process.execPath, [join(installedPackage, "bin", "biobirddna-feedback.mjs"), "--input", feedbackInput, "--expected-integrity", expectedIntegrity, "--output", feedbackOutput], { cwd: consumer, env: childEnvironment });
    assert.equal(JSON.parse(await readFile(feedbackOutput, "utf8")).proposal.status, "pending_human_review");
    evidence.npm.cli_probes.feedback = true;
    evidence.status = "passed";
    evidence.phase = "complete";
  } finally {
    await rm(workspace, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
    await rm(childTemp, { recursive: true, force: true, maxRetries: 5, retryDelay: 200 });
    evidence.cleanup.passed = !(await exists(workspace)) && !(await exists(childTemp));
    const serialized = JSON.stringify(evidence);
    evidence.log_bytes = Buffer.byteLength(serialized);
    evidence.credential_free_log = credentialFree(serialized) && evidence.log_bytes <= 128 * 1024;
    await writeEvidence(evidencePath, evidence);
  }
  assert.equal(evidence.cleanup.passed, true, "temporary native installer state must be removed");
  assert.equal(evidence.credential_free_log, true, "reported evidence must remain bounded and credential-free");
  process.stdout.write("native-installer-smoke: passed\n");
}

main().catch(async (error) => {
  const evidence = {
    status: "failed",
    phase: "unavailable",
    error: safeErrorCode(error),
    credential_free_log: true,
    log_bytes: 0,
  };
  try {
    if (!evidencePath || !(await exists(evidencePath))) {
      await writeEvidence(evidencePath, evidence);
    }
  } catch {
    // Evidence writing is best effort; never disclose child-process output.
  }
  process.stderr.write("native-installer-smoke: failed\n");
  process.exitCode = 1;
});
