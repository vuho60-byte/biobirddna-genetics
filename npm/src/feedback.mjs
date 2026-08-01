import { createHash, timingSafeEqual } from "node:crypto";

const PACKAGE_NAME = "@biobirddna/genetics";
const PACKAGE_VERSION = "0.2.0";
const MAX_INPUT_BYTES = 1024 * 1024;
const TEXT_ENCODER = new TextEncoder();
const OPERATIONS = new Set(["cross", "infer", "pigment", "rank-mates", "normalize-records", "simulate-recovery"]);
const ROOT_FIELDS = new Set(["schema_version", "run_id", "operation", "engine_version", "input_sha256", "output_sha256", "package", "result", "human_review_status"]);
const PACKAGE_FIELDS = new Set(["name", "version", "integrity_sha256"]);
const GRADES = new Set(["G0", "G1", "G2", "G3", "G4", "G5", "G6", "not_assessed"]);

const fail = (code) => { throw new Error(code); };
const object = (value, code) => value !== null && typeof value === "object" && !Array.isArray(value) ? value : fail(code);
const string = (value, code) => typeof value === "string" && value.trim() ? value.trim() : fail(code);
const rawString = (value, code) => typeof value === "string" && value.length > 0 && value.trim() === value ? value : fail(code);
const hash = (value, code) => /^[0-9a-f]{64}$/.test(rawString(value, code)) ? value : fail(code);
const canonical = (value) => Array.isArray(value) ? value.map(canonical) : value && typeof value === "object" ? Object.fromEntries(Object.keys(value).sort().map((key) => [key, canonical(value[key])])) : value;
const stableJson = (value) => JSON.stringify(canonical(value));
function freeze(value) { if (value && typeof value === "object" && !Object.isFrozen(value)) { Object.values(value).forEach(freeze); Object.freeze(value); } return value; }
function exactFields(value, fields, prefix) { for (const key of Object.keys(object(value, prefix))) if (!fields.has(key)) fail(`${prefix}_unknown_field:${key}`); }
function jsonValue(value) { if (value === null || typeof value === "string" || typeof value === "boolean") return true; if (typeof value === "number") return Number.isFinite(value); if (Array.isArray(value)) return value.every(jsonValue); return value && typeof value === "object" && Object.values(value).every(jsonValue); }
const tokenOf = (key) => key.replace(/([a-z])([A-Z])/g, "$1_$2").replace(/-/g, "_").toLowerCase();
const compactTokenOf = (key) => tokenOf(key).replace(/_/g, "");
const mutationWord = (value) => typeof value === "string" && /\b(?:apply|approve|approved|write|delete|update|mutate|mutation|patch|modify)\b/i.test(value);
function mutationObject(value) {
  if (!value || typeof value !== "object" || Array.isArray(value)) return false;
  const fields = Object.fromEntries(Object.entries(value).map(([key, child]) => [tokenOf(key), child]));
  if (fields.apply === true || fields.approve === true || fields.approval === true || fields.write || fields.delete || fields.update || fields.mutate || fields.mutation) return true;
  return ["type", "kind", "action", "status"].some((key) => mutationWord(fields[key]));
}
function rejectMutation(value) {
  if (Array.isArray(value)) { value.forEach(rejectMutation); return; }
  if (!value || typeof value !== "object") return;
  for (const [key, child] of Object.entries(value)) {
    const token = tokenOf(key);
    if (token === "target_files" || token === "target_file") { if (!Array.isArray(child) || child.length !== 0) fail("feedback_mutation_request_forbidden"); continue; }
    if (token === "auto_apply") { if (child !== false) fail("feedback_mutation_request_forbidden"); continue; }
    if (new Set(["mutation", "mutate", "mutation_request", "apply_request", "approve_request", "write_request", "delete_request", "update_request"]).has(token)) fail("feedback_mutation_request_forbidden");
    if (new Set(["apply", "approve", "approval"]).has(token) && child === true) fail("feedback_mutation_request_forbidden");
    if (new Set(["write", "delete", "update"]).has(token) && mutationObject({ [token]: child })) fail("feedback_mutation_request_forbidden");
    if (new Set(["action", "requested_action", "request_action", "proposal_status", "request", "proposal", "mutation_intent", "command", "intent"]).has(token) && (mutationWord(child) || mutationObject(child))) fail("feedback_mutation_request_forbidden");
    rejectMutation(child);
  }
}
function hasConfirmedClaim(value) {
  if (typeof value === "string") return /\b(?:confirmed[ _-]gene|gene[ _-]confirmed|causal[ _-](?:gene|variant)|(?:gene|variant)[ _-]causal|(?:gene|variant|locus)[ _-]identified|confirmed[ _-]as[ _-]causal)\b/i.test(value);
  if (Array.isArray(value)) return value.some(hasConfirmedClaim);
  return value && typeof value === "object" && Object.entries(value).some(([key, child]) => {
    const token = tokenOf(key), compact = compactTokenOf(key), fields = child && typeof child === "object" && !Array.isArray(child) ? Object.fromEntries(Object.entries(child).map(([name, item]) => [compactTokenOf(name), item])) : {};
    const structuredClaim = ["gene", "variant", "geneclaim", "variantclaim", "claim"].includes(compact)
      && (["confirmed", "causal", "identified"].some((state) => fields[state] === true)
        || ["status", "state", "type", "kind"].some((field) => typeof fields[field] === "string" && /(?:confirmed|identified|causal)/i.test(fields[field])));
    return /(?:confirmed[ _-]gene|gene[ _-]confirmed|gene[ _-]confirmation|causal[ _-](?:gene|variant)|(?:gene|variant)[ _-]causal)/i.test(key) || (token === "gene_status" && typeof child === "string" && /(?:confirmed|identified|causal)/i.test(child)) || structuredClaim || (token === "claim" && fields.subject === "gene" && typeof fields.status === "string" && /(?:confirmed|identified|causal)/i.test(fields.status)) || hasConfirmedClaim(child);
  });
}
function rejectNumericalScores(value, ancestry = []) {
  if (Array.isArray(value)) { value.forEach((item) => rejectNumericalScores(item, ancestry)); return; }
  if (!value || typeof value !== "object") return;
  for (const [key, child] of Object.entries(value)) {
    const token = compactTokenOf(key);
    if (new Set(["healthgatescore", "auto" + "researchscore", "auto" + "researchhealthgate", "auto" + "researchhealthgatescore", "tai" + "chelscore", "promotionscore"]).has(token)) fail("feedback_numerical_score_forbidden");
    if (token === "score" && ancestry.includes("auto" + "research") && ancestry.includes("healthgate")) fail("feedback_numerical_score_forbidden");
    rejectNumericalScores(child, [...ancestry, token]);
  }
}
function equalHash(left, right) { return timingSafeEqual(TEXT_ENCODER.encode(left), TEXT_ENCODER.encode(right)); }

export function createFeedbackArtifact(envelope, expectedIntegrity) {
  exactFields(envelope, ROOT_FIELDS, "feedback_envelope");
  if (!jsonValue(envelope)) fail("feedback_envelope_non_json_value");
  if (TEXT_ENCODER.encode(stableJson(envelope)).byteLength > MAX_INPUT_BYTES) fail("feedback_input_too_large");
  rejectMutation(envelope);
  rejectNumericalScores(envelope);
  if (envelope.schema_version !== "1.0.0") fail("feedback_schema_version_unsupported");
  const run_id = rawString(envelope.run_id, "feedback_run_id_invalid");
  if (!OPERATIONS.has(envelope.operation)) fail("feedback_operation_unsupported");
  const engine_version = rawString(envelope.engine_version, "feedback_engine_version_mismatch"); if (engine_version !== PACKAGE_VERSION) fail("feedback_engine_version_mismatch");
  const input_sha256 = hash(envelope.input_sha256, "feedback_input_sha256_invalid"), output_sha256 = hash(envelope.output_sha256, "feedback_output_sha256_invalid");
  exactFields(envelope.package, PACKAGE_FIELDS, "feedback_package");
  if (envelope.package.name !== PACKAGE_NAME) fail("feedback_package_name_mismatch");
  if (envelope.package.version !== PACKAGE_VERSION) fail("feedback_package_version_mismatch");
  const package_integrity_sha256 = hash(envelope.package.integrity_sha256, "feedback_package_integrity_invalid"), expected = hash(expectedIntegrity, "feedback_expected_integrity_invalid");
  if (!equalHash(package_integrity_sha256, expected)) fail("feedback_package_integrity_mismatch");
  const result = object(envelope.result, "feedback_result_invalid"), sourceEvidence = object(result.evidence, "feedback_result_evidence_required"), limits = object(result.limits, "feedback_result_limits_required");
  if (limits.decision_support_only !== true) fail("feedback_decision_support_only_required");
  if (limits.gene_status !== "not_assessed") fail("feedback_gene_status_not_assessed_required");
  if (envelope.human_review_status !== "pending") fail("feedback_human_review_status_required");
  const evidence_grade = sourceEvidence.evidence_grade ?? "not_assessed";
  if (!GRADES.has(evidence_grade)) fail("feedback_evidence_grade_invalid");
  if (hasConfirmedClaim(envelope) && (!/^G[5-6]$/.test(evidence_grade))) fail("feedback_confirmed_gene_claim_requires_G5");
  const missing_data = result.missing_data === undefined ? [] : Array.isArray(result.missing_data) ? result.missing_data : fail("feedback_missing_data_invalid");
  const evidence = { evidence_grade, missing_data: canonical(missing_data), decision_support_only: true, human_review_status: "pending", limits: canonical(limits) };
  const signals = Object.keys(sourceEvidence).length === 0 || missing_data.length > 0 ? [{ code: "INSUFFICIENT_DATA", severity: "warning", rationale: "Production evidence is missing or incomplete; no numerical score is emitted." }] : [];
  const source = { package: PACKAGE_NAME, package_version: PACKAGE_VERSION, package_integrity_sha256 };
  const run = { run_id, operation: envelope.operation, input_sha256, output_sha256, engine_version, schema_version: envelope.schema_version };
  const feedback_id = createHash("sha256").update(stableJson({ source, run, evidence, signals }), "utf8").digest("hex");
  return freeze({ schema_version: "biobirddna_feedback_v1", feedback_id, source, run, evidence, signals, proposal: { status: "pending_human_review", auto_apply: false, target_files: [] } });
}
