import assert from "node:assert/strict";
import test from "node:test";

import { createFeedbackArtifact } from "../src/index.mjs";

const hash = (character) => character.repeat(64);
const feedbackEnvelope = () => ({
  schema_version: "1.0.0",
  run_id: "CM-RUN-NODE-HARDENING",
  operation: "cross",
  engine_version: "0.2.0",
  input_sha256: hash("a"),
  output_sha256: hash("b"),
  package: {
    name: "@biobirddna/genetics",
    version: "0.2.0",
    integrity_sha256: hash("c")
  },
  result: {
    evidence: { evidence_grade: "G0" },
    limits: { decision_support_only: true, gene_status: "not_assessed" },
    missing_data: []
  },
  human_review_status: "pending"
});

const withOutput = (output) => {
  const envelope = feedbackEnvelope();
  envelope.result.output = output;
  return envelope;
};

test("feedback rejects disguised mutation intents", () => {
  for (const output of [
    { mutationIntent: { type: "apply" } },
    { command: "delete" },
    { intent: "approve" }
  ]) {
    assert.throws(
      () => createFeedbackArtifact(withOutput(output), hash("c")),
      /feedback_mutation_request_forbidden/,
      JSON.stringify(output)
    );
  }
});

for (const output of [
  { action: "patch" },
  { requestedAction: "modify" }
]) {
  test(`feedback rejects mutation directive ${JSON.stringify(output)}`, () => {
    assert.throws(
      () => createFeedbackArtifact(withOutput(output), hash("c")),
      /feedback_mutation_request_forbidden/
    );
  });
}

test("feedback rejects structured confirmed gene or variant claims below G5", () => {
  for (const output of [
    { geneClaim: { confirmed: true } },
    { variant: { causal: true } },
    { claim: { identified: true } }
  ]) {
    assert.throws(
      () => createFeedbackArtifact(withOutput(output), hash("c")),
      /feedback_confirmed_gene_claim_requires_G5/,
      JSON.stringify(output)
    );
  }
});

for (const output of [
  { finding: "TYR is confirmed as causal" },
  { finding: "locus identified" }
]) {
  test(`feedback rejects unsupported claim ${JSON.stringify(output)} below G5`, () => {
    assert.throws(
      () => createFeedbackArtifact(withOutput(output), hash("c")),
      /feedback_confirmed_gene_claim_requires_G5/
    );
  });
}

test("feedback rejects camel-case and Pascal-case AutoResearch health scores", () => {
  for (const output of [
    { autoResearchScore: 3.2 },
    { AutoResearchHealthGateScore: 3.2 },
    { AutoResearch: { HealthGate: { score: 3.2 } } }
  ]) {
    assert.throws(
      () => createFeedbackArtifact(withOutput(output), hash("c")),
      /feedback_numerical_score_forbidden/,
      JSON.stringify(output)
    );
  }
});

test("feedback still permits legitimate domain scores and fractions", () => {
  assert.doesNotThrow(() => createFeedbackArtifact(withOutput({
    score: 0.8,
    fraction: 0.2,
    domain: {
      phenotypeScore: 0.7,
      observedFraction: 0.5
    }
  }), hash("c")));
});

test("feedback still permits non-directive observations and unresolved candidate findings", () => {
  assert.doesNotThrow(() => createFeedbackArtifact(withOutput({
    observation: "A pale plumage patch was observed.",
    finding: "TYR remains a candidate locus; causality is not confirmed."
  }), hash("c")));
});
