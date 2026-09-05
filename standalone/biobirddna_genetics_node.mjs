#!/usr/bin/env node
/**
 * BioBirdDNA Genetics — dependency-free Node.js cross simulator.
 *
 * Scope: known/uncertain one-locus crosses for autosomal recessive/dominant,
 * incomplete dominance, and avian Z-linked recessive/dominant models.
 * This file does not confirm genes and does not rank real breeding pairs.
 */

import { readFile } from "node:fs/promises";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";

const MODES = new Set([
  "autosomal_recessive",
  "autosomal_dominant",
  "incomplete_dominant",
  "z_linked_recessive",
  "z_linked_dominant",
  "non_genetic",
]);

function fail(code) {
  throw new Error(code);
}

function modelFrom(value) {
  if (!value || typeof value !== "object" || Array.isArray(value)) fail("model_invalid");
  if (!MODES.has(value.mode)) fail("mode_invalid");
  const model = {
    model_id: String(value.model_id ?? value.mode),
    mode: value.mode,
    normal: String(value.normal_allele ?? "N"),
    variant: String(value.variant_allele ?? "m"),
    normal_label: String(value.normal_label ?? "normal"),
    target_label: String(value.target_label ?? "target_colour"),
    intermediate_label: String(value.intermediate_label ?? "intermediate"),
    penetrance: Number(value.penetrance ?? 1),
    female_probability: Number(value.female_probability ?? 0.5),
    baseline_target_probability: Number(value.baseline_target_probability ?? 0.05),
  };
  if (
    !Number.isFinite(model.female_probability) ||
    model.female_probability < 0 ||
    model.female_probability > 1
  ) {
    fail("female_probability_must_be_between_0_and_1");
  }
  // Z-linked offspring sex arises from maternal Z/W segregation, so the simulated
  // sex ratio is fixed at 0.5 and a configurable female_probability cannot be
  // honored. Fail closed instead of silently ignoring the value (mirrors
  // src/biobirddna_genetics/models.py InheritanceModel.__post_init__ and
  // npm/src/index.mjs modelFrom — keep in sync).
  if (model.mode.startsWith("z_linked") && model.female_probability !== 0.5) {
    fail("z_linked_female_probability_must_be_0_5");
  }
  return model;
}

function canonical(model, sex, genotype) {
  if (model.mode === "non_genetic") return "NA";
  if (typeof genotype !== "string" || genotype.length === 0) fail("genotype_required");
  const tokens = genotype.split("/").map((token) => token.trim());
  const allowed = new Set([model.normal, model.variant]);
  const zLinked = model.mode.startsWith("z_linked");
  if (zLinked && sex === "female") {
    if (tokens.length !== 2 || tokens[1] !== "W" || !allowed.has(tokens[0])) {
      fail("z_linked_female_genotype_invalid");
    }
    return `${tokens[0]}/W`;
  }
  if (tokens.length !== 2 || tokens.some((token) => !allowed.has(token))) {
    fail("genotype_invalid");
  }
  const order = new Map([[model.normal, 0], [model.variant, 1]]);
  tokens.sort((left, right) => order.get(left) - order.get(right));
  return tokens.join("/");
}

function gametes(model, sex, genotype) {
  genotype = canonical(model, sex, genotype);
  if (genotype === "NA") return { NA: 1 };
  const tokens = genotype.split("/");
  if (model.mode.startsWith("z_linked") && sex === "female") {
    return { [tokens[0]]: 0.5, W: 0.5 };
  }
  const counts = new Map();
  for (const token of tokens) counts.set(token, (counts.get(token) ?? 0) + 1);
  return Object.fromEntries([...counts].map(([token, count]) => [token, count / 2]));
}

function phenotype(model, sex, genotype) {
  if (model.mode === "non_genetic") {
    const p = model.baseline_target_probability;
    return { [model.target_label]: p, [model.normal_label]: 1 - p };
  }
  genotype = canonical(model, sex, genotype);
  const copies = genotype
    .split("/")
    .filter((token) => token === model.variant).length;
  let affected = false;
  if (model.mode === "autosomal_recessive" || model.mode === "z_linked_recessive") {
    const required = model.mode === "z_linked_recessive" && sex === "female" ? 1 : 2;
    affected = copies >= required;
  } else if (model.mode === "autosomal_dominant" || model.mode === "z_linked_dominant") {
    affected = copies >= 1;
  } else if (model.mode === "incomplete_dominant") {
    if (copies === 1) return { [model.intermediate_label]: 1 };
    affected = copies === 2;
  }
  if (affected) {
    return {
      [model.target_label]: model.penetrance,
      [model.normal_label]: 1 - model.penetrance,
    };
  }
  return { [model.normal_label]: 1 };
}

function add(destination, key, value) {
  destination[key] = (destination[key] ?? 0) + value;
}

export function crossKnown(modelInput, maleInput, femaleInput) {
  const model = modelFrom(modelInput);
  if (model.mode === "non_genetic") {
    const f = model.female_probability;
    const p = model.baseline_target_probability;
    return {
      model_id: model.model_id,
      mode: model.mode,
      genotype_probabilities: {},
      phenotype_probabilities: {
        [`male:${model.target_label}`]: (1 - f) * p,
        [`male:${model.normal_label}`]: (1 - f) * (1 - p),
        [`female:${model.target_label}`]: f * p,
        [`female:${model.normal_label}`]: f * (1 - p),
      },
      limits: { gene_status: "not_assessed", inheritance_status: "baseline_only" },
    };
  }
  const male = canonical(model, "male", maleInput);
  const female = canonical(model, "female", femaleInput);
  const father = gametes(model, "male", male);
  const mother = gametes(model, "female", female);
  const genotypes = {};
  const phenotypes = {};

  for (const [paternal, paternalProbability] of Object.entries(father)) {
    for (const [maternal, maternalProbability] of Object.entries(mother)) {
      const base = paternalProbability * maternalProbability;
      if (model.mode.startsWith("z_linked")) {
        const sex = maternal === "W" ? "female" : "male";
        const child = maternal === "W"
          ? canonical(model, sex, `${paternal}/W`)
          : canonical(model, sex, `${paternal}/${maternal}`);
        add(genotypes, `${sex}:${child}`, base);
        for (const [label, conditional] of Object.entries(phenotype(model, sex, child))) {
          add(phenotypes, `${sex}:${label}`, base * conditional);
        }
      } else {
        const child = canonical(model, "male", `${paternal}/${maternal}`);
        for (const [sex, sexProbability] of [["male", 1 - model.female_probability], ["female", model.female_probability]]) {
          const probability = base * sexProbability;
          add(genotypes, `${sex}:${child}`, probability);
          for (const [label, conditional] of Object.entries(phenotype(model, sex, child))) {
            add(phenotypes, `${sex}:${label}`, probability * conditional);
          }
        }
      }
    }
  }

  return {
    model_id: model.model_id,
    mode: model.mode,
    male_genotype: male,
    female_genotype: female,
    genotype_probabilities: Object.fromEntries(Object.entries(genotypes).sort()),
    phenotype_probabilities: Object.fromEntries(Object.entries(phenotypes).sort()),
    limits: {
      gene_status: "model_assumption_only",
      inheritance_status: "simulation_not_confirmation",
    },
  };
}

function normalizeDistribution(values) {
  if (!values || typeof values !== "object" || Array.isArray(values)) fail("distribution_invalid");
  const entries = Object.entries(values).map(([key, value]) => [key, Number(value)]);
  if (entries.some(([, value]) => !Number.isFinite(value) || value < 0)) fail("distribution_invalid");
  const total = entries.reduce((sum, [, value]) => sum + value, 0);
  if (!(total > 0)) fail("distribution_zero_total");
  return Object.fromEntries(entries.map(([key, value]) => [key, value / total]));
}

export function crossUncertain(modelInput, maleInputs, femaleInputs) {
  const model = modelFrom(modelInput);
  const maleRaw = normalizeDistribution(maleInputs);
  const femaleRaw = normalizeDistribution(femaleInputs);
  const male = {};
  const female = {};
  for (const [genotype, probability] of Object.entries(maleRaw)) {
    const key = canonical(model, "male", genotype);
    male[key] = (male[key] ?? 0) + probability;
  }
  for (const [genotype, probability] of Object.entries(femaleRaw)) {
    const key = canonical(model, "female", genotype);
    female[key] = (female[key] ?? 0) + probability;
  }
  const genotypes = {};
  const phenotypes = {};
  for (const [maleGenotype, maleProbability] of Object.entries(male)) {
    for (const [femaleGenotype, femaleProbability] of Object.entries(female)) {
      const weight = maleProbability * femaleProbability;
      const result = crossKnown(model, maleGenotype, femaleGenotype);
      for (const [key, probability] of Object.entries(result.genotype_probabilities)) {
        add(genotypes, key, weight * probability);
      }
      for (const [key, probability] of Object.entries(result.phenotype_probabilities)) {
        add(phenotypes, key, weight * probability);
      }
    }
  }
  return {
    model_id: model.model_id,
    mode: model.mode,
    male_genotype_posterior: male,
    female_genotype_posterior: female,
    genotype_probabilities: Object.fromEntries(Object.entries(genotypes).sort()),
    phenotype_probabilities: Object.fromEntries(Object.entries(phenotypes).sort()),
    limits: {
      gene_status: "inferred_parent_genotypes",
      inheritance_status: "simulation_not_confirmation",
    },
  };
}

async function main() {
  const [operation, inputPath] = process.argv.slice(2);
  if (operation !== "cross" || !inputPath) {
    console.error("Usage: node biobirddna_genetics_node.mjs cross <input.json>");
    process.exitCode = 2;
    return;
  }
  const request = JSON.parse(await readFile(inputPath, "utf8"));
  const result = request.male_genotypes
    ? crossUncertain(request.model, request.male_genotypes, request.female_genotypes)
    : crossKnown(request.model, request.male_genotype, request.female_genotype);
  process.stdout.write(`${JSON.stringify(result, null, 2)}\n`);
}

if (
  process.argv[1] &&
  import.meta.url === pathToFileURL(resolve(process.argv[1])).href
) {
  main().catch((error) => {
    console.error(JSON.stringify({ ok: false, error: error.message }, null, 2));
    process.exitCode = 1;
  });
}
