import { createHash } from "node:crypto";
import { adaptBreedingRecordError, normalizeBreedingRecords } from "./breeding-records.mjs";
export { createFeedbackArtifact } from "./feedback.mjs";

export const ENGINE_VERSION = "0.4.0";
const MODES = new Set(["autosomal_recessive", "autosomal_dominant", "incomplete_dominant", "z_linked_recessive", "z_linked_dominant", "non_genetic"]);
const SEXES = new Set(["male", "female", "unknown"]);
const STAGES = ["eggs_laid", "fertile_eggs", "hatched", "survived_to_first_molt"];
const DEFAULT_IDS = ["M1_AR", "M2_AD", "M3_INC", "M4_ZR", "M5_ZD", "M0_NONGENETIC"];
const INFERENCE_ROOT_FIELDS = new Set(["target_label", "male_parent_phenotype", "female_parent_phenotype", "models", "model_priors", "offspring"]);
const INFERENCE_OBSERVATION_FIELDS = new Set(["sex", "phenotype", "count"]);
const MODEL_FIELDS = new Set(["model_id", "mode", "normal_allele", "variant_allele", "normal_label", "target_label", "intermediate_label", "penetrance", "observation_error", "female_probability", "baseline_target_probability", "description"]);
const CROSS_ROOT_FIELDS = new Set(["model", "male_genotype", "female_genotype", "male_genotypes", "female_genotypes"]);
const RANK_ROOT_FIELDS = new Set(["pedigree", "males", "females", "hypotheses", "max_offspring_inbreeding", "block_direct_ancestor_and_full_sibling", "weights"]);
const PEDIGREE_FIELDS = new Set(["individual_id", "sex", "sire_id", "dam_id"]);
const PROFILE_FIELDS = new Set(["individual_id", "sex", "genotype_posteriors", "phenotype_label", "health_risk", "fertility_risk", "welfare_blocked", "notes"]);
const HYPOTHESIS_FIELDS = new Set(["model", "weight"]);
const PAIRING_WEIGHT_FIELDS = new Set(["target", "information", "inbreeding", "health", "uncertainty"]);
const MAX_OFFSPRING_OBSERVATION_COUNT = 1_000_000;
const NORMALIZED_MODELS = new WeakSet();

export class BioBirdDNAError extends Error {
  constructor(code) { super(code); this.name = "BioBirdDNAError"; }
}
const fail = (code) => { throw new BioBirdDNAError(code); };
const isObject = (value) => value !== null && typeof value === "object" && !Array.isArray(value);
const object = (value, code) => isObject(value) ? value : fail(`${code}_must_be_object`);
const array = (value, code) => Array.isArray(value) ? value : fail(`${code}_must_be_array`);
const text = (value, code) => typeof value === "string" && value.trim() ? value.trim() : fail(`${code}_required`);
const exactText = (value, code) => typeof value === "string" && value.length > 0 && value.trim() === value ? value : fail(`${code}_required`);
const finite = (value, code, min = -Infinity, max = Infinity) => typeof value === "number" && Number.isFinite(value) && value >= min && value <= max ? value : fail(`${code}_invalid`);
const integer = (value, code, min = 0, max = Number.MAX_SAFE_INTEGER) => Number.isInteger(value) && value >= min && value <= max ? value : fail(`${code}_invalid`);
const sortedObject = (value) => Object.fromEntries(Object.entries(value).sort(([a], [b]) => a.localeCompare(b)));

function modelFrom(value) {
  if (isObject(value) && NORMALIZED_MODELS.has(value)) return value;
  value = object(value, "model");
  rejectUnknown(value, MODEL_FIELDS, "model");
  if (!MODES.has(value.mode)) fail("mode_invalid");
  const optional = (key, fallback) => Object.hasOwn(value, key) ? value[key] : fallback;
  const model = {
    model_id: exactText(value.model_id, "model_id"), mode: value.mode,
    normal: exactText(optional("normal_allele", "N"), "normal_allele"), variant: exactText(optional("variant_allele", "m"), "variant_allele"),
    normal_label: exactText(optional("normal_label", "normal"), "normal_label"), target_label: exactText(optional("target_label", "target_colour"), "target_label"),
    intermediate_label: exactText(optional("intermediate_label", "intermediate"), "intermediate_label"),
    penetrance: finite(optional("penetrance", 1), "penetrance", 0, 1), female_probability: finite(optional("female_probability", 0.5), "female_probability", 0, 1),
    baseline_target_probability: finite(optional("baseline_target_probability", 0.05), "baseline_target_probability", 0, 1), observation_error: finite(optional("observation_error", 0.01), "observation_error", 0, 1)
  };
  if (Object.hasOwn(value, "description") && typeof value.description !== "string") fail("description_invalid");
  if (model.normal === model.variant) fail("alleles_must_differ");
  if (model.normal_label === model.target_label || (model.mode === "incomplete_dominant" && (model.intermediate_label === model.normal_label || model.intermediate_label === model.target_label))) fail("phenotype_labels_must_be_distinct");
  // Z-linked offspring sex arises from maternal Z/W segregation, so the simulated sex ratio is fixed at
  // 0.5 and a configurable female_probability cannot be honored. Fail closed instead of silently ignoring
  // the value (mirrors src/biobirddna_genetics/models.py InheritanceModel.__post_init__ — keep in sync).
  if (model.mode.startsWith("z_linked") && model.female_probability !== 0.5) fail("z_linked_female_probability_must_be_0_5");
  NORMALIZED_MODELS.add(model);
  return model;
}

function canonical(model, sex, genotype) {
  if (model.mode === "non_genetic") return genotype === "NA" ? "NA" : fail("non_genetic_model_requires_NA_genotype");
  if (typeof genotype !== "string" || !genotype.trim()) fail("genotype_required");
  const tokens = genotype.split("/").map((part) => part.trim());
  const allowed = new Set([model.normal, model.variant]);
  if (model.mode.startsWith("z_linked") && sex === "female") {
    if (tokens.length !== 2 || tokens[1] !== "W" || !allowed.has(tokens[0])) fail("z_linked_female_genotype_invalid");
    return `${tokens[0]}/W`;
  }
  if (tokens.length !== 2 || tokens.some((token) => !allowed.has(token))) fail("genotype_invalid");
  tokens.sort((a, b) => (a === model.normal ? -1 : 1) - (b === model.normal ? -1 : 1));
  return tokens.join("/");
}

function gametes(model, sex, genotype) {
  const normalized = canonical(model, sex, genotype);
  if (normalized === "NA") return { NA: 1 };
  const tokens = normalized.split("/");
  if (model.mode.startsWith("z_linked") && sex === "female") return { [tokens[0]]: 0.5, W: 0.5 };
  const out = accumulator();
  for (const allele of tokens) out[allele] = (out[allele] ?? 0) + 0.5;
  return out;
}

function phenotype(model, sex, genotype) {
  if (model.mode === "non_genetic") return { [model.target_label]: model.baseline_target_probability, [model.normal_label]: 1 - model.baseline_target_probability };
  const copies = canonical(model, sex, genotype).split("/").filter((allele) => allele === model.variant).length;
  if (model.mode === "incomplete_dominant" && copies === 1) return { [model.intermediate_label]: 1 };
  const affected = model.mode.endsWith("recessive") ? copies >= (model.mode.startsWith("z_linked") && sex === "female" ? 1 : 2) : copies >= 1;
  return affected ? { [model.target_label]: model.penetrance, [model.normal_label]: 1 - model.penetrance } : { [model.normal_label]: 1 };
}

const accumulator = () => Object.create(null);
const getOwn = (source, key, fallback) => Object.hasOwn(source, key) ? source[key] : fallback;
const add = (target, key, value) => { target[key] = getOwn(target, key, 0) + value; };

export function crossKnown(modelInput, maleInput, femaleInput) {
  const model = modelFrom(modelInput);
  const male = canonical(model, "male", maleInput), female = canonical(model, "female", femaleInput);
  if (model.mode === "non_genetic") {
    const f = model.female_probability, p = model.baseline_target_probability;
    return { model_id: model.model_id, mode: model.mode, male_genotype: male, female_genotype: female, genotype_probabilities: {}, phenotype_probabilities: {
      [`male:${model.target_label}`]: (1 - f) * p, [`male:${model.normal_label}`]: (1 - f) * (1 - p),
      [`female:${model.target_label}`]: f * p, [`female:${model.normal_label}`]: f * (1 - p)
    }, limits: { gene_status: "not_assessed", interpretation: "non_genetic_baseline_model" } };
  }
  const genotypes = accumulator(), phenotypes = accumulator();
  for (const [paternal, pp] of Object.entries(gametes(model, "male", male))) for (const [maternal, mp] of Object.entries(gametes(model, "female", female))) {
    const base = pp * mp;
    if (model.mode.startsWith("z_linked")) {
      const sex = maternal === "W" ? "female" : "male";
      const child = canonical(model, sex, maternal === "W" ? `${paternal}/W` : `${paternal}/${maternal}`);
      add(genotypes, `${sex}:${child}`, base);
      for (const [label, probability] of Object.entries(phenotype(model, sex, child))) add(phenotypes, `${sex}:${label}`, base * probability);
    } else for (const [sex, sp] of [["male", 1 - model.female_probability], ["female", model.female_probability]]) {
      const child = canonical(model, sex, `${paternal}/${maternal}`), probability = base * sp;
      add(genotypes, `${sex}:${child}`, probability);
      for (const [label, conditional] of Object.entries(phenotype(model, sex, child))) add(phenotypes, `${sex}:${label}`, probability * conditional);
    }
  }
  return { model_id: model.model_id, mode: model.mode, male_genotype: male, female_genotype: female, genotype_probabilities: sortedObject(genotypes), phenotype_probabilities: sortedObject(phenotypes), limits: { gene_status: "model_assumption_only", inheritance_status: "simulation_not_confirmation" } };
}

function distribution(value) {
  value = object(value, "distribution");
  const entries = Object.entries(value).map(([key, probability]) => {
    if (typeof probability !== "number" || !Number.isFinite(probability)) fail(Number.isFinite(probability) ? "probability_distribution_invalid" : typeof probability === "number" ? "probability_distribution_non_finite" : "probability_distribution_invalid");
    if (probability < 0) fail("probability_distribution_negative");
    return [key, probability];
  });
  const total = entries.reduce((sum, [, probability]) => sum + probability, 0); if (!(total > 0)) fail("distribution_zero_total");
  return Object.fromEntries(entries.map(([key, probability]) => [key, probability / total]));
}

export function crossUncertain(modelInput, males, females) {
  const model = modelFrom(modelInput), male = accumulator(), female = accumulator();
  for (const [genotype, probability] of Object.entries(distribution(males))) add(male, canonical(model, "male", genotype), probability);
  for (const [genotype, probability] of Object.entries(distribution(females))) add(female, canonical(model, "female", genotype), probability);
  if (model.mode === "non_genetic") return crossKnown(model, "NA", "NA");
  const genotypes = accumulator(), phenotypes = accumulator(), component_crosses = [];
  for (const [mg, mp] of Object.entries(male)) for (const [fg, fp] of Object.entries(female)) {
    const result = crossKnown(model, mg, fg), weight = mp * fp;
    component_crosses.push({ male_genotype: mg, female_genotype: fg, weight });
    for (const [key, probability] of Object.entries(result.genotype_probabilities)) add(genotypes, key, weight * probability);
    for (const [key, probability] of Object.entries(result.phenotype_probabilities)) add(phenotypes, key, weight * probability);
  }
  return { model_id: model.model_id, mode: model.mode, male_genotype_posterior: sortedObject(male), female_genotype_posterior: sortedObject(female), genotype_probabilities: sortedObject(genotypes), phenotype_probabilities: sortedObject(phenotypes), component_crosses, limits: { gene_status: "inferred_parent_genotypes", inheritance_status: "simulation_not_confirmation" } };
}

function defaultModels(target = "target_colour", normal = "normal", intermediate = "intermediate") {
  return [
    ["M1_AR", "autosomal_recessive"], ["M2_AD", "autosomal_dominant"], ["M3_INC", "incomplete_dominant"], ["M4_ZR", "z_linked_recessive"], ["M5_ZD", "z_linked_dominant"], ["M0_NONGENETIC", "non_genetic"]
  ].map(([model_id, mode]) => modelFrom({ model_id, mode, target_label: target, normal_label: normal, intermediate_label: intermediate }));
}

// Cross-runtime numeric primitives. Each helper below is mirrored
// operation-for-operation in src/biobirddna_genetics/inference.py /
// simulation.py (_portable_log, _portable_exp, _neumaier_sum, _logsumexp,
// _log_gamma, _entropy_bits, _quantile, _order_independent_mean) so both
// runtimes produce bit-identical floats by construction — keep in sync.
//
// Why not the platform libm: glibc and V8 disagree in the last bit for
// roughly 3-10% of log/exp inputs (measured on Linux x86-64), so any
// likelihood built on Math.log/math.log can never be exactly reproducible
// across runtimes. The fdlibm algorithms below use only IEEE-754
// correctly-rounded operations (+, -, *, /) plus explicit bit manipulation,
// which both runtimes define identically.
const DOUBLE_WORDS = new DataView(new ArrayBuffer(8));
function doubleWords(value) { DOUBLE_WORDS.setFloat64(0, value, true); return [DOUBLE_WORDS.getUint32(4, true), DOUBLE_WORDS.getUint32(0, true)]; }
function wordsDouble(high, low) { DOUBLE_WORDS.setUint32(4, high >>> 0, true); DOUBLE_WORDS.setUint32(0, low >>> 0, true); return DOUBLE_WORDS.getFloat64(0, true); }
const LOG_LN2_HI = 6.93147180369123816490e-1, LOG_LN2_LO = 1.90821492927058770002e-10;
const LOG_LG1 = 6.666666666666735130e-1, LOG_LG2 = 3.999999999940941908e-1, LOG_LG3 = 2.857142874366239149e-1, LOG_LG4 = 2.222219843214978396e-1, LOG_LG5 = 1.818357216161805012e-1, LOG_LG6 = 1.531383769920937332e-1, LOG_LG7 = 1.479819860511658591e-1;
const TWO54 = 1.8014398509481984e16;
// fdlibm e_log.c (Sun Microsystems), ported verbatim; mirrored in
// src/biobirddna_genetics/inference.py _portable_log.
function portableLog(value) {
  let x = value, [high, low] = doubleWords(x), k = 0;
  if (high < 0x00100000 || high >= 0x80000000) {
    if (((high & 0x7fffffff) | low) === 0) return -Infinity;
    if (high >= 0x80000000) return NaN;
    k -= 54; x *= TWO54; [high, low] = doubleWords(x);
  }
  if (high >= 0x7ff00000) return x + x;
  k += (high >> 20) - 1023;
  high &= 0x000fffff;
  let i = (high + 0x95f64) & 0x100000;
  x = wordsDouble(high | (i ^ 0x3ff00000), low);
  k += i >> 20;
  const f = x - 1.0;
  if ((0x000fffff & (2 + high)) < 3) {
    if (f === 0.0) { if (k === 0) return 0.0; const dk = k; return dk * LOG_LN2_HI + dk * LOG_LN2_LO; }
    const r = f * f * (0.5 - 0.33333333333333333 * f);
    if (k === 0) return f - r;
    const dk = k;
    return dk * LOG_LN2_HI - ((r - dk * LOG_LN2_LO) - f);
  }
  const s = f / (2.0 + f), dk = k, z = s * s;
  i = high - 0x6147a;
  const w = z * z, j = 0x6b851 - high;
  const t1 = w * (LOG_LG2 + w * (LOG_LG4 + w * LOG_LG6));
  const t2 = z * (LOG_LG1 + w * (LOG_LG3 + w * (LOG_LG5 + w * LOG_LG7)));
  i |= j;
  const r = t2 + t1;
  if (i > 0) {
    const hfsq = 0.5 * f * f;
    if (k === 0) return f - (hfsq - s * (hfsq + r));
    return dk * LOG_LN2_HI - ((hfsq - (s * (hfsq + r) + dk * LOG_LN2_LO)) - f);
  }
  if (k === 0) return f - s * (f - r);
  return dk * LOG_LN2_HI - ((s * (f - r) - dk * LOG_LN2_LO) - f);
}
const EXP_P1 = 1.66666666666666019037e-1, EXP_P2 = -2.77777777770155933842e-3, EXP_P3 = 6.61375632143793436117e-5, EXP_P4 = -1.65339022054652515390e-6, EXP_P5 = 4.13813679705723846039e-8;
const EXP_LN2_HI = 6.93147180369123816490e-1, EXP_LN2_LO = 1.90821492927058770002e-10, EXP_INV_LN2 = 1.44269504088896338700e0;
const EXP_O_THRESHOLD = 7.09782712893383973096e2, EXP_U_THRESHOLD = -7.45133219101941108420e2, EXP_HUGE = 1.0e300, EXP_TWOM1000 = 9.33263618503218878990e-302;
// fdlibm e_exp.c (Sun Microsystems), ported verbatim; mirrored in
// src/biobirddna_genetics/inference.py _portable_exp.
function portableExp(value) {
  let x = value;
  const [highWord, low] = doubleWords(x), xsb = (highWord >>> 31) & 1, hx = highWord & 0x7fffffff;
  if (hx >= 0x40862e42) {
    if (hx >= 0x7ff00000) {
      if (((hx & 0xfffff) | low) !== 0) return x + x;
      return xsb === 1 ? 0.0 : x;
    }
    if (x > EXP_O_THRESHOLD) return EXP_HUGE * EXP_HUGE;
    if (x < EXP_U_THRESHOLD) return EXP_TWOM1000 * EXP_TWOM1000;
  }
  let k = 0, hi = 0.0, lo = 0.0;
  if (hx > 0x3fd62e42) {
    if (hx < 0x3ff0a2e4) {
      hi = x - (xsb === 0 ? EXP_LN2_HI : -EXP_LN2_HI);
      lo = xsb === 0 ? EXP_LN2_LO : -EXP_LN2_LO;
      k = 1 - xsb - xsb;
    } else {
      k = Math.trunc(EXP_INV_LN2 * x + (xsb === 0 ? 0.5 : -0.5));
      const t = k;
      hi = x - t * EXP_LN2_HI;
      lo = t * EXP_LN2_LO;
    }
    x = hi - lo;
  } else if (hx < 0x3e300000) {
    if (EXP_HUGE + x > 1.0) return 1.0 + x;
  }
  const t = x * x;
  const c = x - t * (EXP_P1 + t * (EXP_P2 + t * (EXP_P3 + t * (EXP_P4 + t * EXP_P5))));
  if (k === 0) return 1.0 - ((x * c / (c - 2.0)) - x);
  const y = 1.0 - ((lo - (x * c) / (2.0 - c)) - hi);
  const [highY, lowY] = doubleWords(y);
  if (k >= -1021) return wordsDouble(highY + (k << 20), lowY);
  return wordsDouble(highY + ((k + 1000) << 20), lowY) * EXP_TWOM1000;
}
const LN2 = portableLog(2.0), INV_LN2 = 1.4426950408889634;
// Deterministic (not correctly-rounded) log2: same expression in both
// runtimes, so bit-identical by construction.
function portableLog2(value) { return portableLog(value) * INV_LN2; }
function neumaierSum(values) {
  // Only IEEE-exact operations (+, -, abs, comparison), so the result is
  // bit-identical across runtimes for the same input order. Callers sort
  // summands first to also make the result independent of iteration order.
  let total = 0, compensation = 0;
  for (const value of values) {
    const candidate = total + value;
    if (Math.abs(total) >= Math.abs(value)) compensation += (total - candidate) + value;
    else compensation += (value - candidate) + total;
    total = candidate;
  }
  return total + compensation;
}
function logSumExp(values) {
  if (!values.length) return -Infinity;
  const max = Math.max(...values);
  if (!Number.isFinite(max)) return max;
  // Sorting the exponentiated terms makes the sum independent of the
  // caller's model/configuration ordering (order-independence guarantee).
  const terms = values.map((value) => portableExp(value - max)).sort((a, b) => a - b);
  return max + portableLog(neumaierSum(terms));
}
function parentGenotypes(model, sex, observed) {
  if (model.mode === "non_genetic") return ["NA"];
  const options = model.mode.startsWith("z_linked") && sex === "female" ? [`${model.normal}/W`, `${model.variant}/W`] : [`${model.normal}/${model.normal}`, `${model.normal}/${model.variant}`, `${model.variant}/${model.variant}`];
  if (observed === "unknown" || observed === "UNK") return options;
  return options.filter((genotype) => getOwn(phenotype(model, sex, genotype), observed, 0) > 0);
}
function predictions(model, malePhenotype, femalePhenotype) {
  const out = [];
  for (const male of parentGenotypes(model, "male", malePhenotype)) for (const female of parentGenotypes(model, "female", femalePhenotype)) out.push({ male_genotype: male, female_genotype: female, predicted_probabilities: crossKnown(model, male, female).phenotype_probabilities });
  return out;
}
function scoreModel(model, malePhenotype, femalePhenotype, observations) {
  if (model.mode === "non_genetic") {
    const predicted = crossKnown(model, "NA", "NA").phenotype_probabilities;
    return scoreConfigurations(model, [{ male_genotype: "NA", female_genotype: "NA", predicted_probabilities: predicted }], observations);
  }
  return scoreConfigurations(model, predictions(model, malePhenotype, femalePhenotype), observations);
}
function scoreConfigurations(model, configurations, observations) {
  const observedSexes = new Set(observations.filter((item) => item.count > 0).map((item) => item.sex)), observation_mode = !observedSexes.size || (observedSexes.size === 1 && observedSexes.has("unknown")) ? "sex_collapsed" : observedSexes.has("unknown") ? fail("cannot_mix_known_and_unknown_offspring_sex") : "sex_specific";
  const collapse = (probabilities) => { const result = accumulator(); for (const [key, value] of Object.entries(probabilities)) { const divider = key.indexOf(":"); add(result, key.slice(divider + 1), value); } return result; };
  const counts = new Map(); for (const item of observations) { if (!SEXES.has(item.sex) || !text(item.phenotype, "offspring_phenotype") || !Number.isInteger(item.count) || item.count < 0) fail("offspring_observation_invalid"); const key = observation_mode === "sex_collapsed" ? item.phenotype : `${item.sex}:${item.phenotype}`; counts.set(key, (counts.get(key) ?? 0) + item.count); }
  if (!configurations.length) return { model_id: model.model_id, mode: model.mode, log_evidence: null, observation_mode, observed_counts: Object.fromEntries(counts), configurations: [], best_parent_configuration: null, reason: "parent_phenotype_incompatible_with_model" };
  const prior = 1 / configurations.length, evaluated = configurations.map((configuration) => {
    const raw = observation_mode === "sex_collapsed" ? collapse(configuration.predicted_probabilities) : configuration.predicted_probabilities, categories = [...new Set([...Object.keys(raw), ...counts.keys()])]; const smooth = accumulator();
    for (const category of categories) smooth[category] = (1 - model.observation_error) * getOwn(raw, category, 0) + model.observation_error / categories.length;
    const total = [...counts.values()].reduce((sum, count) => sum + count, 0);
    let logLikelihood = logGamma(total + 1) - [...counts.values()].reduce((sum, count) => sum + logGamma(count + 1), 0);
    for (const [key, count] of counts) {
      const probability = getOwn(smooth, key, 0);
      if (probability <= 0 && count > 0) {
        logLikelihood = -Infinity;
        break;
      }
      logLikelihood += count * portableLog(probability);
    }
    return { ...configuration, prior, log_likelihood: logLikelihood, log_joint: portableLog(prior) + logLikelihood, predicted_probabilities: smooth };
  });
  const log_evidence = logSumExp(evaluated.map((item) => item.log_joint));
  for (const configuration of evaluated) {
    configuration.posterior_within_model = model.mode === "non_genetic" ? 1 : Number.isFinite(log_evidence) ? portableExp(configuration.log_joint - log_evidence) : 0;
    if (!Number.isFinite(configuration.log_likelihood)) configuration.log_likelihood = null;
    if (!Number.isFinite(configuration.log_joint)) configuration.log_joint = null;
    if (model.mode === "non_genetic") delete configuration.log_joint;
  }
  evaluated.sort((left, right) => right.posterior_within_model - left.posterior_within_model);
  return { model_id: model.model_id, mode: model.mode, observation_mode, observed_counts: Object.fromEntries(counts), configurations: evaluated, log_evidence: Number.isFinite(log_evidence) ? log_evidence : null, best_parent_configuration: evaluated[0] };
}
// Lanczos approximation — ported verbatim to src/biobirddna_genetics/
// inference.py _log_gamma so both runtimes share one algorithm; keep in sync.
function logGamma(value) {
  const coefficients = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.984369578019572e-6,
    1.5056327351493116e-7
  ];
  if (value < .5) return portableLog(Math.PI) - portableLog(Math.sin(Math.PI * value)) - logGamma(1 - value);
  const shifted = value - 1;
  let series = coefficients[0];
  for (let index = 1; index < coefficients.length; index += 1) series += coefficients[index] / (shifted + index);
  const scale = shifted + 7.5;
  return .5 * portableLog(2 * Math.PI) + (shifted + .5) * portableLog(scale) - scale + portableLog(series);
}

export function inferCandidateModels(modelsInput, maleParentPhenotype, femaleParentPhenotype, observations, priors = null) {
  const models = modelsInput.map(modelFrom); if (!models.length) fail("candidate_models_required"); if (!observations.length) fail("offspring_observations_required");
  const seenModelIds = new Set(); for (const model of models) { if (seenModelIds.has(model.model_id)) fail(`candidate_model_id_duplicate:${model.model_id}`); seenModelIds.add(model.model_id); }
  const m0Count = models.filter((model) => model.model_id === "M0_NONGENETIC" && model.mode === "non_genetic").length;
  if (m0Count !== 1) fail("candidate_models_require_exactly_one_M0_NONGENETIC");
  if (models.length < 2) fail("candidate_models_require_multiple_models");
  // Z-linked candidates always predict a 0.5 offspring sex ratio (maternal ZW segregation), so mixing
  // them with candidates that assume a different female_probability would bias model comparison
  // silently. Fail closed (mirrors src/biobirddna_genetics/inference.py infer_candidate_models).
  if (models.some((model) => model.mode.startsWith("z_linked")) && models.some((model) => model.female_probability !== 0.5)) fail("female_probability_must_be_0_5_when_z_linked_model_present");
  const rawPriors = Object.fromEntries(models.map((model) => [model.model_id, priors ? finite(Number(getOwn(priors, model.model_id, 0)), "model_prior", 0) : 1])); const total = Object.values(rawPriors).reduce((sum, value) => sum + value, 0); if (!(total > 0)) fail("model_priors_zero_total");
  const ranking = models.map((model) => ({ ...scoreModel(model, maleParentPhenotype, femaleParentPhenotype, observations), prior: rawPriors[model.model_id] / total }));
  const normalizer = logSumExp(ranking.map((item) => item.prior && Number.isFinite(item.log_evidence) ? item.log_evidence + portableLog(item.prior) : -Infinity));
  if (!Number.isFinite(normalizer)) for (const item of ranking) item.posterior = 1 / ranking.length;
  else for (const item of ranking) item.posterior = item.prior && Number.isFinite(item.log_evidence) ? portableExp(item.log_evidence + portableLog(item.prior) - normalizer) : 0;
  ranking.sort((a, b) => b.posterior - a.posterior || a.model_id.localeCompare(b.model_id));
  return { ranking, best_model: ranking[0].model_id, best_model_posterior: ranking[0].posterior, warnings: configurationPriorWarnings(ranking), limits: { claim: "candidate_inheritance_model_only", gene_status: "not_assessed", minimum_language: "phu_hop_hon_voi_mo_hinh" } };
}

const CONFIGURATION_PRIOR_WARNING = "model_ranking_may_reflect_uniform_configuration_prior_not_data";
// Statistical caution, not a biological conclusion. Each configuration of a
// model carries the uniform prior 1/n (n = configurations.length), so
//   log_evidence = logsumexp_c(log_likelihood_c) - log(n).
// For the top-ranked model T and any runner-up R the evidence gap decomposes as
//   log_evidence_T - log_evidence_R
//     = [logsumexp(ll_T) - logsumexp(ll_R)] + [log(n_R) - log(n_T)]
// where the second term comes from configuration counts (model complexity)
// alone. Whenever log_evidence_T + log(n_T) <= log_evidence_R + log(n_R), the
// flat (prior-free) configuration likelihood of R is at least that of T, so
// the ranking gap could be explained by the uniform configuration prior
// rather than the offspring data. Mirrors
// src/biobirddna_genetics/inference.py (_configuration_prior_warnings) —
// keep both sides in sync.
function configurationPriorWarnings(ranking) {
  const top = ranking[0], topConfigurations = top.configurations.length;
  if (!topConfigurations || !Number.isFinite(top.log_evidence)) return [];
  const topFlat = top.log_evidence + portableLog(topConfigurations);
  for (const other of ranking.slice(1)) {
    const otherConfigurations = other.configurations.length;
    if (!otherConfigurations || !Number.isFinite(other.log_evidence)) continue;
    if (topFlat <= other.log_evidence + portableLog(otherConfigurations)) return [CONFIGURATION_PRIOR_WARNING];
  }
  return [];
}

function inferenceRequest(request) {
  request = object(request, "inference_request");
  rejectUnknown(request, INFERENCE_ROOT_FIELDS, "inference_request");
  const target = exactText(request.target_label ?? "target_colour", "target_label");
  const male = exactText(request.male_parent_phenotype, "male_parent_phenotype");
  const female = exactText(request.female_parent_phenotype, "female_parent_phenotype");
  const offspring = array(request.offspring, "offspring").map((raw) => {
    const item = object(raw, "offspring_observation");
    rejectUnknown(item, INFERENCE_OBSERVATION_FIELDS, "offspring_observation");
    const sex = item.sex ?? "unknown";
    if (!SEXES.has(sex)) fail("offspring_sex_invalid");
    return {
      sex,
      phenotype: exactText(item.phenotype, "offspring_phenotype"),
      count: integer(item.count, "offspring_count", 0, MAX_OFFSPRING_OBSERVATION_COUNT)
    };
  });
  const models = request.models === undefined ? defaultModels(target) : array(request.models, "models");
  const priors = request.model_priors === undefined ? null : object(request.model_priors, "model_priors");
  return inferCandidateModels(models, male, female, offspring, priors);
}

function pigment(request) {
  const fields = new Set(["white_area_ratio", "brown_shift_score", "grey_shift_score", "dark_expansion_score", "red_patch_change_score", "pattern_retained", "iris_red_confirmed", "flash_red_eye_risk", "progressive_across_molts", "present_in_first_plumage", "localized_after_injury", "bilateral_symmetry", "image_quality"]); for (const key of Object.keys(request)) if (!fields.has(key)) fail(`pigment_request_unknown_field:${key}`);
  const number = (key) => finite(request[key] === undefined ? 0 : request[key], key, 0, 1), nullableBoolean = (key) => request[key] === undefined || request[key] === null || typeof request[key] === "boolean" ? request[key] ?? null : fail(`${key}_invalid`), boolean = (key) => request[key] === undefined || typeof request[key] === "boolean" ? request[key] ?? false : fail(`${key}_invalid`);
  const white = number("white_area_ratio"), brown = number("brown_shift_score"), grey = number("grey_shift_score"), dark = number("dark_expansion_score"), red = number("red_patch_change_score"), retained = nullableBoolean("pattern_retained"), first = nullableBoolean("present_in_first_plumage"), iris = boolean("iris_red_confirmed"), flash = boolean("flash_red_eye_risk"), progressive = boolean("progressive_across_molts"), injury = boolean("localized_after_injury"), quality = request.image_quality === undefined ? "Q0" : request.image_quality;
  if (!new Set(["Q0", "Q1", "Q2", "Q3", "Q4"]).has(quality)) fail("image_quality_invalid"); if (request.bilateral_symmetry !== undefined && request.bilateral_symmetry !== null) finite(request.bilateral_symmetry, "bilateral_symmetry", 0, 1);
  const labels = ["MEL-LEU", "MEL-PGR", "MEL-ALB", "MEL-INO", "MEL-BRN", "MEL-DIL", "MEL-MEL", "CAR-ALT", "LOC-ACQ", "IMG-ART", "BASE-NORM", "UNK"], limits = { measurement_domain: "visible_phenotype_only", biochemical_pigment_status: "unconfirmed", gene_status: "not_assessed", inheritance_status: "not_assessed", required_language: "kieu_hinh_goi_y" };
  if (quality === "Q0" || quality === "Q1") return { primary_hypothesis: "UNK", alternative_hypotheses: [], scores: Object.fromEntries(labels.map((label) => [label, label === "UNK" ? 1 : 0])), rationale: { UNK: ["image_quality_too_low_for_measurement"] }, needs_recapture: true, limits };
  const scores = {}, rationale = {}; const plus = (label, value, reason) => { scores[label] = (scores[label] ?? 0.05) + value; (rationale[label] ??= []).push(reason); };
  if (flash) plus("IMG-ART", 1.4, "flash_red_eye_risk"); if (injury) plus("LOC-ACQ", 1.8, "localized_after_injury"); if (progressive) plus("MEL-PGR", 2, "white_or_pale_area_progresses_across_molts");
  if (white >= .95) { if (iris && !flash) plus("MEL-ALB", 2.4, "near_complete_white_with_confirmed_red_iris"); else if (retained) plus("MEL-INO", 1.5, "very_pale_with_residual_pattern"); else plus("MEL-LEU", 1, "near_complete_white_without_red_iris_confirmation"); }
  if (white >= .02 && white < .95) { plus("MEL-LEU", .8 + 1.2 * white, "partial_white_plumage"); if (first === false) plus("MEL-PGR", .7, "not_present_in_first_plumage"); }
  if (brown > .15) { plus("MEL-BRN", .5 + 2 * brown, "black_to_brown_visible_shift"); if (retained) plus("MEL-BRN", .4, "pattern_retained"); }
  if (grey > .15) { plus("MEL-DIL", .5 + 2 * grey, "black_to_grey_or_silver_visible_shift"); if (retained) plus("MEL-DIL", .4, "pattern_retained"); }
  if (dark > .15) plus("MEL-MEL", .5 + 2 * dark, "dark_area_expansion"); if (red > .15) plus("CAR-ALT", .5 + 2 * red, "red_or_orange_domain_changed"); if (Math.max(white, brown, grey, dark, red) < .05 && !injury) plus("BASE-NORM", 1.5, "no_strong_abnormal_visible_signal");
  if (quality === "Q2") plus("UNK", .5, "image_quality_supports_visual_screening_only");
  const normalized = Object.fromEntries(labels.map((label) => [label, Math.max(0, scores[label] ?? 0)])), total = Object.values(normalized).reduce((sum, value) => sum + value, 0), ordered = labels.map((label) => [label, total > 0 ? normalized[label] / total : label === "UNK" ? 1 : 0]).sort((left, right) => right[1] - left[1]);
  // Only labels with score > 0 are real candidates: the scoring rules only produce exactly 0 when no rule fired for a label. Strict compare is intentional — genuinely tiny positive scores are kept. Mirrors src/biobirddna_genetics/pigment.py.
  return { primary_hypothesis: ordered[0][0], alternative_hypotheses: ordered.slice(1, 4).filter(([, score]) => score > 0).map(([label]) => label), scores: Object.fromEntries(ordered), rationale, needs_recapture: false, limits };
}

function normalizeRecords(request) { try { return normalizeBreedingRecords(request); } catch (error) { throw adaptBreedingRecordError(error, BioBirdDNAError); } }

function rankEntropy(probabilities) { const terms = Object.values(probabilities).filter((value) => value > 0).map((value) => value * portableLog2(value)).sort((a, b) => a - b); return -neumaierSum(terms); }
function rankDistribution(value) { return distribution(value); }
function rankUncertainty(value) { const normalized = rankDistribution(value); return Object.keys(normalized).length <= 1 ? 0 : rankEntropy(normalized) / portableLog2(Object.keys(normalized).length); }
function rankProfile(value, expectedSex) {
  value = object(value, "genetic_profile");
  rejectUnknown(value, PROFILE_FIELDS, "genetic_profile");
  const individual_id = exactText(value.individual_id, "individual_id"), sex = value.sex;
  if (!SEXES.has(sex)) fail("profile_sex_invalid");
  if (sex !== expectedSex) fail(`profile_not_${expectedSex}:${individual_id}`);
  const rawPosteriors = object(value.genotype_posteriors, "genotype_posteriors"), genotype_posteriors = accumulator();
  for (const [modelId, values] of Object.entries(rawPosteriors)) genotype_posteriors[exactText(modelId, "model_id")] = distribution(values);
  if (Object.hasOwn(value, "phenotype_label") && typeof value.phenotype_label !== "string") fail("phenotype_label_invalid");
  if (Object.hasOwn(value, "notes") && (!Array.isArray(value.notes) || value.notes.some((item) => typeof item !== "string"))) fail("notes_invalid");
  const welfare_blocked = value.welfare_blocked === undefined ? false : typeof value.welfare_blocked === "boolean" ? value.welfare_blocked : fail("welfare_blocked_invalid");
  return { individual_id, sex, genotype_posteriors, health_risk: finite(value.health_risk ?? 0, "health_risk", 0, 1), fertility_risk: finite(value.fertility_risk ?? 0, "fertility_risk", 0, 1), welfare_blocked };
}
function pedigreeState(entries) {
  entries = array(entries, "pedigree"); if (!entries.length) fail("pedigree_individuals_required"); const people = new Map();
  for (const raw of entries) {
    const value = object(raw, "pedigree_individual");
    rejectUnknown(value, PEDIGREE_FIELDS, "pedigree_individual");
    const id = exactText(value.individual_id, "individual_id"), sex = value.sex ?? "unknown";
    if (!SEXES.has(sex)) fail("pedigree_sex_invalid");
    const parent = (name) => value[name] === undefined || value[name] === null ? null : exactText(value[name], name);
    const sire_id = parent("sire_id"), dam_id = parent("dam_id");
    if (people.has(id)) fail("pedigree_duplicate_individual_id");
    if (sire_id === id || dam_id === id) fail("individual_cannot_be_own_parent");
    people.set(id, { id, sire_id, dam_id });
  }
  const indegree = new Map([...people.keys()].map((id) => [id, 0])), children = new Map([...people.keys()].map((id) => [id, []]));
  for (const person of people.values()) for (const parent of [person.sire_id, person.dam_id]) if (parent && people.has(parent)) { indegree.set(person.id, indegree.get(person.id) + 1); children.get(parent).push(person.id); }
  const queue = [...people.keys()].filter((id) => indegree.get(id) === 0), order = [];
  while (queue.length) { const id = queue.shift(); order.push(id); for (const child of children.get(id)) { indegree.set(child, indegree.get(child) - 1); if (indegree.get(child) === 0) queue.push(child); } }
  if (order.length !== people.size) fail("pedigree_cycle_detected"); const index = new Map(order.map((id, i) => [id, i])), matrix = Array.from({ length: order.length }, () => Array(order.length).fill(0));
  for (let i = 0; i < order.length; i += 1) { const person = people.get(order[i]), sire = person.sire_id && index.has(person.sire_id) ? index.get(person.sire_id) : null, dam = person.dam_id && index.has(person.dam_id) ? index.get(person.dam_id) : null; for (let j = 0; j < i; j += 1) matrix[i][j] = matrix[j][i] = .5 * ((sire === null ? 0 : matrix[sire][j]) + (dam === null ? 0 : matrix[dam][j])); matrix[i][i] = sire !== null && dam !== null ? 1 + .5 * matrix[sire][dam] : 1; }
  const require = (id) => { if (!index.has(id)) fail(`pedigree_individual_unknown:${id}`); return index.get(id); };
  const ancestors = (id) => { require(id); const result = new Set(), stack = [id]; while (stack.length) { const current = people.get(stack.pop()); for (const parent of [current.sire_id, current.dam_id]) if (parent && people.has(parent) && !result.has(parent)) { result.add(parent); stack.push(parent); } } return result; };
  return { offspringInbreeding: (male, female) => .5 * matrix[require(male)][require(female)], hardRelation: (male, female) => { if (ancestors(female).has(male) || ancestors(male).has(female)) return "direct_ancestor_pair"; const first = people.get(male), second = people.get(female); return male !== female && first?.sire_id && first?.dam_id && first.sire_id === second?.sire_id && first.dam_id === second?.dam_id ? "full_sibling_pair" : null; } };
}
function rankPredictions(male, female, hypotheses) {
  const predictions = []; let target = 0, uncertainty = 0;
  for (const hypothesis of hypotheses) { const mg = getOwn(male.genotype_posteriors, hypothesis.model.model_id, null), fg = getOwn(female.genotype_posteriors, hypothesis.model.model_id, null); if (!mg || !fg) continue; const probabilities = crossUncertain(hypothesis.model, mg, fg).phenotype_probabilities; predictions.push([hypothesis.weight, probabilities]); target += hypothesis.weight * Object.entries(probabilities).filter(([key]) => key.endsWith(`:${hypothesis.model.target_label}`)).reduce((sum, [, value]) => sum + value, 0); uncertainty += hypothesis.weight * .5 * (rankUncertainty(mg) + rankUncertainty(fg)); }
  const total = predictions.reduce((sum, [weight]) => sum + weight, 0); return total > 0 ? { predictions: predictions.map(([weight, values]) => [weight / total, values]), target: target / total, uncertainty: uncertainty / total } : { predictions: [], target: 0, uncertainty: 1 };
}
function rankInformation(predictions) { if (!predictions.length) return 0; const categories = [...new Set(predictions.flatMap(([, probabilities]) => Object.keys(probabilities)))].sort(), mixture = accumulator(), conditional = predictions.reduce((sum, [weight, probabilities]) => { const total = categories.reduce((sub, key) => sub + getOwn(probabilities, key, 0), 0); if (!(total > 0)) return sum; const normalized = Object.fromEntries(categories.map((key) => [key, getOwn(probabilities, key, 0) / total])); for (const [key, probability] of Object.entries(normalized)) add(mixture, key, weight * probability); return sum + weight * rankEntropy(normalized); }, 0); return Math.max(0, rankEntropy(mixture) - conditional); }
function rankMates(request) {
  request = object(request, "rank_mates_request");
  rejectUnknown(request, RANK_ROOT_FIELDS, "rank_mates_request");
  const weightsInput = request.weights === undefined ? {} : object(request.weights, "weights");
  rejectUnknown(weightsInput, PAIRING_WEIGHT_FIELDS, "weights");
  const weights = { target: finite(weightsInput.target ?? 1, "target_weight", 0), information: finite(weightsInput.information ?? .7, "information_weight", 0), inbreeding: finite(weightsInput.inbreeding ?? 2, "inbreeding_weight", 0), health: finite(weightsInput.health ?? 1.5, "health_weight", 0), uncertainty: finite(weightsInput.uncertainty ?? .4, "uncertainty_weight", 0) };
  const maximum = request.max_offspring_inbreeding === undefined ? .0625 : finite(request.max_offspring_inbreeding, "max_offspring_inbreeding", 0, 1);
  const hardBlock = request.block_direct_ancestor_and_full_sibling === undefined ? true : typeof request.block_direct_ancestor_and_full_sibling === "boolean" ? request.block_direct_ancestor_and_full_sibling : fail("block_direct_ancestor_and_full_sibling_invalid");
  const pedigree = pedigreeState(request.pedigree);
  const hypothesesRaw = array(request.hypotheses, "hypotheses"); if (!hypothesesRaw.length) fail("model_hypotheses_required"); const hypotheses = hypothesesRaw.map((raw) => { const item = object(raw, "model_hypothesis"); rejectUnknown(item, HYPOTHESIS_FIELDS, "model_hypothesis"); return { model: modelFrom(item.model), weight: finite(item.weight, "model_weight", 0) }; }), totalWeight = hypotheses.reduce((sum, item) => sum + item.weight, 0); if (!(totalWeight > 0)) fail("model_hypothesis_weights_zero"); for (const hypothesis of hypotheses) hypothesis.weight /= totalWeight;
  const males = array(request.males, "males").map((item) => rankProfile(item, "male")), females = array(request.females, "females").map((item) => rankProfile(item, "female")), pair_ranking = [];
  for (const male of males) for (const female of females) { const block_reasons = []; if (male.welfare_blocked || female.welfare_blocked) block_reasons.push("welfare_blocked"); if (male.individual_id === female.individual_id) block_reasons.push("same_individual_pair"); const relation = pedigree.hardRelation(male.individual_id, female.individual_id); if (hardBlock && relation) block_reasons.push(relation); const offspring_inbreeding = pedigree.offspringInbreeding(male.individual_id, female.individual_id); if (offspring_inbreeding > maximum) block_reasons.push("offspring_inbreeding_above_configured_limit"); const prediction = rankPredictions(male, female, hypotheses); if (!prediction.predictions.length) block_reasons.push("missing_genotype_posteriors"); const expected_information_gain_bits = rankInformation(prediction.predictions), health_risk = Math.max(male.health_risk, female.health_risk, male.fertility_risk, female.fertility_risk), score = weights.target * prediction.target + weights.information * expected_information_gain_bits - weights.inbreeding * offspring_inbreeding - weights.health * health_risk - weights.uncertainty * prediction.uncertainty; if (!Number.isFinite(score)) fail(`pair_score_non_finite:${male.individual_id}:${female.individual_id}`); pair_ranking.push({ male_id: male.individual_id, female_id: female.individual_id, eligible: !block_reasons.length, block_reasons, score, target_probability: prediction.target, expected_information_gain_bits, offspring_inbreeding, health_risk, genotype_uncertainty: prediction.uncertainty }); }
  pair_ranking.sort((left, right) => Number(right.eligible) - Number(left.eligible) || right.score - left.score); return { pair_ranking, configuration: { max_offspring_inbreeding: maximum, weights }, limits: { purpose: "decision_support_not_automatic_breeding_order", health_and_welfare_review_required: true, gene_status: "model_based" } };
}

function seeded(base, ...parts) { const mask = (1n << 53n) - 1n, prefix = createHash("sha256").update([base, ...parts].join("|")).digest("hex").slice(0, 16); return Number(BigInt(`0x${prefix}`) & mask); }
function rng(seed) { const modulus = 1n << 53n, mask = modulus - 1n; let state = BigInt(seed) & mask; return () => { state = (state * 6364136223846793005n + 1442695040888963407n) & mask; return Number(state) / Number(modulus); }; }
const UTF8_ENCODER = new TextEncoder();
function compareUtf8(left, right) { const leftBytes = UTF8_ENCODER.encode(left), rightBytes = UTF8_ENCODER.encode(right), length = Math.min(leftBytes.length, rightBytes.length); for (let index = 0; index < length; index += 1) if (leftBytes[index] !== rightBytes[index]) return leftBytes[index] - rightBytes[index]; return leftBytes.length - rightBytes.length; }
function sample(probabilities, total, random) { const entries = Object.entries(probabilities).sort(([a], [b]) => compareUtf8(a, b)), totalWeight = entries.reduce((sum, [, probability]) => sum + Number(probability), 0); if (!(totalWeight > 0) || entries.some(([, probability]) => !Number.isFinite(Number(probability)) || Number(probability) < 0)) fail("sampling_probabilities_invalid"); const normalized = entries.map(([key, probability]) => [key, Number(probability) / totalWeight]), counts = Object.fromEntries(normalized.map(([key]) => [key, 0])); for (let i = 0; i < total; i += 1) { const draw = random(); let sum = 0; for (let index = 0; index < normalized.length; index += 1) { const [key, probability] = normalized[index]; sum += probability; if (draw < sum || index === normalized.length - 1) { counts[key] += 1; break; } } } return counts; }
export function __sampleForTest(probabilities, draws) { let index = 0; return sample(probabilities, draws.length, () => draws[index++]); }
// Test-only hook for the cross-runtime numeric primitive parity sweep
// (tests/test_npm_parity.py + npm/test/numeric-parity.test.mjs).
export const __numericPrimitivesForTest = { portableLog, portableExp, portableLog2, neumaierSum, logSumExp, logGamma, entropy, quantile, orderIndependentMean };
function observations(counts) { return Object.entries(counts).filter(([, count]) => count > 0).map(([key, count]) => { const divider = key.indexOf(":"); return { sex: key.slice(0, divider), phenotype: key.slice(divider + 1), count }; }); }
const SIMULATION_ROOT_FIELDS = new Set(["schema_version", "base_seed", "goodness_of_fit_replicates", "misspecification_alpha", "experiments"]);
const SIMULATION_EXPERIMENT_FIELDS = new Set(["experiment_id", "generating_model", "male_genotype", "female_genotype", "male_parent_phenotype", "female_parent_phenotype", "brood_sizes", "replicates", "candidate_model_ids"]);
const SIMULATION_MODEL_FIELDS = new Set(["model_id", "mode", "normal_allele", "variant_allele", "normal_label", "target_label", "intermediate_label", "penetrance", "observation_error", "female_probability", "baseline_target_probability", "description"]);
function rejectUnknown(value, allowed, prefix) { for (const key of Object.keys(object(value, prefix))) if (!allowed.has(key)) fail(`${prefix}_unknown_field:${key}`); }
function rejectNonFiniteNumbers(value, seen = new WeakSet()) {
  if (typeof value === "number" && !Number.isFinite(value)) fail("input_json_non_finite:number");
  if (!value || typeof value !== "object") return;
  if (seen.has(value)) fail("input_json_cycle_forbidden");
  seen.add(value);
  for (const child of Array.isArray(value) ? value : Object.values(value)) rejectNonFiniteNumbers(child, seen);
  seen.delete(value);
}
function modelSnapshot(model) { return { mode: model.mode, normal_allele: model.normal, variant_allele: model.variant, normal_label: model.normal_label, target_label: model.target_label, intermediate_label: model.intermediate_label, penetrance: model.penetrance, observation_error: model.observation_error, female_probability: model.female_probability, baseline_target_probability: model.baseline_target_probability }; }
function sameSnapshot(left, right) { return Object.entries(modelSnapshot(left)).every(([key, value]) => value === modelSnapshot(right)[key]); }
// Same interpolation expression as src/biobirddna_genetics/simulation.py
// _quantile (a * (1 - w) + b * w) so both runtimes agree bit-for-bit.
function quantile(values, probability) { if (!values.length) return null; const ordered = [...values].sort((a, b) => a - b), position = (ordered.length - 1) * probability, low = Math.floor(position), high = Math.ceil(position); if (low === high) return ordered[low]; const weight = position - low; return ordered[low] * (1 - weight) + ordered[high] * weight; }
// Mirrors src/biobirddna_genetics/simulation.py _entropy_bits: Math.log2
// (not Math.log(v)/Math.log(2), which differs in the last bit for ~1/3 of
// inputs) plus sorted Neumaier accumulation — keep both sides in sync.
function entropy(posterior) { const terms = Object.values(posterior).filter((value) => value > 0).map((value) => value * portableLog2(value)).sort((a, b) => a - b); return -neumaierSum(terms); }
// Mirrors src/biobirddna_genetics/simulation.py _order_independent_mean.
function orderIndependentMean(values) { return neumaierSum([...values].sort((a, b) => a - b)) / values.length; }
function wilson(successes, total) { const z = 1.959963984540054, p = successes / total, denominator = 1 + z ** 2 / total, center = (p + z ** 2 / (2 * total)) / denominator, radius = z * Math.sqrt(p * (1 - p) / total + z ** 2 / (4 * total ** 2)) / denominator; return { low: Math.max(0, center - radius), high: Math.min(1, center + radius) }; }
function modelFrequencies(results) { const top = accumulator(), unique = accumulator(); for (const result of results) { for (const id of result.best_models) add(top, id, 1); if (result.best_models.length === 1) add(unique, result.best_models[0], 1); } return Object.fromEntries(Object.keys(top).sort().map((id) => [id, { top_set_count: top[id], top_set_rate: top[id] / results.length, unique_top1_count: getOwn(unique, id, 0), unique_top1_rate: getOwn(unique, id, 0) / results.length }])); }
function modelSetFrequencies(results) { const counts = accumulator(); for (const result of results) { const key = result.best_models.join("|"); add(counts, key, 1); } return Object.fromEntries(Object.keys(counts).sort().map((key) => [key, { count: counts[key], rate: counts[key] / results.length }])); }
function discriminationBins(truePosteriors, recovered) { const bins = Array.from({ length: 5 }, () => []); for (let index = 0; index < truePosteriors.length; index += 1) bins[Math.min(4, Math.floor(truePosteriors[index] * 5))].push([truePosteriors[index], recovered[index]]); return bins.flatMap((items, index) => items.length ? [{ lower: index / 5, upper: (index + 1) / 5, count: items.length, mean_predicted: orderIndependentMean(items.map((item) => item[0])), observed_top1_rate: items.filter((item) => item[1]).length / items.length }] : []); }
function smoothDistribution(probabilities, error) { const keys = Object.keys(probabilities); return Object.fromEntries(keys.map((key) => [key, (1 - error) * probabilities[key] + error / keys.length])); }
export function __devianceForTest(counts, probabilities) { const total = Object.values(counts).reduce((sum, value) => sum + value, 0); let statistic = 0; for (const [key, count] of Object.entries(counts)) { const expected = total * getOwn(probabilities, key, 0); if (count > 0 && expected <= 0) return Infinity; if (count > 0) statistic += 2 * count * portableLog(count / expected); } return Math.max(0, statistic); }
function modelDeviance(model, male, female, counts) {
  const score = scoreModel(model, male, female, observations(counts));
  const configuration = score.best_parent_configuration;
  if (!configuration) return { statistic: Infinity, predicted: null };
  const predicted = configuration.predicted_probabilities, total = Object.values(counts).reduce((sum, value) => sum + value, 0);
  return { statistic: __devianceForTest(counts, predicted), predicted };
}
function goodnessOfFit(models, male, female, counts, bootstrap, alpha, baseSeed, experimentId, broodSize) {
  const candidates = models.map((model) => {
    const observed = modelDeviance(model, male, female, counts);
    if (!observed.predicted) return { model_id: model.model_id, deviance: null, bootstrap_p_value: 0, status: "parent_phenotype_incompatible" };
    let extreme = 0;
    const random = rng(seeded(baseSeed, experimentId, broodSize, model.model_id, "goodness_of_fit"));
    for (let index = 0; index < bootstrap; index += 1) {
      const bootstrapCounts = sample(observed.predicted, Object.values(counts).reduce((sum, value) => sum + value, 0), random);
      if (modelDeviance(model, male, female, bootstrapCounts).statistic >= observed.statistic) extreme += 1;
    }
    return { model_id: model.model_id, deviance: Number.isFinite(observed.statistic) ? observed.statistic : null, bootstrap_p_value: (extreme + 1) / (bootstrap + 1), status: "evaluated" };
  });
  return { method: "parametric_bootstrap_multinomial_deviance", bootstrap_replicates: bootstrap, alpha, candidates, all_candidates_below_alpha: candidates.every((candidate) => candidate.bootstrap_p_value < alpha) };
}
function simulationConfig(request) {
  rejectUnknown(request, SIMULATION_ROOT_FIELDS, "simulation_request");
  if (request.schema_version !== "1.0.0") fail("simulation_schema_version_unsupported");
  const base_seed = integer(request.base_seed, "simulation_base_seed", 0, 2 ** 53 - 1);
  const alpha = finite(request.misspecification_alpha, "misspecification_alpha", Number.EPSILON, 0.25);
  const bootstrap = integer(request.goodness_of_fit_replicates, "goodness_of_fit_replicates", 10, 10000);
  if (!(bootstrap + 1 > 1 / alpha)) fail(`goodness_of_fit_replicates_insufficient_for_alpha:${bootstrap}:${alpha}`);
  const ids = new Set(), experiments = array(request.experiments, "simulation_experiments").map((raw) => {
    rejectUnknown(raw, SIMULATION_EXPERIMENT_FIELDS, "simulation_experiment");
    const experiment_id = exactText(raw.experiment_id, "simulation_experiment_id"); if (ids.has(experiment_id)) fail(`simulation_experiment_id_duplicate:${experiment_id}`); ids.add(experiment_id);
    rejectUnknown(raw.generating_model, SIMULATION_MODEL_FIELDS, "simulation_generating_model");
    if (!("model_id" in raw.generating_model)) fail("simulation_generating_model_model_id_required");
    if ("description" in raw.generating_model && (typeof raw.generating_model.description !== "string" || raw.generating_model.description.trim() !== raw.generating_model.description)) fail("simulation_generating_model_description_invalid");
    for (const field of ["model_id", "mode", "normal_allele", "variant_allele", "normal_label", "target_label", "intermediate_label"]) if (field in raw.generating_model) exactText(raw.generating_model[field], `simulation_generating_model_${field}`);
    const generating = modelFrom(raw.generating_model), candidateIds = array(raw.candidate_model_ids, "simulation_candidate_model_ids").map((id) => exactText(id, "simulation_candidate_model_id"));
    if (new Set(candidateIds).size !== candidateIds.length || candidateIds.filter((id) => id === "M0_NONGENETIC").length !== 1) fail(`simulation_candidate_set_requires_exactly_one_M0_NONGENETIC:${experiment_id}`);
    if (candidateIds.length < 2) fail(`simulation_candidate_set_requires_multiple_models:${experiment_id}`);
    if (candidateIds.some((id) => !DEFAULT_IDS.includes(id))) fail(`simulation_candidate_model_unknown:${experiment_id}`);
    const canonicalModels = defaultModels(generating.target_label, generating.normal_label, generating.intermediate_label), canonicalById = new Map(canonicalModels.map((model) => [model.model_id, model])), candidates = candidateIds.map((id) => canonicalById.get(id));
    const knownId = canonicalModels.find((model) => model.model_id === generating.model_id); if (knownId && !sameSnapshot(generating, knownId)) fail(`simulation_generating_model_snapshot_mismatch:${experiment_id}:${generating.model_id}`);
    const male_genotype = exactText(raw.male_genotype, "simulation_male_genotype"), female_genotype = exactText(raw.female_genotype, "simulation_female_genotype"), male_parent_phenotype = exactText(raw.male_parent_phenotype, "simulation_male_parent_phenotype"), female_parent_phenotype = exactText(raw.female_parent_phenotype, "simulation_female_parent_phenotype");
    if (generating.mode === "non_genetic") {
      if (male_genotype !== "NA") fail(`simulation_non_genetic_genotype_must_be_NA:${experiment_id}:male`);
      if (female_genotype !== "NA") fail(`simulation_non_genetic_genotype_must_be_NA:${experiment_id}:female`);
    }
    for (const [sex, genotype, observed] of [["male", male_genotype, male_parent_phenotype], ["female", female_genotype, female_parent_phenotype]]) {
      const parentDistribution = phenotype(generating, sex, genotype);
      if (observed !== "unknown" && observed !== "UNK" && !(getOwn(parentDistribution, observed, 0) > 0)) fail(`simulation_parent_phenotype_genotype_incompatible:${experiment_id}:${sex}`);
    }
    const brood_sizes = array(raw.brood_sizes, "simulation_brood_sizes").map((value) => integer(value, "simulation_brood_size", 1, 100000)); if (!brood_sizes.length || new Set(brood_sizes).size !== brood_sizes.length) fail(`simulation_brood_sizes_invalid:${experiment_id}`); brood_sizes.sort((a, b) => a - b);
    return { experiment_id, generating, candidateIds, candidates, male_genotype, female_genotype, male_parent_phenotype, female_parent_phenotype, brood_sizes, replicates: integer(raw.replicates, "simulation_replicates", 1, 10000) };
  });
  if (!experiments.length) fail("simulation_experiments_required");
  const runs = experiments.reduce((sum, item) => sum + item.brood_sizes.length * item.replicates, 0), offspring = experiments.reduce((sum, item) => sum + item.brood_sizes.reduce((n, size) => n + size, 0) * item.replicates, 0), draws = experiments.reduce((sum, item) => sum + item.brood_sizes.reduce((n, size) => n + size, 0) * item.replicates * bootstrap * item.candidates.length, 0), refits = experiments.reduce((sum, item) => sum + item.brood_sizes.length * bootstrap * item.candidates.length, 0);
  if (runs > 10000 || offspring > 5000000) fail("simulation_work_budget_exceeded"); if (draws > 5000000 || refits > 250000) fail("simulation_bootstrap_work_budget_exceeded");
  return { base_seed, alpha, bootstrap, experiments };
}
function calibrationByDesign(experiments) {
  const groups = new Map();
  for (const experiment of experiments) for (const block of experiment.by_brood_size) { const key = `${block.brood_size}|${experiment.candidate_model_ids.join(",")}`, group = groups.get(key) ?? { brood_size: block.brood_size, candidate_model_ids: experiment.candidate_model_ids, results: [] }; group.results.push(...block.replicate_results.map((result) => ({ ...result, true_model: experiment.generating_model_id }))); groups.set(key, group); }
  return [...groups.values()].map((group) => {
    const counts = Object.fromEntries([...new Set(group.results.map((result) => result.true_model))].sort().map((id) => [id, group.results.filter((result) => result.true_model === id).length])), trueModels = Object.keys(counts);
    const reason = trueModels.length < 2 ? "requires_at_least_two_generating_models" : trueModels.some((id) => !group.candidate_model_ids.includes(id)) ? "generating_model_outside_candidate_set" : null;
    const models = reason ? [] : group.candidate_model_ids.map((model_id) => ({ model_id, bins: discriminationBins(group.results.map((result) => result.model_posteriors[model_id]), group.results.map((result) => result.true_model === model_id)).map(({ mean_predicted, observed_top1_rate, ...bin }) => ({ ...bin, mean_posterior: mean_predicted, observed_true_model_rate: observed_top1_rate })) }));
    return { brood_size: group.brood_size, candidate_model_ids: group.candidate_model_ids, method: "empirical_reliability_under_declared_simulation_design", status: reason ? "not_estimable" : "estimated", reason, design_model_frequencies: counts, models, limits: { design_dependent: true, not_external_field_calibration: true } };
  });
}
function simulateRecovery(request) {
  const config = simulationConfig(request), output = [];
  for (const experiment of config.experiments) {
    const probabilities = smoothDistribution(crossKnown(experiment.generating, experiment.male_genotype, experiment.female_genotype).phenotype_probabilities, experiment.generating.observation_error), by_brood_size = [];
    for (const brood_size of experiment.brood_sizes) {
      const results = [];
      for (let index = 0; index < experiment.replicates; index += 1) {
        const seed = seeded(config.base_seed, experiment.experiment_id, brood_size, index), counts = sample(probabilities, brood_size, rng(seed)), inference = inferCandidateModels(experiment.candidates, experiment.male_parent_phenotype, experiment.female_parent_phenotype, observations(counts)), posterior = Object.fromEntries(inference.ranking.map((item) => [item.model_id, item.posterior]),), maximum = Math.max(...Object.values(posterior)), best_models = Object.keys(posterior).filter((id) => Math.abs(posterior[id] - maximum) <= 1e-12).sort();
        const logEvidence = inference.ranking.map((item) => item.log_evidence).filter(Number.isFinite), top_two_log2_bayes_factor = logEvidence.length >= 2 ? (logEvidence[0] - logEvidence[1]) / LN2 : null, multiclass_brier = experiment.candidateIds.includes(experiment.generating.model_id) ? neumaierSum(Object.entries(posterior).map(([id, value]) => (value - (id === experiment.generating.model_id ? 1 : 0)) ** 2).sort((a, b) => a - b)) : null;
        results.push({ replicate_index: index, seed, observed_counts: counts, best_model: best_models.length === 1 ? best_models[0] : null, best_models, top1_tied: best_models.length > 1, best_model_posterior: maximum, true_model_posterior: getOwn(posterior, experiment.generating.model_id, null), model_posteriors: sortedObject(posterior), model_entropy_bits: entropy(posterior), top_two_log2_bayes_factor, multiclass_brier });
      }
      const trueInSet = experiment.candidateIds.includes(experiment.generating.model_id), uniqueHits = results.filter((result) => result.best_model === experiment.generating.model_id).length, setHits = results.filter((result) => result.best_models.includes(experiment.generating.model_id)).length, truePosteriors = results.map((result) => result.true_model_posterior).filter((value) => value !== null), pooled = Object.fromEntries(Object.keys(probabilities).map((key) => [key, results.reduce((sum, result) => sum + result.observed_counts[key], 0)])), gof = goodnessOfFit(experiment.candidates, experiment.male_parent_phenotype, experiment.female_parent_phenotype, pooled, config.bootstrap, config.alpha, config.base_seed, experiment.experiment_id, brood_size), misspecified = gof.all_candidates_below_alpha;
      by_brood_size.push({ brood_size, replicates: results.length, small_sample_warning: brood_size < 30, warnings: [...(brood_size < 30 ? ["small_sample_high_uncertainty"] : []), ...(misspecified ? ["candidate_model_set_misspecified_suspected"] : [])], recovery: { true_model_in_candidate_set: trueInSet, top1_rate: trueInSet ? uniqueHits / results.length : null, top_set_contains_true_rate: trueInSet ? setHits / results.length : null, tie_policy: "top1_rate_requires_unique_best_model", wilson_95: trueInSet ? wilson(uniqueHits, results.length) : null, mean_true_model_posterior: truePosteriors.length ? orderIndependentMean(truePosteriors) : null, multiclass_brier_mean: trueInSet ? orderIndependentMean(results.map((result) => result.multiclass_brier)) : null }, uncertainty: { true_model_posterior_quantiles: { p05: quantile(truePosteriors, .05), p50: quantile(truePosteriors, .5), p95: quantile(truePosteriors, .95) }, mean_model_entropy_bits: orderIndependentMean(results.map((result) => result.model_entropy_bits)), top_two_log2_bayes_factor_quantiles: { p05: quantile(results.map((result) => result.top_two_log2_bayes_factor).filter((value) => value !== null), .05), p50: quantile(results.map((result) => result.top_two_log2_bayes_factor).filter((value) => value !== null), .5), p95: quantile(results.map((result) => result.top_two_log2_bayes_factor).filter((value) => value !== null), .95) } }, best_model_frequencies: modelFrequencies(results), best_model_set_frequencies: modelSetFrequencies(results), posterior_discrimination_bins: trueInSet ? discriminationBins(truePosteriors, results.map((result) => result.best_models.length === 1 && result.best_models[0] === experiment.generating.model_id)) : [], goodness_of_fit: gof, replicate_results: results });
    }
    output.push({ experiment_id: experiment.experiment_id, generating_model_id: experiment.generating.model_id, candidate_model_ids: experiment.candidateIds, by_brood_size });
  }
  return { schema_version: "1.0.0", base_seed: config.base_seed, experiments: output, calibration_by_design: calibrationByDesign(output), evidence: { source: "synthetic_simulation_only", evidence_grade: "G0", gene_status: "not_assessed" }, limits: { synthetic_only: true, decision_support_only: true, posterior_scope: "relative_to_candidate_model_set", gene_status: "not_assessed", inheritance_status: "simulation_recovery_not_confirmation", genotype_specific_survival: "not_modeled" } };
}

// Evidence boundary — mirrors src/biobirddna_genetics/evidence.py and the
// operation claim ceilings in src/biobirddna_genetics/service.py; keep all
// three maps in sync (references/scientific-boundaries.md sections 1 and 4).
const EVIDENCE_ORDER = { G0: 0, G1: 1, G2: 2, G3: 3, G4: 4, G5: 5, G6: 6 };
const CLAIM_MINIMUM_GRADE = {
  observation_record: "G0",
  phenotype_hypothesis: "G0",
  simulate_cross: "G0",
  candidate_inheritance_model: "G1",
  supported_inheritance_pattern: "G3",
  strong_inheritance_model: "G4",
  molecular_locus_supported: "G5",
  mechanism_strongly_confirmed: "G6",
};
const OPERATION_CLAIM_CEILINGS = {
  cross: ["simulate_cross", "G0"],
  infer: ["candidate_inheritance_model", "G2"],
  pigment: ["phenotype_hypothesis", "G0"],
  "rank-mates": ["simulate_cross", "G0"],
  "normalize-records": ["observation_record", "G0"],
  "simulate-recovery": ["simulate_cross", "G0"],
};
const claimAllowed = (claim, grade) => Object.hasOwn(CLAIM_MINIMUM_GRADE, claim) && EVIDENCE_ORDER[grade] >= EVIDENCE_ORDER[CLAIM_MINIMUM_GRADE[claim]];
function safeClaimLanguage(claim, grade) {
  if (claimAllowed(claim, grade)) return claim;
  if (EVIDENCE_ORDER[grade] <= EVIDENCE_ORDER.G2) return "exploratory_candidate_model";
  if (EVIDENCE_ORDER[grade] === EVIDENCE_ORDER.G3) return "segregation_pattern_supported";
  if (EVIDENCE_ORDER[grade] === EVIDENCE_ORDER.G4) return "strong_pedigree_model_without_molecular_confirmation";
  return "insufficient_evidence";
}
function applyClaimBoundary(operation, result) {
  if (!Object.hasOwn(OPERATION_CLAIM_CEILINGS, operation)) fail(`operation_claim_ceiling_missing:${operation}`);
  const [claim, grade] = OPERATION_CLAIM_CEILINGS[operation];
  const language = safeClaimLanguage(claim, grade);
  if (!claimAllowed(language, grade)) fail(`claim_not_allowed:${language}:requires_${getOwn(CLAIM_MINIMUM_GRADE, language, "unsupported")}:has_${grade}`);
  if (!isObject(result.limits)) fail(`result_limits_missing:${operation}`);
  result.limits.claim_language = language;
  return result;
}

export function dispatch(operation, request) {
  return applyClaimBoundary(operation, dispatchOperation(operation, request));
}

function dispatchOperation(operation, request) {
  rejectNonFiniteNumbers(request);
  switch (operation) {
    case "cross": {
      request = object(request, "cross_request");
      rejectUnknown(request, CROSS_ROOT_FIELDS, "cross_request");
      const knownMode = request.male_genotype !== undefined && request.male_genotype !== null || request.female_genotype !== undefined && request.female_genotype !== null;
      const uncertainMode = request.male_genotypes !== undefined && request.male_genotypes !== null || request.female_genotypes !== undefined && request.female_genotypes !== null;
      if (knownMode && uncertainMode) fail("cross_parent_input_modes_mutually_exclusive");
      if (uncertainMode) {
        if (request.male_genotypes === undefined || request.male_genotypes === null || request.female_genotypes === undefined || request.female_genotypes === null) fail("both_parent_genotype_distributions_required");
        return crossUncertain(request.model, request.male_genotypes, request.female_genotypes);
      }
      if (!knownMode || request.male_genotype === undefined || request.male_genotype === null || request.female_genotype === undefined || request.female_genotype === null) fail("male_and_female_genotype_required");
      return crossKnown(request.model, request.male_genotype, request.female_genotype);
    }
    case "infer": return inferenceRequest(request);
    case "pigment": return pigment(object(request, "pigment_request"));
    case "rank-mates": return rankMates(object(request, "mate_ranking_request"));
    case "normalize-records": return normalizeRecords(request);
    case "simulate-recovery": return simulateRecovery(object(request, "simulation_request"));
    default: fail(`operation_unknown:${operation}`);
  }
}
