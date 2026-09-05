-- BioBirdDNA Genetics — PostgreSQL schema draft
-- Research architecture only. Review consent, privacy, legal source, and RLS before use.

create extension if not exists pgcrypto;

create table if not exists cases (
  id uuid primary key default gen_random_uuid(),
  case_code text not null unique,
  species_name text not null default 'Pycnonotus jocosus',
  status text not null check (status in (
    'received','identity_validated','phenotype_reviewed','genetics_data_pending',
    'genetics_inference_ready','human_review_required','approved_for_research_cross',
    'needs_more_data','blocked','archived'
  )),
  consent_scope jsonb not null default '{}'::jsonb,
  created_at timestamptz not null default now(),
  updated_at timestamptz not null default now()
);

create table if not exists individuals (
  id uuid primary key default gen_random_uuid(),
  individual_code text not null unique,
  case_id uuid references cases(id),
  sex text not null default 'unknown' check (sex in ('male','female','unknown')),
  sex_method text,
  hatch_date date,
  phenotype_label text,
  phenotype_evidence text,
  image_quality text,
  welfare_blocked boolean not null default false,
  created_at timestamptz not null default now()
);

create table if not exists pedigree_links (
  individual_id uuid primary key references individuals(id) on delete cascade,
  sire_id uuid references individuals(id),
  dam_id uuid references individuals(id),
  parentage_method text,
  parentage_confidence numeric check (parentage_confidence between 0 and 1),
  check (individual_id is distinct from sire_id),
  check (individual_id is distinct from dam_id)
);

create table if not exists phenotype_observations (
  id uuid primary key default gen_random_uuid(),
  individual_id uuid not null references individuals(id) on delete cascade,
  observed_at timestamptz not null,
  molt_stage text,
  ontology_version text not null,
  primary_hypothesis text,
  alternative_hypotheses jsonb not null default '[]'::jsonb,
  feature_values jsonb not null default '{}'::jsonb,
  evidence_level text,
  confidence_level text,
  image_quality text,
  reviewer_id text,
  source_hash text,
  created_at timestamptz not null default now()
);

create table if not exists pairings (
  id uuid primary key default gen_random_uuid(),
  pairing_code text not null unique,
  male_id uuid not null references individuals(id),
  female_id uuid not null references individuals(id),
  paired_at date,
  purpose text not null check (purpose in ('research','target','observation')),
  welfare_review_status text not null default 'pending',
  created_at timestamptz not null default now(),
  check (male_id <> female_id)
);

create table if not exists clutches (
  id uuid primary key default gen_random_uuid(),
  clutch_code text not null unique,
  pairing_id uuid not null references pairings(id) on delete cascade,
  laid_at date,
  egg_count integer check (egg_count >= 0),
  fertile_count integer check (fertile_count >= 0),
  hatched_count integer check (hatched_count >= 0),
  survived_first_molt_count integer check (survived_first_molt_count >= 0),
  notes text,
  created_at timestamptz not null default now()
);

create table if not exists offspring_observations (
  id uuid primary key default gen_random_uuid(),
  clutch_id uuid not null references clutches(id) on delete cascade,
  offspring_id uuid references individuals(id),
  stage text not null check (stage in ('egg','embryo','hatch','fledging','first_molt','adult')),
  sex text not null default 'unknown' check (sex in ('male','female','unknown')),
  sex_method text,
  phenotype_label text not null,
  phenotype_probability jsonb,
  survival_status text,
  observed_at timestamptz,
  created_at timestamptz not null default now()
);

create table if not exists inheritance_models (
  id uuid primary key default gen_random_uuid(),
  model_id text not null,
  version text not null,
  status text not null check (status in ('draft','unit_tested','scientifically_reviewed','active','deprecated')),
  model_definition jsonb not null,
  definition_hash text not null,
  created_at timestamptz not null default now(),
  unique (model_id, version)
);

create table if not exists model_runs (
  id uuid primary key default gen_random_uuid(),
  run_code text not null unique,
  operation text not null check (operation in ('cross','infer','pigment','rank-mates')),
  schema_version text not null,
  engine_version text not null,
  input_hash text not null,
  request_snapshot jsonb not null,
  result_snapshot jsonb not null,
  evidence_grade text,
  human_review_status text not null default 'pending',
  created_at timestamptz not null default now()
);

create table if not exists run_models (
  run_id uuid not null references model_runs(id) on delete cascade,
  model_pk uuid not null references inheritance_models(id),
  weight numeric check (weight between 0 and 1),
  model_snapshot jsonb not null,
  primary key (run_id, model_pk)
);

create table if not exists genotype_posteriors (
  run_id uuid not null references model_runs(id) on delete cascade,
  individual_id uuid not null references individuals(id) on delete cascade,
  model_id text not null,
  genotype_distribution jsonb not null,
  entropy numeric,
  primary key (run_id, individual_id, model_id)
);

create table if not exists reviews (
  id uuid primary key default gen_random_uuid(),
  run_id uuid not null references model_runs(id) on delete cascade,
  reviewer_role text not null,
  reviewer_id text not null,
  decision text not null check (decision in ('approved','needs_more_data','blocked','rejected')),
  rationale text not null,
  created_at timestamptz not null default now()
);

create table if not exists audit_events (
  id uuid primary key default gen_random_uuid(),
  entity_type text not null,
  entity_id text not null,
  event_type text not null,
  actor_id text,
  previous_hash text,
  new_hash text,
  payload jsonb not null default '{}'::jsonb,
  created_at timestamptz not null default now()
);

create index if not exists idx_phenotype_individual_time
  on phenotype_observations(individual_id, observed_at);
create index if not exists idx_pairings_parents
  on pairings(male_id, female_id);
create index if not exists idx_model_runs_input_hash
  on model_runs(input_hash);
