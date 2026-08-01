# BioBirdDNA Genetics Skill v0.2.0

Bá»™ khung skill vÃ  code máº«u cho dá»± Ã¡n **Biologicalbird Community / BioBirdDNA**, táº­p trung vÃ o:

- mÃ´ phá»ng cÃ´ng thá»©c ghÃ©p mÃ u lÃ´ng theo mÃ´ hÃ¬nh di truyá»n;
- suy luáº­n mÃ´ hÃ¬nh di truyá»n á»©ng viÃªn tá»« pháº£ há»‡ vÃ  toÃ n bá»™ Ä‘á»i con;
- xá»­ lÃ½ kiá»ƒu gene cha máº¹ chÆ°a cháº¯c cháº¯n báº±ng phÃ¢n bá»‘ xÃ¡c suáº¥t;
- chuáº©n hÃ³a pairing, clutch, má»i trá»©ng/con vÃ  timeline phenotype mÃ  khÃ´ng Ä‘á»“ng nháº¥t missing vá»›i zero;
- Ä‘o recovery, uncertainty, calibration theo design vÃ  absolute fit trÃªn dá»¯ liá»‡u synthetic;
- tÃ­nh quan há»‡ vÃ  há»‡ sá»‘ cáº­n huyáº¿t Ä‘á»i con;
- xáº¿p háº¡ng cáº·p ghÃ©p theo má»¥c tiÃªu nghiÃªn cá»©u, khÃ´ng chá»‰ theo mÃ u;
- sÃ ng lá»c kiá»ƒu hÃ¬nh sáº¯c tá»‘ T02 báº±ng rule-based code;
- giá»¯ ranh giá»›i rÃµ giá»¯a áº£nh, phenotype, mÃ´ hÃ¬nh di truyá»n vÃ  gene xÃ¡c nháº­n.

## ThÃ nh pháº§n chÃ­nh

```text
BioBirdDNA_Genetics_Skill_v0.2.0/
â”œâ”€â”€ SKILL.md
â”œâ”€â”€ agents/openai.yaml
â”œâ”€â”€ references/
â”œâ”€â”€ planning/
â”œâ”€â”€ knowledge/T02_knowledge_base_v1.0.json
â”œâ”€â”€ assets/schemas/
â”œâ”€â”€ assets/templates/
â”œâ”€â”€ examples/input/
â”œâ”€â”€ examples/output/
â”œâ”€â”€ scripts/
â”œâ”€â”€ src/biobirddna_genetics/
â”œâ”€â”€ npm/
â”œâ”€â”€ standalone/
â””â”€â”€ tests/
```

## TÃ i liá»‡u kiáº¿n trÃºc cho model/coder

- `BioBirdDNA_Genetics_Skill_Blueprint_v0.1.0.md`: báº£n há»£p nháº¥t gáº§n 10.000 tá»« gá»“m skill contract, kiáº¿n trÃºc lá»›p, data contract, cÃ´ng thá»©c, API/database draft, roadmap vÃ  handoff.
- `planning/MODEL_CODER_HANDOFF.md`: prompt, read order vÃ  quy táº¯c patch dÃ nh cho coding model tiáº¿p theo.
- `planning/backlog.json`: backlog cÃ³ phase, priority, acceptance criteria vÃ  dependency.
- `planning/adr/`: cÃ¡c quyáº¿t Ä‘á»‹nh kiáº¿n trÃºc khÃ´ng Ä‘Æ°á»£c thay Ä‘á»•i Ã¢m tháº§m.

## Bá»‘n lá»±a chá»n code

### Báº£n modular

DÃ¹ng Ä‘á»ƒ phÃ¡t triá»ƒn thÃ nh skill production:

```text
src/biobirddna_genetics/
```

CÃ¡c module Ä‘Ã£ cÃ³:

- `cross.py`: Punnett engine cho genotype cháº¯c cháº¯n hoáº·c posterior;
- `inference.py`: likelihood/Bayes thÄƒm dÃ² giá»¯a nhiá»u mÃ´ hÃ¬nh;
- `pedigree.py`: ma tráº­n quan há»‡ vÃ  cáº­n huyáº¿t;
- `optimizer.py`: xáº¿p háº¡ng cáº·p vÃ  expected information gain;
- `pigment.py`: xáº¿p háº¡ng phenotype T02, khÃ´ng xÃ¡c nháº­n gene;
- `breeding_records.py`: chuáº©n hÃ³a há»“ sÆ¡ sinh sáº£n theo giai Ä‘oáº¡n vÃ  Ä‘á»‘i soÃ¡t sá»‘ lÆ°á»£ng;
- `simulation.py`: harness recovery/calibration synthetic cÃ³ seed vÃ  cá»•ng misspecification;
- `evidence.py`: cháº·n claim vÆ°á»£t evidence;
- `service.py`, `cli.py`: orchestration vÃ  CLI.

### Báº£n má»™t file

```text
standalone/biobirddna_genetics_standalone.py
```

Báº£n nÃ y lÃ  snapshot tiá»‡n dá»¥ng cho `cross`, `infer` vÃ  `pigment` báº±ng Python
standard library. Python modular má»›i lÃ  implementation chuáº©n cho Ä‘á»§ sÃ¡u operation
v0.2.0.

### Báº£n Node.js má»™t file

```text
standalone/biobirddna_genetics_node.mjs
```

Báº£n Node.js má»™t file váº«n lÃ  cross simulator tá»‘i giáº£n. Package ESM trong `npm/`
khÃ´ng thay tháº¿ file standalone nÃ y.

### Package Node ESM

```text
npm/
```

Package khÃ´ng cÃ³ runtime dependency vÃ  triá»ƒn khai Ä‘á»§ sÃ¡u operation theo cÃ¹ng
contract v0.2.0 vá»›i Python modular. NÃ³ cÃ²n cÃ³ feedback adapter review-only; test
packed-install xÃ¡c minh cáº£ hai CLI hoáº¡t Ä‘á»™ng ngoÃ i cÃ¢y source. Simulation dÃ¹ng
cÃ¹ng seed/PRNG vÃ  cÃ³ exact cross-runtime parity; inference khÃ³a count quan sÃ¡t á»Ÿ
`1_000_000` Ä‘á»ƒ giá»¯ work budget vÃ  numerical parity.

## YÃªu cáº§u

- Python 3.11 trá»Ÿ lÃªn.
- Node.js 20 trá»Ÿ lÃªn cho package npm/Node tÃ¹y chá»n.
- KhÃ´ng báº¯t buá»™c thÆ° viá»‡n ngoÃ i.

## CÃ i package tÃ¹y chá»n

KhÃ´ng cáº§n cÃ i Ä‘á»ƒ cháº¡y cÃ¡c script. Khi muá»‘n dÃ¹ng entry point package:

```bash
python -m pip install --no-build-isolation -e .
biobirddna-genetics cross --input examples/input/cross_autosomal_recessive.json --pretty
```

`--no-build-isolation` há»¯u Ã­ch trong mÃ´i trÆ°á»ng offline Ä‘Ã£ cÃ³ sáºµn setuptools.

## Package npm/Node public release candidate

YÃªu cáº§u Node.js 20 trá»Ÿ lÃªn:

```bash
cd npm
npm test
npm pack --dry-run
```

Package lÃ  `@biobirddna/genetics`, `UNLICENSED`, khÃ´ng cÃ³ runtime dependency vÃ 
cÃ³ metadata cho public registry scope. ÄÃ¢y chá»‰ lÃ  release candidate Ä‘Ã£ Ä‘Æ°á»£c kiá»ƒm
tarball: khÃ´ng cÃ³ npm login/publish, GitHub push hoáº·c xÃ¡c nháº­n registry trong cÃ¢y
skill nÃ y. Chá»‰ sau hÃ nh Ä‘á»™ng external Ä‘Æ°á»£c owner á»§y quyá»n má»›i xÃ¡c minh báº±ng
`npm view` package/version/integrity. Hai CLI:

```text
biobirddna-genetics
biobirddna-feedback
```

Sau khi táº¡o tarball Ä‘Ã£ kiá»ƒm chá»©ng, cÃ³ thá»ƒ cÃ i local:

```bash
npm install --global ./biobirddna-genetics-0.2.0.tgz
biobirddna-genetics cross --input ./cross.json --output ./result.json
biobirddna-feedback --input ./review-envelope.json \
  --expected-integrity <sha256> --output ./feedback.json
```

CLI feedback chá»‰ táº¡o artifact review/proposal báº¥t biáº¿n cho quy trÃ¬nh AutoResearch.
NÃ³ khÃ´ng gá»i experiment promotion, khÃ´ng tá»± approve/apply vÃ  khÃ´ng dÃ¹ng health
gate 3.0 cá»§a Page lÃ m Ä‘iá»ƒm BioBirdDNA. Input lÃ  normalized production-review
wrapper, khÃ´ng pháº£i raw operation result. Mutation intent, claim gene/variant
Ä‘Æ°á»£c xÃ¡c nháº­n dÆ°á»›i G5 vÃ  AutoResearch health-score trÃ¡ hÃ¬nh Ä‘á»u fail closed; xem
`references/autoresearch-feedback.md`.

## OAuth MCP vÃ  plugin ChatGPT/Codex cÃ¡ nhÃ¢n

GEN-015 thÃªm adapter cloud stateless trong dedicated Supabase project. Endpoint
chá»‰ expose Ä‘Ãºng sÃ¡u compute tool vÃ  má»™t feedback tool, táº¥t cáº£ read-only. OAuth
2.1 kiá»ƒm Ä‘Ãºng MCP audience vÃ  allowed user trÆ°á»›c dispatch; adapter khÃ´ng táº¡o
table, bucket, queue, model-run persistence hoáº·c quyá»n welfare override.

ChatGPT Web dÃ¹ng private custom MCP App Ä‘Ã£ káº¿t ná»‘i OAuth. Universal plugin cá»¥c bá»™
cho ChatGPT desktop Work mode/Codex chá»‰ tham chiáº¿u Ä‘Ãºng registered app ID vÃ  bundle
skill Ä‘Ã£ kiá»ƒm checksum. Token, password, app ID tháº­t vÃ  personal marketplace tree
khÃ´ng thuá»™c release source. Public authenticated review váº«n lÃ  cá»•ng external cá»§a
GEN-016, khÃ´ng pháº£i tráº¡ng thÃ¡i Ä‘Ã£ hoÃ n táº¥t.

Node updater/watcher cháº¡y full source/cloud validation, phÃ¢n loáº¡i sensitive diff,
yÃªu cáº§u security pin Ä‘Ã£ duyá»‡t, deploy immutable snapshot, verify live vÃ  rollback
náº¿u tháº¥t báº¡i. Watcher khÃ´ng gá»i hoáº·c mutate AutoResearch.

## Kiá»ƒm chá»©ng vÃ  cÃ i skill global

Cháº¡y gate Ä‘áº§y Ä‘á»§ trÃªn Ä‘Ãºng cÃ¢y file chuáº©n bá»‹ cÃ i:

```bash
python -B scripts/validate_package.py
python -B scripts/verify_examples.py
python -B scripts/verify_package.py
```

TrÃªn Windows, cÃ i theo cÆ¡ cháº¿ staging, kiá»ƒm checksum rá»“i má»›i atomic rename:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/install_skill.ps1 `
  -DestinationRoot "$env:USERPROFILE\.codex\skills"
```

TrÃªn Linux/macOS:

```bash
./scripts/install_skill.sh "${CODEX_HOME:-$HOME/.codex}/skills"
```

Installer giá»¯ install lock, tá»« chá»‘i ghi Ä‘Ã¨/Ä‘ua race vá»›i target Ä‘Ã£ tá»“n táº¡i,
kiá»ƒm checksum trÆ°á»›c vÃ  sau final rename, vÃ  khÃ´ng Ä‘Æ°a `__pycache__`, `.pyc`
hoáº·c `.pyo` vÃ o báº£n cÃ i. TrÃªn Windows cÃ³ locale khÃ´ng máº·c Ä‘á»‹nh UTF-8, cháº¡y
validator metadata cá»§a Codex báº±ng `python -X utf8 .../quick_validate.py <skill-root>`.
Shell installer fail-loud náº¿u Linux/filesystem khÃ´ng há»— trá»£
`renameat2(RENAME_NOREPLACE)` hoáº·c macOS volume khÃ´ng há»— trá»£
`renamex_np(RENAME_EXCL)`; khÃ´ng fallback sang rename cÃ³ thá»ƒ ghi Ä‘Ã¨.

## Cháº¡y nhanh

### Windows

```bat
run_tests.bat
run_demo.bat
```

### Linux/macOS

```bash
./run_tests.sh
./run_demo.sh
```

### Cháº¡y tá»«ng chá»©c nÄƒng

```bash
python -B scripts/biobirddna_cli.py cross \
  --input examples/input/cross_autosomal_recessive.json \
  --output examples/output/cross_autosomal_recessive.json

python -B scripts/biobirddna_cli.py infer \
  --input examples/input/infer_z_linked_case.json \
  --output examples/output/infer_z_linked_case.json

python -B scripts/biobirddna_cli.py pigment \
  --input examples/input/pigment_progressive_greying.json \
  --output examples/output/pigment_progressive_greying.json

python -B scripts/biobirddna_cli.py rank-mates \
  --input examples/input/mate_ranking_case.json \
  --output examples/output/mate_ranking_case.json

python -B scripts/biobirddna_cli.py normalize-records \
  --input examples/input/normalize_breeding_records.json \
  --output examples/output/normalize_breeding_records.json

python -B scripts/biobirddna_cli.py simulate-recovery \
  --input examples/input/simulate_recovery.json \
  --output examples/output/simulate_recovery.json
```

## VÃ­ dá»¥ cÃ´ng thá»©c

### Láº·n trÃªn NST thÆ°á»ng

```text
N/m Ã— m/m
```

Ká»³ vá»ng mÃ´ hÃ¬nh:

```text
50% N/m â€” phenotype bÃ¬nh thÆ°á»ng nhÆ°ng cÃ³ thá»ƒ mang allele
50% m/m â€” phenotype má»¥c tiÃªu náº¿u penetrance = 1
```

### Láº·n liÃªn káº¿t Z á»Ÿ chim

```text
Trá»‘ng m/m Ã— mÃ¡i N/W
```

Ká»³ vá»ng mÃ´ hÃ¬nh:

```text
50% con trá»‘ng N/m â€” phenotype bÃ¬nh thÆ°á»ng, mang allele
50% con mÃ¡i m/W â€” phenotype má»¥c tiÃªu
```

ÄÃ¢y chá»‰ lÃ  káº¿t quáº£ cá»§a **mÃ´ hÃ¬nh kÃ½ hiá»‡u**. `N` vÃ  `m` khÃ´ng pháº£i tÃªn gene chÃ o mÃ o Ä‘Ã£ xÃ¡c nháº­n.

## Kiáº¿n trÃºc Ä‘á»c theo thá»© tá»±

Cho model/coder má»›i:

1. `SKILL.md`
2. `references/scientific-boundaries.md`
3. `references/architecture.md`
4. `references/data-contracts.md`
5. `references/algorithms.md`
6. `planning/MODEL_CODER_HANDOFF.md`
7. `planning/backlog.json`
8. code trong `src/`
9. test trong `tests/`

## Pháº¡m vi v0.2.0

ÄÃ£ cháº¡y Ä‘Æ°á»£c:

- má»™t locus, hai allele;
- autosomal recessive/dominant/incomplete dominant;
- avian Z-linked recessive/dominant;
- penetrance Ä‘Æ¡n giáº£n;
- lá»—i cháº¥m phenotype Ä‘Æ¡n giáº£n;
- posterior kiá»ƒu gene cha máº¹;
- model comparison thÄƒm dÃ²;
- pedigree relationship matrix;
- mate ranking Ä‘a má»¥c tiÃªu;
- rule-based phenotype screen.

ChÆ°a production-ready cho:

- hai locus/epistasis;
- polygenic trait;
- survival model theo genotype;
- missing-not-at-random;
- genomic relationship matrix tá»« SNP;
- image segmentation/calibration trá»±c tiáº¿p;
- spectrophotometry/Raman/HPLC/LC-MS;
- database vÃ  API multi-user;
- clinical or invasive sampling workflow.

CÃ¡c pháº§n nÃ y Ä‘Ã£ cÃ³ backlog vÃ  interface má»Ÿ rá»™ng trong `references/implementation-roadmap.md`.

## NguyÃªn táº¯c khoa há»c cá»‘t lÃµi

```text
áº¢nh â†’ mÃ´ táº£ kiá»ƒu hÃ¬nh
Pháº£ há»‡ + Ä‘á»i con â†’ mÃ´ hÃ¬nh di truyá»n á»©ng viÃªn
DNA Ä‘á»“ng phÃ¢n ly â†’ locus/biáº¿n thá»ƒ Ä‘Æ°á»£c há»— trá»£
TÃ¡i láº­p/chá»©c nÄƒng â†’ cÆ¡ cháº¿ xÃ¡c nháº­n máº¡nh
```

Skill khÃ´ng cho phÃ©p nháº£y tháº³ng tá»« áº£nh sang gene.

