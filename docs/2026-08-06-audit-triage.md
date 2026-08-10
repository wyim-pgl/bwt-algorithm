# 밤샘 병렬 감사 — 통합 분류 (2026-08-06 06:1x)

세 에이전트가 원고를 구역별로 나눠 수치 주장을 근거 파일과 대조했다. **합계 약 540건 검증, 약 62건 지적**,
여기에 소프트웨어 버그 10건. §4.4는 이미 수정해서 감사 대상에서 제외했다.

| 에이전트 | 범위 | 검증 / 지적 |
|---|---|---|
| audit-results-human-arab | §3.1 human, §3.2 Arabidopsis | ~230 / 13 |
| audit-results-maize | §3.3 (3A, 3A-b, 3B, 3C) | ~232 / 16 |
| audit-frontmatter-discussion | Abstract, §1, §2, §4.1-4.3, §4.5 | 78 / 25 |
| software-health | 코드·테스트 | 출력 정확성 버그 3건 (F절) |

**중요:** 에이전트는 편집 전 스냅샷을 읽었다. 아래 "이미 해결됨"에 든 것은 재작업 불필요.

---

## A. 오늘 밤 적용 완료 (내가 직접 검증함)

| 항목 | 이전 → 이후 | 검증 방법 |
|---|---|---|
| Col-CEN 어셈블리 크기 (4곳) | 134 Mb → **132 Mb** | FASTA 직접 계수 = 132,081,078 bp |
| §3.1 bp precision 최고점 주장 | "most precise point measured here" → ULTRA 대비로 한정 + tantan 70.22 / TRF 93.16 명시 | `figure_curve_data.csv` |
| §3.3 본문 maize 수치 | 41.01 Mb → **41.31**, 1.50 M → **1.51 M** | `score_maize_3a.scan(remeas_maize.bed)` |
| ULTRA maize 스케일링 | 160배 → **151배**, 0.2 h → **0.23 h** (3곳) | 로그 13:47.86=827.86 s, 34:43:56=125,036 s |
| §1 게놈 규모 | "two orders of magnitude" → "twenty-threefold" | 132 Mb ↔ 3.1 Gb = 23.5× |

## B. 이미 해결됨 — 에이전트가 편집 전 스냅샷을 봄

- §3.2 L228 / §4.3 L535 TRF cap: **둘 다 이미 200 bp**. frontmatter 에이전트가 로그와
  manifest로 **독립적으로 200을 확증**했다 → 내 어젯밤 revert가 옳았음이 이중 확인됨.
- Supplementary S1.1 copy-count 공식, S2 환경변수 블록, Methods 스레드 수, §1 "one index
  can serve", 런타임 12.1 h / 31 min, §4.4 메모리 21.86 GB — 전부 적용 완료.

---

## C. 🔴 판단이 필요한 중대 항목 — 아침에 여기부터

### C1. TRASH·NCRF 비용 셀이 표 간 복사됨 (maize) — 셀 8~9개

ULTRA·TRF에서 고쳤던 결함이 **TRASH와 NCRF에는 남아 있었다.** 실측(로그 `Elapsed` +
`Maximum resident set size`):

| 표 / 행 | 게시값 | 실측 |
|---|---|---|
| 3B TRASH(template) | 59.9 h / 12.78 GB | **101.7 h / 120.6 GiB** (126,500,492 KB) |
| 3C TRASH | 59.9 h / 12.78 GB | **97.8 h / 42.8 GiB** (44,910,788 KB) |
| 3A TRASH(template) runtime | 59.9 h | **65.5 h** (65:30:12) |
| 3B TRASH(de novo) runtime | 59.9 h | **58.8 h** (58:50:42) |
| 3A NCRF | <0.01 h / 0.29 GB | **0.132 h / 17.5 GiB** (18,313,536 KB) |
| 3C tantan runtime | 0.5 h | **0.6 h** (35:10.63) |

**파급:** TRASH 3B가 120 GiB면 §4.4의 "BWTandem's footprint is the largest of the de novo
tools compared"와 Abstract의 같은 취지 문장을 재검토해야 한다. (단 그 행은 template 모드다 —
de novo 행은 12.78 GB로 정확하므로 "de novo 중 최대"는 유지될 수 있다. **확인 필요.**)
또 §3.3.1의 "NCRF was by far the fastest"는 실제 0.132 h vs ULTRA 0.23 h로 1.7배 차이일
뿐이고, 메모리는 0.29 GB가 아니라 **17.5 GB로 BWTandem 다음 2위**다 — 문장이 정반대 인상을
준다.

### C2. Table 3C 각주 ‖ 가 사실이 아님

"TRASH de novo and template modes gave identical output on this class" — maize TRASH de novo
exp3C는 **완료된 적이 없다**(로그가 "Working on sequence 4"에서 끊기고 `Elapsed` 블록 없음,
de novo BED 자체가 존재하지 않음). `_trash_exp3C.bed`가 template 파일과 **md5 동일**이다 —
두 모드가 일치한 게 아니라 같은 파일이다. → Table 3C TRASH 행은 **template 단독**으로 표기해야.

### C3. Table 1a / 1c의 일부 셀이 폐기된 BED 값

manifest는 두 표 모두 `remeas_human.bed`(job 5983793)를 선언하는데:

| 셀 | 게시값 | remeas 기준 |
|---|--:|--:|
| Table 1a BWTandem Adj. Prec. | 79.68 | **79.49** |
| Table 1c BWTandem 행 전체 | 98,624 / 3.45 / 66.00 / 13.48 / 29.28 | **98,339 / 3.44 / 65.97 / 13.44 / 29.24** |
| Table 1a ULTRA | 3,216,708 / 53.65 | **3,216,710 / 53.66** (1b·1d는 이미 이 값) |
| Table 1a TRASH | 5,132 / 95.27 / 99.88 | **5,134 / 95.25 / 99.86** (1e는 이미 99.86) |

ULTRA 차이의 원인이 밝혀졌다: 게시 CSV는 **GCF_000001405.26** BED, 재채점은
**GCA_000001405.15** BED. `WP0-B_TABLE1_FINDINGS.md`가 "unexplained +2"로 남긴 항목의 답.

### C4. Methods의 실행 환경 기술이 모든 행에서 틀림

- "each genome processed as an independent SLURM job allocated **four CPU cores and 64 GB**" —
  실제로는 행마다 다르다. 표에 실린 BWTandem은 `run_remeasure.sbatch` = **2 CPU**.
  4 CPU/64 GB는 초기 `filip/bwt/run.sbatch`뿐이고 **표의 어떤 BWTandem 수치도 거기서 오지 않았다.**
- "single-threaded tools were **parallelized across two cores**" — 사실이 아니다. TRF·tantan은
  1 CPU 단일 프로세스, mreps는 염색체별 순차 루프. **문장 삭제가 맞다.**
- 내가 어젯밤 넣은 "Peak memory is the SLURM cgroup maximum (sacct MaxRSS)"는 **BWTandem 행에만
  참**이다. 경쟁 도구 셀은 GNU time 값과 정확히 일치한다(TRASH 15,298,244 KB=14.59 등).
  → 행별 규칙으로 다시 써야 한다.
- 인간 게놈 accession: 밝힌 것은 `GCF_000001405.26`인데 채점된 경쟁 도구 BED는 전부
  **GCA_000001405.15** 산출물이고 BWTandem은 `hg38_primary.fa`다.

### C5. 근거가 없는 주장

- **"NCRF exceeded 800 GB on human"** — 인간 NCRF 실행 스크립트도 로그도 BED도 없다.
  삭제하거나 "Col-CEN 80.96 GB에서의 외삽"으로 낮춰야 한다.
- **"All tools substantially overpredict CentC"** — 도구별 예측 CentC 양이라는 수치가 디스크
  어디에도 없고, tantan은 배열 커버리지 1.10%로 **과예측하지 않는다**. `pre-wp0b-bak`에서
  이월된 문장이다.
- §4.1 "Candidate families can be identified directly from an assembly" — family 식별 실험이
  이 논문에 없다.
- §4.3 "alignment-based methods must first propose a candidate period" — **자기모순**이다.
  §1이 "TRF uses k-tuple matching followed by wraparound DP"라고 적고 있고, k-tuple 위치
  간격에서 후보 주기를 먼저 얻는 것은 TRF도 한다.
- §1 "over 60 heritable human diseases" — 인용된 Mirkin(2007)은 ~40 수준. "20–30%
  inter-copy divergence", 임상 스크리닝 서술도 무인용.

### C6. 정의·해석 오류

- **§3.3.2 "TRASH is the only tool approaching one call per array (0.05)"** — frag score는
  `1/(배열당 call 수)`이므로 0.05는 **배열당 약 20개**이고 one-call-per-array는 1.00이다.
  실제 배열당 call: TRASH ~20, tantan ~106–224, ULTRA ~303–878, BWTandem ~253–475,
  TRF knob180 **1,705**. "hundreds"도 천 단위를 놓친다.
- **§4.1 "Template-guided tools such as NCRF and TRASH"** — TRASH는 세 게놈 모두 **de novo로
  실행**되었고 Table 2에서 85.03% coverage를 낸다. 모티프를 요구하는 것은 NCRF뿐이다.
- **Abstract "Varying one sensitivity threshold across four runs"** — 실제로는 세 가지가
  바뀐다(P는 pass off, B는 identity 0.76+MAX_P 50, F는 0.72+MIN_COPIES 3, H는 0.72).
  Table 1d 캡션은 이미 이를 밝히고 있어 **Abstract와 캡션이 충돌**한다.
- **Abstract "at the cost of splitting each array into many calls"** — Table 3C에서 BWTandem·
  TRF·ULTRA·tantan 모두 frag 0.00으로 공통이다. BWTandem 고유 비용처럼 읽힌다.
- **§4.2 corroboration 비율 두 개가 뒤바뀜** — "ULTRA or tantan" 순서에 맞추면
  ULTRA 53.5%, tantan 59.7%다(현재 59.7/53.5).
- **§1 "TRF restricted to 2–500 bp"** — `trf <fa> 2 7 7 ...`의 앞 2는 **match weight**이지
  최소 주기가 아니다. Methods 자신이 "TRF has no minimum-period parameter"라고 적고 있다.
  `manifest.tsv`의 `period_range=2-500 / 2-200`도 같은 오독을 전파 중.
- **§1 CEN180 "2%"** — 2.25%는 ULTRA 기본설정의 *centromere coverage*다. 같은 지표
  (monomer recall)로는 ULTRA 0.62%. **지표를 맞추면 주장이 더 강해진다.**

### C7. 비대칭 공개 (BWTandem에 불리)

- **tantan에는 최대 주기 옵션 `-w`(기본 100)가 있다.** 올리지 않고 실행했다. Arabidopsis에서
  ULTRA·mreps는 범위를 올려 재실행해 준 것과 비대칭이고, "above period 100 the field empties"
  (Abstract·§3.1·Table 1c)가 여기 기댄다.
- **maize Table 3A**는 period 1–2,000을 탐색한 BWTandem 15.4 h를 period ≤6만 탐색한 TRF 5.2 h /
  ULTRA 0.23 h와 나란히 놓는다. 캡션에 이 비대칭이 없다.
- **TideHunter**는 실제로 실행됐고(Col-CEN 5개 염색체 합 ~22.1 h) 어느 표에도 없다.
  "seven tools" 문구와 배제 근거를 맞춰야 한다.

---

## D. 미해결 · 재확인 필요

- Table 1c의 **98,624건이 재현되지 않는다** — remeas BED 직접 집계로 `>100`은 98,339,
  `≥100`은 99,142. scorer의 주기 정의 확인 필요. Abstract의 stratum과 §2.2의 "2.5%"가 파생.
- Table 1a "Unique Regions" 863,654 — `compute_unique()`가 **chr4 전량 오라벨된 mreps BED**를
  비교 집합에 넣는다. 염색체마다 다른 필터를 거친 값이며 remeas BED로 재계산된 적이 없다.
- §4.5의 unique-region 특성화는 **BWTandem-F(p≤100)** 기준 875,224건인데 Table 1a의 unique는
  1–2,000 실행의 863,654건이다. 어느 실행 기준인지 명시 필요.
- §4.5 flank 주장은 **n=14**(59개 배열 중)이고 11개가 TR-1이다. §3.3.3은 밝히지만 §4.5는 아니다.
- Table 3C BWTandem: **offset 211.5 → 196.03**, coverage 59.02 → 59.03 (구 BED vs remeas).
- Table 2 메모리 단위 혼용: ULTRA 4.41 → **4.20 GiB**, mreps 0.18 → **0.17 GiB**
  (두 재실행 행만 KB/1e6, 나머지는 KiB). Table 3의 ULTRA 3B 10.15 → **9.68**, 3C 3.54 → **3.38**도 동일.
- **좋은 소식:** `WP2.2b_GAPFILL_ABLATION.md`의 Result 3("현재 코드가 published Table 2를
  재현 못함, CEN180 1,380 → 25")은 **철회 가능**하다. manifest 지정 BED를 지정 scorer로
  채점하면 160,740 / 84.59% / 1,380 / 57.17% / 99.73%로 **완전 일치**한다. ablation이
  defaults로 돌았던 탓이었다.

---

## E. 권고 순서

1. **C1·C2 (maize 비용 셀 + 3C 각주)** — 사실 오류이고 로그로 즉시 확정 가능. §4.4의
   "largest footprint" 문장이 여기 걸린다.
2. **C4 (Methods 실행 환경)** — 재현성 항목이라 심사에서 반드시 걸린다.
3. **C3 (Table 1a/1c 셀)** — 표 간 정합성. ULTRA GCF/GCA 건은 원인이 밝혀졌으니 각주로 처리 가능.
4. **C5 (무근거 주장)** — 삭제 또는 근거 확보. 삭제해도 논지는 유지된다.
5. **C6 (정의 오류)** — frag score 오독과 §4.1 TRASH 분류가 특히 눈에 띈다.
6. **C7 (비대칭 공개)** — 전부 BWTandem에 불리한 방향이라 밝히는 편이 심사에 안전하다.

`resume.md`의 하드룰대로, 위 항목 중 BWTandem에 불리한 것들(C1의 NCRF 메모리, C6의 frag,
C7 전부)도 그대로 남겼다.

---

# F. 소프트웨어 감사 (software-health) — 출력 정확성 버그 3건

**아무것도 고치지 않았다.** 셋 다 출력 바이트를 바꾸는 수정이고, 지금은 원고 수치를
확정하는 중이라 밤사이 무단으로 바꿀 수 없다. 에이전트도 같은 판단을 했다.
전체 스위트 **94 passed**, 두 parity 테스트 모두 **실제로 비교했다**(skip 아님, 34 passed).
네 개 `.so` 전부 최신이고 재빌드 트리거 없음.

## F1. 🔴 HIGH — STRfinder CSV가 구조적으로 깨져 있다

`src/models.py:164`
```python
genotype_struct = f"{motif_len}[{cons}]{complete_copies},{truncated}"
```
콤마가 따옴표로 감싸이지 않아 **데이터 행이 13필드, 헤더는 12필드**다. 5번 열부터 전부
한 칸씩 밀린다. 실측: `synth_mixed.fa --format strfinder` → 헤더 12, 데이터 9행 × 13.
STRfinder 규격이 그 필드에 콤마를 요구하므로 **따옴표를 씌우는 것이 수정**이지 콤마 제거가 아니다.

## F2. 🔴 HIGH — C 정렬 경로가 per-copy 상세를 조용히 비운다

`src/motif_utils.py:509-513`이 `variations` / `copy_sequences` / `error_counts`를 `[]`로
하드코딩한다. 재현:
```
use_c=True  : variations=[]                    copy_sequences=0  error_counts=0
use_c=False : variations=['2:1:A>T','4:1:A>G'] copy_sequences=6  error_counts=6
```
→ STRfinder의 Variations 열이 **표준 빌드에서 항상 "-"**, `libalign_accel`이 없는 빌드에서만
실제 주석이 나온다. **도구의 답이 빌드 산출물 존재 여부에 달려 있다** — `test_align_parity.py`가
막으려던 바로 그것이다.

두 테스트가 못 잡는 이유: `BWT_DISABLE_NATIVE`는 Cython `_accelerators`만 게이트하고
**ctypes `_c_align_lib`은 게이트하지 않는다**(`motif_utils.py:10-14`이 독립적으로 로드).
그래서 parity 테스트의 두 경로가 모두 "-"를 내고 같다고 판정한다.

**최고 레버리지 수정은 필드 복사가 아니라 `_c_align_lib`이 `BWT_DISABLE_NATIVE`를 따르게
만드는 것** — 기존 parity 스위트 전체가 그 라이브러리에도 자동 적용된다.

## F3. 🟠 MEDIUM — satellite gap-fill이 서열 끝을 넘는 좌표를 낸다 (기본 동작)

`src/finder.py:459`가 `n = len(text_arr)`로 **sentinel `$`를 포함**한다. 같은 파일
`:564`의 catch-all은 `effective_length(text_arr)`를 올바르게 쓴다. 재현: 18,810 bp 서열에서
`end=18811`, **overhang 1 bp**.

**원고 관련 확인 필요:** BED 구간이 염색체 길이를 넘으면 bedtools가 malformed로 거부한다.
bp 지표는 `bedtools sort | bedtools merge`를 거치므로, **published Col-CEN 수치가 이
1 bp 오버행의 영향을 받았는지 확인해야 한다.** gap-fill은 기본 켜짐이고 Col-CEN 99.8%
monomer recall 헤드라인을 만든 기능이다.

부수: gap-fill은 catch-all(`:612`)과 달리 motif의 ACGT 검증을 하지 않아 `N`이나 `$`가
motif 필드에 들어갈 수 있다.

## F4. 🟠 MEDIUM — BED 6번 열이 정수가 아니다 (기본 동작)

`src/models.py:14`는 `tier: int`인데 `finder.py:540`은 `tier="satellite"`,
`:618`은 `tier="catchall"`을 넣는다. 실제 출력:
```
syn  5130  18810  CGCGCGCC...  80.0  satellite  0.501  +
```
satellite 패스는 기본 켜짐이므로 **실제 센트로미어 데이터의 기본 출력에 나온다.**
`int(fields[5])`를 하는 스코어링 스크립트는 깨진다.

## F5. 그 밖 (낮은 우선순위, 상세는 에이전트 원문)

- **MEDIUM-5** `motif_utils.py:456-466`: C 텍스트 포인터를 `id(sequence)`로 캐시. 문자열
  주소 재사용 시 다른 서열이 캐시 히트. 4/4 재현되지만 **CLI로는 도달 불가**(모든 tier가
  장수명 속성을 넘김). 주목할 점은 `test_align_parity.py`가 `_KEEP` 리스트로 id 재사용을
  막고 있다는 것 — **테스트의 우회가 결함을 가리고 있다.**
- **LOW-6** C의 256M 셀 가드에 Python 대응물이 없어 motif_len ~16 kb에서 OOM 가능.
- **LOW-7** `motif_utils.py:270,730,803`의 `except Exception: pass` — CLAUDE.md 자신의
  "재현 못하는 fallback은 raise해야 한다" 규칙과 반대.
- **INFO-9** `CATCHALL_MIN_IDENTITY`가 end-to-end로 **단조가 아니다**(0.72→0.68에서 커버리지
  481→477 bp). 패스 단독은 완벽히 단조. 원인은 covered 마스크 경쟁 + `finder.py:604`의
  `seg.mean() > 0.3` 일괄 기각. 폭은 0.8%지만 **published frontier가 그 knob에 대해 국소
  단조임이 보장되지 않는다.**
- **INFO-10** 죽은 C 코드 3개(`batch_process_lcp_candidates`, `extend_mismatch_c`,
  `batch_backward_search`) — 호출자 0.

## F6. 테스트 커버리지 구멍

- 가속기 5개 심볼 중 **`align_unit_to_window`만 C-vs-Python 값 비교가 없다.** 그 Python
  대응물은 손으로 쓴 별도 banded DP(`motif_utils.py:273-424`)로, 7월에 31% 불일치를 낸
  `align_accel.c` 루프와 **같은 코드 형태**다.
- **env knob 34개가 동작 테스트 0.** 특히 `TIER1_FMSCAN_MIN_DENSITY`/`MIN_LLR`(v2.1
  정밀도 선도 op-point)와 `CATCHALL_*` 전 계열(v2.2 82% recall 영역) — **원고 수치가
  기대고 있는 바로 그 knob들**이다. `BWT_INDEX_CACHE`도 `test_index_cache.py`가
  `save_index`/`load_index`를 직접 부를 뿐 **env var 배선(`finder.py:54`)은 미검증**이다.
  (내가 어젯밤 추가한 `test_env_var_docs.py`는 문서↔코드 *목록*을 고정할 뿐 동작은 안 본다 —
  다른 축이다.)
- 에이전트가 테스트 3종을 프로토타입해서 **실제로 돌려봤다**(파일은 추가하지 않음):
  1. `align_unit_to_window` 차등 테스트 — **1 passed, 2000 케이스 전부 일치**(회귀 가드, 버그 없음)
  2. `align_repeat_region` 전체 필드 parity — **지금은 실패하도록 설계됨**(F2의 가드)
  3. `CATCHALL_*` knob 4종 — **4 passed, 107 s**. 단조성 테스트를 end-to-end로 짰다가
     실패한 것이 INFO-9의 발견 경위다.
  코드는 에이전트 원문에 그대로 있고, 산출물은 scratchpad `audit/` 아래에 있다.

## F7. 권고 순서 (에이전트 의견 + 내 판단)

1. **F2를 `_c_align_lib`이 `BWT_DISABLE_NATIVE`를 따르게 하는 방식으로** — 기존 parity
   스위트가 공짜로 확장된다. 그다음 프로토타입 테스트 2를 추가.
2. **F3** — 한 줄(`effective_length`). 단 **출력 바이트가 바뀌므로 원고 수치 확정 후에.**
   그 전에 published Col-CEN BED에 오버행 구간이 실제로 있었는지부터 확인.
3. **F1** — 따옴표 처리. STRfinder 포맷은 논문 표에 쓰이지 않으므로 수치 영향 없음.
4. **F6의 `align_unit_to_window` 차등 테스트** — 이미 통과하므로 지금 추가해도 안전하다.

---

# G. 2026-08-06 새벽 적용 기록 (사용자 결정 반영)

**결정: 1=(a) · 2=약화 · 3=제외 · 4=서브에이전트 · 기계적 수정은 전부 진행.**

## G1. 오버행 실측 — 먼저 확인한 것

| BED | 구간 수 | 염색체 끝 초과 | 최대 |
|---|--:|--:|--:|
| Col-CEN (published & remeas) | 161,187 | **0** | — |
| maize remeas | 2,405,401 | **0** | — |
| **human remeas** | 4,009,261 | **19** | +1 bp |

버그는 실재하나 **영향은 human 19건 × 1 bp**로 수치상 무시할 수준. Col-CEN 헤드라인
(99.8% monomer recall)은 **영향 없음**. 형식 결함이므로 수정은 서브에이전트에 위임.

## G2. 적용 완료

**결정 1(a)** §4.4 — "exceeded 800 GB on human"(근거 없음) 삭제, NCRF 80.96 GB(Arabidopsis)
유지 + **maize TRASH template 120.6 GB** 추가. 앞문장(de novo 중 최대)은 human 기준이라 유지.

**결정 2 (약화)** — CentC 과예측(측정 안 했음 + 밴드가 knob180도 포함 + tantan 예외 명시),
family 식별(→"arrays를 찾는다, family 묶기는 미평가"), §4.3 메커니즘(→TRF도 k-tuple 간격에서
주기를 얻음을 인정하고 차이를 "인덱스가 O(k)에 전체 출현 목록을 준다"로 축소), 질병 수
(over 60 → dozens).

**기계적 수정 전부:**
- maize 비용 셀 10개 — 3A NCRF 0.13/17.47, 3A TRASH(tpl) 65.5, 3B TRASH(dn) 58.8,
  3B TRASH(tpl) **101.7/120.64**, 3C TRASH **97.8/42.83**, 3C tantan 0.6,
  ULTRA 3B 9.68 / 3C 3.38 (GiB 통일)
- Table 3C BWTandem — offset **211.5→196.0**, coverage 59.02→59.03. **독립 재현함**:
  무필터 규칙에서 published BED가 211.53(원고와 일치), remeas가 196.03
- Table 1a — adj prec 79.49, ULTRA 3,216,710/53.66, TRASH 5,134/95.25/99.86
- Table 1c 행 전체 — 98,339 / 3.44 / 65.97 / 13.44 / 29.24, 본문·1b 캡션 동반 수정
- Table 2 메모리 GiB 통일 — ULTRA 4.41→**4.20**, mreps 0.18→**0.17** (sacct 원값 확인)
- Methods — 자원 할당을 행별로 서술, "단일스레드 도구 2코어 병렬화" 삭제,
  sacct/GNU time을 **행별 규칙**으로 명시, accession GCF→GCA + BWTandem은 primary FASTA
- §4.2 corroboration 순서 교정(ULTRA 53.5 / tantan 59.7), §1 TRF "2–500"→500 상한,
  §1 CEN180 2%→**0.6%**(지표 일치), frag score 오독 전면 재서술, §4.1 TRASH 분류

## G3. 남긴 것

- **Table 1e `| BWTandem | 79.68 | 79.99 |`** — 첫 값은 79.49로 고칠 수 있으나 **두 번째
  (leave-one-out)는 재측정 BED로 계산된 적이 없다**(`score_remeas_h_5984654.out`이 명시적으로
  제외: 자기 다른 실행이 corroborator가 되어버림). 한쪽만 고치면 한 행이 두 BED를 섞으므로
  **재채점 필요 항목으로 남김.**
- **결정 3(비대칭 공개) 제외** — tantan `-w`, maize 3A period 비대칭, TideHunter/"seven tools".
  사용자 결정.
- **결정 5** Table 1c의 98,624 vs 98,339 — 위에서 98,339로 맞췄으나 scorer의 주기 정의
  (`>100` vs `>=100`, 각각 98,339 / 99,142)는 여전히 확정 필요.
- C7의 Table 1a "Unique Regions" 863,654 — mreps 오라벨 BED 오염 건, 재계산 안 됨.
