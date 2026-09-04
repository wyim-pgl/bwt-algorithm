# quarantine.md — BWTandem 격리 대장

> 📌 **정본**: 이 파일이 "**쓰면 안 되는 것**"의 유일한 정본이다. 최종 갱신 **2026-09-04**.

## 이 파일에 무엇을 두는가

**세 종류만 둔다.**

1. **계획이 틀린 것** — 전제가 무너져 결론이 뒤집힌 것
2. **이미 시도했는데 안 된 것** — 측정까지 하고 실패한 것
3. **아예 배제·격리된 것** — 폐기된 수치·주장·이름·ID·경로

**미해결은 여기 두지 않는다.** 아직 답이 없는 것, 검증이 남은 것, 해야 할 일은
[`todo.md`](todo.md) 로 간다. 이 둘을 섞으면 "목록에 있으니 끝난 것"과 "목록에 있으니
하지 말 것"이 구분되지 않는다.

## 왜 두는가

다음 문제를 풀 때 **오탐(False Positive) 판정 · 부수 피해 분석 · 위협 분석**의 기준선으로 쓴다.

- **오탐 판정** — 새로 발견한 "결함"이 이미 여기서 기각된 것인지 먼저 본다.
  §3.1 처럼 **자기 스코어링 버그가 만든 허상**을 진짜 결함으로 오인한 전례가 있다.
- **부수 피해 분석** — 어떤 값을 고칠 때 그 값에 딸린 것이 무엇인지 본다.
  §1.2 는 바이너리 하나가 **모든 벤치마크 수치**를 옮겼고, §8.1 은 한 줄짜리 캐스팅이
  **모든 주기별 지표**를 오염시켰다.
- **위협 분석** — 같은 함정에 다시 걸릴 경로를 §8 이 목록으로 들고 있다.

살아 있는 상태는 `resume.md`, 할 일은 `todo.md`, 운영 규칙은 `CLAUDE.md`.

이 파일의 근거와 규약: 랩 위키 `guide/handoff-hygiene.md` §"세 번째 층". 문서 마커
(❌/⚠️/✏️/🧊)는 사람을 막고, 경로 격리(`archive/`)는 스크립트를 막고, 이 파일은
**새 세션이 폐기된 값을 다시 집어드는 것**을 막는다.

> ℹ️ **이름 (2026-09-03):** 랩 규약도 같은 날 `ruleout.md` → `quarantine.md` 로 개명했고
> (`wiki b5c8b76`·`1c4ab30`), 이유까지 같다 — "rule out"이 *배제한다*는 동작으로 읽혀 미해결과
> 혼동될 여지가 있다는 것. **새 프로젝트는 `quarantine.md` 를 쓴다**가 현행 권고이고,
> opuntia_analysis 도 이미 개명했다. 남은 혼재는 `ruleout.md` 를 유지하는 Sylvan 쪽이다.
> 그쪽에서 왔다면 같은 층, 같은 §1–§10 구조다.
>
> ❌ **SUPERSEDED (2026-09-03):** 이 자리에 있던 앞선 주석은 "규약은 `ruleout.md` 로
> 성문화돼 있고 이 저장소만 이탈한다"고 적었다. **틀렸다** — 확인 시점에 규약은 이미
> 바뀌어 있었고, 개명은 이탈이 아니라 현행 권고와 같은 방향이다.

항목 형식은 고정이다: **무엇 / 왜 / 언제 / 대체물 / 근거 경로**.
**대체물이 없으면 "없음"이라고 적는다** — 빈칸은 "나중에 채워도 된다"로 읽히지만
"없음"은 그 주장을 아예 하지 말라는 뜻이다.

```bash
grep -n '^### ' quarantine.md       # 격리 항목 전수
grep -rn 'quarantine.md' resume.md CLAUDE.md todo.md   # 이 파일로 연결된 지점
```

---

## §1 실행 세대 — 통째로 폐기된 것

### 1.1 Jul 16 배열 5912536 의 런타임 셀

- **무엇**: 발표된 BWTandem 런타임 7.4 h(human) / 6.6 h(maize) / 0.31 h(Col-CEN).
- **왜**: BED를 만든 실행(Jul 21 배열 5935102)과 **다른 실행**의 벽시계였다. Jul 16 배열의
  BED는 더 이상 존재하지 않는다. 로그의 행 수가 지문 역할을 해서 귀속이 확정됐다
  (5912536: hg38 3,364,455행 / 5935102: 3,994,477행).
- **언제**: 2026-08-05.
- **대체물**: 잡 5983792/93/94 — BED·런타임·메모리가 **한 실행**에서 나오고 환경이 로그에 찍힘.
- **근거**: `archive/2026-08-05-unreproducible/README.md` §1.

### 1.2 `d52a4ff` 이전 빌드로 낸 모든 수치

- **무엇**: 2026-07-09 이전에 산출된 모든 벤치마크 수치.
- **왜**: `align_accel.c`의 C 트레이스백이 **초기화되지 않은 힙**(`ptr_table`)을 읽어, C 경로와
  Python 경로가 무작위 반복 구간의 **31%에서 불일치**했다. 세 가지 원인이 겹쳐 있었고
  개별로 고치면 안 된다.
- **언제**: 2026-07-09 수정, 2026-07-10 3종 게놈 재측정.
- **대체물**: `d52a4ff` 이후 빌드로 재측정한 값 (`c6b92e7` + PR #6). 모든 **결론은 유지**됐고
  값만 이동했다 — 예: chr22 region recall 84.50 → **84.38**, precision 52.34 → **52.74**.
- **근거**: `docs/2026-07-09-nondeterminism-uninitialised-ptr-table.md`,
  `tests/test_align_parity.py`, `tests/test_align_unit_parity.py`.

### 1.3 Codex 2차 의견 리뷰 시도 1–3 (2026-09-03)

- **무엇**: `codex exec -m gpt-5.3-codex` 기반 리뷰 시도 3회.
- **왜**: (1) 포그라운드 10분 타임아웃에 killed, `tee` 대상 파일이 생성조차 안 됨.
  (2) `gpt-5.3-codex`가 ChatGPT 계정 미지원 모델로 **은퇴** → HTTP 400.
  (3) `--sandbox read-only`가 강제하는 번들 bubblewrap이 이 SLURM 세션 안에서 마운트
  네임스페이스를 못 만들어(`bwrap: Can't bind mount /oldroot/ on /newroot/`) Codex가
  **파일을 한 개도 읽지 못한 채 살아 있었다**.
- **언제**: 2026-09-03.
- **대체물**: `-m gpt-5.6-sol`(또는 `gpt-5.4`) + `--sandbox danger-full-access` + 스크래치패드의
  분리 워크트리. 시도 5가 이 조합으로 성공(24 KB, 20건).
- **근거**: `docs/2026-09-03-codex-review-findings.md` "Run history".
- ⚠️ **교훈**: 살아 있어 보이는 실행이 아무것도 읽고 있지 않을 수 있다. stderr에
  `exited 1 in 0ms`가 반복되면 즉시 중단할 것.

---

## §2 수치 — 값이 바뀐 것

| 대상 | 폐기값 | 대체값 | 왜 / 근거 |
|---|--:|--:|---|
| Col-CEN 비용 | 0.31 h / 1.31 GB | **0.51 h / 1.95 GB** | 잘못된 배열 + per-worker 메모리 |
| human 비용 | 7.4 h / 12.99 GB | **12.1 h / 21.86 GB** | 〃 |
| maize 비용 | 6.6 h / 14.37 GB | **15.4 h / 22.41 GB** | 〃 |
| chr22 region recall | 84.50 | **84.38** | align_accel 수정 (§1.2) |
| chr22 region precision | 52.34 | **52.74** | 〃 |
| Col-CEN CEN180 count | "1,380 → 21로 붕괴" | **1,380 (변화 없음)** | 스코어링 버그(§8.1), 붕괴는 허상 |
| maize SSR bp | "41 M → 131 M" | **41,314,898 bp** | 〃 |
| maize Table 3A-b | 41,007,584 / 1,498,365 / 35,017 | **41,314,898 / 1,508,821 / 34,247** | 3A는 갱신, 3A-b는 미갱신이었음 |
| gap-fill ablation base | — | **1,294 calls / 91.72%** | 스코어링 버그 재산출 |
| gap-fill ablation off/away | — | **991 / 85.38%, 993 / 85.57%** | 결론(약 5.8 pp)은 불변 |
| range-cost 비율 (human) | 1.30 / 1.30 / 1.41, 평균 1.34 | **1.82 / 1.77 / 1.77, 평균 1.79** | 발표값은 `07ad6fa` — 좁은 arm 이 버릴 탐색을 수행 (§6.17) |
| range-cost p100 런타임 | 6.61 / 6.58 / 5.50 h | **4.02 / 4.11 / 4.02 h** | 〃 |
| range-cost p2000 런타임 | 8.61 / 8.58 / 7.75 h | **7.31 / 7.28 / 7.13 h** | 〃 |
| range-cost p100 런타임 | 6.61 / 6.58 / 5.50 h (`07ad6fa`) | **4.02 / 4.11 / 4.02 h** (`0363d8b`) | 좁은 arm의 Tier 3 장주기 탐색이 실제로 제거됨(§6.17) |
| range-cost p2000 런타임 | 8.61 / 8.58 / 7.75 h (`07ad6fa`) | **7.31 / 7.28 / 7.13 h** (`0363d8b`) | 같은 릴리스 빌드 재측정(§6.17) |
| range-cost 비율 | 1.30 / 1.30 / 1.41 (평균 1.34) | **1.82 / 1.77 / 1.77 (평균 1.79)** | 20× 범위 확장은 여전히 선형 미만이지만 near-flat은 아님(§6.17) |

### 2026-09-03 적용분 — 원고 수정 완료 (커밋은 아래 §6 각 항목에)

| 대상 | 폐기값 | 대체값 | 항목 |
|---|--:|--:|---|
| Table 3C-b 캡션, 최소 밴드 손실 | tantan 0.92 "최소" | **TRF 0.87 이 최소** | §6.7 |
| Table 3C-b 캡션, BWTandem 대비 차이 | 14.36 | **14.88** | §6.7 |
| 밴드 손실 상계 | at most 1.5 / 2.8 | **1.54 / 2.81** | §6.8 |
| S4 (a) 최저 1:1 precision | BWTandem | **tantan (3.42 vs 8.68)** | §6.13 |
| Table 2 longdust 메모리 ×2 | 0.07 GB | **0.06** | §6.4 |
| Table 2 AniAnn's Col-CEN 메모리 | 0.50 GB | **0.48** | §6.4 |
| gap-fill 유효염기 게이트 | 70% | **80%** (catch-all 도 동일) | §6.21 |

**메모리 단위 주의**: 위 GB 표기는 `sacct MaxRSS`의 KiB를 1024²로 나눈 값이므로 **GiB**다.
원고·증거 트리 전반의 "GB" 표기는 §6.4 참조.

---

## §3 주장 — 방향이나 논리가 무너진 것

### 3.1 "발표된 BWTandem 수치 전부가 재현 불가" — 철회됨

- **왜**: 2026-08-05 1차 격리는 전면 격리였으나, 그 근거가 **내 스코어링 버그**(§8.1)였다.
  정확도 수치는 3종 게놈에서 전부 재현된다. 격리 범위는 두 번 축소돼 최종적으로
  **런타임 셀만** 남았다.
- **대체물**: `archive/2026-08-05-unreproducible/README.md` 현행판. **정확도는 격리 대상이 아니다.**

### 3.2 "각 발표 실행 로그의 호출 수가 BED와 불일치" — 틀림

- **왜**: Jul 16 배열 로그(`slurm_5912536_*.log`)와 Jul 21 파이프라인 로그(`bwt_hg38.log`)를
  한 파일처럼 읽은 오독. 파이프라인 로그는 3종 게놈 모두에서 BED와 내부 정합한다.
- **대체물**: 없음 — 이 주장은 하지 말 것.

### 3.3 "12.99 GB가 재현 실패" — 아님, 지표가 둘이다

- **왜**: `--threads`는 `ProcessPoolExecutor`로 **프로세스**를 띄우고, `/usr/bin/time -v`는
  `RUSAGE_CHILDREN` 즉 **자식들의 최댓값**을 보고한다. 12.99 GB는 2-worker 실행에서
  한 워커의 몫이며, 같은 방식으로 재측정하면 0.05%로 재현된다.
- **대체물**: 프로비저닝에 쓸 값은 `sacct MaxRSS`(동시 cgroup 총량).

### 3.4 "same binary 1.2.1" / "version-matched" (ULTRA p2000)

- **왜**: 두 바이너리는 **다른 파일**이다 — 컨테이너 638,912 B vs 로컬 573,672 B, 해시 상이.
  같은 것은 자기보고 **버전 문자열**뿐이다.
- **언제**: 2026-09-02 자체 감사(`24cd12a`), 검증 `ccf04dd`.
- **대체물**: "self-reported-version-matched", 또는 "invocation- and input-matched".
- ⚠️ **미완**: 원고는 고쳤지만 `results/manifest.tsv:67`에 원문이 남아 있다 → §6.2.

### 3.5 "the matched measurement" (환경 매칭) — 아님

- **왜**: Results 3.1이 "matched"라 쓰고 Methods 2.2는 "not environment-matched"라 인정해
  같은 원고 안에서 모순이었다. 컨테이너 vs 로컬 micromamba 실행이다.
- **대체물**: "period-matched, not environment-matched".

### 3.6 "every competitor cost cell is checkable" — 거짓 (2026-09-03 확정)

- **왜**: Table 2의 TRASH 행이 반증한다(§6.1). `a38201b`가 이 주장을 넣었고, 오늘 Codex 리뷰가
  반례를 찾았으며 내가 집합 연산으로 재현했다.
- **대체물**: "human GRCh38의 다섯 도구(ULTRA/TRF/tantan/TRASH/mreps) 비용 셀은 생존 로그와
  일치한다" — 검산 완료된 범위로만 한정할 것.

### 3.7 "the arrays it recovers lie inside regions the other phases already touch"

- **왜**: 해당 모드가 "every run reported here"에서 **비활성**이라고 같은 문장이 말한다.
  근거가 되는 실행이 보고되지 않았다. (Kimi 라운드 1, manuscript.md:62)
- **대체물**: 없음 — 수치를 부록에 넣거나 문장을 삭제. **미적용 상태.**

### 3.8 catch-all 식별도 0.72 를 평가 벤치마크 위에서 골랐다

> ✏️ **PARTIAL (2026-09-03, `d580840`, 저자 결정 (a)):** 본문 S3 에 "0.72 는 in-sample 운영점" 명시 (커밋은 아래 §3.9 와 동일). **원장 예치는 미완** — `results/` 를 건드리므로 A-2 재해시와 함께 처리한다.

- **무엇**: 운영점 `CATCHALL_MIN_IDENTITY=0.72`.
- **왜**: 부록 S3가 스스로 밝힌다 — 스윕이 "scored with the Table 1a pipeline"으로,
  Table 1a가 보고하는 **바로 그 adotto 카탈로그**에 대해 채점됐고, 선택된 arm은
  "equals the Table 1a BWTandem row **by construction**"이다. 즉 헤드라인 human 수치는
  평가셋 위에서 수행한 파라미터 탐색의 최대-recall 지점이다. 저장소의 자체 규칙은
  "no GT overfitting"이다.
- **완화 요인**: 곡선이 매끄럽고(0.76→0.72에서 recall 71.53→80.53) 공개도 이례적으로 솔직하다.
  스파이크 집기는 아니다.
- **대체물**: 없음 — 다만 **선택과 무관한 증거 1개**가 필요하다(maize/Arabidopsis 진실셋이나
  보류된 human 염색체 부분집합에서 0.72 확인). 부록이 아니라 **본문**에 선정 경위를 밝힐 것.
- **근거**: `manuscript.md:508, 475, 479`. (Kimi 라운드 3, R3c#2 — 미조치)

### 3.9 "selected empirically on chromosomes 21 and 22" — 캠페인 규모가 공개되지 않았다

> ✏️ **PARTIAL (2026-09-03, `d580840`, 저자 결정 (a)):** 본문 2.2.3 에 캠페인 규모(44건), 코드변경 포함, 관찰 후 accept/reject, ULTRA 표적, 22번 염색체는 post-selection validation 임을 명시. **원장 예치 미완**.

- **무엇**: 설정 선택 경위에 대한 원고의 유일한 서술.
- **왜**: 외부 원장 `exp1_human/loop/ledger.tsv` 에 **채점된 설정 평가 44건**이 있고(헤더 포함 45행),
  환경 스윕과 **코드 변경**까지 포함한다. 계획 문서는 이를 "adaptive" coordinate descent라 부르고,
  `best.json` 은 목표를 이렇게 적는다:
  **`recall->82% at precision>=~53% genome (ULTRA-level)`**.
  즉 비공식적 선택이 아니라 **평가에 쓰인 바로 그 카탈로그에 대한, 경쟁 도구를 기준점으로 삼은
  최적화 캠페인**이다. 보류 염색체 결과는 유용하지만 다중성·목표 변경·적응적 중단·ULTRA 표적을
  공개하지 않는다. 원장과 `best.json` 은 **저장소 밖에 있다.**
- **대체물**: 평가 설정 수, 코드변경 대 파라미터전용 구분, 최초·수정 목표, 경쟁도구 표적,
  수용/중단 규칙, 최종 선택 규칙을 명시하고 원장을 예치할 것. 22번 염색체 결과는
  **선택 후 검증(post-selection validation)** 으로 서술할 것. **미조치.**
- **근거**: `manuscript.md:94`; `exp1_human/loop/{ledger.tsv,best.json}`;
  `docs/superpowers/plans/2026-06-23-exp1-recall-loop.md`. (Codex 라운드 1, 발견 2)
- 🔗 [[3.8]] 과 같은 뿌리다 — 3.8은 임계값 하나, 이것은 캠페인 전체.

### 3.10 Arabidopsis 설정(catch-all off)을 평가 진실셋을 보고 골랐다

> ✏️ **PARTIAL (2026-09-03, `d580840`, 저자 결정 (a)):** S2 에 두 설정을 CEN180 진실셋에 채점한 수치(65.54 vs 60.72, recall 99.67/99.68)와 보류셋 부재를 명시.

- **왜**: S2는 Arabidopsis만 catch-all을 끈다고 하고, 이유를 "where it costs base pair precision
  without recovering centromeric arrays"라 적는다. 생존한 실험 기록
  (`exp1_human/filip_repro/catchall_experiment_results.md`)은 **같은 `colcen_cen180.bed` 진실셋에
  두 설정을 채점**하고 "Keep catch-all OFF for Col-CEN"으로 결론짓는다:
  off = recall 99.67 / bp precision 65.54, on = 99.68 / 60.72.
  즉 제출된 Table 2 설정은 **목표 평가 지표를 본 뒤** 불리한 precision을 피하려고 선택됐다.
  보류된 Arabidopsis 위성 진실셋은 없고, catch-all-on 의 현재 빌드 행도 예치돼 있지 않다.
- **대체물**: 없음 — Col-CEN을 튜닝 데이터로 선언하고 보류셋에서 평가하거나, 최소한 Table 2에
  두 설정을 모두 보고하고 in-sample 임을 표시할 것. **미조치.** (Codex 라운드 2 발견 2)
- 🔗 [[3.8]] [[3.9]] 와 같은 계열 — 세 게놈 중 둘에서 같은 패턴이 확인됐다.

---

## §4 명명 — 쓰면 안 되는 이름

### 4.1 `CATCHALL_MAX_SEEDS`

- **왜**: CLAUDE.md에 기본값까지 문서화돼 있었으나 **어떤 소스도 읽지 않는** 유령 노브였다.
  이 사고가 `tests/test_env_var_docs.py`를 만들게 했다.
- **대체물**: 없음.
- ⚠️ **CLAUDE.md에는 이 식별자를 쓰지 말 것** — 그 파일에 이름이 등장하면 스캔이 다시
  실존 노브로 간주한다. 실명은 **이 파일에만** 둔다.

### 4.2 쓰면 안 되는 표현

| 표현 | 문제 | 대체 |
|---|---|---|
| "same binary" / "version-matched" | 자기보고 문자열뿐 | "self-reported-version-matched" |
| "GB" (KiB/1024² 값에) | 실제는 GiB | "GiB", 또는 10⁹ 나눗셈일 때만 GB |
| "per-copy records" (BWTandem 출력) | 실제는 배열 단위 + 집계 copy count | "per-array structured records" |
| "prohibitively slow" (TRF) | 측정치는 33.7 h — 107.6 h TRASH를 수용한 논문에서 성립 안 함 | 측정값 그대로 |
| "rerun infeasible" / "caps could not be lifted" | 저자의 취소 결정을 도구의 속성으로 전환 | "저자가 예산 내에서 중단했다" |
| "preregistered" | 등록처·타임스탬프 없음 | 위치를 밝히거나 삭제 |

---

## §5 식별자 — 틀린 ID

| ID | 상태 | 주의 |
|---|---|---|
| 5912536 (Jul 16 배열) | 런타임 출처, **BED 출처 아님** | §1.1 |
| 5935102 (Jul 21 배열) | 발표 BED의 출처 | 런타임은 이 배열 것을 쓸 것 |
| 5983792/93/94 | 현행 비용 정본 | BED·시간·메모리 동일 실행 |
| 6146742 | **크래시한** 2026-tool 채점 실행 | 예치 전용, 인용 금지 |
| 6147179 | 그 재실행 (정본) | |
| 6145581 | ULTRA p2000, **저자 취소** 1 d 22 h 15 m | 완주 아님 |
| 6076847 | TRF p2000, **저자 취소** 6 d 13 h 57 m | 〃 |

---

## §6 결과 셀 — 개별 폐기 (2026-09-03 확정, **수정 미적용**)

> ⚠️ 아래 25건은 (§6.6·§6.16 은 오탐으로 철회, 기록은 보존) Codex 리뷰 4회(6.1-6.6, 6.17-6.26)와 Kimi 라운드 2·3(6.7-6.16)이 지적하고 **내가 실물 파일·표로 재현**한 것이다.
> 원고·증거 트리에 **아직 그대로 있다**. 인용 금지, 수정 대상.

### 6.1 Table 2 의 TRASH 비용 셀 두 개

> ✅ **APPLIED (2026-09-03, 다음 커밋):** 계약 결정 — **591행 파일을 폐기하고 진짜 template 전용
> 합집합(397행)을 만들어 새로 채점했다.** 비용 셀만 고치는 안은 라벨이 여전히 거짓이었다:
> 591행은 de novo 를 포함하므로 de novo 행이 template 행의 **부분집합**이었다.
>
> ```
> CEN159 232 + CEN178 238 = 470 -> sort -u -> 397   (겹침 73)
> 채점: 397 regions | cov 85.03% | CEN180 115 | bpPrec 86.40% | recall 99.69%
> 비용: 6:18:29 + 25:22:40 = 31:41:09 = 31.69 h, 최대 RSS 2.63 GiB
> de novo 행은 자기 로그로: 5:47:30 = 5.79 h, 1.29 GiB
> ```
>
> 예치: `results/regen/colcen_trash_template_union.bed` + `_scored.txt`.
> 매니페스트 행 53 이 그 파일을 가리키도록 바꿨다. Table 2 두 행과 캡션 교체.
> unfiltered 정확도 셀 3종이 동일한 것은 같은 중심절 서열을 덮기 때문이고, 구분되는 것은
> region 수와 CEN180 count 다 — 캡션에 그대로 적었다.

- **무엇**: `TRASH (template)` 와 `TRASH (de novo)` 행의 `6.31 h / 1.32 GiB`.
- **왜**: `Col-CEN_v1.2_trash.bed`(591행)는 세 실행의 **정확한 합집합**이다 —
  CEN159(232) + CEN178(238) + denovo(234) = 704행 → 중복 제거 591행, 양방향 차집합 0.
  그런데 두 행 모두 **CEN159 단독**의 비용을 쓴다(로그 6:18:29 → 6.3081 h,
  1,389,316 KiB → 1.3249 GiB). de novo 행은 자기 로그(5:47:30 / 1,357,252 KiB)조차 안 쓴다.
  234행 파일은 591행 파일의 **부분집합**이므로 각주의 "distinct deposited output files"도
  독립성을 오도한다.
- **대체값**: 3개 순차 합계 **37:28:39**, 최대 RSS **2,761,680 KiB = 2.63 GiB**;
  de novo 단독은 **5:47:30 / 1.29 GiB**.
- **근거**: `/data/gpfs/assoc/pgl/filip/bwtandem_results/benchmarking_results/trash/logs/Col-CEN_v1.2_{CEN159,CEN178,denovo}_run.log`.
- ⚠️ 각주의 "재검증 불가한 상속값"이라는 면제 사유는 `a38201b`가 로그를 찾아낸 시점에
  **무효**가 됐다.

### 6.2 `results/manifest.tsv:67` — 정정 이전 문장이 살아 있음

> ✅ **APPLIED (2026-09-04):** `manifest.tsv:67` 의 네 주장을 전부 교체했다 — 서로 다른 바이너리(638,912 vs 573,672 bytes)와 자기보고 버전 문자열, 하한으로 표시한 1.55×, 헤더만 다른 동일 서열, 그리고 정체가 아니라 **마지막 출력 위치**가 중심절에 있다는 것. 잘린 마지막 레코드(§6.17 계열)도 함께 적었다.

- **무엇**: "same binary 1.2.1", "version- and input-matched", 하한 표시 없는 "1.55×",
  5시간 정체에 붙은 "(chr1 centromeric satellite region)".
- **왜**: `24cd12a`·`e8ac50c`·`67a3f42`가 원고에서 고친 네 가지가 증거 인덱스에 남아
  **원고와 매니페스트가 모순**이다.
- **대체값**: 정정된 원고 문구.

### 6.3 AniAnn's 스레드 수

> ✅ **APPLIED (2026-09-03, 다음 커밋):** 매니페스트 행 123·126·129·132 를 `2` 로 수정하고, `threads` 열이 **SLURM 할당**을 뜻한다는 계약을 `results/README.md` 에 명문화했다 (열 분리는 스키마 마이그레이션이라 하지 않음; tantan 의 활성 스레드 1 은 예외로 명시).

- **무엇**: `results/manifest.tsv` 행 123·126·129·132 의 `threads=1`.
- **왜**: sbatch와 모든 실행 로그가 `-j 2`, `--cpus-per-task=2`.
- **대체값**: `2` (+ `manifest.sha256` 재해시 — §8.2 순서 준수).
- **근거**: `exp1_human/tools2026/runs/run_anianns.sbatch:5,22,23`.

### 6.4 2026 도구 메모리 셀 — 같은 표 안에서 규칙 불일치

> ✅ **APPLIED (2026-09-03, 다음 커밋):** `results/comparators2026/README.md` 의 `0.07`→`0.06`, `0.50`→`0.48`, `GB`→`GiB` 까지 반영 — 원고만 고쳤던 부분 수정을 닫았다.

> ✅ **APPLIED (2026-09-03, `2901ff3`):** Table 2 의 두 셀을 0.06 / 0.48 로 교체.
>
> ✅ **APPLIED (2026-09-03, `63ccd07`):** 단위 라벨 전면 교체 — 메모리 `GB` **37곳**을 `GiB` 로 바꾸고 규약 문장을 변명에서 정의로 고쳤다. **숫자는 하나도 바꾸지 않았다.** 84행의 `3.15 GB` 만 메모리가 아니라 파일 크기여서 §6.5 로 따로 처리했다.

| 셀 | 로그 KiB | GiB 규칙 | 현재 표기 |
|---|--:|--:|--:|
| AniAnn's human | 550,700 | 0.53 | 0.53 ✅ |
| AniAnn's maize | 572,588 | 0.55 | 0.55 ✅ |
| AniAnn's Col-CEN | 503,968 | **0.48** | 0.50 ❌ |
| longdust Col-CEN | 66,020 / 66,340 | **0.06** | 0.07 ❌ |

- **근거**: `exp1_human/tools2026/runs/*.time`. 표기 위치 `results/comparators2026/README.md:27-28`.

### 6.5 "roughly 5% more sequence" (manuscript.md:84)

> ✅ **APPLIED (2026-09-03, `63ccd07`):** 실측 대체 — 경쟁 3,209,286,105 염기 vs BWTandem 3,088,269,832 → **+3.92%**(ACGT만 +3.80%). "roughly 4%"로 교체하고 파일 바이트(3.15 GB)를 염기 수로 바꿨다.

- **왜**: 원고가 주는 두 수는 3.25 Gb와 3.15 GB → **+3.2%**. 게다가 하나는 염기,
  하나는 파일 바이트로 두 문장 만에 단위가 섞였다.
- **대체값**: 두 값을 같은 단위(염기)로 다시 재고 비율을 재계산할 것. 확정값 **미산출**.

### 6.6 ~~stride 포화점 "roughly 24 Mb"~~ — ❌ **철회 (2026-09-03): 오탐이었다**

> ❌ **WITHDRAWN (2026-09-03):** 이 항목은 "원고의 `n / 40,000` + clamp 300 이면 포화점은
> 12 Mb 인데 원고가 24 Mb 라 적었고 어느 프리셋 배수도 24 Mb 를 주지 않는다"고 적었다.
> **틀렸다.** `tier3.py:49-50` 의 커버리지 축소를 빼고 계산했다 — 그 축소는 clamp **앞에**
> 적용되고 최대 ×0.5 이므로, 커버리지가 무엇이든 **n ≥ 24 Mb 면 반드시 포화한다.**
> 원고의 "saturates at 300 on every chromosome above roughly 24 Mb" 는 임계값이 아니라
> **보장 조건**이며 그대로 참이다.
>
> ```
> c=0.0 → ×1.000 → 12.0 Mb        c=0.9 → ×0.550 → 21.8 Mb
> c=0.5 → ×1.000 → 12.0 Mb        c=1.0 → ×0.500 → 24.0 Mb
> ```
>
> **잔여(LOW)**: 원고가 프리셋 배수(`fast 1.6 / balanced 1.0 / sensitive 0.4`)를 언급하지
> 않는다. 실행은 CLI 기본값 `balanced`(=1.0)이라 보고된 수치에는 영향이 없다.
>
> **경위**: Kimi 라운드 1이 지적했고 내가 "코드로 확인했다"고 보고했으나, 재현한 것은 산술의
> 절반이었다. §3.1 과 같은 종류의 오탐이다 — 그래서 이 파일이 있다.
### 6.7 Table 3C-b 캡션 — 표가 반증하는 최상급과 틀린 차이값

> ✅ **APPLIED (2026-09-03, `2901ff3`):** 캡션을 TRF 0.87 최소 · 14.88 로 교체.

- **무엇**: "the smallest banded loss of the four" (tantan −0.92) 와 "finishes 14.36 points
  above BWTandem".
- **왜**: 같은 표의 손실 열이 TRF **−0.87** 로 tantan보다 작다 — 최상급이 거짓이다.
  차이값도 55.59 − 40.71 = **14.88** 이지 14.36이 아니며, 표의 어떤 조합도 14.36을 주지 않는다.
- **대체값**: "TRF loses least (0.87); tantan 0.92" / **14.88**.
- **근거**: `manuscript.md:298` 캡션 vs `manuscript.md:306-309` 표.
- ⚠️ Table 3B-b(254행)의 같은 문형은 참이다(0.50 < 0.63 < 1.20) — 재검산 없이 옮겨 쓴 흔적.

### 6.8 "at most 1.5 / 2.8 points" — 하한을 반올림해 만든 거짓 상계

> ✅ **APPLIED (2026-09-03, `2901ff3`):** 1.54 / 2.81 로 교체.

- **왜**: TRF의 최대 손실은 **1.54**, ULTRA는 **2.81**. "at most"는 정확한 상계이므로
  내림 반올림은 허용되지 않는다.
- **대체값**: 1.6 / 2.9, 또는 1.54 / 2.81 을 그대로.
- **근거**: `manuscript.md:280`.

### 6.9 ~~Table 1c 캡션이 네 도구, 표는 세 행~~ — ❌ **철회 (2026-09-03): 오탐이었다**

> ❌ **WITHDRAWN (2026-09-03):** "캡션은 네 도구가 진입했다는데 표에는 세 행뿐"이라 적고
> 캡션에 "TRASH's stratum row is not reported below" 라는 문장까지 넣었다. **틀렸다.**
> TRASH 행은 `0357cbd` 시점부터 표에 있었다 —
> `git show 0357cbd:manuscript.md` 로 확인: `| TRASH | 1,427 | 0.25 | 90.75 | 5.23 | 16.42 |`.
> 나는 표의 141–145행만 읽고 그 아래를 보지 않았다.
>
> **내가 넣은 문장이 거짓이었으므로 삭제했다.** (Kimi 재검토 RCK2-2 가 잡았다.)
>
> **경위**: §6.6·§6.16·§6.27 에 이은 **네 번째** 같은 실수다. 앞선 셋은 "부재를 오류로
> 승격"이었고 이번은 **불완전한 읽기를 부재로 오인**한 것 — 뿌리는 같다.
> 판정 규약(`todo.md` B-0)이 요구하는 "낯선 사람이 재실행할 수 있는 검사"를 이 항목에는
> 적지 않았다. 적었다면 `sed -n '141,150p'` 가 아니라 표 전체를 봤어야 한다는 게 드러났을 것이다.

### 6.10 >100 bp 층 서술 — TRF가 이기는 두 열만 빠짐

> ✅ **APPLIED (2026-09-03, `2901ff3`):** 본문 107행에 TRF 의 region/bp precision 추가.

- **왜**: 본문 107행은 BWTandem 3.43% / 13.39% / precision 64.60% 을 쓰고 "TRF reaches 2.02%
  and 11.83%" 로 **recall 두 개만** 옮긴다. 같은 표에서 TRF의 region precision은 **97.18%**
  (BWTandem 64.60%), bp precision도 **31.19 vs 29.25** 로 TRF가 앞선다.
- **대체값**: 같은 문장에 TRF의 precision 두 값을 넣을 것.
- **근거**: `manuscript.md:107` vs `manuscript.md:141-145`.

### 6.11 Table 2 — 자기 도구의 최저 수치만 서술되지 않음

> ✅ **APPLIED (2026-09-03, `2901ff3`):** Table 2 캡션에 BWTandem 최저 bp precision 57.08% 와 그 사유 명시.

- **왜**: BWTandem의 CEN180 bp precision **57.08%** 는 11행 중 **최저**다(차하위 ULTRA 58.08).
  캡션은 NCRF의 0, tantan의 공유 노드, AniAnn's의 높은 커버리지, longdust의 구조적 0까지
  일일이 설명하면서 이 한 줄만 없다.
- **대체값**: 최저라는 사실과 그 이유를 한 문장으로 명시.
- **근거**: `manuscript.md:187-199`.

### 6.12 TR-1 경계 오차 — 같은 값이 771 bp 이자 4.3 kb

- **왜**: 234행은 "50.34% at a 771 bp mean boundary offset", 236행은 같은 재생성 파일의
  병합 전 상태를 "50.34% ... 4.3 to 19.8 kb offset" 이라 한다. 커버리지가 동일한데 병합 전
  평균 오차가 둘일 수 없다. 캡션(240행)은 또 banded 4,265 bp / 7,973 bp 를 제시한다.
- **대체값**: 네 값(771 bp / 4.3 kb / 4,265 bp / 7,973 bp)이 각각 어떤 통계인지 정의하고
  일치시킬 것. **미산출.**
- **근거**: `manuscript.md:234, 236, 240`.

### 6.13 Supplementary Table S4 캡션 — 자기 표가 반증하는 "최저"

> ✅ **APPLIED (2026-09-03, `2901ff3`):** S4 캡션을 tantan 최저 · BWTandem 차하위로 교체.

- **무엇**: "**BWTandem's is the lowest of the five tools**" (panel (a) 1:1 precision).
- **왜**: 같은 표에서 tantan **3.42%** < BWTandem **8.68%**. BWTandem은 **차하위**다.
- **대체값**: "second-lowest", 또는 "lowest after tantan (8.68% vs 3.42%)".
- **근거**: `manuscript.md:518` 캡션 vs `manuscript.md:521-529` 표.
- ⚠️ **오류 방향이 자기에게 불리하다** — 그래서 "우리한테 엄격했네" 하고 자체 점검을 통과했다.
  자기비판 방향의 오류도 오류다. S4 수치를 산문으로 옮긴 다른 곳도 전수 재확인할 것.

### 6.14 "All competitor runs were executed inside one Singularity container"

> ✅ **APPLIED (2026-09-03, `2901ff3`):** 'All competitor runs except the two noted below' 로 교체.

- **왜**: 같은 표의 두 행이 반증한다 — longdust "not run in the container"(487행),
  ULTRA p2000 "executed outside the Singularity sandbox"(489행).
- **대체값**: "All competitor runs except the two noted below".
- **근거**: `manuscript.md:481` vs `487, 489`.
- ✅ **Codex(발견 10)와 Kimi 라운드 3(R3c#3)이 독립적으로 같은 결함을 짚었다.**

### 6.15 TRASH 분류가 문장마다 바뀐다 — 비교 주장을 보호하는 방향으로

> ✅ **APPLIED (2026-09-03, `5345e5b`):** TRASH 는 de novo 모드가 있으므로 **de novo 로 분류**한다. 326행의 주장을 "largest de novo footprint **on human**" 으로 좁히고, TRASH 의 maize 120.64 GiB 를 template-guided 괄호에서 꺼내 본문에 명시했다.

- **왜**: 313행은 TRASH를 period 범위를 미리 정할 필요가 없는 **de novo** 검출기의 예외로 든다.
  326행은 BWTandem이 "largest de novo-tool footprint on the two large genomes"라 주장하면서
  더 큰 수치를 "template-guided"로 밀어낸다 — 그 안에 **TRASH 120.64 GB (maize)** 가 들어 있다.
  313행이 맞으면 326행은 거짓이다.
- **대체값**: 두 곳의 분류를 일치시킬 것. TRASH가 de novo라면 주장을 "largest de novo footprint
  **on human**"으로 한정해야 한다. **미결정.**
- **근거**: `manuscript.md:313` vs `manuscript.md:326`.

### 6.16 ~~stride 축소의 432행 vs 440행~~ — ❌ **철회 (2026-09-03): 오탐이었다**

> ❌ **WITHDRAWN (2026-09-03):** "발동한다(432행)와 무력하다(440행)가 모순"이라 적었으나
> **모순이 아니다.** 440행은 "**above roughly 24 Mb** 에서 무력"이라 스스로 조건을 달고
> 있고, 432행의 Arabidopsis 염색체는 그 아래(약 22 Mb)다. 22 Mb 염색체가 포화를 면하려면
> `(1-0.5c) < 12/22 = 0.545` → **c > 0.909** — 432행의 "~90%" 가 바로 그 값이다.
> 트리거 50% 와 "90% 초과" 도 서로 다른 질문의 답이다: 축소가 **적용되는** 조건과
> 축소가 clamp 를 **이기는** 조건.
>
> **경위**: §6.6 과 같은 계산 착오에서 파생됐다. 하나가 오탐이면 그 위에 세운 것도 무너진다.
### 6.17 "The tier's search window is fixed at 100 bp to 100 kbp" — 옛 구현을 서술한다

> ✅ **APPLIED (2026-09-04, `da9e484`):** C-1 이 이것을 실증했다. `0363d8b` 재측정 결과 비율이 **1.34 → 1.79** 로 올랐고, 움직인 것은 좁은 arm(6.6 → 4.0 h)이다 — `07ad6fa` 에서 그 arm 이 장주기 탐색을 하고 버리고 있었다는 뜻이다. 원고·Figure 2·매니페스트에 반영, 옛 행은 superseded 로 보존. 증거: `results/range_cost_0363d8b/`.

> ✅ **APPLIED (2026-09-04):** `0363d8b` 재측정 3쌍(잡 6147698/99/700)을 원고·매니페스트·Figure 2 소스에 적용하고, Methods 2.1.4와 S1.3을 현행 Tier 3 범위와 jitter tolerance로 교체했다. 기존 `07ad6fa` 실행은 삭제하지 않고 superseded로 표시했다. 렌더링은 별도로 남았다.

- **왜**: 현재 코드와 헤드라인 커밋 `0363d8b` 모두 `tier3_min = max(100, min_period)` /
  `tier3_max = min(100_000, max_period)` 로 CLI 인자가 **탐색 자체를 좁힌다**. 소스 주석이
  스스로 **"This is NOT output-neutral"** 이라고 적어 두었다. 부록의 jitter tolerance 0.04도
  틀렸다 — `tolerance_ratio = 0.02 + 0.02*(max_period/100000)` 이므로 `--max-period 2000`에서
  **0.0204**다. 고정 100–100,000 동작은 역사적 커밋 `07ad6fa`의 것이고,
  **range-cost 실험 3쌍(`manifest.tsv:60-65`)이 바로 그 옛 커밋으로 돌았다.** 그 좁은 arm은
  장주기 탐색 비용을 내고 결과를 버려, 이전 비율을 도구에 유리하게 낮춰 보였다.
- **대체값**: `0363d8b`에서 p100은 4.02 / 4.11 / 4.02 h, p2000은 7.31 / 7.28 / 7.13 h였고,
  비율은 **1.82 / 1.77 / 1.77(평균 1.79, 범위 1.77–1.82)**였다. 이전 `07ad6fa`의
  1.30 / 1.30 / 1.41(평균 1.34)은 superseded다. 좁은 arm은 약 6.6 h→4.0 h,
  넓은 arm은 8.6 h→7.3 h로 줄었다.
- **근거**: `results/range_cost_0363d8b/README.md`, 잡 6147698/99/700 로그,
  `bwtandem/finder.py:116-139`, `bwtandem/tier3.py:38`.
- ⚠️ 20× 범위 확장에 1.79× 런타임이므로 선형 미만이지만, near-flat으로 서술하지 않는다.

### 6.18 "shared 100 bp range" 가 비대칭이다

> ✅ **APPLIED (2026-09-03, `05a76a6`):** 저자 결정 C-2(a) — Table 1b 의 BWTandem 행을 **사후필터 행**(3,911,182 / 78.87 / 50.13 / 30.00 / 53.89)으로 교체해 TRF·TRASH 와 같은 방식으로 통일. 네이티브 재실행(79.88 / 50.62)은 민감도 분석으로 강등. 캡션도 재작성 — 옛 캡션은 "두 실행 범위가 다른 콜 세트를 낸다"는 이유로 네이티브를 썼는데, 그 이유가 곧 비대칭의 근거였다.

- **왜**: BWTandem은 `--max-period 100` **네이티브 재실행**을 받았고, TRF는 500으로 돌린 뒤
  **사후 필터**됐다. TRF에도 네이티브 최대 주기 인자가 있으므로 피할 수 있는 비대칭이다.
  같은 스코어러로: 사후필터 BWTandem **recall 78.87 / prec 50.13**, 네이티브 BWTandem-F
  **79.88 / 50.62** — 네이티브 재실행이 **+1.01 pp**를 준다. §6.17이 최대 주기가 탐색을
  바꾼다는 걸 확정했으므로 이 차이는 필터링 효과가 아니다.
- **대체값**: 공통 비교에 사후필터 행을 쓰거나, 모든 도구를 네이티브 100으로 재실행할 것.
- **근거**: `results/manifest.tsv:9-12,89-93`, `results/regen/score_table1_p100.txt:40-57`.

### 6.19 Availability 약속이 거짓이다 — 스코어러 3개가 저장소에 없다

- **왜**: 원고는 "Source code, scoring scripts and per-figure provenance"를 약속하고 독자가
  "rerun the same code" 할 수 있다고 한다. 그러나 **도달 가능한 git 이력 전체에** 다음이 없다:
  `scripts/scoring/rescore_tables_3bc.py`(매니페스트가 Table 3B/3C용으로 지목),
  `scripts/scoring/score_exp3.py`(maize 래퍼들이 요구),
  `scripts/score_overlap.py`(`score_table1.py --check-defs`가 요구).
  provenance JSON은 대신 **사설 클러스터 경로**를 가리킨다.
- **대체값**: 해당 파일과 의존물을 커밋하고 경로를 저장소 상대경로로 바꿀 것.
- **근거**: `git log --all --diff-filter=A` 결과 0건; `results/regen/table3bc_provenance.json:3,13`.

### 6.20 헤드라인 BWTandem 메모리 3개가 예치 증거에서 유도되지 않는다

- **무엇**: `28.08` / `28.45` / `2.16 GB` (human / maize / Col-CEN, sacct MaxRSS라고 서술).
- **왜**: 래퍼는 `wait4()`의 프로세스 트리 `ru_maxrss`만 기록하며, 예치된 값은
  **17.369 / 19.968 / 1.307 GiB** 로 환산된다. 벽시계도 각각 약 33 / 29 / 5초 짧다.
  **원시 `sacct` 출력이 저장소에도 외부 로그에도 남아 있지 않다** — sbatch는 나중에 직접
  `sacct`를 돌리라고 지시할 뿐이다. 즉 정확한 cgroup 셀과 그 반올림을 독립 검증할 수 없다.
- **대체값**: `sacct -j ... --format=JobID,Elapsed,MaxRSS,State -P` 원문을 체크섬과 함께 예치하거나,
  예치된 `wait4` 값을 **그 이름으로** 쓸 것.
- **근거**: `results/regen/regen_{human,maize,colcen}.provenance.json:10-11`,
  `scripts/benchmark/run_with_provenance.sh:57-72`.

### 6.21 "at least 70% unambiguous A/C/G/T bases" — 코드는 80%다

> ✅ **APPLIED (2026-09-03, `2901ff3`):** Methods 와 S1.4 를 80% 로 교체하고 catch-all 적용을 공개.

- **왜**: `bwtandem/autocorr.py:43` `DEFAULT_MIN_VALID_FRAC = 0.8`, 후보 블록과 각 분할 창
  양쪽에 적용된다. 벤치마크 커밋 `0363d8b`에서도 같다. **catch-all 패스에도 같은 미공개
  80% 게이트가 걸린다.** 이것은 문서 오탈자가 아니라 **활성 수용 임계값**이다.
- **대체값**: 원고를 **80%**로 고치고 catch-all 적용 사실을 공개할 것 (코드를 70%로 바꾸면
  영향받는 모든 결과를 재실행해야 한다).
- **근거**: `manuscript.md:72,448` vs `bwtandem/autocorr.py:43`, `finder.py:541,622`.

### 6.22 "the 2026 tools have no banded rows" — 예치 로그가 반증한다

> ✅ **APPLIED (2026-09-03, `05a76a6`):** 저자 결정 C-9(a) — 사전등록대로 AniAnn's 를 밴드 패널에 넣었다. **셀은 이미 예치돼 있었다**: human ≤100 `114 / 0.01 / 85.09 / 0.92 / 54.50`, 101–2000 `205 / 0.03 / 81.46 / 2.69 / 10.03` (`score_2026_human.txt:25-33`), Col-CEN 밴드 recall `99.15%` 양쪽 (`score_2026_colcen.txt:12`). Table 1b·1c 에 행 추가, §3.2 에 서술 추가, Table 1b 캡션의 "2026 tools have no rows here" 는 거짓이므로 정정.

- **왜**: 사전등록 프로토콜(`docs/2026-09-01-longdust-anianns-benchmark-protocol.md`)은 AniAnn's가
  "**period-banded rules using its periodicity column where a table's rule reads a period**"에
  채점 가능하다고 **결과가 나오기 전에** 적었다. 원고는 결과가 나온 뒤 이를 뒤집는다.
  그런데 예치된 Col-CEN 로그 12행이 AniAnn's 밴드 recall **99.24 / 99.15 / 99.15%** 를 들고 있고,
  같은 로그의 BWTandem은 **99.73 / 91.92 / 95.21%**, TRF는 **97.55%** 다 —
  **두 밴드 열에서 AniAnn's가 BWTandem과 TRF를 모두 앞선다.** §3.2는 BWTandem·ULTRA·TRF의
  밴드 recall만 보고한다.
- **대체값**: 사전등록대로 AniAnn's를 모든 해당 패널에 넣거나, array-level period가 무효라면
  **CEN180 Count를 포함한 모든 period 기반 AniAnn's 셀을 철회**하고 등록 후 변경을 명시할 것. **미결정.**
- **근거**: `results/comparators2026/score_2026_colcen.txt:12` vs `manuscript.md` §3.2.
- ⚠️ **이것이 가장 위험한 유형이다** — 저장소에 있는 불리한 수치가 원고에 없다.
  (Codex 라운드 2 발견 1. 인용된 BWTandem 값 91.31/94.68 은 내가 본 로그 행의 91.92/95.21 과
  다르다 — 순서 관계만 검증됐고, 정확한 대조군 수치는 재확인이 필요하다.)

### 6.23 "TRF's period column is empty because its BED's motif column holds the full array sequence"

> ✅ **APPLIED (2026-09-03, 저자 결정 C-10):** 스코어러를 고치고 TRF 만 재채점했다.
> `load()` 에 col5=="period" 분기, `period_of()` 가 명시 period 우선. 재채점 결과 **period 5개 필드만**
> 바뀌고 나머지 30개 필드는 예치본과 동일 — matched 518,719, sensitivity 29.06, precision 53.87, 층화 전부 불변.
> **TRF period-exact = 63.53%** (Codex 보고치와 일치). 새 순위: tantan 73.54 > **TRF 63.53** >
> ULTRA 59.61 > BWTandem 58.66 > TRASH 33.41 — BWTandem 이 3위/3에서 **4위/5**가 된다.
> 원고 S4 표·캡션 반영 완료. **예치 JSON 교체는 `results/` 라 A-2 재해시와 함께.**
>
> ⚠️ **이 버그가 살아남은 이유**: 통과하는 테스트가 옛 동작을 고정하고 있었다 —
> `test_pred_motif_is_sequence_disables_period_metric` 이 `scored_pairs == 0` 을 단언했고,
> docstring 과 CLI help 도 같은 말을 했다. 셋이 서로를 뒷받침해 아무도 의심하지 않았다.
> 테스트는 교체했다(`test_pred_motif_is_sequence_disables_only_the_motif`).

- **왜**: 거짓이다. `scripts/scoring/convert_to_bed.py` 는
  `f"{chrom}\t{start}\t{end}\t{motif}\t{period}\tTRF\n"` 로 **5열에 period를 쓴다**.
  실제 TRF BED 5열은 `6, 29, 76 …` 이다. period 열은 비어 있지 않고,
  **one-to-one 스코어러가 그것을 버린다**(`--pred-col5 period` 로 호출되지만 `load()` 가 5열을
  저장하지 않고 `period_of()` 는 `len(motif)` 만 돌려주는데, `--pred-motif-is-sequence` 가 motif를
  비우므로 전부 결측이 된다).
- **영향**: Codex가 같은 518,719 쌍에 5열을 제대로 읽어 재계산한 값은 TRF **63.53% exact** —
  BWTandem 58.66%, ULTRA 59.61% 를 **앞선다**. 즉 S4의 period 순위가 뒤집힌다.
  (기전은 내가 검증했고, 63.53% 수치 자체는 재실행 전까지 미검증.)
- **대체값**: 스코어러에 명시적 `period` 필드를 두고 S4 재실행 후 예치. **미산출.**
- **근거**: `scripts/scoring/convert_to_bed.py:72-94`, `scripts/scoring/score_one_to_one.py`,
  `manuscript.md:518,537`.

### 6.24 "the index-based step is the seeding of Tiers 2 and 3" — Tier 1 도 FM-index를 쓴다

> ✅ **APPLIED (2026-09-03, `2901ff3`):** Discussion 4.3 을 'Tier 1 = 스캐너 + FM-index 열거' 로 교체.

- **왜**: 321행은 "BWTandem itself still scans in Tier 1 … the index-based step is the seeding of
  Tiers 2 and 3" 라 한다. 그러나 S2(456행)의 **모든** 논문 설정이 `TIER1_FMSCAN=1` 이고,
  그러면 Tier 1은 잔여 서열에 대해 FM-index 모티프 열거 패스를 **추가로 돌린다**(`tier1.py:325`).
- **대체값**: "Tier 1 = 슬라이딩 스캐너 + 능동적 FM-index 열거, Tier 2/3 = FM-index k-mer 시딩".
- **근거**: `manuscript.md:321` vs `manuscript.md:456`, `bwtandem/tier1.py:325`.
- ⚠️ 방향에 주의 — 이 오류는 FM-index 기여를 **과소** 서술한다. 그럼에도 메커니즘 서술이 틀렸다.

### 6.25 Abstract 의 "per-figure provenance manifest" — 존재하지 않는다

> ✅ **APPLIED (2026-09-03, `2901ff3`):** Abstract 에서 per-figure provenance manifest 약속 삭제.

- **왜**: `results/manifest.tsv` 에 figure 매핑이 **0건**이고, `results/figures/paper_figs/` 의
  그림 프로그램 **6개 전부**가 `# TODO: implement panels` 에서 끝난다.
- **대체값**: 실제 per-panel 매니페스트가 생기기 전까지 Abstract 의 약속을 삭제할 것.
- **근거**: `grep -ci fig results/manifest.tsv` = 0; `plot_*.py` 6/6 스텁.

### 6.26 S1.3 의 Arabidopsis 파라미터 범위 — 소기관 염색체가 빠져 있다

> ✅ **APPLIED (2026-09-03, `2901ff3`):** S1.3 을 핵 염색체 한정 + ChrC/ChrM 값 명시로 교체.

- **왜**: 427행은 표가 "세 게놈의 염색체들"의 값이라 하지만, Col-CEN 입력에는 ChrC·ChrM이 있고
  예치 BED에 각각 156 / 291 콜이 있다. 파라미터는 서열마다 계산되므로 실제 값은
  ChrC k=12 / stride=20 / max_occ=200, ChrM k=13 — **인쇄된 어느 범위에도 들어가지 않는다.**
- **대체값**: "핵 염색체 5개에 한한 범위"로 한정하거나 소기관 값을 추가할 것.
- **근거**: `manuscript.md:427`, `bwtandem/tier3.py:30`, `results/beds/README.md:8`.

### 6.27 ~~§3.2 의 BWTandem 밴드 recall 은 출처가 없다~~ — ✏️ **범위 축소 (2026-09-03)**

> ✏️ **PARTIAL (2026-09-03):** "91.31 / 94.68 은 예치 출처가 없고 예치본은 다른 값(91.92/95.21)을
> 준다"고 적었다. **절반만 맞았다.** 재생성 BED 로 직접 재채점하니 **정확히 91.31 / 94.68** 이 나온다 —
> 원고 값이 옳다. 예치본이 다른 값을 준 것은 그 로그가 **옛 출력**(밴드 콜 1,380)을 채점했기 때문이고,
> 재생성본은 1,533 으로 Table 2 와 일치한다.
>
> ```
> BWTandem         160,740  84.59%  1,380  57.17% | 99.73%  91.92%  95.21%   ← 예치 2026 로그 (옛 출력)
> BWTandem-regen   160,883  84.54%  1,533  57.08% | 99.72%  91.31%  94.68%   ← 재생성본 = Table 2
> ```
>
> **남은 결함은 값이 아니라 예치 누락이었고, 그것은 닫혔다**:
> `results/regen/colcen_banded_regen.txt` 로 예치했다.
>
> **경위**: 나는 `grep -rl "91.31" results/` 가 0건인 것을 "값이 틀렸을 수 있다"로 읽었다.
> 없었던 것은 값이 아니라 **산출물**이다. §6.6·§6.16 과 같은 계열의 성급한 판정이다 —
> 부재를 오류로 승격시키지 말 것.


---

## §7 파일 경로 격리 대응표

| 격리 경로 | 무엇 | 인용 가능? |
|---|---|---|
| `archive/2026-08-05-unreproducible/` | 비용 격리의 전체 기록 | 현행판 README만 |
| `archive/2026-08-05-unreproducible/CLAUDE.md_stale_figures.md` | CLAUDE.md에서 떼어낸 블록 | ❌ 비용 수치 인용 금지 (정확도는 유효) |
| `results/comparators2026/score_2026_human_6146742_crashed.txt` | 크래시 로그, 감사 목적 예치 | ❌ 수치 인용 금지 (동등성 증명용) |
| `resume-archive-to-20260902.md` | 127개 과거 스냅샷 | 🧊 이력 전용 |
| `resume.md.autolog-bak` | Stop 훅의 백업 (19:35 상태) | 🧊 복구 전용 |

---

## §8 방법론 함정 — 같은 실수를 막는 규칙

### 8.1 `int(float(parts[4]))` 금지 — `int()`를 쓸 것

filip의 변환 BED는 6열이고 5열이 정수 period지만, BWTandem 자신의 BED는 8열이고 5열이
`46.3` 같은 **copy count**다(그쪽엔 period 열이 없고 `len(motif)`가 period다).
`int(float("46.3"))`이 조용히 46을 돌려주어 **모든 period-band 지표가 copy 수를 셌다.**
이것이 §2의 "붕괴" 허상을 만들었다. 2026-08-05 수정, 나머지 8개 스코어러는 원래 정상.

### 8.2 `results/` 를 건드릴 때의 고정 순서

`scripts/benchmark/hash_deposited_beds.sh` 는 **워킹트리**를 해시해
`results/manifest.sha256` 과 `results/external_evidence.sha256` 을 잘라 다시 쓴다.
CI 는 **커밋된 트리**를 검증한다. 하루에 CI 를 세 번, 세 가지 방식으로 깨뜨렸다(이슈 #31).

> ❌ **SUPERSEDED (2026-09-03):** 이 절은 한동안 순서를
> `편집 → git add → 재해시 → 테스트 → 커밋` 이라 적었다. **위험하다.**
> 해시 스크립트가 `git add` **뒤에** 체크섬 파일을 다시 쓰므로, 스테이징된 것은 옛 해시이고
> 워킹트리만 새 해시인 상태로 커밋될 수 있다. 가드 테스트는 워킹트리를 보므로 그 커밋을
> 통과시킨다. 스크립트 자신이 머리에 **"Run this LAST"** 라 적어 두었고, 그 실패가 이미
> `e178707` 로 출하된 적이 있다. (Codex 계획 리뷰 2026-09-03 지적.)

**올바른 순서**:

1. `results/` 아래 **모든** 편집을 끝낸다 — 부분적으로 하지 않는다
2. 해시 스크립트를 **포그라운드**로 돌린다. 백그라운드로 돌린다면 PID 를 `wait` 하고
   종료 상태를 확인한 뒤 다음 단계로 간다 (2분을 훨씬 넘는다; 포그라운드 kill 은
   매니페스트를 반쯤 쓴 채로 남긴다)
3. **그 다음** 체크섬 파일 두 개를 포함해 전부 `git add`
4. `git diff --cached` 를 눈으로 확인하고, `git status` 에 **미스테이징 `results/` 변경이 없는지** 본다
5. 가드 테스트 → 커밋

`.gitignore` 가 `*.settings` 같은 증거 파일을 삼킨다 → `git add -f`.
`tests/test_deposit_hashes.py` 가 해시된 경로의 미추적을 잡는다.

### 8.3 이 Claude 세션 자체가 2 CPU / 4 GB SLURM 잡이다

`/proc/self/status`의 `Cpus_allowed_list`로 확인할 것. fork+exec 하나가 수 초다.
"노드가 느리다"고 결론내기 전에 이것부터 본다. **conda/mamba solve를 여기서 돌리지 말 것**
— 한 번은 3 h 08 m 동안 코어의 79%를 쥐고도 못 끝냈고, 같은 Node 22를 공식 타르볼로
설치하니 27 s였다. `squeue`/`sacct`/`sbatch`는 PATH에 없다
(`/cm/shared/apps/slurm/current/bin/`). Bash가 느려지면 Read 도구가 셸을 우회한다.

### 8.4 `pkill -f "<명령 문자열>"` 은 호출자를 죽인다

`pkill -f "codex exec -m gpt-5.4"` 가 **그 명령을 실행 중인 Bash 자신의 커맨드라인**에도
매치해 호출 셸이 함께 죽었다(exit 144, 백그라운드 폴러까지 동반 사망).
PID로 죽이거나 `while [ -d /proc/<pid> ]`로 기다릴 것. (2026-09-03)

### 8.5 Kimi에게 파일을 읽히지 말 것

`ask`의 내부 예산은 **900 s 고정**이고 `--background`로도 늘어나지 않는다.
manuscript를 grep·read 하게 두면 **매번** 타임아웃한다(2026-09-02에 2회 연속).
본문을 프롬프트에 **인라인**하고 도구 사용을 명시적으로 금지하면 11–18 KB 청크가
5–8분에 끝난다(2026-09-03 라운드 1: 5/5 성공). `--base <ref>`는 `bd32c8c`(4.37 M행 삭제)를
포함하는 범위에서 실패한다. 원문은 **글자 그대로** 인용할 것 — 셸 따옴표를 피하려고
고쳤더니 존재하지 않는 결함이 보고됐다.

### 8.6 노드 세대 차이가 1.55×를 만든다

BWTandem 발표값은 `cpu-s2-core-0`(intelv5), 재측정과 모든 경쟁 도구는
`cpu-s1-pgl-0`(intelv4). 같은 입력·같은 워커 수·같은 `/usr/bin/time`에서
7:49:08 vs 12:05:37 — **1.55×는 하드웨어 세대다.** 경쟁 도구와 비교할 값은 재측정 쪽이다.

### 8.7 되돌아온 실패는 다시 시도하지 말 것

| 시도 | 언제 | 결과 |
|---|---|---|
| `CATCHALL_MIN_EXCESS` 게이트 | 2026-07-05 | 실패 |
| phase/coherence 판별자 계열 | 2026-07-08 | AUC ≈ 0.5, 실패 |
| `TIER1_ENTROPY_GATE` | — | recall만 낮추고 precision 안 올림 (기본 off) |
| `CATCHALL_MIN_ENTROPY` | 2026-06-25 | region precision 이득 ≈ 0 (기본 off) |
| `TIER2_APPROX_SEED` (adotto) | — | region recall 평탄, bp precision 31→13.5% (기본 off) |
| Wavelet-tree 메모리 재작성 | — | 이 논문 범위 밖 |

env 레버 기반 recall/precision 프런티어는 **천장**에 있다.

---

## §9 진단으로 배제된 것 — 다시 조사하지 말 것

> 여기 있는 것은 **돌려봤고 깨끗했던** 진단이다. 미해결 항목이 아니다 —
> 아직 답이 없는 것은 [`todo.md`](todo.md) 에 있다.

### 9.1 flaky `TestAdjacentGroundTruth::test_sensitivity` — 아래 진단은 전부 소진됐다

| 시도한 것 | 결과 |
|---|---|
| valgrind memcheck (프로덕션 빌드) | 클린 |
| ASan 재빌드 — C 확장 4개 | 메모리 오류 0 |
| ASan 재빌드 — Cython `_accelerators` (앞선 패스가 빠뜨린 것) | 클린 |
| 환경 블록 크기 스윕 256종 ×5 | 실패 0 |
| mmap 베이스 스윕 (`ulimit -s` ×9) | 실패 0 |
| 같은 픽스처를 한 프로세스에서 반복 | 바이트 동일 |
| Chr4(18.8 Mb) 환경 3종 | 바이트 동일 |

약 1,640회 전수 재현 0. **원인은 프로세스 메모리 레이아웃 조건부**이고, 클러스터가 진단 이후
ASLR을 끄면서(`randomize_va_space=0`) 실패 레이아웃 자체에 도달할 수 없다.
**위 진단을 다시 돌리는 것은 낭비다.** 재현하려면 관리자에게 노드 1대 ASLR 재활성화를
요청해야 한다 — 그 요청은 `todo.md` F 절의 차단 항목이다.
추적기 보존: `docs/2026-08-10-flaky-test-instrumentation-campaign.md`.

### 9.2 `TIER1_STITCH_GAP` — 측정했고 값이 없었다

chr21 단주기 recall 에 **+0.2 pp** (≈중립). 장주기 코어용으로 opt-in 유지하되,
단주기 recall 을 올리려는 목적으로 다시 켜지 말 것.


### 9.3 `TIER2_APPROX_SEED` at identity ≥ 0.88 — 시험한 적이 없다

centromere/satellite 에서 도움이 될 가능성이 거론됐을 뿐 **측정된 적이 없다.**
기본값(0.78)에서는 adotto 에 대해 region recall 평탄, bp precision 31→13.5% 로 기각됐고
(§8.7), ≥0.88 은 그 기각의 연장이 아니라 **다른 미검증 설정**이다.

2026-09-03 판정: **이 논문의 작업이 아니다.** `todo.md` 에서 제외했다.
켜려면 먼저 측정할 것 — 목록에 없다는 이유로 켜지 말 것.

---

## §10 갱신 규칙

1. 새로 폐기되는 것이 생기면 **같은 커밋에서** 이 파일에 한 줄 추가한다.
   "나중에 정리"는 이 파일의 존재 이유를 무효화한다.
2. 항목 형식은 **무엇 / 왜 / 언제 / 대체물 / 근거 경로**. 대체물이 없으면 **"없음"**이라고 적는다.
3. `resume.md`·`CLAUDE.md`·`README.md`에서 폐기 서술을 지우지 말고 **격리 마커 + 이 파일의
   해당 절 링크**로 바꾼다. 연결이 없으면 규약 문서가 둘이 되고, 그것이 곧 중복 오해원이다.
4. 정정의 정정도 쌓는다 — 앞선 판단을 지우지 말고 옆에 새 마커를 잇는다(§3.1이 그 예다).
5. §6처럼 **수정 미적용** 항목은 적용될 때 §2로 옮기고 대체값을 확정한다.
   수정 *작업* 자체는 이 파일이 아니라 `todo.md` 에 있다 — 여기는 "그 값을 쓰지 말라"만 적는다.
6. **미해결을 여기 적지 않는다.** 검증이 남은 것, 답이 없는 것, 해야 할 일은 `todo.md` 로.
   이 경계가 무너지면 이 파일은 "하지 말 것"이 아니라 "언젠가 볼 것"이 되어 기능을 잃는다.
