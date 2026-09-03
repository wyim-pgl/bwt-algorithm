# ruleout.md — BWTandem 폐기 대장

> 📌 **정본**: 이 파일이 "쓰면 안 되는 것"의 유일한 정본이다. 폐기된 수치·주장·이름·ID를
> 인용하기 전에 여기부터 grep 한다. 최종 갱신 **2026-09-03**.
> 살아 있는 상태·다음 작업은 `resume.md`, 운영 규칙은 `CLAUDE.md`.

이 파일의 근거와 규약: 랩 위키 `guide/handoff-hygiene.md` §"세 번째 층 — 프로젝트 단위
`ruleout.md`". 문서 마커(❌/⚠️/✏️/🧊)는 사람을 막고, 경로 격리(`archive/`)는 스크립트를
막고, 이 파일은 **새 세션이 폐기된 값을 다시 집어드는 것**을 막는다.

항목 형식은 고정이다: **무엇 / 왜 / 언제 / 대체물 / 근거 경로**.
**대체물이 없으면 "없음"이라고 적는다** — 빈칸은 "나중에 채워도 된다"로 읽히지만
"없음"은 그 주장을 아예 하지 말라는 뜻이다.

```bash
grep -n '^### ' ruleout.md          # 폐기 항목 전수
grep -rn 'ruleout.md' resume.md CLAUDE.md   # 이 파일로 연결된 격리 지점
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

> ⚠️ 아래 12건은 Codex 리뷰(6.1-6.6)와 Kimi 라운드 2(6.7-6.12)가 지적하고 **내가 실물 파일·표로 재현**한 것이다.
> 원고·증거 트리에 **아직 그대로 있다**. 인용 금지, 수정 대상.

### 6.1 Table 2 의 TRASH 비용 셀 두 개

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

- **무엇**: "same binary 1.2.1", "version- and input-matched", 하한 표시 없는 "1.55×",
  5시간 정체에 붙은 "(chr1 centromeric satellite region)".
- **왜**: `24cd12a`·`e8ac50c`·`67a3f42`가 원고에서 고친 네 가지가 증거 인덱스에 남아
  **원고와 매니페스트가 모순**이다.
- **대체값**: 정정된 원고 문구.

### 6.3 AniAnn's 스레드 수

- **무엇**: `results/manifest.tsv` 행 123·126·129·132 의 `threads=1`.
- **왜**: sbatch와 모든 실행 로그가 `-j 2`, `--cpus-per-task=2`.
- **대체값**: `2` (+ `manifest.sha256` 재해시 — §8.2 순서 준수).
- **근거**: `exp1_human/tools2026/runs/run_anianns.sbatch:5,22,23`.

### 6.4 2026 도구 메모리 셀 — 같은 표 안에서 규칙 불일치

| 셀 | 로그 KiB | GiB 규칙 | 현재 표기 |
|---|--:|--:|--:|
| AniAnn's human | 550,700 | 0.53 | 0.53 ✅ |
| AniAnn's maize | 572,588 | 0.55 | 0.55 ✅ |
| AniAnn's Col-CEN | 503,968 | **0.48** | 0.50 ❌ |
| longdust Col-CEN | 66,020 / 66,340 | **0.06** | 0.07 ❌ |

- **근거**: `exp1_human/tools2026/runs/*.time`. 표기 위치 `results/comparators2026/README.md:27-28`.

### 6.5 "roughly 5% more sequence" (manuscript.md:84)

- **왜**: 원고가 주는 두 수는 3.25 Gb와 3.15 GB → **+3.2%**. 게다가 하나는 염기,
  하나는 파일 바이트로 두 문장 만에 단위가 섞였다.
- **대체값**: 두 값을 같은 단위(염기)로 다시 재고 비율을 재계산할 것. 확정값 **미산출**.

### 6.6 stride 포화점 "roughly 24 Mb" (manuscript.md:66)

- **왜**: 원고의 공식 `n / 40,000` + clamp 300 → 포화점은 **12 Mb**다. 24 Mb가 되려면 배수
  0.5가 필요한데 실제 프리셋 배수는 `fast 1.6 / balanced 1.0 / sensitive 0.4`
  (`bwtandem/tier3.py:23-25,31,63`) — **어느 것도 24 Mb를 주지 않는다.**
  원고가 프리셋 배수를 언급조차 하지 않아 재현이 불가능하다.
- **대체값**: 프리셋 배수를 포함한 공식과 그로부터 유도한 포화점. **미산출.**

### 6.7 Table 3C-b 캡션 — 표가 반증하는 최상급과 틀린 차이값

- **무엇**: "the smallest banded loss of the four" (tantan −0.92) 와 "finishes 14.36 points
  above BWTandem".
- **왜**: 같은 표의 손실 열이 TRF **−0.87** 로 tantan보다 작다 — 최상급이 거짓이다.
  차이값도 55.59 − 40.71 = **14.88** 이지 14.36이 아니며, 표의 어떤 조합도 14.36을 주지 않는다.
- **대체값**: "TRF loses least (0.87); tantan 0.92" / **14.88**.
- **근거**: `manuscript.md:298` 캡션 vs `manuscript.md:306-309` 표.
- ⚠️ Table 3B-b(254행)의 같은 문형은 참이다(0.50 < 0.63 < 1.20) — 재검산 없이 옮겨 쓴 흔적.

### 6.8 "at most 1.5 / 2.8 points" — 하한을 반올림해 만든 거짓 상계

- **왜**: TRF의 최대 손실은 **1.54**, ULTRA는 **2.81**. "at most"는 정확한 상계이므로
  내림 반올림은 허용되지 않는다.
- **대체값**: 1.6 / 2.9, 또는 1.54 / 2.81 을 그대로.
- **근거**: `manuscript.md:280`.

### 6.9 Table 1c — 캡션이 네 도구, 표는 세 행 (TRASH 누락)

- **왜**: 캡션은 "the stratum BWTandem, TRF, TRASH and the re-run tantan entered" 라 하지만
  표에는 BWTandem·TRF·tantan 세 행뿐이다. 경쟁 도구 한 행이 조용히 빠져 있다.
- **대체값**: TRASH 행을 넣거나, 캡션에서 빼고 제외 사유를 밝힐 것. **미산출.**
- **근거**: `manuscript.md:139` vs `manuscript.md:141-145`.

### 6.10 >100 bp 층 서술 — TRF가 이기는 두 열만 빠짐

- **왜**: 본문 107행은 BWTandem 3.43% / 13.39% / precision 64.60% 을 쓰고 "TRF reaches 2.02%
  and 11.83%" 로 **recall 두 개만** 옮긴다. 같은 표에서 TRF의 region precision은 **97.18%**
  (BWTandem 64.60%), bp precision도 **31.19 vs 29.25** 로 TRF가 앞선다.
- **대체값**: 같은 문장에 TRF의 precision 두 값을 넣을 것.
- **근거**: `manuscript.md:107` vs `manuscript.md:141-145`.

### 6.11 Table 2 — 자기 도구의 최저 수치만 서술되지 않음

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

`scripts/benchmark/hash_deposited_beds.sh`는 **워킹트리**를 해시해
`results/manifest.sha256`을 잘라 다시 쓰는데, CI는 **커밋된 트리**를 검증한다.
하루에 CI를 세 번, 세 가지 방식으로 깨뜨렸다(이슈 #31).

1. 편집 완료 → `git add` → **그 다음** 재해시 → 테스트 → 커밋
2. 실행 중에는 `results/` 아래를 절대 건드리지 않는다
3. 2분을 훨씬 넘는다(외부 파일 약 75개, 일부 수 GB): `run_in_background` + 600 s 타임아웃.
   포그라운드 kill은 매니페스트를 반쯤 쓴 채로 남긴다
4. `.gitignore`가 `*.settings` 같은 증거 파일을 삼킨다 → `git add -f`.
   `tests/test_deposit_hashes.py`가 해시된 경로의 미추적을 잡는다

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

## §9 미해결 — 폐기도 채택도 아닌 것 (인용 금지)

> 목록에 없으니 써도 된다"로 읽히는 것을 막기 위해 명시적으로 격리한다.

### 9.1 flaky `TestAdjacentGroundTruth::test_sensitivity`

memcheck·ASan(Cython 포함) 전부 클린, 환경 블록·mmap 베이스 스윕 약 1,640회 재현 0.
클러스터가 진단 이후 ASLR을 끄면서(`randomize_va_space=0`) **실패 레이아웃에 도달할 수
없는 상태**다. 원인 미확정 — "고쳤다"고 쓰지 말 것. 재현하려면 관리자에게 노드 1대
ASLR 재활성화를 요청해야 한다. 추적기는 `docs/2026-08-10-flaky-test-instrumentation-campaign.md`.

### 9.2 `TIER1_STITCH_GAP`

chr21 단주기 recall에 ≈중립(+0.2 pp). 폐기도 채택도 아님 — 장주기 코어용으로 opt-in 유지.

### 9.3 `TIER2_APPROX_SEED` at identity ≥ 0.88

centromere/satellite에서 도움이 될 수 있으나 **미시험**. 시험 없이 켜지 말 것.

### 9.4 2026-09-03 원고 교차 검증의 미확정 발견

Kimi 라운드 1 5건 + Codex 20건 중 **§6에 올린 6건만 내가 실물로 재현**했다.
나머지(예: TideHunter 제외 사유 불일치, 2,000 bp 상한 논지 문제, S1/S2 상수 개수 불일치,
merge gap 10p→100p 불연속)는 **미검증**이다. 라운드 2·3이 진행 중이며,
검증 전에는 원고 수정 근거로 쓰지 말 것.
출력: `docs/2026-09-03-codex-review-findings.md`, 스크래치패드 `msreview/out_*.md`.

---

## §10 갱신 규칙

1. 새로 폐기되는 것이 생기면 **같은 커밋에서** 이 파일에 한 줄 추가한다.
   "나중에 정리"는 이 파일의 존재 이유를 무효화한다.
2. 항목 형식은 **무엇 / 왜 / 언제 / 대체물 / 근거 경로**. 대체물이 없으면 **"없음"**이라고 적는다.
3. `resume.md`·`CLAUDE.md`·`README.md`에서 폐기 서술을 지우지 말고 **격리 마커 + 이 파일의
   해당 절 링크**로 바꾼다. 연결이 없으면 규약 문서가 둘이 되고, 그것이 곧 중복 오해원이다.
4. 정정의 정정도 쌓는다 — 앞선 판단을 지우지 말고 옆에 새 마커를 잇는다(§3.1이 그 예다).
5. §6처럼 **수정 미적용** 항목은 적용될 때 §2로 옮기고 대체값을 확정한다.
