# resume.md — bwtandem 핸드오프 노트

**최종 업데이트:** 2026-08-27
**감사 대상 트리:** 이 체크아웃, HEAD `7a80ed7` (**2026-04-01**) — 푸시됨. 미커밋인 것은 워킹 트리뿐.
**주의:** 활성 원고 브랜치는 `perf/exp1-human-sensitivity` `e0acd56` (2026-08-19), **80커밋 앞섬** (`/data/gpfs/assoc/pgl/devel/bwt-algorithm`)
**현재 단계:** 벤치마크 재생성 캠페인 실행 중 (SLURM 6141841/6141842, 핀 `0363d8b`) — 상세는 바로 아래 최신 스냅샷
**읽는 법:** 이 파일은 스냅샷이 위로 쌓인다. **최신 스냅샷 하나만 현재 상태**이고, 격리 배너(🧊) 아래는 전부 과거 기록이다. 마커 규약: wiki `guide/handoff-hygiene.md`

---

## ▶️ 현재 상태 스냅샷 (2026-08-31 09:23 — 여기서 다시 시작)

### origin 오설정 해결 (2026-08-27 감사 미해결 항목)

remote 토폴로지를 devel 클론과 동일한 관례로 정리:
- `origin` → **wyim-pgl/bwt-algorithm** (실제 트래커, https)
- `framazan` → framazan/bwtandem 별도 remote로 보존

처리 내역: `git remote rename origin framazan` 후 새 `origin` 추가, 양쪽 `fetch --prune` 실행. 검증: HEAD `7a80ed7`이 origin/main의 ancestor (behind 54 — 감사용 동결 상태이므로 pull 하지 않음, 워킹 트리 무접촉). `main` 업스트림을 `origin/main`으로 재지정. resume.md 하단 감사 항목에 ❌ SUPERSEDED 마커 부착.

### 다음 작업
- Issue #11 판독 (400건 blinded sheet + dot plot, 저자 2인 직접 판독)
- reciprocal-overlap 민감도 재계산, #22 archive 정책, 태그·릴리스

## ▶️ 이전 스냅샷 (2026-08-31 09:40)

### origin 오설정 해결 (2026-08-27 감사 항목)

이 체크아웃의 remote 토폴로지를 devel 클론과 동일 관례로 정리:
- `origin` → **`wyim-pgl/bwt-algorithm`** (https, 실제 트래커·푸시 대상)
- `framazan` → `framazan/bwtandem` 별도 remote로 보존 (의도적 협업 저장소)
- 양쪽 fetch --prune 완료, `main` 업스트림을 `origin/main`으로 재지정
- 검증: HEAD `7a80ed7`이 origin/main(`6a8ae40`)의 ancestor 확인 (behind 54 — 이 트리는 감사용 동결 상태이므로 pull 하지 않음, 워킹 트리 무접촉)

resume.md 하단 감사 기록의 해당 항목에 ❌ SUPERSEDED 마커 부착.

### 다음 작업
- Issue #11 판독 (400건, 저자 2인) — 남은 최우선 저자 병목
- reciprocal-overlap 민감도 재계산, #22 archive 정책, 태그·릴리스

## ▶️ 이전 스냅샷 (2026-08-31 08:53)

### 캠페인 9/9 최종 착지 · 원고 커밋 `8e284bd` 푸시 완료

**p100·idsweep 쌍 채점 완료** — 6145300(p100 4종) 및 6145311(idsweep 5종) 종료 후 원고의 pending 마커 전량 실측값으로 교체. Table 1b(F 79.88/50.62/32.55/55.51), Table 1d(P/B/F/H 54.69→81.60% regRecall + runtime), Supplementary S3(62.59→81.46%), Figure 1 CSV+PNG/PDF 완성.

### 검증 게이트 3중 통과

1. 두 채점 보고서 baseline 행이 v2 채점표와 소수점 일치 (같은 코드·입력 확정)
2. idsweep 0.72 arm이 deposited `regen_human.bed` 콜 단위 재현 (adjPrec 79.45까지 전 패널 일치)
3. S2 특성화가 구빌드 값과 사실상 동일 (unique 875,084 vs 구 875,224; P run "pass와 무관" 문장 재검증)

### 증거 배포 · manifest 재구축

provenance 9건 + 채점 2건 + 특성화 2건을 `results/regen/`에 배포. manifest.tsv에 10행 추가(잡·노드·sacct MaxRSS·SHA 정합), manifest.sha256 전체 재계산 후 **43/43 검증 OK**. #7에 완결 코멘트 게시.

### 원고 최종 정정

Introduction·Section 3.1·Discussion·Methods·Data Availability의 구 잡 번호(6129408/6129410)와 pending 문장 전부 제거. Intro line 42 및 Section 3.1 line 107 수정 완료 후 `perf/exp1-human-sensitivity` 브랜치에 푸시.

### 함정 기록

frozen 스코어러는 스크립트 위치 기준으로 GT 경로를 잘못 풀므로 `score_table1_regen.py` 래퍼(`--bwt-bed`) 필수. `--adj both`는 무효 — `published,loo` 콤마 리스트로 실행해야 함.

### 다음 작업

- Issue #11 판독(400건 blinded sheet + dot plot, 저자 2인 직접 판독)
- reciprocal-overlap 민감도 재계산(Methods withheld 항목 1건)
- #22 archive 정책·태그·릴리스(저자 결정)

## ▶️ 이전 스냅샷 (2026-08-31 09:20)

### p100/idsweep 캠페인 완결 → 원고 착지 (`8e284bd` push)

**9/9 CLEAN 착지** (마지막 idsweep 0.68: 4,221,750콜, 7h49m). 채점 2건(6145300 p100, 6145311 idsweep) 완료 후 원고의 6129408/6129410 pending 마커 전량 교체:
- **Table 1b**: native matched-range F 행 (3,993,151 / 79.88 / 50.62 / 32.55 / 55.51)
- **Table 1d**: P/B/F/H 4 operating point + published-rule support rate (77.12/78.28/79.85/78.08) + SLURM elapsed (4스레드 명기)
- **Supplementary S2**: F run unique 875,084건(21.9%) 특성화 — `analyze_unique_regions.py` 재실행, catalog 밖 91.3%는 P run에서도 동일(pass 무관 문장 재검증)
- **Supplementary S3**: 5-arm sweep 62.59→81.46% region recall
- **Figure 1**: CSV 채움 + fail-closed 렌더러로 PNG/PDF 생성 (렌더러도 커밋)
- Intro "Three properties", Section 3.1/Discussion/Methods/DA 문면 정합

**검증 게이트 (모두 통과):** ① 두 채점 보고서의 baseline 행 전부 v2 표와 소수점 일치 ② idsweep 0.72 arm이 deposited `regen_human.bed`를 콜 단위 재현 (adjPrec 79.45까지 전 패널 일치) ③ S2 특성화가 구빌드 값과 사실상 동일 (875,084 vs 875,224)

**증거 배포:** `results/regen/`에 provenance 9건+채점 2건+특성화 2건, manifest.tsv +10행, manifest.sha256 재구축 43/43 OK. #7에 완결 코멘트 2건.

**함정 기록 (재발 방지):** frozen 스코어러는 스크립트 위치로 GT 경로를 유도하므로 반드시 `score_table1_regen.py` 래퍼로 실행 (`--bwt-bed`); `--adj both`는 무효 — `published,loo` 콤마 리스트여야 함.

### 다음 작업
- Issue #11 판독 (400건, 저자 2인) — 남은 유일한 저자 병목
- reciprocal-overlap 민감도 분석 재계산 (Methods에 withheld로 명시된 잔여 항목)
- #22 archive 정책, 태그·릴리스 (저자 결정)

## ▶️ 이전 스냅샷 (2026-08-31 08:21)

### 캠페인 9/9 전체 착지 완료 · 병렬 채점 제출

**idsweep 0.68 (6143151_4) CLEAN 완료** — 4,221,750콜, 7h49m, 피크 17.4GB. 핀 `0363d8b` 9개 task 전부 provenance CLEAN 확정. 두 차례 전멸 이후 무손실 착지.

현재 **2건 채점 병렬 실행**:
- **6145300**: p100 4종 (published+loo 4패널, Table 1b/1d용)
- **6145311**: idsweep 5종 (off/0.80/0.76/0.72/0.68, Supplementary S3용) — 방금 제출

idsweep 콜 수 곡선 검증: off 2.77M → 0.80 3.30M → 0.76 3.32M → 0.72 4.01M → 0.68 4.22M (identity 완화에 따라 단조증가, 정상).

### 다음 작업

- 6145300·6145311 완료 대기 → 각 결과 검증
- Table 1b/1d 최종 값 정리 · Supplementary S3 채우기
- 캠페인 완결 후 원고·위키 일괄 갱신

## ▶️ 이전 스냅샷 (2026-08-31 08:13)

### 채점 1차 완료 · 기준 행 검증 통과 → p100 4종 결과 확정

**잡 6145292 완료** (exit 0) — 내장 검증 통과: 기준 행(BWTandem, ULTRA, tantan 경쟁군, 101–2000 stratum)이 v2 채점표와 **소수점까지 일치**. 같은 코드·같은 입력에서 p100 4종만 추가된 것 확인.

**p100 4 variant 최종 결과** (Table 1b/1d 채울 값):
- P (catch-all off): 2,085,734 calls, 54.69 regRecall, 62.16 regPrec
- B (id 0.76+cap50): 3,410,417 calls, 71.89 regRecall, 53.28 regPrec  
- F (id 0.72, ≥3 copies): 3,993,151 calls, 79.88 regRecall, 50.62 regPrec
- H (id 0.72): 4,355,741 calls, 81.60 regRecall, 48.44 regPrec

matched-range(≤100) 대비: ULTRA 81.62/53.66, tantan 78.00/57.66 — H가 ULTRA와 동급 recall, F가 균형점.

**조정 규칙 오류 수정** — `--adj both` 대신 `--adj published,loo` 필요 (published 규칙은 ULTRA+tantan corroborator). 잡 **6145300** 재제출(~25분 예상). Table 1d Support rate 열 채우기 위함.

idsweep 0.68 계속 실행 중.

### 다음 작업

- 잡 6145300 완료 대기 → Support rate 값 확인 → Table 1b/1d 최종 채우기
- idsweep 0.68 완료 → 5종 일괄 채점(S3) 제출
- 캠페인 완결 후 원고·위키 일괄 갱신

## ▶️ 이전 스냅샷 (2026-08-31 08:06)

### 낡은 벤치마크 산출물 경로 격리 완료 (`RETIRED_DO_NOT_USE/`)

handoff-hygiene 규약(경로 격리)에 따라 stale 체크아웃의 낡은 산출물을 `RETIRED_DO_NOT_USE/`로 이동 (~1.4GB, 삭제 아님):
- `results/` — #7 CRITICAL stale corpus (bwtandem_Chr4_v3.bed의 all-N 오탐)
- `ath_cen/` — #8 좌표 혼용, 수치 사용 불가 판정
- `oropetium/` — #18 캠페인 제외 (저자 결정), FASTA만 보존
- `TRASH/` — 구 도구+산출물, 기준선은 filip 고정 BED
- `figures/`, `slurm-*.out` 7건 — 구 산출물 기반

안전 확인:
- 경로 참조 0건, 깨진 심링크 0건 (이 체크아웃 + exp1_human 양쪽)
- `RETIRED_DO_NOT_USE/README.md`에 항목별 격리 사유·현행 정본 포인터 기록
- 추적 파일 41건 D 표시 — 의도된 격리, `git restore` 금지 (README·resume.md 명시)
- 활성 증거(devel 워크트리)·v1 채점표는 미접촉

**진행 중:** 채점 잡 6145292 (`score_table1_regen.py` 래퍼, v2 동일 설정) · idsweep 0.68 (6143151_4) RUNNING

### 다음 작업

- 6145292 완료 → 기준 행이 v2 표와 일치 검증 → Table 1b/1d 값·잡 번호 업데이트
- idsweep 0.68 착지 → 5종 일괄 채점(S3) 제출
- 캠페인 완결 후 위키 일괄 갱신

## ▶️ 이전 스냅샷 (2026-08-31 08:20)

### 낡은 벤치마크 산출물 경로 격리 완료 (`RETIRED_DO_NOT_USE/`)

handoff-hygiene 규약(경로 격리)에 따라 이 stale 체크아웃의 superseded 산출물을 `RETIRED_DO_NOT_USE/`로 이동 (~1.4GB, 삭제 아님):
- `results/` (bwtandem_Chr4_v3.bed의 all-N 오탐 = #7 stale corpus 포함), `ath_cen/` (#8 좌표 혼용), `oropetium/` (#18 캠페인 제외), `TRASH/` (구 도구+산출물), `figures/`, `slurm-*.out` 7건
- 격리 사유·현행 정본 포인터는 `RETIRED_DO_NOT_USE/README.md`에 기록
- 깨진 심링크 검사: 이 체크아웃·exp1_human 모두 **0건** (규약 필수 단계)
- 추적 파일 41건이 `git status`에 D로 표시됨 — **의도된 격리, `git restore` 금지** (README에 명시)
- 활성 증거(devel/exp1_human, devel/bwt-algorithm)는 미접촉. regen의 v1 채점표는 "기록용 보존" 결정이 있어 격리 제외

**진행 중:** score_p100 채점 잡 **6145292** (래퍼 `score_table1_regen.py` 방식 — 1차 6145291은 frozen 스코어러 직접 호출로 GT 경로 오해석 즉사) · idsweep 0.68 (6143151_4) RUNNING

### 다음 작업
- 6145292 완료 → 기준 행이 v2 표와 일치하는지 검증 → Table 1b/1d 채우기
- idsweep 0.68 착지 → 5종 일괄 채점 (S3)
- 캠페인 완결 시 위키 일괄 갱신

## ▶️ 이전 스냅샷 (2026-08-31 07:54)

### 6145291 실패(GT 경로 오류) → 6145292 재제출

**첫 시도 실패:** 잡 6145291이 즉시 실패. 원인: worktree의 `score_table1.py`가 구버전 스코어러를 직접 호출하면서 GT 경로(SOURCES 기준)가 worktree 디렉토리 기준으로 잘못 풀림. 반면 원고의 Table 1b/1d는 `0c3c982`에서 커밋된 regen-safe 래퍼(`score_table1_regen.py`, `--bwt-bed` 오버라이드)로 생성됨.

**재제출:** 올바른 래퍼(`score_table1_regen.py`)로 sbatch **6145292** 제출 (v2와 동일 설정 + p100 4종 extra). 기준 행들이 v2 표와 정확히 일치하는지가 자체 검증됨.

**검증:** 6145292 완료 후 결과 확인 → matched-range·operating-point 값 채우기.

**진행 중:** idsweep 0.68 착지 (6143151_4, 여전히 RUNNING).

### 다음 작업

- 6145292 완료 대기 → 결과 검증 → Table 1b·1d 값 + 잡 번호 업데이트
- idsweep 0.68 완료 → 5종 일괄 채점(S3) 제출
- 캠페인 완결 후 위키·원고 일괄 갱신

## ▶️ 이전 스냅샷 (2026-08-31 07:53)

### 0.72 재현성 확정 + #25 수정·close + Table 1b/1d 채점 제출

**재현성 검증 완료:** idsweep 0.72 arm의 BED(`bwt_human_idsweep_0.72_p2000_0363.bed`)가 deposited `regen_human.bed`와 바이트 수 동일(211,949,715B) + 정렬 후 SHA-256 완전 일치. 독립 재실행이 원고 human 재생성을 콜 단위로 100% 재현.

**#25 이슈 처리:** 추적된 `__pycache__/*.pyc` 2개(dirty-gate 경쟁 근본 원인)를 `git rm --cached`로 추적 해제, 커밋 `77ba9fe`를 `perf/exp1-human-sensitivity`에 push. 전체 테스트 **179 passed** 확인 후 issue close.

**#7 코멘트:** 8/9 착지 표(콜 수·시간·메모리) + 0.72 재현성 증거 게시.

**Table 1b/1d 채점 제출:** SLURM **6145291** (workdir `score_p100_work`, `--adj both`). p100 4종을 `--extra` 행으로 v2 채점표 생성, 기존 corroborator 관례 유지. 완료 후 matched-range·operating-point curve 값 채우기 예정.

**주의:** 원고의 Table 1b/1d 캡션이 죽은 잡 번호(6129408) 인용 중 — 값 교체 시 6141841/6143150으로 고칠 예정.

### 다음 작업

- idsweep 0.68 착지 대기(7.4h 경과, RUNNING)
- 6145291 완료 → Table 값 + 잡 번호 업데이트
- 0.68 착지 → 5종 일괄 채점(Supplementary S3용) 제출
- 캠페인 완결 후 위키 일괄 갱신

## ▶️ 이전 스냅샷 (2026-08-31)

### 8/9 착지 CLEAN + 0.72 재현성 확정 + #25 처리 + Table 1b/1d 채점 가동

**착지 (전부 provenance CLEAN, 핀 `0363d8b`):** p100 P 2,085,734 / B 3,410,417 / F 3,993,151 / H 4,355,741 · idsweep off 2,770,930 / 0.80 3,304,372 / 0.76 3,322,038 / **0.72 4,014,108**. 잔여: idsweep 0.68 (6143151_4) RUNNING.

**재현성 확정:** idsweep 0.72 BED가 deposited `regen_human.bed`와 **바이트 수 동일(211,949,715) + 정렬 후 SHA-256 완전 일치**(`10c3aef0…`) — 행 순서만 스레드 수(2→4) 차이. 독립 재실행이 deposited human 재생성을 콜 단위 100% 재현.

**이슈 처리:**
- **#25 신규 등록→수정→close**: 추적된 `__pycache__/*.pyc` 2개가 dirty-gate 경쟁의 근본 원인 → `git rm --cached` 커밋 `77ba9fe` push (perf/exp1-human-sensitivity, 179 passed 확인 후). worktree의 skip-worktree 플래그는 잔여 arm 실행 중이라 유지.
- **#7 코멘트**: 착지 표 + 0.72 재현성 증거 게시.

**Table 1b/1d 채점 가동:** `score_table1.py`(worktree본 = 브랜치 tip과 diff 0)에 p100 4종을 `--extra`로 넣어 sbatch **6145291** 제출 (--adj both, workdir `score_p100_work`, 출력 `score_table1_p100.txt`). extra 행은 corroborator로 쓰이지 않아 기존 관례 유지.

### 다음 작업
- 6145291 채점 완료 → Table 1b(P matched-range)·1d(P/B/F/H) 값 채우기 (runtime은 provenance elapsed_wall_s)
- idsweep 0.68 착지 → 5종 일괄 채점 (Supplementary S3 / identity sweep) — 같은 방식 --extra 5종
- 위키 갱신은 캠페인 완결 시 일괄

## ▶️ 이전 스냅샷 (2026-08-30 13:42)

### 2차 즉사 (8/9 task) — dirty-gate 경쟁 → skip-worktree로 해소, 재제출 완료

git 수정(gitshim)은 유효했으나 **8개 task가 다시 0~1초 즉사** (이번엔 git이 아니라 다른 함정). 원인: `6141841_0`(P variant)이 먼저 cpu-50에서 실행되며 Python import가 **저장소에 추적된 `src/c_extensions/__pycache__/*.pyc` 2개를 재생성** → worktree 더티 → 몇 초 늦게 도착한 8개가 클린 트리 게이트 검사에서 침묵 사망. 9개 task가 worktree 1개를 공유하는 한 필연적 경쟁 (8/28의 `pip wheel` dirty-gate 사건과 같은 뿌리).

**수정:** 두 `.pyc`에 `git update-index --skip-worktree` 설정 — 커밋·내용 변경 없음, 핀 `0363d8b` 유지, 소스 `.py`는 계속 게이트 보호. 추적된 `.c` 4개는 이번 잡이 in-tree 빌드를 안 하므로 안전 확인.

**현재 상태:**
- `6141841_0` (P) — **RUNNING** (cpu-50)
- `6143150_[1-3]` (B/F/H) + `6143151_[0-4]` (idsweep) — 재제출 완료, 12h, PENDING
- 모니터 교체 완료 — 9개 task 전체 감시 중

### 다음 작업

- 6141841_0 / 6143150 / 6143151 착지 모니터링
- 캠페인 후 추적된 `.pyc` 2개 `git rm --cached` + `.gitignore` 제거

## ▶️ 이전 스냅샷 (2026-08-30 14:15)

### 2차 즉사 (8/9 task) — 추적 .pyc dirty-gate 경쟁 → skip-worktree로 해소, 재제출 완료

git 수정(6141841/42)은 유효했으나 **cpu-48의 8개 task가 다시 0~1초 즉사** (err에 오류 메시지 없음 = `set -e`의 `test` 게이트 침묵 사망). 유일하게 먼저 게이트를 통과한 `6141841_0`(P variant)이 cpu-50에서 실행되며 Python import가 **추적된 `src/c_extensions/__pycache__/*.pyc` 2개를 재생성** → 트리 더티 → 몇 초 늦게 게이트에 도달한 8개가 클린 트리 검사에서 죽음. 9개 task가 한 worktree를 공유하는 한 필연적 경쟁 (2026-08-28 `pip wheel` dirty-gate 사건과 같은 뿌리: 추적된 빌드 생성물).

**수정:** worktree에서 두 .pyc에 `git update-index --skip-worktree` 설정 — 커밋·내용 변경 없음, 핀 `0363d8b` 유지, 소스 `.py`는 계속 게이트 보호. 게이트 클린 확인 후 실패분만 재제출: **6143150_[1-3] (B/F/H) / 6143151_[0-4] (idsweep), 12h 확인, PENDING.** `6141841_0`(P)은 계속 RUNNING (cpu-50, 29분+). 추적된 `.c` 4개는 이번 잡이 in-tree 빌드를 안 하므로 안전 확인.

**후속 과제 (캠페인 종료·커밋 금제 해제 후):** 추적된 `__pycache__/*.pyc` 2개 저장소에서 제거(git rm --cached + .gitignore).

### 다음 작업
- 6141841_0 / 6143150 / 6143151 착지 모니터링 (모니터 교체 가동) → Table 1b/1d, Figure 1, Supplementary S3
- Issue #11 판독 수렴
- 캠페인 후 .pyc 추적 해제 커밋

## ▶️ 이전 스냅샷 (2026-08-30 11:54)

### git worktree 미인식 원인 확정 및 수정·재제출 완료

**진단 결과:** 컴퓨트 노드 `/usr/bin/git` = 1.8.3.1 (el7 스톡, worktree `commondir` 미지원). 로그인 노드 git 2.47과 차이.

**수정:** `maize_rescore/gitshim/` 디렉토리에 `~/micromamba/bin/git` (2.49) 심링크 생성. 두 sbatch(`run_human_p100_0363`, `run_human_idsweep_0363`) 머리에 `PATH="$M/gitshim:..."` 고정. `run_with_provenance.sh` 내부 bare `git` 호출 4곳 모두 커버 확인.

**2차 함정 대응:** 1차 재제출(6141839/40)이 남아 있던 `SBATCH_TIMELIMIT=13-23:00:00` 환경변수에 잡혀 13일 walltime으로 진입. 즉시 취소하고 `SBATCH_*` 전부 unset 후 재제출.

**재제출 결과:** 6141841 (regen_p100, P/B/F/H ×4) / 6141842 (regen_idsweep ×5) — 12:00:00 TimeLimit 확인, PENDING.

### 다음 작업

- 6141841/6141842 착지 모니터링 (5분 간격, 이벤트 알림) → Table 1b/1d, Figure 1, Supplementary S3 채우기
- Issue #11 판독 수렴
- Pronghorn 위키에 gotcha 2건 기록 (컴퓨트 노드 git 1.8.3.1 + worktree, SBATCH_* 환경변수)

## ▶️ 이전 스냅샷 (2026-08-30 11:55)

### worktree git 사망 원인 확정 → 수정·재제출 완료 (6141841/6141842)

**프로브 6140000 (cpu-48) 결과 — 원인 확정:** 컴퓨트 노드 `/usr/bin/git` = **1.8.3.1** (el7 스톡, worktree `commondir` 미지원). 로그인 노드 `/usr/bin/git`은 2.47이라 사전 테스트는 통과됨 — 노드 이미지 전용 함정. 죽은 6129408/6129410은 micromamba 없는 PATH로 제출돼 bare `git`이 1.8.3.1로 풀림. 최신 git(2.49)로는 같은 노드·같은 worktree에서 정상 동작 확인.

**수정:** `maize_rescore/gitshim/git → ~/micromamba/bin/git` 심링크 디렉토리 생성, 두 sbatch(`run_human_p100_0363`, `run_human_idsweep_0363`) 머리에서 PATH 앞에 고정 — 스크립트 본문 + `run_with_provenance.sh` 내부 bare git 4곳 모두 커버. git만 담아 python 등은 안 가림.

**재제출 중 2차 함정 재발:** 1차 재제출(6141839/40)이 `SBATCH_TIMELIMIT=13-23:00:00` 셸 환경변수에 잡혀 13일 walltime — 즉시 취소 후 `SBATCH_*` 전부 unset하고 재제출. **6141841 (regen_p100 P/B/F/H ×4) / 6141842 (regen_idsweep ×5), 12:00:00 확인, PENDING.** 게이트(핀 `0363d8b`·클린 트리)는 제출 전 사전 검증 통과.

산출물 오염 없음 (0초 사망이라 out 파일 미생성). worktree·입력·환경변수 설정 전부 기존 그대로.

### 다음 작업
- 6141841/6141842 착지 모니터링 (모니터 가동) → Table 1b/1d, Figure 1, Supplementary S3 채우기
- Issue #11 판독 수렴
- pronghorn 위키에 gotcha 2건 기록: 컴퓨트 노드 git 1.8.3.1 + worktree 함정, SBATCH_TIMELIMIT 재발

## ▶️ 이전 스냅샷 (2026-08-30 09:42)

### SLURM 6129408/6129410 (native p100 + idsweep 5) 즉사 — git worktree 미인식

9개 task 전부 오늘 새벽 03:25 시작 후 0~1초 만에 즉사. 로그 오류: `fatal: Not a git repository: /data/gpfs/assoc/pgl/devel/bwt-algorithm/.git/worktrees/worktree` — run_with_provenance.sh 내 `git rev-parse HEAD` 첫 호출에서 `set -e`로 죽음.

**원인:** 9개 task 모두 cpu-48 노드 배정. 이 에러는 구식 git(< 2.5)이 worktree의 `commondir` 레이아웃을 모를 때 전형적. CentOS 7 스톡 git 1.8.3.1이 컴퓨트 노드에 남아있을 가능성 높음(로그인 노드는 2.47로 업그레이드됨). 이전 캠페인(6110900~6124640)은 일반 저장소에서 실행해 구식 git도 처리 가능했음.

**검증:** worktree 자체는 멀쩡함 — 로그인 노드에서 방금 `git rev-parse HEAD` → `0363d8b` 정상, `.git/worktrees/worktree/` 메타 무손상. 프로브 잡 6140000 (cpu-48 git 버전 확인)을 넣어 대기 중.

**해결:** sbatch 머리에 GPFS 최신 git을 PATH 앞 고정 (예: `gitshim/git → ~/micromamba/bin/git` 심링크). 그 후 두 잡 재제출.

### 다음 작업
- 프로브 6140000 완료 대기 (git 버전 확인)
- PATH 수정 → gitshim 준비 후 sbatch 2개 재제출
- 착지 모니터링 재개 (Table 1b/1d, Figure 1, Supplementary S3 채우기)

## ▶️ 이전 스냅샷 (2026-08-29 12:40)

### Codex 2차 교차 검토 완료 및 push (`935f94e`)

Codex 최종 검토 2라운드 통과. ① 1차: Figure 1 pending 누락·새 수치 근거 미배포 → `33f3ab6`으로 수정. ② 2차: results/README.md 낡은 provenance 서술 2건·원고 잡 역할 오기 → `935f94e`으로 수정. 최종 4개 커밋(`14c55f8`, `917a3d8`, `33f3ab6`, `935f94e`) push 완료.

**검증된 것:** deposited BED 3종 provenance SHA와 바이트 일치, manifest 28/28 해시 통과, 원고·README 정합성 통과.

**대기 중:** Slurm 6129408/6129410 (native p100 P/B/F/H 4개, identity sweep 5개) 착지 → Table 1b/1d, Figure 1, Supplementary S3 채우기. Issue #11 판독(400건, 저자 2인) 진행 중.

### 다음 작업
- 6129408/6129410 착지 모니터링
- 착지 후 표·Figure 갱신
- Issue #11 판독 수렴

## ▶️ 이전 스냅샷 (2026-08-29 10:56)

### maize 단일 관례 재채점 완료 및 Table 1–3 교체 (`14c55f8`)

Codex가 원고 Table 3B/3C의 진정한 채점 관례를 재구현(원시 호출 기반, gap-merge 아님)해 재채점 스크립트를 완성했습니다. 각 표마다 "구 BED로 원고값이 재현되는가"를 먼저 검증한 뒤 교체했습니다.

**검증 결과 확정:**
- Table 1a/1c: 경쟁 도구 4행이 재생성 채점표와 소수점까지 일치 → 동일 관례 확인
- Table 2: colcen 게이트 통과 (구 행 재현) → REGEN 채택
- Table 3B/3C: 구 BED로 **80.32 / 50.49 / 275bp / 925bp / 58.55% 재현** (경쟁 5행도 전부 일치)

**주요 발견 — CentC 점수 뒤집힘:** 재생성 후 BWTandem 58.55%는 TRASH-template 58.68%에 **0.13점 뒤집습니다.** "최고 커버리지", "0.35점 리드" 주장을 초록·Section 3.3.3·line 94에서 3자 동률 서술로 정정, Table 3C 재정렬 완료. TR-1 offset도 771 bp (TRF 885 아래) 수정.

**측정 관례 정정:** 원고 Methods에서 "*Peak memory는 sacct MaxRSS*" 명시(GNU time 과소계상)했으므로 rusage 값 17.37 GB를 sacct 값 **28.08 GB**로 정정.

### 다음 작업
- Table 1b 방식(native p100 vs v2 post-hoc) 결정·재실행 — Codex 회신 대기
- Table 3B-b 산문 수치 교체 (knob180/TR-1/banded TRF 리드값)
- adjusted precision 개명 (lines 107, 113–123, 146–167, 303)
- Methods 수정 — gap-fill/tier5/환경 버전/Col-CEN 라벨
- Data Availability + 재생성 BED·manifest 배포

## ▶️ 이전 스냅샷 (2026-08-29 10:50)

### colcen 재채점 완료 및 측정 관례 함정 발견·수정

colcen validate 게이트 통과(구 행 84.59/1,380/57.17/99.73 정확히 재현). 중요 발견: 원고 line 82가 측정 관례를 명시—**"Peak memory for the BWTandem rows is the SLURM cgroup maximum (sacct MaxRSS)"**. GNU time은 멀티워커 실행 과소계상하므로 의도적으로 sacct 기준. 제가 처음 입력한 rusage 값(17.37 GB)을 관례에 맞는 sacct 값(28.08 GB)으로 정정했습니다.

**표 교체 완료 (모두 검증):**
- **Table 1a** — `4,014,108 | 80.53 | 50.50 | 79.45 | 43.34 | 42.74 | 893,480 | 12.6h | 28.08GB`. 경쟁 4행 소수점까지 일치, Unique Regions는 구 BED로 원고값 893,547 재현 후 신규 값 산출
- **Table 1c** — `102,926 | 3.43 | 64.60 | 13.39 | 29.25` 재현 확인
- **Table 2** — `160,883 | 84.54 | 1,533 | 57.08 | 99.72 | 0.67h | 2.16GB` 확정

**Codex 재채점 진행 중**: Table 3B/3C는 gap-merge 방식 채점기 아닌 원시 호출 기반 채점기 사용. 재채점 스크립트(검산 포함) 완성, 제가 실행 중.

### 다음 작업
- maize 표 재채점 결과 → Table 3B/3C 교체
- Table 1b 방식(native p100 vs v2 post-hoc) 결정·재실행
- adjusted precision 개명, Methods·초록·Data Availability 갱신

## ▶️ 이전 스냅샷 (2026-08-29 10:48)

### Codex 백그라운드 위임 + Table 1a/1c 행 교체 완료

Codex에게 4건 위임: ① Table 3B/3C 생성 스크립트·채점 관례 특정, ② TRASH 포함 전 행을 단일 관례로 재채점해 완성본 생산, ③ Table 1b 방식(native p100 vs v2 post-hoc) 결정·재실행, ④ "신규 BWT + 구 경쟁 도구" 혼합 표/문장 전수 점검. 원고 편집은 분리.

**적용 완료 (검증)**:
- **Table 1a** — BWTandem 행 교체: `4,014,108 | 80.53 | 50.50 | 79.45 | 43.34 | **42.74**`. 경쟁 도구 4행이 v2 채점표와 일치함을 선행 확인 후 시행.
- **Unique Regions 893,480** — 재계산·확정. 구 BED로 구 값(893,547) 정확히 재현해 계산 관례가 원고와 동일함을 검증.
- **Table 1c** — BWTandem 행 교체: `102,926 | 3.43 | 64.60 | 13.39 | 29.25`.
- ⚠️ **Runtime/Memory** 교체함 (21.86GB → 17.37GB): 구 값은 SLURM job accounting(COW 중복 부풀림 가능), 신규는 rusage 실측. 경쟁 도구 측정 방식 확인 필요.

**진행 중**: colcen 채점 재실행(Table 2용, `--validate` 게이트 포함). Codex 회신 대기.

### 다음 작업
- colcen 채점 완료 → Table 2 교체
- Codex 회신 → Table 1b/3B/3C 처리 + 혼합 표 전수 검토
- Runtime/Memory 측정 관례 확인 → 21.86GB 유지 또는 17.37GB 확정

## ▶️ 이전 스냅샷 (2026-08-29 07:58)

### 패키징 결함 발견·수정 완료 (`cd37bd7`)

#14 CI 실패의 원인 확정: setuptools가 src-layout으로 오인해 `src` 패키지가 휠에서 누락, 최상위 모듈 오염(main, models, coverage). 진입점 실패 + 배포판 덮어쓰기 지뢰.

**수정 완료:** `packages = ["src", "src.c_extensions"]` 명시, 로컬 검증 통과. 덤: license deprecation(2027년 중단 예고) 정리, build residue `.gitignore` 추가. CI 감시 중.

**저자 판단 필요:** setup.py의 "릴리스 때 `src` → `bwtandem` 변경" 주석이 provenance 재현성과 충돌. 재생성 3건 명령이 `python -m src.main ...`으로 박제. 선택: (a) `src` 유지(재현성 보존) 또는 (b) 변경 + provenance 주석 병기. 태그 직전 결정.

**기타:** #11 판독 0/400건, bwtandem SLURM 잡 없음.

### 다음 작업
- CI 완료 대기 (cd37bd7)
- 패키지명 선택사항 검토 → 태그 직전 결정
- #11 판독 진행도 확인

## ▶️ 이전 스냅샷 (2026-08-29 00:41)

### tantan-w2000 v2 채점표 재실행 완료 및 패키지 메타데이터 커밋

**v2 채점표** — tantan `-w 2000` 포함 재실행 완료 (exit 0). 결과 `score_table1_regen_v2.txt`에 저장, STRATUM 101–2000 행 포함 확인.

**패키지 메타데이터 커밋** (#14):
- LICENSE (MIT, 작성자 결정 2026-08-28)
- CITATION.cff
- pyproject.toml (Cython 빌드 설정)
- setup.py (wheels에 _accelerators 동봉)
- .github/workflows/ (CI)

**푸시 완료**: `0363d8b..0c3c982` (원본 `perf/exp1-human-sensitivity`).

### 다음 작업
- v2 채점 결과 검증 (STRATUM 101–2000 항목, tantan overlap 재확인)
- maize triage 6128318 완료 모니터
- 원고 수정 항목 정리 (tier5 원인 명시, 1bp 규칙 표기)

## ▶️ 이전 스냅샷 (2026-08-29 00:21)

### Codex 적대적 검토 완료 및 v2 채점표 재실행 (2026-08-29)

Codex 보고서(`ADVERSARIAL_REVIEW.md`) 핵심 6건 판정:
1. **tantan w2000 누락** (필수 수정) — period >100 37,700건 미포함 → EXTRA 행으로 v2 재실행 필수
2. **ADJUSTED PRECISION** — "catalog-or-cross-tool overlap support rate"로 표기 강등, precision 해석 금지
3. **tier5 bp −77M** — #12 N==N 결함 제거의 실제 효과, 30Mb급 병적 콜 202개 소멸
4. **period cap 취약** — 현재 오염 없음, 준비물 커밋 시 schema 명시 개선 검토
5. **1bp 규칙 관대함** — 원고 methods 표기 필수
6. **BWT-only 809,886** ✅ — #11 감사 모집단 유효 **확정**

**v2 채점표 재실행 진행 중** (`--extra tantan-w2000:...` 플래그, workdir v2, ~13분 예상). 병행: maize triage(6128318) 감지기 가동.

### 다음 작업
- v2 채점표 + maize triage 완료 모니터
- v2 표 검증 → maize 채점 → 준비물 커밋(LICENSE, CITATION.cff, pyproject.toml)
- 원고 수정 항목 정리(명칭 강등·tier5 원인 명시·1bp 규칙 표기)

## ▶️ 이전 스냅샷 (2026-08-29 00:13)

### 세 게놈 재생성 캠페인 최종 착지 (2026-08-28 21:30)

**maize 착지 CLEAN** (잡 6124640, exit 0): 2,406,800콜, 15.8h, 피크 20GB. Acceptance gate **전 항목 PASS** (스키마 2,406,800행 검증 완료). colcen·human·maize **세 게놈 모두 게이트 통과**.

**human triage 신호** (잡 6127752, 완료): tier 1–3 변화는 소폭(+353/+50/+27콜)이지만, **catch-all bp 106M→29M 급감 감지**. #12 분절 수정의 실제 효과인지 집계 버그인지 판별 필요.

**maize triage 제출됨** (잡 6128318, 32G/4h): 완료 감시기 가동 중.

### Codex 적대적 검토 시작 (새 스레드)

5개 항목: (1) recall/precision 구현 재확인 (2) corroboration 순환 논증 (3) tantan `-w 100` 함정 (4) catch-all bp 급감 판별 (5) BWT-only 81만 건 정합성.

### 다음 작업

- maize triage 6128318 완료 대기
- Codex 검토 결과 수령
- catch-all bp 급감 원인 판별 (집계 vs #12 실제 효과)
- Codex/maize 결과 도착 후 준비물(LICENSE 등) 일괄 커밋

## ▶️ 이전 스냅샷 (2026-08-28 20:39)

### 전체 레인 상태 점검 완료 (2026-08-28 20:40)

**Codex human 채점 — 완료됨 (13:14, exit 0)**
- regRecall 80.53% / regPrec 50.50% (GT 178만 리전 기준)
- 보정 정밀도 79.5–79.8% — #11 감사 중인 BWTandem-only 81만 건이 이 41% 미확인 콜의 진위 결정
- 스트래터지 `.strat.bed` 0바이트 정상 (ULTRA/tantan은 기간 101–2000 콜 없음)

**Acceptance gate 직접 실행: 전항목 PASS**
- 4,014,108행 (SLURM/provenance/카운트/SHA/스키마) 클린
- gate 패스 → triage 잡 6127752 제출 완료 (32G/4h)

**#11 감사 판정 예제 5장 생성 완료** (`audit11/examples.html`)
- EX1–2: SUPPORTED (사다리형 대각선 유지) / EX3: UNSUPPORTED (무작위) / EX4: UNSURE (30% 변이 경계) / EX5: UNSUPPORTED (<3 copies)
- EX2는 원고에서 보던 12% SNP 범벅 사례와 직접 비교용

**Maize 6124640: 정상 진행 중** (9번째 염색체 Tier 2, 로그 20:30까지)

### 다음 작업

- triage 6127752 착지 모니터 (예상 21시 전후)
- maize 6124640 완료 모니터 
- 두 잡 완료 후 준비물 커밋 (LICENSE, CITATION.cff, pyproject.toml)
- #11 최종 판독 패키지 제시 (블라인드 시트 + v2 도플롯 400장)

## ▶️ 이전 스냅샷 (2026-08-28 16:20)

### #11 감사 패키지 v2 플롯 생성 시작

**v2 플롯 스크립트 작성 및 실행** — `make_dotplots_v2.py` 작성하여 flank를 max(5×period, 100)bp(상한 500bp)로 확장, 콜 경계를 빨간 선으로 표시. 시트의 좌표·period만 사용하므로 블라인딩 유지 확인.

**레이아웃 강화**:
- v2 갤러리 상단: "SNP가 흩어져도 배경 위로 사다리형 대각선이 이어지면 SUPPORTED" 판정 기준 명시
- unit-shift 패널: 무작위 기대치(0.25) 점선 추가 → flank 배경과 비교 기준선 제공

**백그라운드 렌더링 진행 중** — 400장 × 4 패널/장이라 계산량 커서 병렬 실행. 완료되면 자동 알림.

**생성 예정 파일**:
- `audit11/dotplots_v2/` (400장)
- `audit11/review_v2.html` (v2 갤러리)

원본 `dotplots/`, `review.html` 유지 → 판독은 v2로 진행, methods 기록.

### 다음 작업

- v2 렌더링 완료 알림 대기
- 완료 후 판독 패키지 최종 확인 및 사용자 제시

## ▶️ 이전 스냅샷 (2026-08-28 13:55)

### 병렬 4레인 동시 시동 완료 (2026-08-28 14:10)

4개 독립 축 동시 가동:
- **레인 1**: maize 재생성 (잡 6124640) — 실행 중
- **레인 2**: human 채점 재실행 (Codex `task-mtddlbta-7h49mn`) — 진단 진행 중
- **레인 3** (신규): #11 감사 패키지 빌드 — 백그라운드 실행 시작
- **레인 4** (신규): T3 arm 귀속 분석 — 에이전트 시작

### #11 감사 패키지 파이프라인 (레인 3)

입력 확인 완료: hg38_primary.fa(.fai), 도구 BED 5종 고정 기준선, regen_human.bed 최종본(CLEAN).

파이프라인 시작: `run_audit_prep.sh` 백그라운드 실행 (`/data/gpfs/assoc/pgl/devel/exp1_human/audit11_prep/run.log`).

단계: (1) `bedtools intersect -v` → BWTandem-only 콜 추출 (2) seed 20260827로 400건 층화 시트+봉인 키 생성 (3) dot plot 400장 렌더링(self + unit-shift identity). matplotlib/numpy 환경 확인됨.

### T3 arm 귀속 분석 (레인 4)

6111318(narrow) / 6111319(wide) 산출물 검증: tolerance_ratio(0.0400→0.0204) 귀속이 맞는지 직접 확인, 원고 쓸 문장 생성.

### 다음 작업

- 레인 3 완료 → 판독 패키지(400건 시트+dot plot) 제시
- 레인 4 완료 → T3 귀속 문장 확인
- maize/human 착지 → 재채점 workflow 진행
- 전체 완료 후 untracked 준비물 커밋 해제(LICENSE, CITATION.cff, pyproject.toml)

## ▶️ 이전 스냅샷 (2026-08-28 13:51)

### 저자 결정: MIT 라이선스 확정

**LICENSE 선택**: MIT (copyright 2026 Won C. Yim) — 이슈 #14 결정 완료.

**반영 완료**:
- `LICENSE` — MIT 전문, untracked 준비됨
- `CITATION.cff` — `license: MIT` 활성화
- `pyproject.toml` — `license = {file = "LICENSE"}` 주석 해제

**상태**: 세 파일 모두 untracked. 추적 트리 클린 유지 → dirty-gate 안전, 핀(`0363d8b`) 검증 일치. maize 잡 완료 후 준비물 일괄 커밋 때 함께 병합.

**이슈 #14**: [코멘트 등록](https://github.com/wyim-pgl/bwt-algorithm/issues/14#issuecomment-5457629996) — 저자 몫 병목 2건 중 1건 해결.

### 다음 작업

1. maize 잡 6124640 착지 모니터
2. human 채점(Codex 재진단) 완료 → acceptance_gate.sh 실행
3. 두 잡 완료 후 준비물 일괄 커밋(`LICENSE`, `CITATION.cff`, `pyproject.toml`, `pyproject.lock`)
4. #11 판독 패키지 완성(400건 시트 + dot plot)

## ▶️ 이전 스냅샷 (2026-08-28 13:11)

### 사용자 확약 메모리 기록 (#11)

메모리 파일 추가: `issue11-audit-user-commitment.md` — 400건 블라인드 시트 + dot plot 준비물이 완성되면 사용자가 직접 판독하겠다는 확약 기록. 논문 최우선 목표와 연결.

### 준비 워크플로우 계획 확정

Codex 진단 완료 후 자동 진행 순서:
1. `bedtools intersect -v` → BWTandem-only BED 생성
2. `sample_specificity_audit.py` (seed 고정) → 400건 층화 시트 + 봉인 키
3. dot plot 이미지 사전 렌더링까지 → 판독 준비 패키지 완성

### 현재 진행 상황

- **Codex 백그라운드**: `task-mtddlbta-7h49mn` 진행 중 (human 6110901 채점 재진단)
- **maize 6124640**: PENDING 상태 모니터 중

### 다음 작업

1. Codex 진단 결과 대기 → intersect 실행
2. maize 착지 후 동일 workflow 태우기
3. 두 잡 모두 완료 후 untracked 준비물 커밋 해제

## ▶️ 이전 스냅샷 (2026-08-28 13:03)

### 랩위키 Codex 1차 진단 가이드 신설 완료

**`guide/codex-first-triage.md` 신설** (`a900d7e`, push 완료) — 잡 실패 또는 오탐 신호 발생 시 **상태 변경 전에 Codex를 1차 진단자로 쓰는 절차** 문서화. `/codex:rescue` 세션 경로, CLI (`codex exec`/`codex resume`) 사용법, 필수 컨텍스트 체크리스트(핀 커밋·캠페인 제약·로그 경로), Codex 권고도 실측 검증 후 채택하는 규칙(실제 사례: Sylvan MaxRSS 반려), 그리고 실제 오탐 4건 표 추가(`SBATCH_TIMELIMIT` 미스터리, `.extern` OOM artifact 등). index.md·log.md 등재 완료.

### 논문 최우선 목표 메모리 기록

프로젝트 최우선: **acceptance까지의 critical path** 문서화. (1) 재생성 캠페인 완결 (colcen ✅ 판정 완료 / human/maize 진행 중) → (2) 재채점 + T3 arm 귀속 → (3) 준비물 커밋(pyproject/CITATION/gate/스크립트) → (4) 저자 결정 2건(LICENSE 선택, human call 블라인드 감사). 병목: #14·#11은 코드 대체 불가 사항.

### 다음 작업

1. Codex 백그라운드 태스크(`task-mtddlbta-7h49mn`) 결과 대기 — human 재실행 진단 중
2. maize 잡 6124640 착지 모니터 → 동일 acceptance gate·재채점 workflow
3. colcen/human/maize 세 잡 모두 완료 후 untracked 준비물 커밋 해제

## ▶️ 이전 스냅샷 (2026-08-28 00:10)

### maize 재제출 및 dirty-gate 사건 복구

**6110902 maize FAILED**: dirty-gate 게이트 작동 — `pip wheel .` 빌드 중 추적 파일 `src/_accelerators.c`(cython 생성물)이 재생성되면서 트리가 더티 상태. `run_with_provenance.sh`의 23초 차단이 예상대로 작동했고, 오염된 실행 방지 완료. 

**복구**:
- `git restore src/_accelerators.c` + `build/` 제거 → 트리 클린 복원
- HEAD = 핀(`0363d8b`) 검증 일치
- maize **재제출 6124640** (62.5G/30h, `SBATCH_TIMELIMIT` unset 확인)

**교훈**: 캠페인 기간 중 "untracked 파일만 안전"에 **"추적 생성물(.c)을 건드리는 빌드 명령 금지"** 규칙 추가 필요 (`pip wheel`, in-place cython 빌드 해당). 다음 wiki 갱신 때 기록.

### 현재 실행 상황

- **6110900 colcen**: ✅ COMPLETED (피크 1.95G, 31분) — acceptance gate 통과 대기
- **6110901 human**: RUNNING (~계속 진행)
- **6124640 maize**: PENDING (backfill 대기)

### 다음 작업

1. colcen 결과 → acceptance_gate.sh 실행 → triage → score_colcen.py --validate
2. 6124640 maize 착지 후 동일 workflow
3. 두 잡 완료 후 untracked 준비물 커밋 해제

## ▶️ 이전 스냅샷 (2026-08-27 13:55)

### Col-CEN 재생성 착지·채점 완료 — headline 수치 안정

- **잡 6110900 COMPLETED** (40m11s, MaxRSS 1.37GB, 161,330 calls). acceptance gate 전 항목 PASS (게이트 스크립트의 `A && B || C` 버그 1건 수정 후).
- **채점**: `--validate` 게이트 통과(구 BED 84.59/1,380/99.73 재현) + 신 BED: **coverage 84.54%(−0.05pp) / CEN180 band 1,533(+153) / recall 99.72%(−0.01pp)**. banded count만 실질 변화 — satellite 분할(649 whole-block→801 runs)이 원인. Col-CEN은 N 0bp라 모호염기 수정 효과 없음. Tier3 +65콜은 tolerance_ratio — T3 arm(6111318/9)이 귀속 예정.
- **이슈 트래커 정리 완료**: 20건 전부 코멘트, **9건 close**(#8/#10/#12/#13/#18/#19/#20/#21/#23), 2건 재제목(#9/#17), colcen 결과를 #5/#7에 게시.
- **untracked 준비물** (human/maize 시작 후 커밋): `score_one_to_one.py`(#15, fixture 검증), `sample_specificity_audit.py`(#11 블라인드 샘플러), `pyproject.toml`+`CITATION.cff`+CI(#14, LICENSE만 저자 대기), `results/comparator_baselines.md`(#16), `acceptance_gate.sh`, `t3_attrib_colcen.sbatch`.

### 다음 작업
1. human(6110901)/maize(6110902) 착지 → 동일 게이트→triage(32G sbatch)→재채점, T3 arm 결과 귀속
2. 두 잡이 **시작**하면 커밋 제약 해제 → untracked 준비물 커밋, #17 README 정합
3. LICENSE 선택(저자) → #14 마감, 태그

## ▶️ 이전 스냅샷 (2026-08-27 13:17)

### 코덱스 상담 결과 반영 및 게이트·Tier 3 arm 구축

코덱스 판정: 실행 중인 잡은 그대로 두고, 착지 후 도구와 스코러에서 문제 해결. triage가 게이트가 아님(대조·검증 없음) 발견 → **`acceptance_gate.sh` 신작** — 4중 합의 강제 (SLURM+provenance+로그 계산+SHA 재계산+스키마). 커밋 제약(핀)과 더티 제약(래퍼) 딜레마는 untracked 파일로만 준비하면 양쪽 검사 모두 통과함을 확인.

Q4 미스터리 해결: 13d23h의 원인은 셸 환경 `SBATCH_TIMELIMIT=13-23:00:00` 하드코딩 — SLURM이 `SBATCH_*` 환경변수를 지시자보다 우선. unset 후 T3 arm 제출 시 6:00:00 정확히 걸림 (실증).

### Tier 3 귀속 통제 arm 제출 (6111318/6111319)

full 재실행 대신 Tier-3 단독 arm 2개 제출 (colcen, 16G/6h):
- **6111318** (narrow): `--max-period 2000` (기존)
- **6111319** (wide): `--max-period 100000` + 사후 필터
같은 커밋 고정, tolerance_ratio 효과만 분리.

### 현재 실행 상황

- **6110900 colcen**: RUNNING (~12분 경과, ~30분 내 완료 예상)
- **6110901 human / 6110902 maize**: PENDING (62.5G/30h)
- **6111318/6111319 T3 arm**: PENDING (backfill 대기)

### 다음 작업

1. colcen 착지 즉시 acceptance gate → triage → `score_colcen.py --validate` 실행  
2. human/maize 잡 착지 후 (스코러 배선 커밋 금제 풀림) 동일 workflow  
3. pronghorn 가이드에 `SBATCH_TIMELIMIT` gotcha 기록

## ▶️ 이전 스냅샷 (2026-08-27 13:06)

### 슬룸 잡 메모리/타임리밋 재조정 및 모니터링

세 잡 확인 시 과대 요청 발견. `scontrol`로 실측 기준 재조정:
- **6110900**(colcen): 120G/13d23h → **8G/4h** (피크 1.95G, 31분)
- **6110901**(human): 120G/13d23h → **62.5G/30h** (피크 24.8G, 16h)
- **6110902**(maize): 120G/13d23h → **62.5G/30h** (피크 22.4G, 15.5h)

파티션 혼잡(대기 183/실행 63) 속에서 colcen 8G/4h로 backfill 진입 가능. 부수 발견: 로그인 노드 `scontrol` 전부 `libreadline.so.6` 부재로 고장, conda shim + `LD_LIBRARY_PATH` 우회 적용. MinMemoryNode는 MB 단위만 수락(8G 아닌 8000).

### 다음 작업

1. 모니터 가동 중 — colcen 착지 감지 시 `triage_regen.py` + `--validate` 실행
2. human/maize 착지 후 동일 triage → 1a/1c/3A/3B/3C 재채점
3. pronghorn 가이드에 MinMemoryNode MB 단위 gotcha 기록

## ▶️ 이전 스냅샷 (2026-08-27 12:19)

### 핸드오프 마커 규약 성문화 및 resume.md 적용

위키 `guide/handoff-hygiene.md` 신설. bwtandem의 실제 오독 2건(재검증 후에도 남은 "15건 미확인", 푸시 후에도 남은 "미커밋")을 계기로 스냅샷 스택형 핸드오프 노트의 superseded 주장 표시 규약 문서화. 마커 4개(`❌ SUPERSEDED` / `⚠️ CAUTION` / `✏️ PARTIAL` / `🧊 ARCHIVED BELOW`) 정의, 모두 날짜+포인터 필수. CLAUDE.md 규칙 연결, `index.md`·`log.md` 등재 (fcd3ccc 푸시, 링크 검사 0 broken).

resume.md 적용: 머리말 갱신 (현재 단계 캠페인 실행 중), 읽는 법 추가, 격리 배너 도입, 위험 지점 10곳에 마커 적용. 재생성 잡 3개(6110900–2) 여전히 PENDING(Priority).

### 다음 작업

1. 잡 완료 대기(모니터 중) → colcen 착지 즉시 `triage_regen.py` 실행  
2. human/maize 착지 후 triage + 신 BED 재채점  
3. #14 (LICENSE) — 저자 결정 대기

## ▶️ 이전 스냅샷 (2026-08-27 11:58)

### 코덱스 권고 채택 완료 및 캠페인 커밋 핀 갱신

코덱스 세션 결과 5가지 결정 수신. 채택 방침:
- **#16** (비교 도구): 재실행 없음 — 기존 BED 고정 기준선, 탐색 범위 병기
- **#11** (BWTandem-only 감사): 400개 블라인드 층화 수동 감사는 저자(2인) 직접 실행 필요
- **스코어링**: 구 BED로 `--validate` 게이트 → 신 BED를 별도 라벨로 재채점 (1a/1c/2/3A/3B/3C)
- **태깅**: full hash로 충분, #14(LICENSE) 이후
- **triage**: `triage_regen.py` Chr4 스모크 검증 ✅ (tier 4만 −1콜/−9,173 bp — 수동 실측과 일치)

triage 스크립트 커밋(`0363d8b`) 후 HEAD가 핀(`d0417f5`)을 지나간 사고 발견 → 잡 시작 시 자기거부할 상황이었으나, 잡이 아직 PENDING인 동안 핀을 갱신(`0363d8b`). 검증: `git diff d0417f5..0363d8b -- src/`가 0줄이므로 검출기 계보 동일.

### 다음 작업

1. 세 잡(6110900/6110901/6110902) 완료 대기 (모니터 가동)
2. colcen 착지 즉시 `triage_regen.py` + `score_colcen.py --validate` 실행
3. human/maize 착지 후 동일 triage → 1a/1c/3A/3B/3C 재채점

---

> 🧊 **ARCHIVED BELOW (2026-08-27).** 이 선 아래는 전부 **과거 스냅샷과 감사 당시 기록**이다.
> 현재 상태는 상단 최신 스냅샷과 wiki `projects/bwt-tandem-repeat-finder/` 만 신뢰할 것.
> 아래 문장들은 작성 시점에는 참이었으나 이후 뒤집힌 것이 많다 — 특히 위험한 곳에는
> 개별 마커(❌/⚠️/✏️)를 달아 두었다. 마커 없는 문장도 **인용 전 최신 스냅샷과 대조**할 것.

## ▶️ 이전 스냅샷 (2026-08-27 12:10)

### 벤치마크 재생성 캠페인 가동 (#3/#4/#5)

- SLURM 잡 제출: **6110900**(colcen, ~31분) / **6110901**(human, ~12-16h) / **6110902**(maize, ~15h). 파티션 `cpu-s1-pgl-0`, t2/p1-2000, v2.2 gate 노브(원 remeasure와 동일), catch-all 0.72/3 (human·maize), OFF (colcen).
- 캠페인 커밋 핀: `/data/gpfs/assoc/pgl/devel/exp1_human/regen/EXPECT_COMMIT` = `0363d8b` (HEAD 불일치 시 잡이 스스로 거부). 출력: `regen/regen_<genome>.bed` + `.provenance.json` (full SHA-256, 커널 실측).
- 신규 커밋: `d0417f5`(캠페인 sbatch), `0363d8b`(triage_regen.py — Chr4 쌍에서 tier4만 −1콜/−9,173bp로 수동 실측과 일치 확인). 모두 푸시됨.
- **Oropetium(#18) 캠페인 제외** — 저자 결정.

### 코덱스 권고 (세션 01a04491-fe10-7c10-9c87-42a152632fc1) — 채택 방침
1. **#16**: 이번 라운드 비교 도구 재실행 없음. 기존 BED = 고정 기준선, 탐색 범위 병기. (human TRF 2000bp 재실행은 6.6일 미완으로 실행 불가 확인)
2. **#11**: human BWTandem-only 콜 400개 블라인드 층화 수동 감사 (period 4계층×100, dot plot+unit-shift, 2인, Wilson CI) — **저자 실행 필요**
3. 스코어링: `score_colcen.py --validate`를 구 BED로 게이트 → 신 BED 별도 라벨 채점, 1a/1c/2/3A/3B/3C 재채점. 신 BED numeric tier 포맷은 스코어러 호환 확인됨
4. 태깅: full hash로 충분, 태그는 #14(LICENSE) 후 — 실행 비차단
5. triage 순서: 스키마 → N 효과 → satellite → catch-all → Tier3 band → 재채점 (`triage_regen.py`가 1-5 구현)

### 다음 작업
1. 잡 완료 대기 (모니터 가동 중) → colcen 착지 즉시 `triage_regen.py` + `score_colcen.py --validate`(구 BED) → 신 BED 채점
2. human/maize 착지 후 동일 triage → 1a/1c/3A/3B/3C 재채점
3. #14 (LICENSE/CITATION/pyproject) — 라이선스 저자 결정 대기

## ▶️ 이전 스냅샷 (2026-08-27 11:55)

### 재생성 캠페인 시작 (#3/#4/#5)

3개 잡 제출 완료: Job 6110900(colcen), 6110901(human), 6110902(maize), 모두 큐 대기 중.

**구성:** 원 v2.2 remeasurement 설정 (gate, threads 2, period 1-2000), catch-all id 0.72, MIN_COPIES 3(human/maize)/OFF(colcen), run_with_provenance.sh 경유, 캠페인 커밋 검증(`EXPECT_COMMIT=d0417f5`), 모든 출력 SHA-256 기록, 입력 FASTA 3개 확인(3.0G/2.1G/129M).

**병행:** 모니터 가동(완료/실패 시 알림), 코덱스 상의 중(#16/#11 등 4개 결정 대기), Oropetium(#18) 제외 기록.

### 다음 작업

1. Col-CEN 결과 검증 (예상 ~31분 후)
2. 코덱스 권고 수신 → 스코어링 정책 결정
3. 전체 재생성 완료 후 인수 검증

## ▶️ 이전 스냅샷 (2026-08-27 11:45)

### 세 과제 완료·푸시 (`1cb7f48`)

네이티브 경로 N 처리 (`0f9f92a`), #12 구조 수정 (`c375456`), 벤치마크 준비 (`1cb7f48`) — 모두 실행 및 검증 완료.

**1️⃣ 네이티브 경로 N 처리:** 2-bit 패킹이 N과 `$`를 A 코드(0)로 매핑해 packed/raw 발산이 발생 — `N==A`가 매치로 계산. 규칙 적용: "모호 염기는 자기 자신과도 매치하지 않는다." Cython(_total_mismatches, extend, scan_unit), C(tier1 run 확장), Python fallback 3개 계층 수정. 신규 테스트 18건, 수정 전 빌드에서 7건 실패 확인.

**2️⃣ #12 구조 수정:** satellite 블록이 가장자리 3개 윈도우만으로 전체(최대 100 kb)를 주장하던 과잉. `_segment_periodic_block()`이 O(n) cumsum으로 증거 기반 분할 수행. **Chr4 실측: 비주기 실서열 7,675 bp trim** (N 제외) — `fe3f9cd`의 N 전용 수정이 보지 못한 과잉 주장을 네이티브 수정이 잡음.

**3️⃣ 벤치마크 준비:** environment.yml/core.lock.yml 고정(3.11.15), `score_colcen.py` --validate 게이트 실행 가능(84.59/1,380/99.73 재현), `run_with_provenance.sh` 더티 거부·rusage·SHA-256, `results/manifest.sha256` 신설. **179 passed**.

### 다음 작업

1. 벤치마크 재생성 실행 — 저자 결정 4건(코퍼스, 도구 파라미터 정책, specificity truth, 라이선스) 필요
2. #14 (CI/pyproject.toml)
3. 이슈 정리: #21 close, #19/#20/#23/#12 반영 확인

## ▶️ 이전 스냅샷 (2026-08-27 11:45)

### 총 7커밋 푸시 완료 (`perf/exp1-human-sensitivity` @ `1cb7f48`)

2차 배치 (코덱스 적대적 검토 반영 + 3대 과제):
- `3d3a665`: 코덱스 검토가 확인한 4결함 수정 — `_cumsum_windows` off-by-one(마지막 window 누락), `_str_autocorr_identity`의 N==N, `{"all","typo"}` 우회, `--tiers all,tier1` 회귀. Tier 3 "출력 중립" 주장 철회(tolerance_ratio 0.0400→0.0204 실측). Venn·stress 회귀 테스트 추가.
- `0f9f92a`: **네이티브 경로 N 처리** — 2-bit 패킹이 N/$를 A 코드로 매핑해 packed vs raw 발산(A*40+N*40이 packed에서 한 콜), tier1_scan.c 런 확장, hamming/extend/anchor 전부. 규칙: 모호 염기는 자기 자신과도 불일치. 18테스트, 수정 전 빌드에서 7건 실패 확인.
- `c375456`: **#12 구조 수정** — satellite fill이 3개 샘플 윈도우로 블록 전체(최대 100kb)를 주장하던 것을 O(n) 증거 기반 분할로 교체. Chr4 실측: 비주기 실서열 7,675bp trim(Chr4:3,990,030), 콜 수 4,717 유지.
- `1cb7f48`: **#7/#13 준비** — environment.yml을 실측 환경(3.11.15)으로 고정+osx-64 제거, `environment.core.lock.yml` 신설, `score_colcen.py` fail-closed + `--validate` 게이트(실데이터로 84.59/1,380/99.73 재현 확인, exit 0), `run_with_provenance.sh`(더티 트리 거부, rusage 실측, full SHA-256), `results/manifest.sha256`.

**전체 스위트 179 passed.** 코덱스 세션: `01a04460-2f8d-7a21-9512-6b36243e4901`.

### 다음 작업
1. 벤치마크 재생성 실행 — `run_with_provenance.sh`로 클린 커밋에서 GRCh38/Mo17/ColCEN (#3/#4/#5). **저자 결정 필요:** 코퍼스 선택, 비교 도구 파라미터 정책(#16), specificity truth(#11), Oropetium 존치(#18).
2. #14 (LICENSE/CITATION/pyproject/CI) — 라이선스 선택은 저자 결정.
3. 이슈 트래커 정리: #21 close 가능, #19/#20/#23은 이번 커밋들로 해소(검증 포함), #12는 구조 수정까지 반영됨 — 단 3-window 지명 방식 자체는 유지되므로 코덱스 재확인 권장.

## ▶️ 이전 스냅샷 (2026-08-27 11:10)

### 코덱스 적대적 검토 완료, 4개 결함 발견·수정·테스트 추가

코덱스 리뷰(`task-mtbtu7bl-9t7og1`) 결과 수신. 4개 실제 결함 확인·수정:

1. **`_cumsum_windows` off-by-one**: `m = eq.size - window` → 마지막 full window 누락. 실측 28 vs 29, `eq.size==window`일 때 `None` 반환. catch-all·Tier 2 양쪽 영향.
2. **`_str_autocorr_identity`가 여전히 N==N 세기**: `refine_repeat`의 primitive-period 보정과 satellite period 탐색 2곳에 도달. `fe3f9cd`의 제목이 범위 과장.
3. **`{"all","typo"}` 우회**: `all`에서 early return → 집합 순회 순서에 따라 전 tier 조용히 켜짐. #19 미해결.
4. **`--tiers all,tier1` 거부**: `34a66fd`에서 추가한 회귀. "unknown tier(s): all; choose from ... or all" 자가당착 메시지.

**과장 바로잡음**: Tier 3 "출력 중립" 주장은 `max_period` 한계에 의존. `tolerance_ratio = 0.02 + 0.02 * (max_period/100000)` → 기본값 2000에서 0.0400 → 0.0204. Chr4 바이트 동일은 그 한 설정에서만 성립.

**테스트 추가**: `test_venn_binning.py`(7건), `test_stress_crash_accounting.py`(2건). 각 그룹을 수정 되돌린 상태에서 검증—Venn 4건, stress 1건, autocorr 4건 실패 확인. **159 passed.**

**커밋**: `3d3a665` (fix: act on adversarial review...) 푸시 완료. `wyim-pgl/bwt-algorithm` 원격 동기.

### 다음 작업

1. 남은 네이티브 경로 N 처리 (#12의 구조적 문제: `extend_with_mismatches`·`anchor_scan_boundaries`·Tier 1 스캐너)
2. 벤치마크 재생성 (GRCh38·Mo17·Arabidopsis) — Tier 3 파라미터·Venn binning·`--mask soft/both` 영향 반영 필수
3. 남은 이슈 (#7·#9·#11·#13·#14·#15·#16·#17·#18·#22) 범위 확정 및 처리

## ▶️ 이전 스냅샷 (2026-08-27 11:01)

### 3개 커밋 푸시 완료, 코덱스 적대적 검토 대기

`perf/exp1-human-sensitivity` 브랜치에 3개 커밋 fast-forward 푸시 완료:
- `a163ad9`: TRF 좌표·period 인자·멀티프로세스 fail-closed·무효 tier 거부·Tier 3 period 연결
- `34a66fd`: #19(혼합 tier)·#20(크래시 케이스)·#21(fixture staticmethod)·#23(Venn 경계) 고정
- `fe3f9cd`: **#12 CRITICAL** — 자기상관 N-신호 버그 전수정 (4개 소비자 모두)

원격(`wyim-pgl/bwt-algorithm`) HEAD: `fe3f9cd`로 업데이트됨. 로컬/원격 동기.

`archive/` 폴더는 untracked 그대로 유지 (#22 보관 정책 미결정).

코덱스 검토(`task-mtbtu7bl-9t7og1`)는 실행 중:
- catch-all 분모 유지 시 `CATCHALL_MIN_IDENTITY` 의미 변경 여부 검증
- 0.8 기본값 방어 가능성
- N-신호 처리 경로 잔존 여부 (`motif_utils._str_autocorr_identity` 등)

### 다음 작업

1. 코덱스 적대적 검토 결과 수신 대기
2. 필요시 후속 수정 커밋 (검토에서 이슈 발견 시)
3. 벤치마크 재생성 (GRCh38·Mo17·Arabidopsis) — `fe3f9cd` 로직 반영 필수
4. 남은 이슈 처리 (#7·#9·#11·#13·#14·#15·#16·#17·#18·#22)

## ▶️ 이전 스냅샷 (2026-08-27 10:59)

### 3개 커밋 완료, 코덱스 검토 진행 중

포팅 후 3개 커밋 추가:
- `a163ad9`: TRF 좌표·period 인자·멀티프로세스 fail-closed·무효 tier 거부·Tier 3 period 연결
- `34a66fd`: #19(혼합 tier)·#20(크래시 케이스)·#21(fixture staticmethod)·#23(Venn 경계) 고정
- `fe3f9cd`: **#12 CRITICAL** — 자기상관이 `N==N` 증거로 세던 것, 4개 소비자 전부 고정

테스트 추가: 회귀 테스트 17건 (결함 주입 검증), 143 passed.

**실측 변화:** `fe3f9cd`가 Chr4에서 1개 콜 제거 (72.1%가 N인 갭, 3,955,856–3,957,243 지점). 실제 갭이 많은 GRCh38·Mo17에서는 훨씬 클 것으로 예상 → **원고 벤치마크 재생성 필수.**

미결정: #7·#9·#11·#13·#14·#15·#16·#17·#18·#22 (재범위·저자 판단·Oropetium).

### 다음 작업

1. 코덱스 적대적 검토 결과 대기
2. 벤치마크 재생성 (GRCh38·Mo17·Arabidopsis)
3. 남은 이슈 처리 (재범위 조정 또는 폐기 결정)

## 1. 지금까지 한 일

### 1차 감사 (REVIEW.md, untracked)
Codex 전체 저장소 감사. 세션 `01a041b5-a458-7fc1-9e86-7435f4eeba42`.
판정: **논문 제출 준비 안 됨.** CRITICAL 2건 + HIGH 6건 + MEDIUM/LOW 다수.
이 감사에서 일부 수정이 워킹 트리에 적용됨 (아래 3절).

### 2차 재검토 (Codex, 세션 `01a043d9-eb67-7892-a315-e6b367d75fec`)
1차 감사의 각 블로커가 실제로 해결됐는지 현재 코드 대비 재검증 + 신규 이슈 발굴.
판정: **여전히 제출 준비 안 됨.**

`codex resume 01a043d9-eb67-7892-a315-e6b367d75fec` 로 이어받을 수 있음.

---

## 2. 확정된 사실 — 저장소 상태

### git 계보 — 2026-08-27 정정

앞서 이 노트에 적었던 두 진술은 **틀렸음.** 확인 결과:

- **`7a80ed7`은 푸시돼 있음.** `wyim-pgl/bwt-algorithm`에 존재하며 `main`(`6a8ae40`)과 원고 브랜치 `e0acd56` 양쪽의 ancestor. 미커밋인 것은 워킹 트리뿐.
- **`framazan/bwtandem`은 무관한 남의 저장소가 아님.** 의도적으로 연동된 협업 저장소. `/data/gpfs/assoc/pgl/devel/bwt-algorithm`에 `framazan` remote가 명시 등록돼 있고, `wyim-pgl/bwt-algorithm`에 `backup/main-pre-framazan-sync`, `merge-framazan-updates` 브랜치가 있음. 다만 `b7a5689`와 `7a80ed7` 사이 merge-base는 없음(별개 계보).

여전히 유효한 문제:
- ❌ **SUPERSEDED (2026-08-31):** 아래 origin 오설정은 **해결됨** — origin을 `wyim-pgl/bwt-algorithm`으로 교체, `framazan`은 devel 클론과 동일하게 별도 remote로 보존, 양쪽 fetch·main 업스트림 재지정 완료. 최신 스냅샷 참조.
- ~~**이 체크아웃의 `origin`이 `framazan/bwtandem`을 가리키는 건 부적절.**~~ 계정 `wyim-pgl`은 그쪽에 READ 권한뿐이고 푸시 대상이 아님. 로컬 `origin/main` ref도 `7a80ed7`로 stale (실제 remote HEAD는 `b7a5689`).
- 실제 프로젝트 트래커는 **`wyim-pgl/bwt-algorithm`** (public, ADMIN). 기존 이슈 #3 GRCh38 / #4 Maize Mo17 / #5 Col-CEN.

### 감사 대상이 5개월 낡음 (중대)

두 감사 모두 **2026-04-01 스냅샷**(`7a80ed7`)을 봤음. 활성 원고 브랜치 `perf/exp1-human-sensitivity`(`e0acd56`, 2026-08-19)는 **80커밋 앞섬** — 최근 커밋 메시지에 `align_accel` 버그 수정, de-novo sensitivity overhaul, 3종 재측정 등이 보임.

> ❌ **SUPERSEDED (2026-08-27):** 아래 "나머지 15건 미확인"은 3차 코덱스 재검토로
> **전건 검증 완료** — 결과는 #10만 CLOSE, #12/#13/#15는 "해결됨" 판정이 반박됨.
> 게이트도 해소됨. 처분표는 최신 스냅샷·wiki update 2026-08-27c 참조.

→ **등록된 17개 이슈 중 상당수가 최신 브랜치에서 이미 해결됐을 수 있음.**
직접 확인된 것: **#14**(LICENSE/CITATION.cff/pyproject/.github 없음)와 **#20**(`except Exception: … continue`)은 `e0acd56`에도 **여전히 적용됨.** 나머지 15건은 **미확인.**
이 항목이 `PLAN.md` 실행 전체를 게이트함.

### 워킹 트리

> ⚠️ **CAUTION (2026-08-27):** 아래는 **stale 체크아웃**(`/data/…/data/bwtandem`)의 목록이다.
> 여기 있던 수정들은 원고 브랜치(`devel/bwt-algorithm`)에 `a163ad9`부터 별도 커밋으로
> **이식·확장되어 푸시됐다.** 이 체크아웃의 워킹 트리 자체는 여전히 미커밋 상태로 남아 있다.

수정됨(미커밋): `CLAUDE.md`, `README.md`, `environment.yml`, `scripts/venn_compare.py`, `src/finder.py`, `src/main.py`, `src/models.py`, `src/c_extensions/__pycache__/build.cpython-311.pyc`
Untracked: `REVIEW.md`, `TRASH.md`, `TRASH/`, `docs/`, `oropetium/`, `results/bwt_chr4_v3.{bed,log}`, `tests/test_pipeline_constraints.py`

**1·2차 감사의 모든 수정과 신규 테스트는 아직 커밋되지 않았음.** `HEAD`에 없음.

---

## 3. 적용된 수정 (미커밋, 2차 재검토에서 검증)

> ❌ **SUPERSEDED (2026-08-27):** "미커밋"은 이 stale 체크아웃 얘기다. 아래 표의 수정은
> 원고 브랜치에 **커밋·푸시 완료** (`a163ad9`, `34a66fd`, `fe3f9cd` 등 — 최신 스냅샷 참조).
> 특히 ColCEN 행의 "PARTIAL — 내부 불일치 유발"은 대상 파일 자체가 원고 계보에 없고,
> 후속 스코어러 `score_colcen.py`가 `--validate` 게이트로 교체·검증됐다.

| 항목 | 위치 | 재검증 결과 |
|---|---|---|
| TRF `.dat` 시작좌표 1-based | `src/models.py:85` | **FIXED** |
| 멀티프로세스 실패 시 부분 출력 차단 | `src/main.py:167,186-196` | **FIXED** (회귀 테스트 없음) |
| TRASH.md 수치 정정 (32,509 calls / ~26분 / ~5×) | `TRASH.md` | **FIXED** (문서만; 벤치마크 자체는 stale) |
| stress test 문서 드리프트 (20kb / 107TP·0FP·1FN) | `README.md:495` | **FIXED** |
| satellite ambiguous-base 제외 + 80% valid 임계 | `src/finder.py:459-504` | 구현됨, 단 결과 재생성 안 됨 |
| Tier 3 CLI period 구간 연결 | `src/finder.py` | 구현됨 |
| ColCEN GFF3 좌표 변환 | `ath_cen/bench_results/analyze_centromere.py:27,177` | **PARTIAL — 오히려 내부 불일치 유발** |

---

## 4. 미해결 블로커 (심각도순)

> ✏️ **PARTIAL (2026-08-27):** 이 목록은 감사 시점 기준. 이후 처리 현황 —
> 1(벤치마크 stale)·2(ColCEN)는 **재생성 캠페인 + `--validate` 게이트로 진행 중**,
> 5(satellite whole-block)는 `c375456`으로 **구조 수정됨**, 6의 environment.yml은
> `1cb7f48`로 **교정됨**, 7(stress test)은 `34a66fd`로 **수정됨**.
> 3(specificity)과 6의 LICENSE/CI(#14)만 실질 잔존. "2차에서 새로 발견" 3건도 전부 수정됨.

1. **CRITICAL — 벤치마크 코퍼스 stale.** `results/bwtandem_Chr4_v3.bed:1`에 수정 전의 all-`N` satellite 오탐(`Chr4:0–1005`)이 그대로 있음. ColCEN/Oropetium/TRASH/telomere/Venn 전부 고정 커밋에서 재생성 필요.
2. **CRITICAL — ColCEN coverage 여전히 신뢰 불가.** half-open 변환 후에도 `analyze_centromere.py`의 41/127/144/148행이 inclusive 산술(`+1`, `<=`)을 씀. TRF는 1-based 미변환(65-80행), mreps도 원본 유지(84-101행) → 도구 간 좌표계 불일치. README:631의 98.2%는 사용 불가.
3. **HIGH — 실제 게놈 specificity 미측정.** 합성 truth와 annotation coverage/overlap만 있고 precision/FP rate 없음.
4. **HIGH — 핵심 근거와 신규 테스트가 버전관리 밖.** fresh clone으로 재현 불가.
5. **HIGH — satellite whole-block 추론 미검증.** 3개 window 샘플로 최대 100kb 라벨링 (`finder.py:462,471-477,517-527`).
6. **HIGH — 릴리스 재현성 미비.** 의존성 unpinned, `environment.yml`이 Python 3.13 요구 + 환경명 `bwt`, 개인 절대경로 잔존, LICENSE/CITATION.cff/pyproject/CI/release 메타데이터 전무.
7. **HIGH — 검증 실패 처리 부실.** `tests/test_random_stress.py:139`가 예외를 삼키고 `continue` → 크래시 케이스가 분모에서 사라진 채 통과 가능.

### 2차에서 새로 발견 (1차 REVIEW.md에 없던 것)
- **MEDIUM** — 유효/무효 tier 이름 혼합 시 조용히 통과 (`finder.py:540`). `{"tier1","typo"}`가 경고 없이 실행됨.
- **MEDIUM** — satellite ambiguous-base 필터가 gap merging 경로에는 미적용 (`finder.py:324`).
- **LOW** — Venn binning이 `end`가 bin 경계에 딱 떨어질 때 bin을 하나 더 셈 (`scripts/venn_compare.py:84`).

### 검증 불가

> ❌ **SUPERSEDED (2026-08-27):** 이후 실측 환경(conda `bwtandem`, Python 3.11.15)에서
> 원고 브랜치 스위트를 반복 실행 — **115 → 179 passed.** CANNOT VERIFY 상태는 해소됨.

2차 Codex 셸에 `pytest`/`conda`가 없어 "34/34 통과"를 독립 재확인하지 못함 → **CANNOT VERIFY**.

---

## 5. 이슈 등록 완료 — `wyim-pgl/bwt-algorithm`

Codex가 초안 17건을 작성했으나 `gh` 접근이 샌드박스에 막혀, Claude 서브에이전트가 등록을 수행함.
FIXED 항목은 제외. 실패 0건, 드롭된 라벨 0건. 기존 #3/#4/#5는 그대로 둠.

| 이슈 | 심각도 | 제목 |
|---|---|---|
| #7 | CRITICAL | satellite 수정 이후 원고 벤치마크 재생성 |
| #8 | CRITICAL | ColCEN coverage 좌표 처리를 half-open으로 일관화 |
| #9 | HIGH | 원고 근거를 버전관리/아카이브로 |
| #10 | HIGH | Cython 없을 때 명확히 실패하거나 완전 동작 |
| #11 | HIGH | 실제 게놈 specificity·FP 검증 추가 |
| #12 | HIGH | satellite whole-block 채우기·병합 검증 및 제한 (finder.py:324 gap merging 포함) |
| #13 | HIGH | 벤치마크 의존성·명령·입력·실행 메타데이터 freeze |
| #14 | HIGH | license·citation·패키지 메타데이터·릴리스·CI |
| #20 | HIGH | stress test 케이스 실패 시 테스트가 실패하도록 |
| #15 | MEDIUM | 1:1 및 경계 인식 벤치마크 매칭 |
| #16 | MEDIUM | 외부 도구 탐색 범위·파라미터 정합 |
| #17 | MEDIUM | 근거 없는 비교·인과 주장 제거/보강 |
| #18 | MEDIUM | Oropetium 검증을 재현 가능한 스크립트로 |
| #19 | MEDIUM | 유효/무효 tier 혼합 선택 거부 |
| #21 | LOW | pytest 10 대비 class-scoped fixture 수정 |
| #22 | LOW | 저장소 위생·벤치마크 산출물 보관 정책 |
| #23 | LOW | Venn bin 경계 half-open 처리 수정 |

기존 이슈: #3 Human GRCh38 / #4 Maize Mo17 satellite / #5 Arabidopsis Col-CEN.

> ⚠️ **CAUTION (2026-08-27):** 아래 `PLAN.md` Phase 표는 **stale 트리 기준**이라
> 상당 부분 obsolete — Phase 0(저장소 위상)·Phase 2 #8/#10 대상·Phase 3-4의
> Chr4/Oropetium 재실행은 3차 재검토에서 폐기/재범위됐다. Oropetium은 캠페인 제외
> 확정(#18). 현행 실행 계획은 최신 스냅샷의 "다음 작업"이다.

**완료:** Codex가 `PLAN.md` (465줄) 작성. 6단계 Phase 구성.

| Phase | 내용 | 담당 이슈 | Exit criterion |
|---|---|---|---|
| 0 | 권위 있는 저장소·근거 경계 확립 | #9, #22 | 올바른 대상 저장소, 감사 수정 커밋, 산출물 정책 |
| 1 | 과학적 주장·평가 계약 확정 | #11, #16, #17 | 검증·비교 프로토콜 승인 |
| 2 | 검출기·좌표·테스트 동작 수정 | #8, #10, #12, #19, #20, #21, #23 | 수정 병합 + 경계/회귀 테스트 통과 |
| 3 | 평가 harness 구축 + RC freeze | #13, #15, #18 | 환경 락, 독립 검증된 테스트, 불변 RC |
| 4 | 고정 벤치마크 실행 | #7 (+ #3, #4, #5) | 모든 원고 결과를 RC에서 provenance와 함께 재생성 |
| 5 | 패키징·아카이브·릴리스 승인 | #14 | 태그·아카이브·CI 검증·fresh clone 재현 |

**핵심 의존 순서:** Phase 3 종료 전에는 전체 게놈 벤치마크를 돌리지 말 것 — 그 이전 실행은 진단용이며 원고에 쓸 수 없음. #8은 ColCEN 재실행 **전에** 수정·경계 테스트 완료 필수.

**저자 판단 필요 (REQUIRES AUTHOR DECISION) 6건:**
1. git remote 토폴로지 (origin 교체 / upstream 보존 / fresh clone)
2. 실제 게놈 specificity truth 소스 (수동 큐레이션 / held-out annotation / 시뮬레이션 / 층화 샘플)
3. 비교 도구 파라미터 정책 (각 도구 기본값 / 정합 범위 / 양쪽 모두)
4. Cython을 필수로 둘지, 완전한 fallback을 제공할지
5. satellite scanner를 core/default로 둘지 optional mode로 둘지
6. 프로젝트 라이선스 선택

---

## 6. 다음 액션

> ❌ **SUPERSEDED (2026-08-27):** 아래 1–6은 전부 실행됐거나 재범위됐다
> (2·3·4·6 완료, 5는 캠페인 실행 중, 1은 PLAN.md obsolete 판정).
> 현행 액션은 최신 스냅샷의 "다음 작업".

1. `PLAN.md` 확인 후 단계별 실행 착수.
2. `origin` remote 오설정 정리 — 제3자 저장소를 가리키는 상태 해소, 로컬 HEAD `7a80ed7` 푸시 경로 확보.
3. 워킹 트리의 수정 + `tests/test_pipeline_constraints.py`를 release candidate로 **커밋** (#9).
4. ColCEN 좌표 산술 전면 수정 + 경계 케이스 테스트 (#8) — 재실행 **전에** 완료해야 함.
5. 고정 커밋에서 전체 벤치마크 재실행 (#7, #13) — 명령/커밋해시/환경락/wall time/peak mem/CPU/체크섬 기록.
6. "34/34 통과"를 문서화된 환경에서 독립 재현.
