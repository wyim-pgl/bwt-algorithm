# RETIRED_DO_NOT_USE — 격리된 낡은 벤치마크 산출물 (2026-08-31)

이 디렉토리의 모든 파일은 **원고 근거로 인용 금지**. 위키 `guide/handoff-hygiene.md`의
경로 격리 규약에 따라 이동됨. 삭제가 아니라 격리이므로 내용은 그대로 보존돼 있다.

## 왜 격리됐나

이 체크아웃(HEAD `7a80ed7`, 2026-04-01)은 감사 대상 stale 트리다. 아래 산출물은
2026-08-27~29 재생성 캠페인(핀 `0363d8b`, `/data/gpfs/assoc/pgl/devel/exp1_human/regen/`)으로
전부 superseded됐다:

| 항목 | 격리 사유 |
|---|---|
| `results/` | 2026-03-30~31 Chr4/합성 벤치마크 산출물. `bwtandem_Chr4_v3.bed`는 수정 전 all-N satellite 오탐(`Chr4:0–1005`) 포함 — 이슈 #7 CRITICAL에 명시된 stale corpus |
| `ath_cen/` | 구 Col-CEN 벤치마크. `analyze_centromere.py`의 inclusive/half-open 좌표 혼용(#8) — README의 98.2% 수치는 사용 불가 판정 |
| `oropetium/` | #18 캠페인 제외 (저자 결정, 2026-08-27). 게놈 FASTA 포함 — 재사용 시 여기서 꺼내되 BED 산출물은 인용 금지 |
| `TRASH/` | 구 TRASH 도구 설치본 + colcen 실행 산출물 (2026-04-01). 캠페인 기준선은 filip 디렉토리의 고정 BED |
| `figures/` | 구 산출물 기반 그림 (2026-03-31) |
| `slurm-*.out` | 2026-03-31 실행 로그 |

## 현행 정본

- 재생성 BED + provenance: `/data/gpfs/assoc/pgl/devel/exp1_human/regen/` (colcen/human/maize, 전부 gate PASS)
- 채점표: `score_table1_regen_v2.txt` (v2가 정본, v1은 기록용 보존)
- 원고·수치: `/data/gpfs/assoc/pgl/devel/bwt-algorithm` `perf/exp1-human-sensitivity` 브랜치의 `manuscript.md` + `results/manifest.tsv`
- 상태 추적: 이 체크아웃의 `resume.md` 최신 스냅샷 + 위키 `projects/bwt-tandem-repeat-finder/`

주의: `results/` 등 일부는 git 추적 파일이라 `git status`에 D(deleted)로 나타난다.
이 트리는 원고 트리가 아니므로 복원(`git restore`)하지 말 것 — 격리가 의도다.
