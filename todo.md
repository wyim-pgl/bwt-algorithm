# todo.md — BWTandem 남은 작업

> 📌 **정본**: 이 파일이 "앞으로 할 일"의 유일한 정본이다 — **미해결은 전부 여기 있다.**
> 이미 틀렸거나 시도해서 안 됐거나 배제된 것은 [`quarantine.md`](quarantine.md),
> 살아 있는 상태는 `resume.md`, 운영 규칙은 `CLAUDE.md`.
> 최종 갱신 **2026-09-03**.

## 표기 규약

| 마크 | 뜻 |
|---|---|
| `- [ ]` | 미착수 |
| `- [~]` | 진행 중 |
| `- [x]` | **완료** — 뒤에 `✅ (YYYY-MM-DD, 근거)` 를 붙인다 |
| `- [!]` | 차단됨 — 뒤에 `⛔ 무엇에 막혔는지` 를 붙인다 |
| `- [?]` | 저자 결정 대기 — 내가 임의로 정할 수 없는 것 |

**완료 표시 규칙**: 체크만 바꾸지 말고 **근거**(커밋 해시, 잡 ID, 파일 경로)를 같은 줄에 남긴다.
근거 없는 완료는 다음 세션이 다시 조사한다.

```bash
grep -n '^- \[ \]' todo.md    # 미착수
grep -n '^- \[!\]' todo.md    # 차단
grep -n '^- \[?\]' todo.md    # 저자 결정 대기
```

---

## A. 원고 결함 수정 — 내가 고칠 수 있는 것

`quarantine.md` §6 의 확정 결함 중 값·문구 교체로 끝나는 것. §6.17/C-1은
`0363d8b` 재측정값으로 2026-09-04에 적용했다.

### A-1. 원고만 수정 (재해시 불필요, CI 위험 없음)

- [x] **§6.7** Table 3C-b 캡션: "smallest banded loss" → TRF 0.87 이 최소. 차이값 14.36 → **14.88**
- [x] **§6.8** "at most 1.5 / 2.8 points" → 실제 최대 **1.54 / 2.81**
- [x] **§6.9** Table 1c: 캡션은 네 도구, 표는 세 행 — TRASH 행 추가하거나 캡션에서 제외 사유 명시
- [x] **§6.10** 본문 107행에 TRF 의 precision 두 값 추가 (97.18% vs 64.60%, 31.19 vs 29.25)
- [x] **§6.11** Table 2 캡션에 BWTandem 의 CEN180 bp precision 최저(57.08%) 사실 한 문장
- [x] **§6.13** S4 캡션 "lowest of the five tools" → **second-lowest** (tantan 3.42% < 8.68%)
- [x] **§6.14** "All competitor runs ... inside one Singularity container" → 예외 2건 명시
- [x] **§6.21** "at least 70%" → **80%**, catch-all 에도 적용됨을 공개
- [x] ~~§6.6 stride 공식~~ · ~~§6.16 stride 축소 모순~~ — **오탐으로 철회** ✅ (2026-09-03, `quarantine.md` §6.6·§6.16)
- [x] **§6.5** "roughly 5% more sequence" → 두 값을 같은 단위(염기)로 재계산 후 교체
- [x] **§6.4** 2026 도구 메모리 셀 2개: AniAnn's Col-CEN 0.50 → **0.48**, longdust 0.07 → **0.06**
- [ ] **§6.12** TR-1 경계 오차 네 값(771 bp / 4.3 kb / 4,265 / 7,973)의 통계 정의 후 일치
- [x] **§6.24** Discussion 4.3 — Tier 1 도 FM-index 열거를 돌린다(`TIER1_FMSCAN=1`). 메커니즘 서술 정정
- [x] **§6.25** Abstract 의 "per-figure provenance manifest" 약속 삭제 (매핑 0건, 그림 프로그램 6/6 스텁)
- [x] **§6.26** S1.3 파라미터 범위를 "핵 염색체 5개"로 한정하거나 ChrC/ChrM 값 추가
- [x] **Codex R3-10** Discussion 4.2 에 게놈×패스 행렬 명시 — 어블레이션이 있는 것과 없는 것을 구분
- [x] **Codex R3-12** Limitations 에 native-path flake 명시 ✅ (2026-09-03) §4.4 에 네 번째 한계로 추가 — 소진된 진단과 재현 불가 사유까지
- [x] **Codex R2-9** 밴드 손실 반올림 ✅ (2026-09-03) 예치본 기준 전수 정정 — 3C-b 표 −0.92→−0.91, 캡션 0.92→0.91, 3B-b 캡션 0.50→0.51, ULTRA 상계 2.81→**2.82**(내 앞선 수정도 뺄셈이라 틀렸다), §4.4 의 같은 주장도
- [x] **Codex R2-8(정정)** 전체 게놈 실행은 3.11.14/2.3.1 이 맞다 — 다른 것은 p100 패널·식별도 스윕이며,
      S2 에 그 별도 환경을 공개했다. 매니페스트 스코어러 해시 2건은 `results/` 라 A-2 로 이월
- [x] 원고 전체 **GiB/GB 라벨 37곳 교체** ✅ (2026-09-03, 숫자 불변)

### A-2. `results/` 수정 — ⚠️ 재해시 필수

**순서 엄수** (`quarantine.md` §8.2 — 2026-09-03 정정): 모든 편집 완료 → **재해시(포그라운드, 또는 wait 로 완료 확인)** →
그 다음 체크섬 2개 포함 `git add` → `git diff --cached` 확인 + 미스테이징 results 변경 없음 → 테스트 → 커밋.
**옛 순서(`git add` 를 먼저)는 위험하다** — 스테이징엔 옛 해시, 워킹트리엔 새 해시인 커밋이 나온다.
재해시 중 `results/` 를 건드리지 않는다. 아래 넷을 **한 번에** 처리해 재해시를 1회로 끝낸다.

- [ ] **§6.1** Table 2 TRASH 비용 셀 — template 행은 3실행 합집합임을 반영(합계 37:28:39 / 2.63 GiB),
      de novo 행은 자기 로그값(5:47:30 / 1.29 GiB)으로. 캡션의 "재검증 불가" 면제 문구 삭제
- [ ] **§6.2** `results/manifest.tsv:67` — 정정 전 문장 4개를 원고 현행 문구로 교체
- [ ] **§6.3** AniAnn's 스레드 `1` → `2` (매니페스트 행 123·126·129·132)
- [ ] **Codex 발견 11** tantan 스레드/명령: 실제 활성 스레드 1, 확장 명령(`-w 500/2000/200`) 명시.
      `allocated_cpus` 와 도구 스레드를 열로 분리
- [ ] **Codex 발견 14** 매니페스트의 `PENDING` 행 5개(8, 16–19)를 superseded 로 표시하거나 분리
- [x] ~~재해시 후 pytest~~ — **항목 아님**: A-2 절차의 마지막 단계이지 독립 항목이 아니다
- [x] **Codex 라운드 2** — Results ↔ 예치 산출물, 9건 ✅ (2026-09-03, `msreview/codex_out_R2.txt`)
- [x] **Codex 라운드 3** — Discussion 메커니즘 ↔ 코드, 14건 ✅ (2026-09-03, `msreview/codex_out_R3.txt`)
- [x] R2·R3 의 최중대 6건 실물 검증 후 §6 에 추가 ✅ (2026-09-03, §6.22–6.26 · §3.10)
- [ ] R2·R3 의 나머지 17건 검증 또는 명시적 기각
- [ ] Codex 첫 리뷰의 **미검증 14건** 처리 — 검증하거나 명시적으로 기각
- [ ] Kimi 라운드 1–3 의 MEDIUM/LOW 중 미검증분 선별 검증
- [x] ~~`TIER2_APPROX_SEED` at identity ≥ 0.88~~ — **제외** (2026-09-03): 미시험 연구 방향이지
      논문 작업이 아니다. `quarantine.md` §9.3 으로 이관
- [x] **C-1. Tier 3 탐색 창 — (b) `0363d8b` 재실행 적용** ✅ (2026-09-04)
      SLURM `6147698`·`6147699`·`6147700` 세 쌍의 비율은 **1.82 / 1.77 / 1.77
      (평균 1.79, 범위 1.77–1.82)**. p100은 6.6 h→4.0 h, p2000은 8.6 h→7.3 h로
      줄었다. `07ad6fa`에서 p100 arm이 장주기 탐색을 수행하고 결과를 버려 기존
      1.30–1.41(평균 1.34)이 너무 유리하게 낮았다. 원고·manifest·Figure 2 소스·격리 대장에 반영.
- [x] **C-2. matched-range — (a) 사후필터로 통일** ✅ (2026-09-03) Table 1b 행·캡션·본문 교체, 네이티브 재실행(79.88/50.62)은 민감도 분석으로 강등
- [~] **C-3. 튜닝 캠페인 — (a) 최소 공개 채택** ✅ 본문 반영 (2026-09-03, `d580840`) · ⏳ 원장 예치는 A-2 재해시와 함께
- [x] **C-4. catch-all 0.72 — (a) 최소 공개 채택** ✅ (2026-09-03, `d580840`) S3 에 in-sample 명시
- [x] **C-9. AniAnn's 밴드 셀 — (a) 사전등록 준수** ✅ (2026-09-03) 예치돼 있던 셀을 Table 1b·1c·§3.2 에 반영, 캡션의 거짓 배제 사유 정정
- [x] **C-10. TRF period 열 복구 + S4 재채점 — 저자 지시대로 수행** ✅ (2026-09-03) 스코어러 2곳 수정, 테스트 교체(39 passed), TRF period-exact **63.53%**, 원고 S4 반영. ⏳ 예치 JSON 교체는 A-2 재해시와 함께
- [x] **C-11. Arabidopsis catch-all off — (a) 최소 공개 채택** ✅ (2026-09-03, `d580840`) S2 에 65.54 vs 60.72 와 보류셋 부재 명시
- [x] **C-5. TRASH = de novo** ✅ (2026-09-03, 저자 결정: TRASH 에 de novo 모드가 있다) — 326행을 human 한정으로 좁히고 TRASH 를 de novo 로 재분류. §6.15
- [x] **C-6. 투고처 — Bioinformatics Application Note** ✅ (2026-09-03, 저자 결정)
- [x] **C-7. 제목 — 현행 유지** ✅ (2026-09-03) `BWTandem: FM-index seeding for wide-period-range tandem repeat detection in assembled genomes` — 저자가 지정한 문구가 원고 현행과 동일. **#27(FM-index 인과 프레이밍 지적)은 저자 유지 결정으로 종결**
- [ ] **C-8. Abstract 재구성 — 전체 이슈 정리 후로 연기** (저자 결정 2026-09-03)
---

## D. 재실행·재측정이 필요한 것

- [x] **§6.20** `sacct -j ... --format=JobID,Elapsed,MaxRSS,State -P` 원문을 체크섬과 함께 예치. ✅ (2026-09-03, `e4ae632`) 34개 잡 sacct 예치, 헤드라인 3개 정확히 일치
- [x] **§6.27** Col-CEN 밴드 recall 을 재생성 BED 로 재채점 — 원고의 91.31/94.68 은 예치 출처가 없다 ✅ (2026-09-03, `e4ae632`) **오탐이었다** — 재생성본이 91.31/94.68 을 그대로 준다. 산출물만 없었고 이제 예치됨
- [x] **TRF S4 JSON 교체** — `results/one_to_one/one_to_one_trf_annot_r50.json` 을 재채점본으로. ✅ (2026-09-03, `e4ae632`)
- [x] **원장 예치** — `exp1_human/loop/{ledger.tsv,best.json}` 를 `results/` 에 예치 (C-3 (a) 의 나머지 절반). ✅ (2026-09-03, `e4ae632`) `results/tuning_ledger/` + README
- [x] **§6.19** 누락 스코어러 3개 커밋 — `rescore_tables_3bc.py`, `score_exp3.py`, `score_overlap.py`. ✅ (2026-09-03, `e4ae632`) 3개 커밋 + 클러스터 경로를 env 오버라이드로
- [x] **Codex 발견 9** CEN180 identity 층화를 재생성 BED 로 재채점, 다섯째 층(99–100%, 100 monomers) 포함 ✅ (2026-09-03, `e4ae632`) 재생성본으로 재채점, 다섯째 층 추가
- [x] **Codex 발견 8** §2.2.4 채점 경로 정확화 ✅ (2026-09-03) "no second scoring implementation exists" 는 거짓 — human 은 EXTRA·argv 도 덮고, maize 는 `score_maize_postmerge.py`(좌표 전용 진단)를 쓰지 정본 `rescore_tables_3bc.py` 가 아니다
- [x] **#30** 경쟁 도구 GNU-time 로그 예치 + 해시 (현재 저장소 밖) ✅ (2026-09-03, `e4ae632`) 40개 예치, 5개 도구 전부 재현 확인
- [x] **Codex 발견 10** CEN180 진실셋 생성 과정(BLAST 버전·명령·원시 68,840 hit) 복구 또는 ✅ (2026-09-03, `e4ae632`) **원시 hit 68,840개가 살아 있었다** — 예치하고 원고 정정
- [ ] **#14** [HIGH] license, citation, package metadata, release versioning, CI
- [x] **#26** 2026-tool 벤치마크 반영 ✅ (2026-09-03) C-9(a) 로 Table 1b·1c·§3.2 에 진입, 본문 언급 32곳 · **GitHub 닫힘**
- [x] **#27** 제목의 FM-index 인과 프레이밍 — **저자 유지 결정으로 종결** ✅ (2026-09-03, C-7)
- [x] ~~#28~~ — **C-8 에 흡수**: Abstract 재구성 안에서 함께 처리된다
- [x] ~~#29~~ — **C-1 에 흡수**: 재실행 결과가 이 주장의 운명을 정한다
- [x] **#30** 경쟁 도구 GNU-time 로그 예치 ✅ (2026-09-03, `e4ae632`) 40개 + 5도구 재현 확인 · **GitHub 닫힘**
- [x] **#31** deposit-hash 워크플로 ✅ 완화됨 — `tests/test_deposit_hashes.py` 가 오늘 두 번 잡아냈고,
      순서는 `quarantine.md` §8.2 에 성문화. 스크립트 자체 재작성은 하지 않는다 · **GitHub 닫힘**
- [x] **#32** 동등성 주장에 검사 첨부 ✅ `24cd12a`·`ccf04dd`·`e8ac50c` 가 두 건 모두 처리.
      재발 방지는 `quarantine.md` 의 "대체물 없으면 없음이라 적는다" 규약이 맡는다 · **GitHub 닫힘**
- [~] **#33** 2차 의견 리뷰 — Kimi 3라운드·Codex 4라운드 **전부 완료**. 남은 것은 미검증분 처리(B절)

---

## F. 기타

> ⚠️ **SLURM 제출 함정 (2026-09-03)**: 이 환경은 `SBATCH_PARTITION`·`SBATCH_ACCOUNT`·**`SBATCH_TIMELIMIT`**
> 를 export 해 두고 있어 `--time` 이 무시된다(24h 요청이 13-23:00:00 이 됐다). `resume` 의 회피법
> `SBATCH_X= sbatch` 는 `SBATCH_TIMELIMIT` 에서는 *빈 시간 스펙* 오류를 낸다. **`env -u` 로 지울 것**:
> `env -u SBATCH_PARTITION -u SBATCH_ACCOUNT -u SBATCH_TIMELIMIT sbatch --partition=cpu-s2-core-0 --account=cpu-s2-pgl-0 ...`
> s1 은 자기 6일짜리 annotation 잡으로 늘 막혀 있다 — 이걸로 시작 예정이 9/5 에서 즉시로 바뀌었다.

- [x] **origin 푸시** ✅ (2026-09-03) `67f11b1..f57fcb0`, 21 커밋
- [x] 위키의 "미커밋" ⚠️ 정리 + 개명 반영 ✅ (2026-09-03, wiki `50e863d`)
- [!] **flaky `TestAdjacentGroundTruth::test_sensitivity`** ⛔ 클러스터 ASLR 이 꺼져 있어
      실패 레이아웃에 도달 불가. 관리자에게 노드 1대 ASLR 재활성화 요청 필요 (`quarantine.md` §9.1)
- [~] **Filip 그림 6개 — 전달 완료 (2026-09-03), 통합됨**
      6개 전부 구현·렌더됨(`results/figures/paper_figs/rendered/`, PNG+PDF). 스텁 0개.
      Fig 1 은 C-2 결정(사후필터 arm)을 독립적으로 같은 방향으로 그렸다.
      **Fig 2 소스와 CSV는 C-1 재측정으로 교체 완료** — 패널 A 범위 1.77–1.82,
      평균 1.79. 렌더링은 이 작업 후 별도로 남았다. 상세는 그 디렉터리의 README.
- [x] `quarantine.md` 신설 (§1–§10, 폐기 대장) ✅ (2026-09-03, `0357cbd`)
- [x] resume.md 격리 — 📌 정본 포인터 + ❌ MOVED 배너 4곳, 앵커 0 broken ✅ (2026-09-03, `0357cbd`)
- [x] Codex 2차 의견 리뷰 전사본 예치 ✅ (2026-09-03, `0357cbd` / `docs/2026-09-03-codex-review-findings.md`)
- [x] Kimi 원고 교차검증 라운드 1–3 (11청크, 타임아웃 0) ✅ (2026-09-03, `msreview/out_R*.md`)
- [x] Codex 원고 라운드 1 (코드 ↔ 원고 대조) ✅ (2026-09-03, `msreview/codex_out_R1.txt`)
- [x] 확정 결함 21건을 `quarantine.md` §6 에 기록 ✅ (2026-09-03, 전건 실물 파일 재현)
- [x] 랩 위키 반영 ✅ (2026-09-03, wiki `eecf950`)
