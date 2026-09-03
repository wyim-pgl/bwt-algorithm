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

`quarantine.md` §6 의 확정 결함 중 값·문구 교체로 끝나는 것. **§6.17 처리가 먼저다**(C-1 참조) —
그것이 A-2 의 서술 방향을 결정한다.

### A-1. 원고만 수정 (재해시 불필요, CI 위험 없음)

- [ ] **§6.7** Table 3C-b 캡션: "smallest banded loss" → TRF 0.87 이 최소. 차이값 14.36 → **14.88**
- [ ] **§6.8** "at most 1.5 / 2.8 points" → 실제 최대 **1.54 / 2.81**
- [ ] **§6.9** Table 1c: 캡션은 네 도구, 표는 세 행 — TRASH 행 추가하거나 캡션에서 제외 사유 명시
- [ ] **§6.10** 본문 107행에 TRF 의 precision 두 값 추가 (97.18% vs 64.60%, 31.19 vs 29.25)
- [ ] **§6.11** Table 2 캡션에 BWTandem 의 CEN180 bp precision 최저(57.08%) 사실 한 문장
- [ ] **§6.13** S4 캡션 "lowest of the five tools" → **second-lowest** (tantan 3.42% < 8.68%)
- [ ] **§6.14** "All competitor runs ... inside one Singularity container" → 예외 2건 명시
- [ ] **§6.16** stride 축소 조건: 432행과 440행 중 틀린 쪽 정정 (약 24 Mb 미만에서만 활성)
- [ ] **§6.21** "at least 70%" → **80%**, catch-all 에도 적용됨을 공개
- [ ] **§6.5** "roughly 5% more sequence" → 두 값을 같은 단위(염기)로 재계산 후 교체
- [ ] **§6.6** stride 공식에 프리셋 배수 포함, 포화점 재유도
- [ ] **§6.4** 2026 도구 메모리 셀 2개: AniAnn's Col-CEN 0.50 → **0.48**, longdust 0.07 → **0.06**
- [ ] **§6.12** TR-1 경계 오차 네 값(771 bp / 4.3 kb / 4,265 / 7,973)의 통계 정의 후 일치
- [ ] **§6.24** Discussion 4.3 — Tier 1 도 FM-index 열거를 돌린다(`TIER1_FMSCAN=1`). 메커니즘 서술 정정
- [ ] **§6.25** Abstract 의 "per-figure provenance manifest" 약속 삭제 (매핑 0건, 그림 프로그램 6/6 스텁)
- [ ] **§6.26** S1.3 파라미터 범위를 "핵 염색체 5개"로 한정하거나 ChrC/ChrM 값 추가
- [ ] **Codex R3-10** Discussion 4.2 에 세 게놈 × 두 패스 전체 행렬 명시
      (gap-fill 은 세 게놈 전부, catch-all 은 human·maize 만)
- [ ] **Codex R3-12** Limitations 에 native-path 미해결 flake 를 명시
- [ ] **Codex R2-9** 밴드 손실 4건이 예치본과 0.01 pp 어긋남 — 반올림 규약 통일
- [ ] **Codex R2-8** S2 의 Python/numpy 버전(3.11.14 / 2.3.1)이 예치 provenance(3.11.15 / 2.4.6)와 불일치,
      매니페스트 스코어러 해시 2건도 실제와 다름
- [ ] 위 수정 후 원고 전체에서 **GiB/GB 표기 일괄 점검** (§6.4·Codex 발견 12 — 숫자는 그대로, 라벨만)

### A-2. `results/` 수정 — ⚠️ 재해시 필수

**순서 엄수** (`quarantine.md` §8.2): 편집 → `git add` → 재해시(백그라운드 600초) → 테스트 → 커밋.
재해시 중 `results/` 를 건드리지 않는다. 아래 넷을 **한 번에** 처리해 재해시를 1회로 끝낸다.

- [ ] **§6.1** Table 2 TRASH 비용 셀 — template 행은 3실행 합집합임을 반영(합계 37:28:39 / 2.63 GiB),
      de novo 행은 자기 로그값(5:47:30 / 1.29 GiB)으로. 캡션의 "재검증 불가" 면제 문구 삭제
- [ ] **§6.2** `results/manifest.tsv:67` — 정정 전 문장 4개를 원고 현행 문구로 교체
- [ ] **§6.3** AniAnn's 스레드 `1` → `2` (매니페스트 행 123·126·129·132)
- [ ] **Codex 발견 11** tantan 스레드/명령: 실제 활성 스레드 1, 확장 명령(`-w 500/2000/200`) 명시.
      `allocated_cpus` 와 도구 스레드를 열로 분리
- [ ] **Codex 발견 14** 매니페스트의 `PENDING` 행 5개(8, 16–19)를 superseded 로 표시하거나 분리
- [ ] 재해시 후 `pytest tests/test_deposit_hashes.py tests/test_env_var_docs.py`

---

## B. 남은 검증 — 목록을 닫기 전에

- [x] **Codex 라운드 2** — Results ↔ 예치 산출물, 9건 ✅ (2026-09-03, `msreview/codex_out_R2.txt`)
- [x] **Codex 라운드 3** — Discussion 메커니즘 ↔ 코드, 14건 ✅ (2026-09-03, `msreview/codex_out_R3.txt`)
- [x] R2·R3 의 최중대 6건 실물 검증 후 §6 에 추가 ✅ (2026-09-03, §6.22–6.26 · §3.10)
- [ ] R2·R3 의 나머지 17건 검증 또는 명시적 기각
- [ ] Codex 첫 리뷰의 **미검증 14건** 처리 — 검증하거나 명시적으로 기각
- [ ] Kimi 라운드 1–3 의 MEDIUM/LOW 중 미검증분 선별 검증
- [ ] `TIER2_APPROX_SEED` at identity ≥ 0.88 — centromere/satellite 에 도움될 수 있으나 **미시험**.
      시험 없이 켜지 말 것 (이 논문 범위 밖일 수 있음)

---

## C. 저자 결정이 필요한 것 — 내가 정할 수 없다

- [?] **C-1. §6.17 — Tier 3 탐색 창 주장을 어떻게 할 것인가.**
      원고는 "the command-line maximum period filters its report rather than bounding its search"
      라 하지만 코드는 탐색을 좁히고 주석이 스스로 "NOT output-neutral" 이라 적는다.
      range-cost 실험 3쌍은 옛 커밋 `07ad6fa` 로 돌았다.
      **선택지**: (a) 두 코드 시대를 분리 서술하고 현재 구현의 범위 비민감성 주장을 철회,
      (b) 100 vs 2,000 실험을 `0363d8b`/현재 코드로 재실행. **논문 중심 주장이 걸려 있다.**
- [?] **C-2. §6.18 — matched-range 비대칭 처리.**
      BWTandem 네이티브 재실행(+1.01 pp) vs TRF 사후 필터.
      **선택지**: (a) 공통 비교를 사후필터 행으로 통일, (b) 모든 도구를 네이티브 100 으로 재실행
- [?] **C-3. §3.9 — 튜닝 캠페인 공개 범위.**
      원장 44건·ULTRA 표적을 어디까지 밝힐지, 22번 염색체 결과를 "선택 후 검증"으로 재서술할지
- [?] **C-4. §3.8 — catch-all 0.72 선정의 선택-무관 증거.**
      maize/Arabidopsis 진실셋 또는 보류 human 염색체에서 재확인할지
- [?] **C-9. §6.22 — AniAnn's 밴드 셀을 사전등록대로 넣을 것인가, period 기반 셀을 전부 철회할 것인가.**
      사전등록 프로토콜은 채점 가능하다고 적었고, 예치 로그의 AniAnn's 밴드 recall(99.15%)이
      BWTandem·TRF를 모두 앞선다. **저장소에 있는 불리한 수치가 원고에 없는 상태다.**
- [?] **C-10. §6.23 — TRF period 열 복구 후 S4 재실행 여부.**
      스코어러가 TRF 의 period 열을 버린다. 제대로 읽으면 TRF 의 period 일치도가
      BWTandem·ULTRA 를 앞선다는 재계산이 있다(미검증). 재실행하면 S4 순위가 바뀐다.
- [?] **C-11. §3.10 — Arabidopsis catch-all off 선택을 in-sample 로 표시할 것인가, 재설계할 것인가.**
      선택 근거가 평가 진실셋에 대한 bp precision 비교(65.54 vs 60.72)였다.
- [?] **C-5. §6.15 — TRASH 를 de novo 로 볼 것인가 template-guided 로 볼 것인가.**
      313행과 326행이 반대로 적고 있고, 분류에 따라 "largest de novo footprint" 주장이 바뀐다
- [?] **C-6. 투고처** — GigaScience full paper(재현성 심사) vs Bioinformatics Original.
      분석: `docs/2026-09-01-venue-and-positioning-analysis.md`
- [?] **C-7. 제목** — 후보 5개. #27 은 현행 제목이 부록이 반증하는 메커니즘을 함의한다고 지적
- [?] **C-8. Abstract 재구성** — #26·#28, #27/#28 의 증거 반영

---

## D. 재실행·재측정이 필요한 것

- [ ] **§6.20** `sacct -j ... --format=JobID,Elapsed,MaxRSS,State -P` 원문을 체크섬과 함께 예치.
      불가하면 예치된 `wait4` 값(17.37 / 19.97 / 1.31 GiB)을 그 이름으로 쓰도록 원고 수정
- [ ] **§6.19** 누락 스코어러 3개 커밋 — `rescore_tables_3bc.py`, `score_exp3.py`, `score_overlap.py`.
      사설 클러스터 경로를 저장소 상대경로로 교체 (**Availability 약속이 현재 거짓**)
- [ ] **Codex 발견 9** CEN180 identity 층화를 재생성 BED 로 재채점, 다섯째 층(99–100%, 100 monomers) 포함
- [ ] **Codex 발견 8** 2026 도구를 Table 3B/3C 정본 스코어러로 재채점하고 출력 예치
- [ ] **#30** 경쟁 도구 GNU-time 로그 예치 + 해시 (현재 저장소 밖)
- [ ] **Codex 발견 10** CEN180 진실셋 생성 과정(BLAST 버전·명령·원시 68,840 hit) 복구 또는
      "재현 불가한 legacy 편의 진실셋"으로 명시

---

## E. 열린 이슈 (GitHub)

- [ ] **#14** [HIGH] license, citation, package metadata, release versioning, CI
- [ ] **#26** [MEDIUM] 2026-tool 벤치마크를 남은 서술에 반영
- [ ] **#27** [HIGH] 제목·Abstract 의 FM-index 인과 프레이밍 → C-7
- [ ] **#28** [HIGH] Abstract: caveat 이 우리 수치에만 붙고 경쟁 도구엔 안 붙음 → C-8
- [ ] **#29** [MEDIUM] range-cost "sublinearity" 주장 → C-1 과 연동
- [ ] **#30** [HIGH] 경쟁 도구 GNU-time 로그 예치 → D
- [ ] **#31** [HIGH] deposit-hash 워크플로 (재해시 순서) → A-2 에서 실증됨
- [ ] **#32** [MEDIUM] 증거 트리의 동등성 주장에 검사 첨부
- [~] **#33** [MEDIUM] 2차 의견 리뷰 — Kimi 완료, Codex 라운드 1 완료, R2·R3 진행 중

---

## F. 기타

- [ ] `quarantine.md` §6.7–6.21, §3.8–3.9 추가분 커밋 (현재 워킹트리에만 있음)
- [ ] origin 푸시 (로컬이 1 커밋 앞섬)
- [ ] 위키의 "미커밋" ⚠️ 정리(`0357cbd` 로 해소) **및 `ruleout.md` → `quarantine.md` 개명 반영**
- [!] **flaky `TestAdjacentGroundTruth::test_sensitivity`** ⛔ 클러스터 ASLR 이 꺼져 있어
      실패 레이아웃에 도달 불가. 관리자에게 노드 1대 ASLR 재활성화 요청 필요 (`quarantine.md` §9.1)
- [ ] Filip 그림 6개 — 우리 손 밖. 브리프 `results/figures/paper_figs/HANDOFF_FILIP.md`

---

## 완료 기록

- [x] `quarantine.md` 신설 (§1–§10, 폐기 대장) ✅ (2026-09-03, `0357cbd`)
- [x] resume.md 격리 — 📌 정본 포인터 + ❌ MOVED 배너 4곳, 앵커 0 broken ✅ (2026-09-03, `0357cbd`)
- [x] Codex 2차 의견 리뷰 전사본 예치 ✅ (2026-09-03, `0357cbd` / `docs/2026-09-03-codex-review-findings.md`)
- [x] Kimi 원고 교차검증 라운드 1–3 (11청크, 타임아웃 0) ✅ (2026-09-03, `msreview/out_R*.md`)
- [x] Codex 원고 라운드 1 (코드 ↔ 원고 대조) ✅ (2026-09-03, `msreview/codex_out_R1.txt`)
- [x] 확정 결함 21건을 `quarantine.md` §6 에 기록 ✅ (2026-09-03, 전건 실물 파일 재현)
- [x] 랩 위키 반영 ✅ (2026-09-03, wiki `eecf950`)
