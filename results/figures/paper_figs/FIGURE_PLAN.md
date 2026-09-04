> **2026-09-02 부록**: longdust·AniAnn's 벤치마크 편입(Methods 2.2.4)에 따른
> 패널 델타는 `HANDOFF_FILIP.md`에 반영돼 있다 — Fig 1A/1B에 `competitor2026`
> 점 3개(AniAnn's 0.10%는 축 위에 유지+각주), Fig 2C 비용 점 2개(FASTA 스코프
> 각주), Fig 4A/4B에 2026 행(4B의 longdust 0.00 세 개는 생략 금지), Fig 1C/4C는
> 부재 사유 캡션. 금지 도판 D의 근거였던 "모든 평가 도구가 모든 어레이 검출"은
> 이제 원 비교군 한정 서술로 바뀌었지만, 그 도판은 여전히 금지다(체리피킹이라는
> 본질은 불변).

## 1. Original Paper 트랙 우선순위 도판 목록

권장 구성은 본문 Figure 1–4와 보충 Figure S1–S2입니다. 핵심 메시지는 Figure 1과 2만으로도 완결되며, Figure 3–4는 한계와 다중 생물종 검증을 담당합니다.

### Figure 1. 공유 탐색 범위에서의 정확도–정밀도 절충

- **우선순위:** 1
- **형식:** 기존 Figure 1 재설계, 3패널
- **패널 구성:**
  - **A:** region precision–recall 산점도. BWTandem P/B/F/H를 선으로 연결하고 ULTRA, tantan, TRF, TRASH를 개별 점으로 표시.
  - **B:** bp precision–recall 산점도. 동일한 색·기호 체계 유지.
  - **C:** overlap rule별 matched-range region recall을 선 그래프로 표시. x축은 one-base, reciprocal 0.25, reciprocal 0.50; ULTRA, BWTandem, tantan만 강조하고 TRF/TRASH는 흐리게 표시하거나 보충자료로 이동.
- **정확한 데이터 출처:**
  - A–B: `results/figures/figure_curve_data.csv`
    - 행: BWTandem P/B/F/H 및 competitor ULTRA/tantan/TRF/TRASH.
    - 열: `region_recall`, `region_precision`, `bp_recall`, `bp_precision`.
    - 핵심 점: ULTRA 81.62/53.66 region recall/precision; BWTandem-F 79.88/50.62; H 81.60/48.44.
  - C:
    - `results/regen/recip_none.txt`, `MATCHED RANGE`의 `regRecall%`
    - `results/regen/recip_0.25.txt`, 같은 블록과 열
    - `results/regen/recip_0.50.txt`, 같은 블록과 열
    - 검증된 ULTRA/BWTandem 값: 81.62/78.87, 72.67/67.31, 31.05/17.61.
- **방어하는 주장:** BWTandem은 공유 범위 정확도에서 선두가 아니지만 경쟁 가능한 operating region에 있으며, ULTRA의 region-recall 1위는 overlap rule 변화에도 유지된다는 정직한 포지셔닝.
- **설계 주의:** A–B의 F/H는 native `--max-period 100` 실행인 반면 C의 BWTandem은 full-range 결과를 ≤100으로 사후 필터링한 값이다. 따라서 C에서 P/B/F/H 연결선을 재사용하면 안 되고 “regenerated full-range output, post-hoc ≤100”으로 별도 표기해야 한다.
- **seaborn 난이도:** **중간.** 산점도 자체는 쉽지만, 점 라벨 충돌 회피와 서로 다른 BWTandem 실행 계열을 시각적으로 구분해야 한다.

### Figure 2. 요청 최대 period에 따른 계산비용과 인간 전장 실행 비용

- **우선순위:** 2
- **형식:** 3패널
- **패널 구성:**
  - **A:** BWTandem 인간 paired runs. 각 replicate를 100→2,000 bp로 잇는 paired slope plot; y축 runtime, x축 maximum period. 세 쌍을 각각 연결하고 1.82×, 1.77×, 1.77×를 주석으로 표시.
  - **B:** maize에서 ULTRA와 TRF의 maximum period–runtime 관계. x축과 y축 모두 로그 척도 권장. ULTRA와 TRF는 각각 선으로 연결하되 “cross-tool mechanism comparison이 아닌 observed scaling”이라고 명시.
  - **C:** 인간 실행의 core-hours와 memory를 나란히 보여주는 두 개의 작은 horizontal dot plot. BWTandem, ULTRA, TRF만 포함하고 각 점 옆에 탐색 상한을 `2,000`, `100`, `500 bp`로 직접 표기.
- **정확한 데이터 출처:**
  - A: `results/manifest.tsv`, `range-rep` 행:
    - p100-r1 04:01:05 → p2000-r1 07:18:21
    - p100-r2 04:06:38 → p2000-r2 07:16:41
    - p100-r3 04:01:22 → p2000-r3 07:07:45
    - 원고 §3.1의 비율: 1.82, 1.77, 1.77(평균 1.79).
  - B:
    - ULTRA: manuscript Table 3A row, line 203, 0.23 h at period 6; Table 3C row, line 278, 12.9 h at period 200; Table 3B row, line 233, 34.7 h at period 500.
    - TRF: Table 3A line 202, 5.2 h at period 6; Table 3C line 277, 5.5 h at period 200; Table 3B line 232, 5.5 h at period 500.
  - C:
    - Core-hours: abstract line 21 및 Introduction line 42—BWTandem 25.3, ULTRA 59.6, TRF 33.7 core-h.
    - Memory: Table 1a lines 113–116—BWTandem 28.08 GB, ULTRA 1.68 GB, TRF 1.45 GB.
    - Wall clock와 범위: 같은 Table 1a—12.6 h/2,000 bp, 29.8 h/100 bp, 33.7 h/500 bp.
- **방어하는 주장:** BWTandem의 핵심 장점은 “정확도 1위”가 아니라 maximum-period 범위를 20배 늘렸을 때 릴리스 빌드 런타임이 평균 1.79배로 선형 미만 증가했고, 인간 전장 실행에서 ULTRA보다 적은 측정 core-hours를 사용했다는 점.
- **설계 주의:**
  - A의 `0363d8b` narrow arm은 `--max-period` 100에서 Tier 3의 장주기 탐색을 수행하지 않는다. 축 제목은 **maximum period**이며, narrow arm이 그 작업을 수행하고 버렸던 `07ad6fa` 1.30–1.41 비율은 superseded다.
  - B는 maize, A/C는 human이다. 패널 제목에 genome을 크게 명시하고 하나의 공통 회귀선이나 scaling exponent를 계산하지 않는다.
  - C의 실행은 range-matched가 아니며 FASTA scope도 다르다. 이를 캡션과 패널 내 각주에 표시한다.
- **seaborn 난이도:** **중간.** paired plot과 로그축은 간단하지만, 서로 다른 genome·범위·provenance를 오해 없이 분리하는 주석 설계가 중요하다.

### Figure 3. 민감도 설정과 검증되지 않은 고유 call의 비용

- **우선순위:** 3
- **형식:** 2패널
- **패널 구성:**
  - **A:** full-range identity sweep의 region recall과 region precision. x축을 `off`, `0.80`, `0.76`, `0.72`, `0.68` 순서의 범주형으로 두고 두 metric을 별도 선 또는 같은 패널의 명확히 구분된 점선으로 표시.
  - **B:** audit verdict의 period stratum별 100% stacked bar. `SUPPORTED`, `UNSUPPORTED`, `UNSURE`를 1–6, 7–20, 21–100, 101–2,000 bp별로 표시하고 전체 4/346/50을 우측에 별도 요약.
- **정확한 데이터 출처:**
  - A: `results/regen/score_table1_idsweep.txt`, `BASELINE` 블록:
    - `BWT-id-off`: 62.59 recall, 54.72 precision
    - 0.80: 71.36, 53.42
    - 0.76: 71.53, 53.35
    - 0.72: 80.53, 50.50
    - 0.68: 81.46, 48.92
    - 동일 값은 manuscript Supplementary Table S3, lines 485–493에도 있음.
  - B: `results/audit11/aggregate_reviewer2_20260831.txt`
    - 1–6: 3/59/38
    - 7–20: 1/89/10
    - 21–100: 0/99/1
    - 101–2,000: 0/99/1
    - 전체: 4/346/50
    - definitive verdict 기준 Wilson 95% CI 0.4–2.9%; 보수적 전체 supported rate 1.0%.
- **방어하는 주장:** recall-favouring setting은 precision 하락과 연결되며, 다른 도구와 catalog 모두에서 벗어난 call은 수동 audit에서 대부분 지지되지 않았다.
- **설계 주의:** A의 identity-sweep population과 B의 audited full-range unmatched-call population은 동일한 표본 단위가 아니다. 두 패널 사이에 화살표를 그어 인과관계처럼 보이게 하지 말고 “관련된 두 한계”로만 병치한다.
- **seaborn 난이도:** **낮음–중간.** 범주형 선 그래프와 stacked bar는 단순하지만 Wilson CI가 전체 요약에만 해당함을 정확히 표시해야 한다.

### Figure 4. 식물 위성 반복에서의 coverage–boundary–period assignment 절충

- **우선순위:** 4
- **형식:** 3패널
- **패널 구성:**
  - **A:** Arabidopsis Col-CEN의 centromere coverage와 CEN180 monomer recall을 두 개의 aligned dot plot으로 표시. BWTandem, ULTRA, TRF, tantan 중심; template-guided 도구는 다른 기호로 구분.
  - **B:** maize knob180, TR-1, CentC별 unfiltered coverage를 grouped dot plot으로 표시. 막대보다 점이 1–3% 차이를 과장하지 않아 적절하다.
  - **C:** 동일 세 class에서 period-band 적용 전후 coverage를 slope plot으로 표시. BWTandem, TRF, ULTRA, tantan 각각 unfiltered→banded 연결.
- **정확한 데이터 출처:**
  - A: manuscript Table 2, lines 180–189:
    - 열 `Cen. Coverage (%)`, `CEN180 Monomer Recall (%)`
    - 예: BWTandem 84.54/99.72, ULTRA 84.44/99.80, TRF 84.39/97.55, tantan 81.50/99.24.
    - 실행 provenance는 `results/manifest.tsv`의 experiment `2` 행.
  - B:
    - `results/regen/table3bc_replacement.md`, Table 3B/3C replacement의 coverage 열.
    - 또는 manuscript Table 3B-b lines 241–260 및 Table 3C lines 273–279.
    - BWTandem: knob180 79.79, TR-1 50.34, CentC 58.55.
    - ULTRA: 78.71, 36.36, 57.73.
    - TRF: 80.01, 45.77, 58.50.
    - tantan: 71.92, 22.74, 56.51.
  - C:
    - `results/regen/table3bc_replacement.md`, Table 3B-b/3C-b replacement.
    - BWTandem 손실: −15.59, −33.01, −17.84 pp.
    - ULTRA: −1.20, −2.82, −1.89 pp.
    - TRF: −0.63, −1.54, −0.87 pp.
    - tantan: −0.50 또는 파일 정밀값 −0.51, −1.38, −0.91 pp. 원고 반올림값과 파일값이 일부 0.01 pp 다르므로 figure에는 source file 값을 일관되게 사용한다.
- **방어하는 주장:** BWTandem은 식물 위성 반복에서 unfiltered coverage상 leading group에 있지만 일관된 선두는 아니며, period assignment의 불안정성이 실질적인 한계다.
- **seaborn 난이도:** **중간–높음.** class, tool, scoring rule의 세 차원을 동시에 보여주되 coverage와 offset을 한 축에 섞지 않아야 한다.

### Supplementary Figure S1. 염색체 분할에 대한 인간 benchmark 순위 안정성

- **우선순위:** 보충 1
- **패널 구성:**
  - **A:** matched-range region recall
  - **B:** matched-range region precision
  - x축: 22-chromosome subset, chr21/22 subset, all 24 chromosomes; 선: BWTandem, ULTRA, tantan.
- **정확한 데이터 출처:**
  - `results/regen/heldout_22chrom.txt`
  - `results/regen/heldout_select_chr21_22.txt`
  - `results/regen/heldout_all24.txt`
  - 각 파일의 `MATCHED RANGE` 블록, `regRecall%`, `regPrec%`.
  - 예: all24에서 BWTandem 78.87/50.13, ULTRA 81.62/53.66; chr21/22에서 80.03/51.49와 82.14/56.65.
- **방어하는 주장:** 공유 범위에서 ULTRA 우위와 BWTandem의 근접한 성능이 선택된 염색체 분할에서 크게 뒤집히지 않는다.
- **설계 주의:** 파일명이 “heldout”이더라도 설정 선택 과정과 완전히 독립적인 external validation이라는 증거는 이 파일만으로 확인되지 않는다. 제목에서 “held-out validation” 대신 “chromosome-subset sensitivity analysis”를 사용한다.
- **seaborn 난이도:** **낮음.**

### Supplementary Figure S2. Maize satellite coordinate post-merge 절충

- **우선순위:** 보충 2
- **패널 구성:**
  - 각 class(knob180, TR-1, CentC)를 열로 둔 3패널.
  - x축: coordinate-only merge gap 0–10,000 bp, 로그 또는 symlog.
  - 왼쪽 y축 대신 두 개의 세로 facet 사용: coverage와 mean boundary offset. 동일 그래프의 이중 y축은 피한다.
  - BWTandem 중심으로 표시하고, 필요하면 다른 도구는 보조선으로 추가.
- **정확한 데이터 출처:**
  - `results/regen/maize_extra_evidence.json`, `coordinate_postmerge` 객체의 class/tool/gap별 `coverage_percent`, `mean_offset_bp`, `calls_per_array`.
  - 텍스트 검증: `results/regen/score_maize_regen.txt`.
  - 예: BWTandem knob180은 gap 0에서 79.79%/253 bp/286.6 calls per array, gap 2,000에서 90.06%/6,215 bp/37.8, gap 5,000에서 98.38%/130,412 bp/4.3.
- **방어하는 주장:** fragmentation을 후처리로 줄이면 coverage는 증가하지만 boundary fidelity가 급격히 악화되므로 “전체 array를 검출했다”는 식의 단순화가 부적절하다.
- **seaborn 난이도:** **높음.** 로그축, 3개 지표, 여러 gap 수준을 오해 없이 분리해야 한다.

---

## 2. Application Note로 전환할 때 남길 2개 도판

### 1순위: Figure 2 — 범위–비용과 인간 전장 계산비용

이 도판이 방법의 가장 차별적인 가치 제안을 직접 보여준다. BWTandem은 정확도 1위가 아니므로 Application Note에서는 “무엇을 가능하게 하는가”가 중요하며, 100→2,000 bp에서 1.77–1.82배(평균 1.79) 비용 증가와 인간 실행의 25.3 대 59.6 core-hours가 그 답이다. 단, genome과 범위가 일치하지 않는 비교를 패널별로 명확히 분리해야 한다.

### 2순위: Figure 1 — 공유 범위 정확도–정밀도 절충

비용 도판만 남기면 “정확도를 희생해 빨라진 도구”로 보일 수 있다. Figure 1은 ULTRA가 1위임을 숨기지 않으면서 BWTandem-F/H가 실용적인 경쟁 영역에 있다는 점을 보여준다. reciprocal-overlap 패널까지 포함하면 결과 선택성에 대한 방어도 가능하다.

Figure 3의 audit는 중요한 한계지만 표 한 개로 압축 가능하고, Figure 4의 식물 결과는 Application Note의 제한된 지면에서는 보충자료로 이동하는 편이 낫다.

---

## 3. 명시적으로 만들지 말아야 할 도판

### 피해야 할 Figure A. “20× wider and 2.4× faster” 단일 막대그래프

- **위험:** 인간 BWTandem과 ULTRA의 서로 다른 range와 FASTA scope를 하나의 공정한 speed benchmark처럼 보이게 한다.
- **노출되는 반박:** ULTRA는 100 bp cap, BWTandem은 2,000 bp cap이며 경쟁 도구 입력이 약 5% 더 컸고 ULTRA 비용 provenance는 inherited/approximate다.
- **대안:** Figure 2처럼 paired BWTandem scaling, maize scaling, human cost를 분리한다.

### 피해야 할 Figure B. Full-range bp recall만 비교하는 도구 순위 막대그래프

- **위험:** BWTandem 43.34% 대 ULTRA 38.14%를 정확도 우위처럼 보이게 한다.
- **노출되는 반박:** `results/regen/score_table1_p100.txt`의 baseline은 서로 다른 검색 범위를 섞는다. matched range에서는 ULTRA bp recall 38.14%, BWTandem-F 32.55%다.
- **대안:** Figure 1의 matched-range 패널을 사용하고 long-period stratum은 표 또는 보충 inset으로만 제시한다.

### 피해야 할 Figure C. 0.02 pp “tie” 확대 inset

- **위험:** y축을 80–82%로 잘라 H 81.60과 ULTRA 81.62만 강조하면 정밀도 손실과 다른 overlap rule 결과를 숨긴다는 인상을 준다.
- **노출되는 반박:** H의 region precision은 48.44%로 ULTRA 53.66%보다 5.22 pp 낮고, reciprocal 0.50에서는 ULTRA 31.05% 대 BWTandem 17.61%다.
- **대안:** 전체 precision–recall 공간과 overlap-rule sensitivity를 함께 보여준다.

### 피해야 할 Figure D. “BWTandem만 전체 maize satellite array를 검출” genome-browser 예시

- **위험:** cherry-picking이며 저장소 결과와도 맞지 않는다.
- **노출되는 반박:** manuscript lines 222와 263은 모든 평가 도구가 모든 curated knob180, TR-1, CentC array를 검출했다고 보고한다. BWTandem 자체도 array당 수백 call로 파편화되어 있다.
- **대안:** Figure 4 또는 Supplementary Figure S2처럼 전체 class의 coverage–boundary trade-off를 표시한다.

### 피해야 할 Figure E. Unique-vs-shared “분포” violin/box plot

- **위험:** 실제 raw observation-level 분포 파일 없이 요약 통계만으로 violin이나 box distribution을 재구성하게 된다.
- **근거:** `results/regen/s2_F_p100.txt`와 `s2_P_p100.txt`에는 median, quartile, 평균 entropy, period-category 비율만 있고 개별 call 값은 없다.
- **추가 반박 위험:** F의 unique call은 shared보다 평균 entropy가 오히려 높다(1.241 대 1.030 bits). 따라서 “unique calls are lower entropy”라는 단순 시각 서사는 파일 내용과 반대다.
- **대안:** Supplementary Table S2를 유지한다. 시각화가 꼭 필요하면 요약값 bar/dot plot임을 명시하되 novelty 증거로 읽지 않는다.

### 피해야 할 Figure F. Audit supported rate를 모집단 전체의 진실률로 직접 외삽한 pie chart

- **위험:** stratified sample의 단순 4/400을 809,886개 전체에 그대로 적용한 모집단 추정처럼 보인다.
- **노출되는 반박:** 각 period stratum에서 100개씩 뽑은 균등 층화 표본이며 실제 모집단 stratum 크기로 가중된 추정치가 제공되지 않았다. 또한 single-reader audit다.
- **대안:** Figure 3B처럼 stratum별 raw verdict count와 전체 기술 요약만 제시한다.

### 피해야 할 Figure G. Col-CEN에서 BWTandem의 “위성 검출 우위” 막대그래프

- **위험:** 99.72% monomer recall만 확대하면 BWTandem이 최고처럼 보인다.
- **노출되는 반박:** ULTRA는 99.80%, BWTandem의 CEN180 bp precision은 비교 가능한 검출 도구 중 가장 낮은 57.08%이며 tantan은 0.20 h로 더 빠르다(Table 2).
- **대안:** Figure 4A처럼 coverage, recall, tool class를 함께 보여준다.

---

## 4. 모든 제안·회피 도판의 한 문장 캡션 초안

### 본문 및 보충 도판

- **Figure 1:** “GRCh38의 공통 period 범위(≤100 bp)에서 BWTandem operating points는 recall 증가와 precision 감소의 절충을 보이며, ULTRA는 세 가지 region-overlap 규칙 모두에서 가장 높은 region recall을 유지한다.”

- **Figure 2:** “BWTandem의 `0363d8b` 인간 paired runs에서 maximum period를 100에서 2,000 bp로 늘릴 때 runtime이 1.77–1.82배(평균 1.79) 증가한 반면, 별도의 maize 실행에서 ULTRA의 비용은 maximum period와 함께 급증했으며, 인간 cross-tool 비용 비교는 탐색 범위와 FASTA scope가 일치하지 않는다. `07ad6fa`의 좁은 arm은 버릴 장주기 탐색 비용을 내서 이전 비율을 너무 낮췄다.”

- **Figure 3:** “Catch-all identity 기준을 완화하면 full-range region recall이 62.59%에서 81.46%로 증가하고 precision은 54.72%에서 48.92%로 감소하며, 다른 도구와 catalog에 모두 없는 stratified calls의 단일 판독 audit에서는 400개 중 4개만 지지되었다.”

- **Figure 4:** “Arabidopsis와 maize에서 BWTandem의 unfiltered satellite coverage는 leading group에 속하지만 일관된 선두는 아니며, nominal period-band filtering은 BWTandem coverage를 경쟁 도구보다 훨씬 크게 감소시킨다.”

- **Supplementary Figure S1:** “공통 period 범위의 region recall과 precision 순위는 22개 염색체, chr21/22 및 전체 24개 염색체 분석에서 대체로 유지되며, 이 분석은 external validation이 아닌 chromosome-subset sensitivity analysis이다.”

- **Supplementary Figure S2:** “Maize satellite calls의 coordinate-only 병합 간격을 늘리면 BWTandem의 coverage와 contiguity는 증가하지만 mean boundary offset도 급격히 증가하여 단일 연속 array 호출과 정확한 경계 사이의 절충을 드러낸다.”

### 만들지 말아야 할 도판

- **Avoid Figure A:** “서로 다른 maximum period와 FASTA scope에서 측정한 BWTandem과 ULTRA runtime을 하나의 직접 speedup으로 표시하는 것은 공정한 range-matched 비교가 아니다.”

- **Avoid Figure B:** “서로 다른 탐색 범위의 full-output bp recall 순위는 범위 효과와 검출 정확도를 분리하지 못하므로 도구 정확도 순위로 해석할 수 없다.”

- **Avoid Figure C:** “81.60%와 81.62%만 확대하면 BWTandem-H의 낮은 precision과 stricter reciprocal-overlap 규칙에서 확대되는 ULTRA 우위를 가린다.”

- **Avoid Figure D:** “모든 평가 도구가 모든 curated maize array를 검출했고 BWTandem도 고도로 파편화되어 있으므로 BWTandem만 전체 array를 검출한다는 예시 도판은 근거가 없다.”

- **Avoid Figure E:** “개별 call-level 자료 없이 unique/shared 요약 통계로 재구성한 violin 또는 box plot은 존재하지 않는 분포 정보를 암시한다.”

- **Avoid Figure F:** “균등 층화된 single-reader audit의 4/400 결과를 전체 unmatched-call 모집단의 비가중 진실률로 표시해서는 안 된다.”

- **Avoid Figure G:** “Col-CEN monomer recall만 표시하면 ULTRA의 더 높은 recall, BWTandem의 낮은 bp precision, tantan의 더 낮은 비용을 누락하게 된다.”

Codex session ID: 01a0599a-f533-7141-b686-62899f84a08c
Resume in Codex: codex resume 01a0599a-f533-7141-b686-62899f84a08c
