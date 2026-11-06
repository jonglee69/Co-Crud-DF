# co-crud-df — Co 도핑 모의 crud 표면제염 DF 계산기

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/<GITHUB_USER>/co-crud-df/blob/main/DF_calculator.ipynb)

이 계산기는 STS 모재 위 마그네타이트 crud(Co 도핑)를 표면제염한 결과에서 **제염계수(DF)**를 객관적으로 산출한다.

- 제염 후에는 crud가 모두 제거되므로 잔류 Co가 검출되지 않는다. 그렇다고 "오염도 0"이라고 쓸 수는 없다.
- 대신 **Currie(1968) / ISO 11929 검출한계 이론**으로 잔류 Co의 상한을 구한다.
- 결과는 **"DF ≥ 하한값 (95% 신뢰)"** 형태로 보고한다.

방법의 근거와 실험 절차는 [`docs/지침_DF정량평가.md`](docs/지침_DF정량평가.md)를 먼저 읽을 것.

---

## 1. 바로 실행하기

### (a) 브라우저에서 바로 실행 (설치 불필요)
1. 위의 **Open in Colab** 배지를 누른다.
2. 셀을 위에서부터 차례로 실행한다(Shift+Enter).
3. 예제 결과를 확인한 뒤, 자기 CSV를 업로드해 계산한다.

### (b) 내 PC에서 실행
```bash
git clone https://github.com/<GITHUB_USER>/co-crud-df.git
cd co-crud-df
pip install -r requirements.txt
python df_calc.py examples/example_data.csv
```
결과는 `results/` 폴더에 저장된다.

| 파일 | 내용 |
|---|---|
| `report.md` | 최종 DF, 특성 한계, 시편별 판정, 방법 검증, 보고 문장 초안 |
| `results.csv` | 시편별 계산값 (x, u(x), σ0, L_C, L_D, L_Q, 회수율 등) |
| `plot.png` | 시편별 과잉 Co와 L_C·L_D·L_Q 선 |

---

## 2. 입력 데이터 작성

[`templates/template.csv`](templates/template.csv)를 복사해서 실측값을 넣는다. Excel에서 편집한 뒤 **CSV(UTF-8)**로 저장한다.

| 열 | 필수 | 단위 | 설명 |
|---|---|---|---|
| `sample_id` | ✔ | — | 시편 번호 |
| `group` | ✔ | — | A, L, B, C, D, E, F (아래 표) |
| `area_cm2` | ✔ | cm² | 분석 대상 면적(양면 합) |
| `Co_ug`, `Ni_ug` | ✔ | µg | 용해액 **전체**에 든 질량 = ICP 농도(µg/L) × 용액 부피(L) × 희석배수 |
| `u_Co_ug`, `u_Ni_ug` | | µg | 측정 표준불확도. 비우면 상대 2%로 계산 |
| `dissolved_mass_mg` | 권장 | mg | 분석을 위해 녹인 질량(전후 질량차). 층 두께 계산에 쓴다 |
| `crud_mass_mg` | A만 | mg | crud 생성 시 질량 증가. A0 방법 b에 쓴다 |
| `spike_Co_ug` | D만 | µg | 첨가한 Co 표준 질량 |
| `note` | | — | 메모 |

| group | 시편 | 최소 수 | 용도 |
|---|---|---|---|
| **A** | 제염 전 시편, 산화막만 선택 용해 | 3 | 초기 crud Co, A0 (방법 a·b) |
| **L** | 제염액 | 1 이상 | A0 방법 c (물질수지, 확인용) |
| **B** | 제염 시편, 표면 ~1 µm 용해 | 3 | 잔류 Co 판정 |
| **C** | 공정 블랭크 (crud 없는 모재, B와 동일 처리) | **5** | 모재 Co/Ni 비 r, σ0 |
| **D** | B와 같은 시편에 Co 스파이크 | 2 | 회수율 85~115% |
| **E** | 불완전 제염 (양성대조) | 2 | x > L_Q 로 검출 능력 입증 |
| F | 모재 원소재 전체 용해 | — | 참고용 (계산에는 쓰지 않음) |

> ⚠️ **시편 전체를 녹인 데이터는 넣지 말 것.** crud의 Co가 벌크 모재 Co에 희석되어 의미 있는 DF가 나오지 않는다.

---

## 3. 계산 내용

1. **모재 Co/Ni 비**: 블랭크(C)에서 r = Co/Ni를 구한다 (평균 r̄, 표준편차 s_r).
2. **과잉 Co**: x = (Co − r̄·Ni) / A
   - u(x)² = [u(Co)² + r̄²·u(Ni)² + Ni²·s_r²] / A²   (GUM 불확도 전파)
   - 블랭크 자신의 x는 leave-one-out r로 다시 계산한다. 자기 데이터로 만든 r을 쓰면 x가 0에 가깝게 나와 과소평가되기 때문이다.
3. **σ0**: max(순신호 0일 때 전파한 불확도, 블랭크 LOO 산포). 경험값은 배경 Co에 비례하도록 환산해 각 시편에 적용한다.
4. **특성 한계** (α = β = 0.05, 자유도 n_C − 1 인 t 분포)
   - 결정 한계 L_C = k_α·σ0
   - 검출 한계 L_D = (k_α + k_β)·σ0
   - 정량 한계 L_Q = 10·σ0
5. **A0**: 방법 a(기본값), b, c의 평균과 단측 하한 A0,low. 방법 사이의 일치도도 검사한다.
6. **DF**
   - 경우 1, 모든 B가 x < L_C: **DF ≥ A0,low / max(L_D)**
   - 경우 2, 일부가 검출되었지만 L_Q 미만: DF ≥ A0,low / max(x̄의 단측 상한, L_D)
   - 경우 3, 모든 B가 x ≥ L_Q: DF = Ā0 / x̄ ± u
7. **검증**: 블랭크 평균 ≈ 0, 스파이크 회수율, 양성대조 검출, A0 방법 간 일치, 분석 층 두께 ≤ 10 µm.

### 옵션
```bash
python df_calc.py data.csv --out results_run1 \
    --rel-u-co 0.02 --rel-u-ni 0.02 \   # 불확도 열이 비었을 때 쓰는 상대 정밀도
    --w-co-ppm 5000 --u-w-rel 0.05 \    # 방법 b: crud 내 Co 농도와 상대 불확도
    --a0-method a \                     # DF에 쓸 A0 방법 (a|b|c)
    --alpha 0.05 --beta 0.05 --kq 10 \
    --use-z                             # 블랭크가 많을 때(≥ 20) t 대신 z 사용
```

---

## 4. 예제 결과 해석

`examples/example_data.csv`는 다음 가정으로 만든 **가상 데이터**다.
- 모재: Co 1,500 ppm, Ni 9%
- crud: 5 µm, Co 5,000 ppm (A0 ≈ 13 µg/cm²)
- 분석 층: 1 µm

결과는 **DF ≥ 45**다. 지침서의 단순 계산(약 120)보다 작은 이유는 두 가지다.
- 블랭크가 5개뿐이라 t 분포를 쓰게 되어 k가 2.13이다.
- 모재 Co/Ni 비의 산포(약 4%)가 반영되었다.

**하한을 올리는 방법**
- 블랭크 수를 늘린다 (k 감소).
- 같은 히트의 모재를 쓰고, 용해 조건을 일정하게 해서 s_r을 줄인다.
- 분석 층을 더 얇게 만든다.

---

## 5. 참고문헌
- L. A. Currie, *Anal. Chem.* **40** (1968) 586–593.
- ISO 11929-1:2019, Determination of the characteristic limits for measurements of ionizing radiation.
- JCGM 100:2008 (GUM); Eurachem/CITAC *Quantifying Uncertainty in Analytical Measurement* (2012).
