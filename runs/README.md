# 실행 방법

Sentaurus가 설치된 서버에서 Sentaurus Workbench(SWB)로 실행한다. 모든 실행은 같은 코드(`code/A_thermal_budget_fps.txt`)를 쓰고, 파라미터 값만 다르다.

## 절차

1. SWB에서 새 프로젝트(Traditional)를 만들고 도구로 sprocess를 추가한다.
2. sprocess 명령 파일에 `code/A_thermal_budget_fps.txt` 내용 전체를 붙여 넣는다. 코드는 고치지 않는다.
3. 파라미터 5개를 이 순서로 만든다. 이름은 대소문자까지 똑같아야 한다.

   | 파라미터 | 뜻 | 기준값 |
   |---|---|---|
   | `Tspike` | spike RTA 피크 온도 (°C) | 1050 |
   | `Eimp` | BF₂ 주입 에너지 (keV) | 10 |
   | `pai` | Ge PAI 여부 (0 또는 1) | 0 |
   | `cdose` | 탄소 도즈 (cm⁻²), 0이면 없음 | 0 |
   | `Tpost` | 후속 열처리 온도 (°C) | |

4. 실행하려는 폴더의 `params.md`에 적힌 값을 입력한다. 바꾸지 않는 파라미터도 값을 하나씩 넣어야 한다.
5. 전체 노드를 실행한다. 노드 하나에 약 2분이 걸린다.
6. 끝나면 `n<노드번호>_A.csv`를 내려받아 해당 폴더에 넣는다.

## 환경이 같은지 확인하기

기준 레시피(Tspike 1050, Eimp 10, pai 0, cdose 0)에 Tpost 700을 넣고 한 노드를 돌렸을 때 아래 값이 나오면 같은 환경이다.

| 항목 | 값 |
|---|---|
| Xj0 | 79.70 nm |
| Rs0 | 145.8 Ω/□ |
| 60분 뒤 ΔXj | +8.30 nm |
| 240분 뒤 Xj | 103.8 nm |

실행 로그 첫머리의 버전이 T-2022.03인지도 본다.

## 폴더

| 폴더 | 내용 | 상태 |
|---|---|---|
| `partA_1114/` | 진단 8개 + 개선 32개 노드 | 완료 |
| `member1_carbon_dose/` | 탄소 양 5e14, 2e15 | 실행 예정 |
| `member2_base_midT/` | 기준 레시피 675·725·775 °C | 실행 예정 |
| `member3_paic_midT/` | PAI + 탄소 725·775 °C | 실행 예정 |

노드 번호는 프로젝트마다 새로 매겨진다. 서로 다른 실행의 CSV를 한 폴더에 섞으면 같은 이름끼리 덮어쓰이므로 폴더를 나눠 둔다.

## 후처리

```
python code/tcad_post.py winA runs/partA_1114/n*_A.csv --dxj 10 --drs 10 --out A_window.png
python code/tcad_post.py remA runs/partA_1114/n*_A.csv --tpost 750 --base noPAI_noC_E10_spk1050 --out A_remedy_750C.png
```

`remA`는 `--base noPAI_noC_E10_spk1050`을 줘야 기준 레시피가 맞게 잡힌다.

명령 파일 확장자를 `.txt`로 둔 것은 Windows에서 `.cmd`를 더블클릭하면 배치 파일로 실행되기 때문이다.
