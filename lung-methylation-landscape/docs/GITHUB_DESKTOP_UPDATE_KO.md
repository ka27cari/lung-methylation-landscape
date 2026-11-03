# GitHub Desktop: 새 코드로 교체하기

1. GitHub Desktop에서 `lung-methylation-landscape` 저장소를 선택하고 **Repository → Show in Explorer**로 실제 저장소 폴더를 열어.
2. 아직 이전 ZIP을 복사하지 않았다면, 이번 ZIP에서 `lung-methylation-landscape` **안에 있는 파일과 폴더**를 저장소 폴더에 그대로 복사하면 끝임.
3. 이전 ZIP을 복사했다면, 별도로 수정한 내용이 없는지 먼저 확인하고, 저장소 폴더 안의 예전 `pipeline/`와 `r/` 두 폴더를 삭제한 뒤 이번 ZIP의 파일들을 복사해 덮어써. 이전 번호 파일이 남아 있으면 안 됨. **숨김 폴더 `.git/`은 절대로 삭제하지 마.**
4. GitHub Desktop의 **Changes**에서 새 파일과 삭제된 이전 경로를 보고 `Organize publication analysis code`로 커밋하면 됨.
5. **Publish repository**에서 기관 공개 승인이 아직이라면 비공개(Keep this code private)를 선택해.

확인할 파일 구조:

```text
lung-methylation-landscape/
  README.md
  run.py
  pipeline/                 # Python 단계가 한 폴더에 모여 있음
  src/pfl/                  # 수치분석 함수
  r/analysis/               # 코호트 및 EWAS
  r/figures/                # 그림 분석 주제별 R 스크립트
  r/source_data/            # Source Data 내보내기
  tools/                    # 수치자료 도구
  tests/                    # Python 자동 테스트
  config/  data/  docs/  env/
```

**이 ZIP에 `rebuild/`, 개인자료, 원본 RDS는 포함돼 있지 않음.** 본문 Figure Source Data는 저널용 별도 파일임.
