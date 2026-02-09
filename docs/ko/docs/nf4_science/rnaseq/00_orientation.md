# 오리엔테이션

교육 환경에는 이 교육 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 및 데이터가 포함되어 있으므로 별도로 설치할 필요가 없습니다.
다만, 로그인하려면 (무료) 계정이 필요하며, 인터페이스에 익숙해지는 데 몇 분 정도 시간을 할애해야 합니다.

아직 완료하지 않았다면, 더 진행하기 전에 [환경 설정](../../envsetup/) 단기 과정을 먼저 학습하세요.

## 제공되는 자료

이 교육 과정 전반에 걸쳐 `nf4-science/rnaseq/` 디렉토리에서 작업하게 됩니다. 교육 워크스페이스를 열 때 이 디렉토리로 이동해야 합니다.
이 디렉토리에는 필요한 모든 코드 파일, 테스트 데이터 및 보조 파일이 포함되어 있습니다.

이 디렉토리의 내용을 자유롭게 탐색해 보세요. 가장 쉬운 방법은 VSCode 인터페이스의 교육 워크스페이스 왼쪽에 있는 파일 탐색기를 사용하는 것입니다.
또는 `tree` 명령을 사용할 수도 있습니다.
과정 전반에 걸쳐 `tree` 명령의 출력을 사용하여 디렉토리 구조와 내용을 읽기 쉬운 형태로 표시하며, 때로는 명확성을 위해 약간 수정하기도 합니다.

여기서는 두 번째 수준까지 목차를 생성합니다:

```bash
tree . -L 3
```

??? success "디렉토리 내용"

    ```console
    rnaseq
    ├── data
    │   ├── genome.fa
    │   ├── paired-end.csv
    │   ├── reads
    │   │   ├── ENCSR000COQ1_1.fastq.gz
    │   │   ├── ENCSR000COQ1_2.fastq.gz
    │   │   ├── ENCSR000COQ2_1.fastq.gz
    │   │   ├── ENCSR000COQ2_2.fastq.gz
    │   │   ├── ENCSR000COR1_1.fastq.gz
    │   │   ├── ENCSR000COR1_2.fastq.gz
    │   │   ├── ENCSR000COR2_1.fastq.gz
    │   │   ├── ENCSR000COR2_2.fastq.gz
    │   │   ├── ENCSR000CPO1_1.fastq.gz
    │   │   ├── ENCSR000CPO1_2.fastq.gz
    │   │   ├── ENCSR000CPO2_1.fastq.gz
    │   │   └── ENCSR000CPO2_2.fastq.gz
    │   └── single-end.csv
    ├── nextflow.config
    ├── rnaseq.nf
    └── solutions
        ├── modules
        │   ├── fastqc.nf
        │   ├── fastqc_pe.nf
        │   ├── hisat2_align.nf
        │   ├── hisat2_align_pe.nf
        │   ├── multiqc.nf
        │   ├── trim_galore.nf
        │   └── trim_galore_pe.nf
        ├── rnaseq-2.1.nf
        ├── rnaseq-2.2.nf
        ├── rnaseq-2.3.nf
        ├── rnaseq-3.1.nf
        ├── rnaseq-3.2.nf
        └── rnaseq_pe-3.3.nf
    ```

!!!note

    내용이 많아 보여도 걱정하지 마세요. 과정의 각 단계에서 관련된 부분을 다룰 것입니다.
    이것은 단지 전체적인 개요를 제공하기 위한 것입니다.

**시작하기 위해 알아야 할 내용은 다음과 같습니다:**

- **`rnaseq.nf` 파일**은 개발할 워크플로우 스크립트의 개요입니다.

- **`nextflow.config` 파일**은 최소한의 환경 속성을 설정하는 구성 파일입니다. 지금은 무시해도 됩니다.

- **`data` 디렉토리**에는 입력 데이터와 관련 리소스가 포함되어 있습니다:

  - 인간 염색체 20의 작은 영역(hg19/b37)으로 구성된 `genome.fa`라는 _참조 게놈_.
  - 파일 크기를 줄이기 위해 작은 영역으로 부분 추출된 _RNAseq 데이터_가 `reads/` 디렉토리에 있습니다.
  - 일괄 처리를 위한 예제 데이터 파일의 ID와 경로를 나열하는 _CSV 파일_.

- **`solutions` 디렉토리**에는 과정의 각 단계에서 완성된 워크플로우 스크립트와 모듈이 포함되어 있습니다.
  이는 작업을 확인하고 문제를 해결하기 위한 참조 자료로 사용됩니다.
  파일 이름의 번호는 과정의 해당 부분 단계에 해당합니다.

!!!tip

    어떤 이유로든 이 디렉토리에서 벗어났다면, 다음 명령을 실행하여 언제든지 돌아올 수 있습니다:

    ```bash
    cd /workspaces/training/nf4-science/rnaseq
    ```

이제 과정을 시작하려면 이 페이지의 오른쪽 하단 모서리에 있는 화살표를 클릭하세요.
