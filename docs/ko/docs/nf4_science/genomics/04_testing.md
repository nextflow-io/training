# Part 4: 테스트 추가하기

이 교육 과정의 첫 번째 파트에서는 완전히 선형적이며 각 샘플의 데이터를 독립적으로 처리하는 변이 호출 파이프라인을 구축했습니다.

두 번째 파트에서는 GATK를 사용한 공동 변이 호출을 구현하기 위해 채널과 채널 연산자를 사용하는 방법을 보여드렸습니다.

세 번째 파트에서는 파이프라인을 모듈화했습니다.

교육의 이번 파트에서는 Nextflow와 잘 통합되며 파이프라인에 모듈 수준 및 workflow 수준 테스트를 추가하는 것을 간단하게 만들어주는 테스트 프레임워크인 [**nf-test**](https://www.nf-test.com/)를 사용하는 방법을 보여드리겠습니다. 이번 파트를 따라하려면 Part 1, Part 2, Part 3뿐만 아니라 nf-test의 기본 사항과 테스트가 중요한 이유를 다루는 [nf-test side quest](../../side_quests/nf-test.md)를 완료해야 합니다.

---

## 0. 워밍업

!!! note "참고"

    올바른 작업 디렉토리에 있는지 확인하십시오:
    `cd /workspaces/training/nf4-science/genomics`

이 교육 과정의 이전 파트를 완료했다면, 적절한 모듈 디렉토리 구조를 갖춘 작동하는 버전의 genomics 파이프라인을 가지고 있어야 합니다.

??? abstract "디렉토리 내용"

    ```console
    modules/
    ├── gatk
    │   ├── haplotypecaller
    │   │   └── main.nf
    │   └── jointgenotyping
    │       └── main.nf
    └── samtools
        └── index
            └── main.nf
    ```

이 모듈 디렉토리는 필요한 경우 `solutions` 디렉토리에서 찾을 수 있습니다.

Part 3과 동일한 workflow로 시작하겠습니다. 이는 `genomics-4.nf` 파일에 제공되어 있습니다. [nf-test side quest](../../side_quests/nf-test.md)와 마찬가지로, 이 파이프라인의 세 가지 프로세스와 workflow 수준 테스트에 몇 가지 다른 유형의 테스트를 추가할 것입니다.

### 0.1. workflow 실행 확인

테스트를 추가하기 전에 workflow가 예상대로 실행되는지 확인하십시오.

```bash
nextflow run genomics-4.nf -resume
```

이 교육 과정의 처음부터 작업해 오셨다면 이제 매우 익숙해 보일 것입니다.

??? success "명령 출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

    executor >  local (7)
    [18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
    [30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
    [a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
    ```

이전과 마찬가지로 이제 프로젝트 디렉토리 내에 `work` 디렉토리와 `results_genomics` 디렉토리가 생성됩니다. 나중에 테스트에서 이러한 결과를 실제로 사용할 것입니다. 하지만 지금부터는 파이프라인을 테스트하기 위해 `nf-test` 패키지를 사용할 것입니다.

### 0.2. `nf-test` 초기화

[nf-test side quest](../../side_quests/nf-test.md)와 마찬가지로 `nf-test` 패키지를 초기화해야 합니다.

```bash
nf-test init
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

    Project configured. Configuration is stored in nf-test.config
    ```

??? abstract "nf-test.config 내용"

    ```groovy title="nf-test.config"
    config {

        testsDir "tests"
        workDir ".nf-test"
        configFile "tests/nextflow.config"
        profile ""

    }
    ```

또한 설정 파일 스텁이 포함된 `tests` 디렉토리를 생성합니다.

### 핵심 요약

이제 genomics 파이프라인에 대한 테스트를 작성할 준비가 되었습니다.

### 다음 단계는?

프로세스 호출이 성공했고 올바른 출력을 생성했는지 평가하는 기본 테스트를 작성합니다.

---

## 1. 성공 및 일치하는 출력에 대한 프로세스 테스트

BAM 파일에 대한 인덱스 파일을 생성하여 효율적인 임의 접근을 가능하게 하는 `SAMTOOLS_INDEX` 프로세스를 테스트하는 것으로 시작하겠습니다. 이것은 다음과 같은 이유로 좋은 첫 번째 테스트 사례입니다:

1. 단일하고 명확하게 정의된 입력(BAM 파일)을 가집니다
2. 예측 가능한 출력(BAI 인덱스 파일)을 생성합니다
3. 동일한 입력에 대해 출력이 동일해야 합니다

### 1.1. 테스트 파일 스텁 생성

먼저 테스트 파일 스텁을 생성합니다:

```bash
nf-test generate process modules/samtools/index/main.nf
```

??? success "명령 출력"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

이것은 `main.nf`와 동일한 디렉토리에 파일을 생성합니다.
파일 탐색기에서 디렉토리로 이동하여 파일을 열면 다음 코드가 포함되어 있어야 합니다:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

시작 어설션은 [nf-test side quest](../../side_quests/nf-test.md)에서 익숙해야 합니다:

- `assert process.success`는 프로세스가 성공적으로 실행되고 오류 없이 완료될 것으로 예상한다는 것을 나타냅니다.
- `snapshot(process.out).match()`는 실행 결과가 이전 실행에서 얻은 결과와 동일할 것으로 예상한다는 것을 나타냅니다(해당되는 경우).
  이에 대해서는 나중에 더 자세히 설명합니다.

이것을 시작점으로 사용하여 samtools index 프로세스에 대한 올바른 테스트 입력과 해당하는 경우 매개변수를 추가해야 합니다.

### 1.2. 테스트 파일 이동 및 스크립트 경로 업데이트

테스트를 작성하기 전에 파일을 최종 위치로 이동해야 합니다. 각 모듈에 디렉토리를 추가한 이유 중 하나는 이제 각 모듈의 `main.nf` 파일과 함께 위치한 `tests` 디렉토리에 테스트를 배치할 수 있기 때문입니다. 해당 디렉토리를 생성하고 테스트 파일을 그곳으로 이동합니다.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

이제 테스트 파일의 `script` 섹션을 상대 경로로 단순화할 수 있습니다:

=== "변경 후"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "../main.nf"
    process "SAMTOOLS_INDEX"
    ```

=== "변경 전"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"
    ```

이것은 테스트에게 전체 경로를 지정하지 않고도 모듈의 `main.nf` 파일을 찾을 위치를 알려줍니다.

### 1.3. SAMTOOLS_INDEX에 대한 테스트 입력 제공

스텁 파일에는 `samtools index`의 입력에 적합한 실제 테스트 입력으로 교체해야 하는 자리 표시자가 포함되어 있습니다. 적절한 입력은 `data/bam` 디렉토리에서 사용할 수 있는 BAM 파일입니다.

=== "변경 후"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        input[0] = file("${projectDir}/data/bam/reads_son.bam")
        """
    }
    ```

=== "변경 전"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
    process {
        """
        // define inputs of the process here. Example:
        // input[0] = file("test-file.txt")
        """
    }
    ```

### 1.4. 기능에 따라 테스트 이름 지정

이전에 배운 것처럼 테스트 컨텍스트에서 의미가 있는 것으로 테스트 이름을 변경하는 것이 좋은 관행입니다.

=== "변경 후"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should index reads_son.bam correctly") {
    ```

    이것은 임의의 문자열을 사용하므로 원하는 것을 넣을 수 있습니다.
    여기서는 파일 이름과 형식을 참조하도록 선택합니다.

=== "변경 전"

    ```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
    test("Should run without failures") {
    ```

### 1.5. 테스트 실행 및 출력 검사

테스트를 실행합니다:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.717s)
      Snapshots:
        1 created [Should index reads_son.bam correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 7.727s
    ```

이전에 배운 것처럼, 이것은 프로세스의 성공에 대한 기본 어설션을 검증하고 프로세스의 출력을 기반으로 스냅샷 파일을 생성했습니다. `tests/modules/samtools/index/tests/main.nf.test.snap` 파일에서 스냅샷 파일의 내용을 볼 수 있습니다:

```json title="modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.3",
      "nextflow": "25.10.2"
    },
    "timestamp": "2026-01-27T15:09:48.394063389"
  }
}
```

테스트를 다시 실행하면 출력이 스냅샷과 동일하기 때문에 통과하는 것을 볼 수 있습니다:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.938s)


    SUCCESS: Executed 1 tests in 7.987s
    ```

### 1.6. SAMTOOLS_INDEX에 더 많은 테스트 추가

다양한 잠재적 문제를 테스트하기 위해 다양한 입력 파일을 테스트하는 것이 유용할 때가 있습니다. 테스트 데이터의 trio에서 어머니와 아버지의 BAM 파일에 대한 테스트를 추가합니다.

```groovy
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

그런 다음 테스트를 다시 실행할 수 있습니다:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process SAMTOOLS_INDEX

      Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (7.185s)
      Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (6.576s)
      Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (6.31s)
      Snapshots:
        2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


    Snapshot Summary:
      2 created

    SUCCESS: Executed 3 tests in 20.117s
    ```

`--update-snapshot` 매개변수의 효과를 언급하는 경고에 주목하십시오.

!!! note "참고"

    여기서는 파이프라인의 과학적 출력을 시연하기 위해 이전에 사용했던 테스트 데이터를 사용하고 있습니다.
    프로덕션 환경에서 이러한 테스트를 운영할 계획이었다면 테스트 목적으로 더 작은 입력을 생성했을 것입니다.

    일반적으로 프로세스 기능을 평가하는 데 필요하고 충분한 가장 작은 데이터 조각을 사용하여 단위 테스트를 가능한 한 가볍게 유지하는 것이 중요합니다. 그렇지 않으면 총 실행 시간이 심각하게 증가할 수 있습니다.
    실행하는 데 너무 오래 걸리는 테스트 스위트는 신속함을 위해 건너뛸 가능성이 높은 테스트 스위트입니다.

### 핵심 요약

genomics 프로세스에 대한 첫 번째 모듈 테스트를 작성하여 `SAMTOOLS_INDEX`가 다른 BAM 파일에 대해 인덱스 파일을 올바르게 생성하는지 확인했습니다. 테스트 스위트는 다음을 보장합니다:

1. 프로세스가 성공적으로 실행됩니다
2. 인덱스 파일이 생성됩니다
3. 출력이 실행 간에 일관됩니다
4. 프로세스가 모든 샘플 BAM 파일에 대해 작동합니다

### 다음 단계는?

연결된 프로세스를 처리하기 위해 setup 메서드를 사용하여 genomics workflow의 다른 프로세스에 대한 테스트를 작성하는 방법을 배웁니다. 또한 VCF 파일인 출력에 예상되는 변이 호출이 포함되어 있는지 평가합니다.

---

## 2. 연결된 프로세스에 테스트 추가 및 내용 테스트

`GATK_HAPLOTYPECALLER`를 테스트하려면 프로세스에 `SAMTOOLS_INDEX` 출력을 입력으로 제공해야 합니다. `SAMTOOLS_INDEX`를 실행하고 출력을 검색한 다음 workflow의 테스트 데이터와 함께 저장하여 이를 수행할 수 있습니다. 이것은 실제로 세련된 파이프라인에 권장되는 접근 방식이지만, nf-test는 `setup` 메서드를 사용하는 대안적 접근 방식을 제공합니다.

setup 메서드를 사용하면 테스트 설정의 일부로 `SAMTOOLS_INDEX` 프로세스를 트리거한 다음 그 출력을 `GATK_HAPLOTYPECALLER`의 입력으로 사용할 수 있습니다. 이것에는 비용이 있습니다: `GATK_HAPLOTYPECALLER`에 대한 테스트를 실행할 때마다 `SAMTOOLS_INDEX` 프로세스를 실행해야 합니다. 그러나 아마도 우리는 여전히 workflow를 개발 중이고 나중에 변경해야 할 수 있는 테스트 데이터를 미리 생성하고 싶지 않을 수 있습니다. `SAMTOOLS_INDEX` 프로세스도 매우 빠르므로 출력을 미리 생성하고 저장하는 이점이 무시할 만할 수 있습니다. setup 메서드가 작동하는 방식은 다음과 같습니다.

### 2.1. 테스트 파일 생성 및 배치

이전과 마찬가지로 먼저 파일 스텁을 생성합니다:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

??? success "명령 출력"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/haplotypecaller/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/haplotypecaller/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

다음 테스트 스텁을 생성합니다:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 2.2. 테스트 파일 이동 및 스크립트 경로 업데이트

모듈의 `main.nf` 파일과 함께 위치한 테스트 파일을 위한 디렉토리를 생성합니다:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

그리고 테스트 스텁 파일을 그곳으로 이동합니다:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

마지막으로 스크립트 경로를 업데이트하는 것을 잊지 마십시오:

=== "변경 후"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "../main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

=== "변경 전"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
        name "Test Process GATK_HAPLOTYPECALLER"
        script "modules/gatk/haplotypecaller/main.nf"
        process "GATK_HAPLOTYPECALLER"
    ```

### 2.3. setup 메서드를 사용하여 입력 제공

`when` 블록 앞에 `setup` 블록을 삽입하여 원본 입력 파일 중 하나에서 `SAMTOOLS_INDEX` 프로세스 실행을 트리거할 수 있습니다. 또한 이전과 마찬가지로 테스트 이름을 의미 있는 것으로 변경하는 것을 기억하십시오.

=== "변경 후"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
        test("Should call son's haplotype correctly") {

            setup {
                run("SAMTOOLS_INDEX") {
                    script "../../../samtools/index/main.nf"
                    process {
                        """
                        input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                        """
                    }
                }
            }

            when {
    ```

=== "변경 전"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
    test("Should run without failures") {

        when {
    ```

그런 다음 테스트 입력을 지정하는 `when` 블록에서 해당 프로세스의 출력을 참조할 수 있습니다:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

해당 변경 사항을 적용하고 테스트를 다시 실행합니다:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Snapshots:
        1 created [Should call son's haplotype correctly]


    Snapshot Summary:
      1 created

    SUCCESS: Executed 1 tests in 40.555s
    ```

이전처럼 스냅샷 파일도 생성합니다.

### 2.4. 다시 실행하고 실패 관찰

흥미롭게도 정확히 같은 명령을 다시 실행하면 이번에는 테스트가 실패합니다.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? failure "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' FAILED (40.123s)

      java.lang.RuntimeException: Different Snapshot:
      [                                                                                           [
          {                                                                                           {
              "0": [                                                                                      "0": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ],                                                                                          ],
              "1": [                                                                                      "1": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "idx": [                                                                                    "idx": [
                  "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
              ],                                                                                          ],
              "vcf": [                                                                                    "vcf": [
                  "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
              ]                                                                                           ]
          }                                                                                           }
      ]                                                                                           ]

      Nextflow stdout:

      Nextflow stderr:


        Obsolete snapshots can only be checked if all tests of a file are executed successful.


    FAILURE: Executed 1 tests in 40.156s (1 failed)
    ```

오류 메시지는 두 실행의 스냅샷 간에 차이가 있음을 알려줍니다. 구체적으로 VCF 파일에 대한 md5sum 값이 다릅니다.

왜 그럴까요? 간단히 말해서, HaplotypeCaller 도구는 매번 다른(정의상) 타임스탬프를 VCF 헤더에 포함합니다.
결과적으로, 변이 호출 자체의 내용이 동일하더라도 파일이 동일한 md5sum을 가질 것으로 기대할 수 없습니다.

이것을 어떻게 처리할까요?

### 2.5. 내용 어설션 메서드를 사용하여 특정 변이 확인

문제를 해결하는 한 가지 방법은 [다른 종류의 어설션](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions)을 사용하는 것입니다.
이 경우 동일성을 주장하는 대신 특정 내용을 확인할 것입니다.
더 정확히는, 도구가 VCF 파일의 라인을 읽고 특정 라인의 존재를 확인하도록 할 것입니다.

실제로 `then` 블록의 두 번째 어설션을 다음과 같이 교체합니다:

=== "변경 후"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
            then {
                assert process.success
                assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
                assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
            }
    ```

=== "변경 전"

    ```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
    then {
        assert process.success
        assert snapshot(process.out).match()
    }
    ```

여기서는 VCF 출력 파일의 전체 내용을 읽고 내용 일치를 검색하는데, 이는 작은 테스트 파일에서는 괜찮지만 더 큰 파일에서는 하고 싶지 않을 것입니다.
대신 특정 라인을 읽도록 선택할 수 있습니다.

이 접근 방식은 테스트할 '신호'로 사용할 것을 더 신중하게 선택해야 합니다.
긍정적인 측면에서, 분석 도구가 추가 개발을 거치면서 '어려운' 기능(예: 희귀 변이)을 일관되게 식별할 수 있는지 매우 정밀하게 테스트하는 데 사용할 수 있습니다.

### 2.6. 다시 실행하고 성공 관찰

이런 식으로 테스트를 수정하면 테스트를 여러 번 실행할 수 있으며 일관되게 통과할 것입니다.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)


    SUCCESS: Executed 1 tests in 40.555s
    ```

### 2.7. 더 많은 테스트 추가

어머니와 아버지 샘플에 대한 유사한 테스트를 추가합니다:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. 테스트 명령 실행

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_HAPLOTYPECALLER

      Test [c5156c2b] 'Should call son's haplotype correctly' PASSED (40.53s)
      Test [10de94a8] 'Should call mother's haplotype correctly' PASSED (41.47s)
      Test [c0386fc7] 'Should call father's haplotype correctly' PASSED (45.556s)


    SUCCESS: Executed 3 tests in 127.586s
    ```

파이프라인의 두 번째 단계에 대한 기본 테스트 계획이 완료되었습니다. 세 번째이자 마지막 모듈 수준 테스트로 넘어갑니다!

### 핵심 요약

다음을 배웠습니다:

1. 다른 프로세스의 출력에 의존하는 프로세스를 테스트하는 방법
2. VCF 출력 파일에서 특정 genomic 변이를 확인하는 방법
3. 특정 내용을 확인하여 비결정적 출력을 처리하는 방법
4. 여러 샘플에 걸쳐 변이 호출을 테스트하는 방법

### 다음 단계는?

공동 유전자형 분석 단계에 대해 미리 생성된 테스트 데이터를 사용하는 테스트를 작성하는 방법을 배웁니다.

---

## 3. 미리 생성된 테스트 데이터 사용

공동 유전자형 분석 단계의 경우 다른 접근 방식을 사용합니다 - 미리 생성된 테스트 데이터를 사용합니다. 이것은 다음과 같은 경우에 종종 선호됩니다:

1. 여러 의존성을 가진 복잡한 프로세스
2. 실행하는 데 오랜 시간이 걸리는 프로세스
3. 안정적이고 프로덕션 파이프라인의 일부인 프로세스

### 3.1. 테스트 데이터 생성

이 섹션의 시작 부분에서 생성한 결과를 검사합니다:

```bash
tree results_genomics/
```

```console title="결과 디렉토리 내용"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── gvcf
│   ├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
│   ├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
│   ├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
│   ├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
│   ├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
│   └── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx
└── indexed_bam
    ├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
    ├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
    ├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
    ├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
    ├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
    └── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai

2 directories, 14 files
```

공동 유전자형 분석 단계는 haplotype caller 단계에서 생성된 VCF 파일을 인덱스와 함께 입력으로 필요로 합니다. 그래서 우리가 가지고 있는 결과를 `jointgenotyping` 모듈의 테스트 디렉토리로 복사합시다.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp results_genomics/gvcf/*.g.vcf results_genomics/gvcf/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

이제 이 파일들을 공동 유전자형 분석 단계에 대해 작성할 테스트의 입력으로 사용할 수 있습니다.

### 3.2. 테스트 파일 스텁 생성

이전과 마찬가지로 먼저 파일 스텁을 생성합니다:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

??? success "명령 출력"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/modules/gatk/jointgenotyping/main.nf'
    Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/gatk/jointgenotyping/main.nf.test'

    SUCCESS: Generated 1 test files.
    ```

다음 테스트 스텁을 생성합니다:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

### 3.3. 테스트 파일 이동 및 스크립트 경로 업데이트

이번에는 이미 모듈의 `main.nf` 파일과 함께 위치한 테스트 디렉토리가 있으므로 테스트 스텁 파일을 그곳으로 이동할 수 있습니다:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

그리고 스크립트 경로를 업데이트하는 것을 잊지 마십시오:

=== "변경 후"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "../main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

=== "변경 전"

    ```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"
    ```

### 3.4. 입력 제공

프로세스 입력 정의를 기반으로 입력을 채우고 그에 따라 테스트 이름을 변경합니다:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. 내용 어설션 사용

공동 유전자형 분석 단계의 출력은 또 다른 VCF 파일이므로 다시 내용 어설션을 사용할 것입니다.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

출력 파일에서 특정 변이의 내용을 확인함으로써 이 테스트는 다음을 검증합니다:

1. 공동 유전자형 분석 프로세스가 성공적으로 실행됩니다
2. 출력 VCF에 올바른 순서로 세 샘플이 모두 포함됩니다
3. 특정 변이가 다음과 함께 올바르게 호출됩니다:
   - 각 샘플에 대한 정확한 유전자형(아버지는 0/1, 어머니와 아들은 1/1)
   - 올바른 읽기 깊이 및 유전자형 품질
   - 대립 유전자 빈도(AF=0.833)와 같은 집단 수준 통계

전체 파일을 스냅샷하지는 않았지만 특정 변이를 확인함으로써 공동 유전자형 분석 프로세스가 예상대로 작동하고 있다고 확신할 수 있습니다.

### 3.6. 테스트 실행

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Process GATK_JOINTGENOTYPING

      Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (53.827s)


    SUCCESS: Executed 1 tests in 53.837s
    ```

테스트가 통과하여 공동 유전자형 분석 프로세스가 올바르게 다음을 수행함을 확인했습니다:

1. 개별 샘플 VCF를 결합합니다
2. 공동 변이 호출을 수행합니다
3. 실행 간에 일관된 유전자형 호출을 가진 다중 샘플 VCF를 생성합니다

### 핵심 요약

다음을 배웠습니다:

- 이전에 생성된 결과를 테스트 입력으로 사용하는 방법
- 미리 생성된 테스트 데이터를 사용하여 테스트를 작성하는 방법

### 다음 단계는?

전체 변이 호출 파이프라인이 엔드 투 엔드로 작동하는지 확인하기 위해 workflow 수준 테스트를 추가합니다.

---

## 4. workflow 수준 테스트 추가

이제 BAM 파일에서 공동 유전자형까지 전체 변이 호출 파이프라인을 테스트할 것입니다. 이것은 다음을 확인합니다:

1. 모든 프로세스가 함께 올바르게 작동합니다
2. 단계 간에 데이터가 적절하게 흐릅니다
3. 최종 변이 호출이 일관됩니다

### 4.1. workflow 테스트 생성

전체 파이프라인에 대한 테스트 파일을 생성합니다:

```bash
nf-test generate pipeline genomics-4.nf
```

??? success "명령 출력"

    ```console
    Load source file '/workspaces/training/nf4-science/genomics/genomics-4.nf'
    Wrote pipeline test file '/workspaces/training/nf4-science/genomics/tests/genomics-4.nf.test'

    SUCCESS: Generated 1 test files.
    ```

기본 테스트 스텁을 생성합니다:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

이름을 의미 있는 것으로 수정하기만 하면 됩니다(곧 이것이 왜 유용한지 알게 될 것입니다).

=== "변경 후"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run the pipeline without failures") {
    ```

=== "변경 전"

    ```groovy title="tests/genomics-4.nf.test" linenums="6" hl_lines="1"
        test("Should run without failures") {
    ```

!!! note "참고"

    이 경우 테스트 파일은 `nf-test`가 생성한 위치에 그대로 둘 수 있습니다.

### 4.2. 입력 매개변수 지정

여전히 입력을 지정해야 하는데, 이것은 모듈 수준 테스트와 비교하여 workflow 수준에서 약간 다르게 수행됩니다.
프로파일을 지정하는 것을 포함하여 이를 수행하는 몇 가지 방법이 있습니다.
그러나 더 간단한 방법은 `nf-test init`이 원래 `tests` 디렉토리에 생성한 `nextflow.config` 파일에서 `params {}` 블록을 설정하는 것입니다.

```groovy title="tests/nextflow.config" linenums="1"
/*
========================================================================================
    Nextflow config file for running tests
========================================================================================
*/

// Output directory for workflow outputs
outputDir = 'results_genomics'

/*
 * Pipeline parameters
 */

params {
    // Primary input (file of input files, one per line)
    reads_bam = "${projectDir}/data/sample_bams.txt"

    // Accessory files
    reference = "${projectDir}/data/ref/ref.fasta"
    reference_index = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict = "${projectDir}/data/ref/ref.dict"
    intervals = "${projectDir}/data/ref/intervals.bed"

    // Base name for final output file
    cohort_name = "family_trio"
}
```

테스트를 실행하면 `nf-test`는 이 설정 파일을 가져와서 그에 따라 입력을 가져옵니다.

### 4.3. workflow 테스트 실행

```bash
nf-test test tests/genomics-4.nf.test
```

??? success "명령 출력"

    ```console
    🚀 nf-test 0.9.3
    https://www.nf-test.com
    (c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


    Test Workflow genomics-4.nf

      Test [1b4c6936] 'Should run the pipeline without failures' PASSED (171.019s)


    SUCCESS: Executed 1 tests in 171.056s
    ```

테스트가 통과하여 전체 변이 호출 파이프라인이 다음을 수행함을 확인했습니다:

1. 모든 샘플을 성공적으로 처리합니다
2. 모든 단계를 올바르게 연결합니다

### 4.4. 모든 테스트 실행

nf-test에는 한 가지 더 기능이 있습니다. 모든 테스트를 한 번에 실행할 수 있습니다! nf-test가 모든 디렉토리에서 nf-test 파일
