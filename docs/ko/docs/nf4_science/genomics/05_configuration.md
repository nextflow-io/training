# 파트 3: 리소스 프로파일링 및 최적화

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이것은 임시 자리표시자입니다

!!!note "참고"

    이 교육 모듈은 현재 재개발 중입니다.

---

TODO

### 1.1. 리소스 사용률 보고서를 생성하기 위해 워크플로우 실행하기

Nextflow가 보고서를 자동으로 생성하도록 하려면, 명령줄에 `-with-report <파일명>.html`을 추가하기만 하면 됩니다.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-1.html
```

보고서는 html 파일이며, 다운로드하여 브라우저에서 열 수 있습니다. 또한 왼쪽의 파일 탐색기에서 마우스 오른쪽 버튼을 클릭하고 `Show preview`를 클릭하여 VS Code에서 볼 수도 있습니다.

시간을 내어 보고서를 살펴보고 리소스 조정 기회를 찾을 수 있는지 확인하십시오.
할당된 리소스의 백분율로 사용률 결과를 보여주는 탭을 클릭하십시오.
사용 가능한 모든 기능을 설명하는 [문서](https://www.nextflow.io/docs/latest/reports.html)가 있습니다.

<!-- TODO: insert images -->

한 가지 관찰 사항은 `GATK_JOINTGENOTYPING`이 CPU를 매우 많이 사용하는 것으로 보이며, 이는 복잡한 계산을 많이 수행하기 때문에 당연합니다.
따라서 이를 증가시켜 실행 시간을 단축할 수 있는지 시도해 볼 수 있습니다.

그러나 메모리 할당은 과도하게 설정한 것으로 보입니다. 모든 프로세스가 할당된 양의 일부만 사용하고 있습니다.
이를 줄여서 리소스를 절약해야 합니다.

### 1.2. 특정 프로세스에 대한 리소스 할당 조정하기

`withName` 프로세스 선택자를 사용하여 특정 프로세스에 대한 리소스 할당을 지정할 수 있습니다.
프로세스 블록에 단독으로 있을 때의 구문은 다음과 같습니다:

```groovy title="구문"
process {
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

이것을 `nextflow.config` 파일의 기존 프로세스 블록에 추가해 보겠습니다.

```groovy title="nextflow.config" linenums="11"
process {
    // 모든 프로세스에 대한 기본값
    cpus = 2
    memory = 2.GB
    // 특정 프로세스에 대한 할당
    withName: 'GATK_JOINTGENOTYPING' {
        cpus = 4
    }
}
```

이렇게 지정하면 기본 설정이 `GATK_JOINTGENOTYPING` 프로세스를 **제외한** 모든 프로세스에 적용됩니다. `GATK_JOINTGENOTYPING` 프로세스는 훨씬 더 많은 CPU를 받는 특별한 경우입니다.
이것이 효과가 있기를 바랍니다.

### 1.3. 수정된 구성으로 다시 실행하기

수정된 구성과 보고 플래그를 켠 상태로 워크플로우를 다시 실행해 보겠습니다. 하지만 보고서를 구별할 수 있도록 다른 이름을 지정합니다.

```bash
nextflow run main.nf -profile my_laptop -with-report report-config-2.html
```

다시 한번, 실행 시간에서 큰 차이를 느끼지 못할 것입니다. 이는 워크로드가 매우 작고 도구들이 '실제' 작업을 수행하는 것보다 보조 작업에 더 많은 시간을 소비하기 때문입니다.

그러나 두 번째 보고서는 리소스 사용률이 이제 더 균형 잡혔음을 보여줍니다.

<!-- **TODO: screenshots?** -->

보시다시피, 이 접근 방식은 프로세스마다 다른 리소스 요구 사항이 있을 때 유용합니다. 추측이 아닌 실제 데이터를 기반으로 각 프로세스에 설정한 리소스 할당을 적정 크기로 조정할 수 있게 해줍니다.

!!!note "참고"

    이것은 리소스 사용을 최적화하기 위해 할 수 있는 작은 맛보기일 뿐입니다.
    Nextflow 자체에는 리소스 제한으로 인해 실패한 작업을 재시도하는 정말 멋진 [동적 재시도 로직](https://www.nextflow.io/docs/latest/process.html#dynamic-task-resources)이 내장되어 있습니다.
    또한 Seqera Platform은 리소스 할당을 자동으로 최적화하기 위한 AI 기반 도구도 제공합니다.

    이 두 접근 방식은 이 교육 과정의 향후 부분에서 다룰 예정입니다.

그렇긴 하지만, 사용하는 컴퓨팅 실행자와 컴퓨팅 인프라에 따라 할당할 수 있는(또는 할당해야 하는) 것에 몇 가지 제약이 있을 수 있습니다. 예를 들어, 클러스터에서 다른 곳에서 실행할 때는 적용되지 않는 특정 제한 내에서 유지하도록 요구할 수 있습니다.
