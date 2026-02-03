---
title: 환경 옵션
description: Nextflow 교육을 위한 환경 설정 옵션
hide:
  - toc
  - footer
---

# 환경 옵션

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

학습자가 소프트웨어 관리에 시간과 노력을 들이지 않고 Nextflow 학습에 집중할 수 있도록 일관되고 철저히 테스트된 환경을 제공하는 것을 목표로 합니다.
이를 위해 모든 과정을 진행하는 데 필요한 모든 소프트웨어, 코드 파일 및 예제 데이터가 포함된 컨테이너화된 환경을 개발했습니다.

이 컨테이너화된 환경은 Github Codespaces에서 즉시 실행하거나 Devcontainers 확장 기능이 있는 VS Code에서 로컬로 실행할 수 있습니다.

<div class="grid cards" markdown>

- :material-cloud-outline:{ .lg .middle } **Github Codespaces**

  ***

  GitHub Codespaces는 클라우드의 가상 머신에 의해 지원되는, 모든 도구와 데이터가 포함된 사전 구축된 교육 환경을 제공할 수 있는 웹 기반 서비스입니다. Github 계정이 있는 누구나 무료로 이용할 수 있습니다.

  [Github Codespaces 사용하기:material-arrow-right:](01_setup.md){ .md-button .md-button--primary .mt-1 }

- :material-laptop:{ .lg .middle } **로컬 Devcontainers**

  ***

  Devcontainers가 있는 VS Code는 모든 교육 도구가 사전 구성된 로컬 실행 컨테이너화된 개발 환경을 제공합니다. Codespaces와 동일한 사전 구축된 환경을 제공하지만 전적으로 로컬 하드웨어에서 실행됩니다.

  [로컬에서 Devcontainers 사용하기 :material-arrow-right:](03_devcontainer.md){ .md-button .md-button--primary .mt-1 }

</div>

## 수동 설치 안내

위의 옵션이 적합하지 않은 경우 소프트웨어 의존성을 수동으로 설치하고 교육 저장소를 복제하여 자체 로컬 시스템에서 이 환경을 복제할 수 있습니다.

[수동 설치 :material-arrow-right:](02_local.md){ .md-button .md-button--primary .mt-1 }

---

!!! info "Gitpod 지원 중단"

    Nextflow Training은 2025년 2월까지 [Gitpod](https://gitpod.io)를 사용했습니다.
    그러나 Gitpod 제작자들은 [Gitpod Flex](https://www.gitpod.io/blog/introducing-gitpod-flex) 시스템을 위해 무료 기능을 종료하기로 결정했습니다.
    그 이유로 사전 설정 없이 원클릭 개발자 환경을 제공하는 GitHub Codespaces로 전환했습니다.

    Gitpod에 가입한 시기와 서비스 종료 시기에 따라 이전 클라우드 IDE에서 교육을 시작할 수 있지만 앞으로 안정적인 액세스를 보장할 수 없습니다:
    [Gitpod에서 열기](https://gitpod.io/#https://github.com/nextflow-io/training).
