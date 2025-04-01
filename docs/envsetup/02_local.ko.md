# 로컬 설치

어떠한 이유로든 GitHub Codespaces 세션을 **사용할 수 없는 경우**, 모든 환경을 로컬에 설치하여 사용할수 있습니다.

로컬 머신의 운영체제나 설정에 따라 요구 사항이 다를수 있습니다.

## 요구 사항

Nextflow는 POSIX 호환 시스템(Linux, macOS, Windows Subsystem for Linux 등)에서 사용 가능합니다.

**필수 구성 요소**

- Bash
- [Java 11 이상 (최대 21까지 지원)](https://www.oracle.com/technetwork/java/javase/downloads/index.html)
- [Git](https://git-scm.com/)
- [Docker](https://docs.docker.com/get-docker/)

**선택 구성 요소**

- [Singularity](https://github.com/sylabs/singularity) 2.5.x (이상)
- [Conda](https://conda.io/) 4.5 (이상)
- [Graphviz](http://www.graphviz.org/)
- [AWS CLI](https://aws.amazon.com/cli/)
- 사전에 설정된 AWS Batch 작업 실행 환경

## Nextflow 다운로드

터미널에서 다음 명령어를 실행하세요:

```bash
wget -qO- https://get.nextflow.io | bash
```

또는 `curl` 명령어를 사용할수 있습니다:

```bash
curl -s https://get.nextflow.io | bash
```

다음으로, 다운로드된 실행 파일에 실행 권한을 부여하세요:

```bash
chmod +x nextflow
```

마지막으로, `nextflow` 실행 파일이 `$PATH`에 포함되어 있는지 확인하세요. 실행 파일은 `/usr/local/bin`, `/bin/` 등의 경로에 위치할 수 있습니다.

## Docker

Docker Desktop이 컴퓨터에서 실행되고 있는지 확인하세요. 설치하지 않았다면 [이 링크](https://docs.docker.com/get-docker/)를 통해 다운로드할 수 있습니다.

## 교육 자료

[이 링크](https://training.nextflow.io/)에서 교육 자료를 확인할 수 있습니다.

자료를 다운로드하려면 다음 명령어를 실행하세요요:

```bash
git clone https://github.com/nextflow-io/training.git
```

이후, 터미널에서 `cd` 명령어를 사용해 해당 디렉토리로 이동하세요. 기본적으로 `hello-nextflow` 디렉토리로 이동하면 됩니다.

## 설치 확인

다음 명령어를 실행하여 `nextflow`가 올바르게 설치되었는지 확인합니다:

```bash
nextflow info
```

해당 명령어를 실행하면 현재 버전, 시스템 정보, 런타임 정보가 출력됩니다.

!!! question "연습 문제"

    환경이 제대로 작동하는지 확인하려면 다음 명령어를 실행해 보세요:

    ```bash
    nextflow info
    ```

    실행 결과에는 다음과 같은 Nextflow 버전 및 런타임 정보가 출력됩니다 (실제 버전은 다를 수 있습니다):

    ```console
    Version: 23.10.1 build 5891
    Created: 12-01-2024 22:01 UTC
    System: Linux 6.1.75-060175-generic
    Runtime: Groovy 3.0.19 on OpenJDK 64-Bit Server VM 11.0.1-internal+0-adhoc..src
    Encoding: UTF-8 (UTF-8)
    ```
