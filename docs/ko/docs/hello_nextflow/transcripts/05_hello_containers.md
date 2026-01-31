# Part 5: Hello Containers - 스크립트

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "중요 사항"

    이 페이지는 스크립트만 표시합니다. 전체 단계별 지침은 [교육 자료](../05_hello_containers.md)로 돌아가십시오.

    스크립트에 표시된 섹션 번호는 참고용으로만 제공되며 자료의 모든 섹션 번호를 포함하지 않을 수 있습니다.

## 환영합니다

안녕하세요, Hello Nextflow 교육 과정의 Part Five에 오신 것을 환영합니다.

이 장의 제목은 Hello Containers입니다. Nextflow가 Docker 및 Singularity와 같은 도구와 통합되어 소프트웨어 컨테이너를 사용하여 파이프라인 사용자에게 소프트웨어를 제공하는 방법에 대해 이야기하겠습니다.

이는 사람들이 파이프라인을 실행할 때 모든 다양한 도구를 직접 설치할 필요가 없다는 것을 의미합니다. Nextflow가 대신 해줄 것입니다.

컨테이너는 매우 강력한 기술이며 재현성과 사용 편의성에 매우 중요합니다. 먼저 컨테이너 자체에 대한 간략한 소개로 시작하여 몇 가지 docker 명령을 수동으로 실행한 다음, 동일한 컨테이너를 Nextflow 파이프라인에 넣을 것입니다.

좋습니다. 시작하겠습니다.

이전과 마찬가지로 교육 자료를 로드하는 것으로 시작하겠습니다. training.nextflow.io로 이동하십시오. Hello Nextflow, Chapter Five, Hello Containers.

Codespaces 환경으로 이동하여 왼쪽에서 hello containers dot nf를 볼 수 있습니다.

이전과 마찬가지로 이것은 이전 4장에서 완료한 것과 동일한 스크립트이므로 익숙할 것입니다.

입력 파일과 배치 이름을 지정하는 명령줄 매개변수가 있습니다. 세 개의 모듈을 포함하고 있으며, 세 개의 프로세스를 실행하는 워크플로우가 있습니다.

## 0. 워밍업: hello-containers.nf 실행

이 워크플로우를 다시 실행하고 예상한 출력을 생성하는지 다시 확인하십시오. 지금은 실제로 닫고 터미널로 들어가겠습니다.

## 1. 컨테이너를 '수동으로' 사용

이 장을 시작하기 위해 컨테이너 기술에 대해 조금 복습하겠습니다. docker나 singularity 또는 다른 컨테이너 기술에 매우 익숙하다면 이것을 복습으로 간주하거나 완전히 건너뛰어도 됩니다.

Nextflow는 다양한 유형의 컨테이너 기술을 지원합니다. 여기에는 Docker, Singularity, Podman, Shifter, Charliecloud 등이 포함됩니다.

이 교육에서는 Docker에 중점을 둘 것입니다. 이는 code spaces에 사전 설치되어 있으며 특히 자신의 컴퓨터나 노트북에서 개발하는 경우 가장 인기 있는 컨테이너 기술 중 하나입니다.

공유 HPC의 학술 환경에서 작업하는 경우 Singularity를 사용할 수 있고 Docker는 사용할 수 없을 수 있습니다. 괜찮습니다. 모든 개념은 정확히 동일합니다. 몇 가지 수동 명령은 다르지만 Docker를 이해하면 singularity도 이해할 수 있습니다.

실제로 Singularity도 Code Spaces 환경에 설치되어 있습니다. 원한다면 Docker 대신 Singularity를 사용하여 동일한 작업을 시도할 수 있습니다.

좋습니다, 컨테이너 기술이란 무엇입니까? Docker의 아이디어는 원격 소스에서 이미지를 가져올 수 있다는 것입니다. 로컬 머신으로 가져온 다음 해당 이미지를 기반으로 컨테이너를 생성합니다.

이 실행 중인 컨테이너는 컴퓨터에서 실행되는 가상 머신과 비슷합니다. 환경에서 격리되어 있으며 운영 체제 및 사용 가능한 소프트웨어 세트와 함께 사전 패키지되어 있습니다.

## 1.1. 컨테이너 이미지 가져오기

기존 이미지를 가져오는 데 필요한 구문은 "docker pull"입니다. 터미널에 입력하겠지만 이제 작업할 이미지가 필요합니다.

이미지를 직접 빌드할 수 있습니다. Docker Hub 또는 quay.io와 같은 공개 레지스트리에서 찾을 수 있습니다. 그러나 이미지를 빠르게 얻는 정말 좋은 방법은 Seqera Containers를 사용하는 것입니다.

이것은 2024년에 구축한 무료 커뮤니티 서비스로 로그인 없이 사용할 수 있습니다.

seqera.io/containers로 이동하거나 여기 상단에서 containers를 클릭하면 검색 인터페이스가 나타나며 Conda 또는 Python 패키지 인덱스에서 사용 가능한 도구의 이름을 입력할 수 있습니다.

기본적으로 Bioconda 및 Conda Forge 채널을 검색하지만 원하는 경우 Conda 채널 앞에 접두사를 붙일 수 있습니다.

재미를 위해 cowpy를 사용하겠습니다. cowpy를 입력하겠습니다. Python Package Index 및 Conda Forge에서 결과를 제공합니다. 컨테이너에 추가하기 위해 클릭하겠습니다. 원한다면 여기에 여러 패키지를 추가할 수 있습니다. Docker를 선택하고 linux/amd64를 선택한 다음 Get Container를 클릭합니다.

이것은 아직 생성되지 않은 경우 요청 시 이미지를 빌드하고 복사할 수 있는 URL을 제공합니다.

관심이 있으시면 view Build Details를 클릭하면 사용된 conda 환경 파일과 빌드에 대한 전체 빌드 로그 및 보안 스캔 결과를 보여주는 페이지로 이동합니다.

code spaces로 돌아가면 이제 이 컨테이너 이름을 붙여넣고 엔터를 칠 수 있습니다.

Docker는 이제 이 컨테이너 이미지 내의 모든 다른 레이어를 다운로드하고 이제 이 이미지를 사용할 수 있다고 알려줍니다.

## Singularity 이미지 가져오기

singularity를 사용하는 경우 프로세스는 기본적으로 동일합니다. 이미지 패키지를 선택하고 cowpy를 선택합니다. 이제 Docker 대신 Singularity를 선택하고 Get Container를 클릭합니다. oras://를 사용하는 이미지 URL을 제공합니다. 또는 원하는 경우 해당 상자를 선택하여 https://를 사용할 수 있습니다. 해당 URL을 복사합니다. 이제 Code Spaces로 이동합니다. 실제로 이 공간에 Apptainer가 설치되어 있는데, 이는 Singularity와 동일하지만 서로 별칭이 지정되어 있습니다. 따라서 apptainer pull을 수행한 다음 cowpy sif라고 부르겠지만 원하는 대로 부를 수 있습니다. URL을 붙여넣습니다. 그러면 해당 이미지를 다운로드합니다.

ls -lh를 수행하면 cowpy.sif를 볼 수 있습니다.

Singularity는 Docker와 다릅니다. singularity는 모든 이미지를 플랫 파일에 저장하는 반면 Docker는 호스트 머신에 모든 레이어를 별도로 보관하는 레지스트리가 있으며 이 모든 것을 추적하기 위해 실행 중인 데몬이 있습니다.

## 1.2. 컨테이너를 사용하여 cowpy를 일회성 명령으로 실행

좋습니다, Docker로 돌아가겠습니다. 이제 docker run을 수행하여 생성한 이 이미지를 실행해 볼 수 있습니다.

dash dash rm을 수행하겠습니다. 이것은 이미지의 일회성 실행만 수행합니다. 그리고 이미지 URL을 붙여넣겠습니다. 그런 다음 마지막으로 실행하려는 명령으로 완료합니다.

생성한 이미지에는 cowpy가 설치되어 있으므로 cowpy를 시도해 보겠습니다.

됩니다. 명령이 실행되었습니다. cowpy가 로컬에 설치되어 있지 않습니다. 실행하려고 하면 존재하지 않는 것을 볼 수 있습니다. 그러나 이 명령에서는 Docker를 사용하여 실행했으며 이 출력을 올바르게 생성했습니다.

## 1.3. 컨테이너를 사용하여 cowpy를 대화식으로 실행

원한다면 더 나아가 컨테이너를 대화식으로 실행하고 내부를 둘러볼 수 있습니다. 다시 "docker run dash dash rm"을 수행합니다. 이제 대화형 터미널을 원한다고 Docker에 알리는 dash it를 수행하겠습니다. 이미지 URL을 다시 수행하고 이번에는 cowpy 대신 bin bash를 수행하겠습니다. 왜냐하면 실행하려는 명령이 bash이기 때문입니다.

이것은 실행 중인 컨테이너로 이동하며 이제 프롬프트가 변경된 것을 볼 수 있습니다.

LS slash를 수행하면 여기의 디렉토리가 다른 것을 볼 수 있습니다.

GitHub Code Spaces에서 실행 중인 오른쪽에 두 번째 터미널을 열고 LS slash를 수행하면 workspaces 및 temp와 같은 디렉토리가 있는 것을 볼 수 있지만 Docker의 여기는 다릅니다.

따라서 이 환경은 Docker 내에서 완전히 분리되어 있으며 호스트 환경에서 격리되어 있습니다. 이것은 좋은 일입니다. 왜냐하면 이 명령의 실행을 Docker 이미지로 격리하고 다른 호스트 시스템의 다른 사람들 사이에서 재현 가능하게 유지하기 때문입니다.

Docker 이미지 내에서 호스트 시스템의 데이터를 사용하려면 명시적으로 컨테이너에 마운트해야 합니다.

잠시 후에 그렇게 하겠습니다.

## 1.3.2. 원하는 도구 명령 실행

먼저 cowpy를 실행할 수 있는지 확인해 보겠습니다. 다시 명령줄에서 직접 사용할 수 있으며 더 복잡한 작업을 시작하고 인수를 전달할 수 있습니다. Hello containers 그리고 소 대신 tux 펭귄을 수행하겠습니다. 다른 무엇이 있는지 보겠습니다.

cheese를 수행하겠습니다. 훌륭합니다. Dragon과 Cow는 어떻습니까? 꽤 좋습니다.

## 1.3.3. 컨테이너 종료

좋습니다. 이 컨테이너에 데이터가 없기 때문에 더 이상 할 수 없습니다. 그래서 이 실행 중인 이미지를 종료하고 컨테이너에 일부 데이터를 마운트할 수 있는지 확인하겠습니다. control D를 수행하거나 exit를 입력하여 수행할 수 있습니다. 좋습니다, 이제 일반 GitHub code space로 돌아왔습니다.

## 1.3.4. 컨테이너에 데이터 마운트

Docker 컨테이너에 일부 데이터를 마운트하려면 dash V를 사용해야 합니다. 이전 docker 명령을 가져와서 시작 부분으로 이동하여 dash v를 수행하겠습니다. 현재 로컬 작업 디렉토리에 대해 "."를 수행한 다음 콜론을 사용하여 호스트 디렉토리에서 마운트될 위치를 지정하고 slash data를 수행합니다. 따라서 이것은 이 특정 디렉토리를 slash data의 컨테이너에 마운트합니다.

이제 LS slash를 수행하면 data라는 새 디렉토리가 있는 것을 볼 수 있으며 LS data를 수행하면 사이드바에 있는 모든 파일을 볼 수 있습니다. 훌륭합니다.

## 1.3.5. 마운트된 데이터 사용

이제 Docker 이미지 내의 호스트 시스템에 있는 일부 파일을 사용할 수 있습니다. cat data greetings csv라고 말할 수 있습니다. 기억하시다시피 이것은 이전의 다양한 인사말이 포함된 CSV 파일이며 cowpy로 파이프할 수 있습니다. 훌륭합니다. 이제 진전이 있습니다.

좋습니다. Docker를 대화식으로 실행하는 것은 충분합니다. Docker가 무엇인지, 일회성 방식으로 명령을 실행하는 방법과 이미지를 대화식으로 사용하는 방법에 대해 대략적으로 느끼셨기를 바랍니다. singularity를 사용하는 경우 apptainer exec 또는 apptainer run 또는 singularity exec 또는 singularity run과 같은 작업을 수행한다는 점을 제외하고 명령은 모두 매우 유사합니다.

## 2. Nextflow에서 컨테이너 사용

다음으로 Nextflow 워크플로우로 돌아가서 Nextflow 파이프라인 내에서 이 기술을 사용하는 방법을 살펴보겠습니다.

터미널을 닫고 Hello Containers를 다시 열겠습니다.

## 2.1. cowpy 모듈 작성

cowpy 예제를 고수하기 위해 cowpy를 사용하는 워크플로우에 새 프로세스를 생성하겠습니다. 모듈로 이동하여 새 파일을 생성하고 cowpy nf라고 부르겠습니다. 이제 조금 속임수를 쓰고 교육 자료에서 코드를 복사하여 저장을 누르겠습니다. 그리고 살펴보겠습니다.

이것은 간단한 프로세스입니다. 이제 프로세스의 구성 요소가 어떻게 생겼는지 이해하기를 바랍니다. 다시 결과로 이동하는 publishDir이 있습니다. 입력 파일과 character라는 문자열의 두 가지 입력이 있습니다. cowpy 입력 파일의 출력이 있으며 잠시 전에 docker 이미지 내에서 수동으로 실행한 것과 정확히 동일하게 보이는 스크립트가 있습니다: 파일을 인쇄하는 cat, cowpy에 파이프하고 사용하려는 cowpy 문자 유형을 지정하고 여기에서 출력으로 전달하는 출력 파일로 출력합니다.

## 2.2. 워크플로우에 cowpy 추가

좋습니다, 워크플로우로 돌아가서 이 새 프로세스를 가져오겠습니다. 그래서 modules cowpy nf에서 cowpy. 원하는 문자를 지정할 수 있도록 새 매개변수를 생성하겠습니다. 기본적으로 Turkey라고 하겠습니다. 그런 다음 워크플로우 끝에서 이 새 프로세스를 호출하겠습니다.

cowpy. 그리고 Collect Greetings의 출력을 사용하겠습니다. 따라서 collect greetings out, out file입니다. 그런 다음 방금 만든 새 params인 두 번째 인수가 필요합니다. params dot character.

## 2.2.4. 워크플로우를 실행하여 작동하는지 확인

좋습니다, 새 프로세스가 작동하는지 확인하겠습니다. Nextflow run hello containers. 이것은 처음 세 개의 프로세스를 실행한 다음 끝에 cowpy를 실행하려고 시도해야 합니다.

오류가 발생했습니다. 여기서 말하는 것은 cowpy에 오류가 있었고 종료 상태 127이 있었으며 확실히 command sh cowpy command not found입니다.

cowpy에 사용할 수 있는 Docker 이미지가 있다고 Nextflow에 알리지 않았으므로 호스트 시스템에서 실행하려고 시도했으며 호스트 시스템에 cowpy가 설치되어 있지 않아 오류가 발생했습니다.

## 2.3. 컨테이너를 사용하여 실행

따라서 해야 할 일은 사용 가능한 컨테이너가 있다고 Nextflow에 알려야 하는 것입니다. cowpy 프로세스로 이동하여 프로세스 상단에 container라는 새 지시문을 추가하겠습니다.

그런 다음 이미지를 찾아 URL을 복사하고 문자열에 넣습니다.

이것만으로는 충분하지 않습니다. X Flow 파이프라인에는 소프트웨어를 지정하는 여러 가지 방법이 있을 수 있기 때문입니다. 예를 들어 conda conda-forge cowpy도 수행할 수 있습니다. 그리고 Nextflow는 사용하려는 기술을 알아야 합니다.

## 2.3.2. nextflow.config 파일을 통해 Docker 사용 활성화

따라서 Docker가 활성화된 상태로 실행하려면 다음 장에서 더 자세히 다룰 내용인 Nextflow config 파일을 사용하여 약간 앞서나가겠습니다. 이 디렉토리에 Nextflow Config라는 파일이 있으며 여기에 이미 docker.enabled False가 있습니다.

Docker를 활성화하기 위해 True로 변경한 다음 워크플로우를 다시 실행할 수 있습니다.

## 2.3.3. Docker가 활성화된 상태로 워크플로우 실행

Nextflow run hello containers nf 이번에는 cowpy가 성공적으로 실행되었습니다. Results를 살펴보겠습니다. cowpy collected test 그리고 우리의 Turkey가 있습니다. 훌륭합니다.

따라서 백그라운드에서 Nextflow는 해당 프로세스에 사용 가능한 컨테이너가 있다는 것을 알았습니다.

이미지를 가져와서 명령을 실행했습니다.

## 2.3.4. Nextflow가 컨테이너화된 작업을 시작한 방법 검사

궁금하시다면 work 디렉토리를 살펴봄으로써 실제로 수행한 작업을 정확히 볼 수 있습니다. code work를 수행한 다음 해시를 수행하고 command run을 수행합니다. 기억하시다시피 이것은 해당 작업에 대해 실제로 실행되는 파일입니다. NXF launch라는 함수를 찾을 수 있습니다. 그리고 여기에서 Nextflow가 사용한 정확한 docker 명령을 볼 수 있습니다. 이것은 이전에 터미널에서 수동으로 수행한 것과 매우 유사합니다. Docker run. 이 호스트 디렉토리를 컨테이너에 바인딩한 다음 컨테이너 URL을 지정합니다.

따라서 여기에는 마법이 없습니다. Nextflow가 파이프라인에서 컨테이너를 쉽게 지정할 수 있는 방식으로 자동으로 무거운 작업을 수행하고 있을 뿐이며, 워크플로우를 실행하는 다른 사람들이 쉽게 사용할 수 있습니다. 그리고 그 사람들은 더 이상 분석 파이프라인을 실행하기 위해 소프트웨어 관리에 대해 생각할 필요가 없습니다.

매우 매우 간단하고 매우 편리하며 재현성도 뛰어납니다. 전반적으로 좋습니다.

좋습니다, 수고하셨습니다. Chapter Five의 끝입니다. 다음 비디오에서 Hello Nextflow 교육의 마지막 부분인 part six에 참여하여 Nextflow 구성에 대해 더 자세히 이야기하겠습니다.

다음 비디오에서 뵙겠습니다.

[다음 비디오 스크립트 :octicons-arrow-right-24:](06_hello_config.md)
