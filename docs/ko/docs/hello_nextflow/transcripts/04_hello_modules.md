# Part 4: Hello Modules - 대본

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "중요 참고 사항"

    이 페이지는 대본만 표시합니다. 전체 단계별 지침은 [교육 자료](../04_hello_modules.md)로 돌아가십시오.

    대본에 표시된 섹션 번호는 참고용으로만 제공되며 자료의 모든 섹션 번호를 포함하지 않을 수 있습니다.

## 환영합니다

안녕하세요, Hello Nextflow 교육 과정의 Part 4에 오신 것을 환영합니다.

이 장은 Hello Modules라고 하며, Nextflow 코드를 모듈화하는 방법에 대해 이야기하겠습니다. 우리가 할 일은 하나의 workflow 스크립트를 가져와서 별도의 파일로 분리하는 것입니다.

이렇게 하면 workflow가 커질수록 코드를 더 쉽게 탐색하고 유지 관리할 수 있으며, pipeline 간에 모듈을 공유할 수 있게 되어 동일한 도구를 사용하는 여러 pipeline이 있는 경우 해당 process를 한 번만 작성하면 됩니다.

이에 대한 전형적인 예가 nf-core modules 저장소로, 수천 개의 다양한 도구가 바로 사용 가능한 모듈로 제공되어 workflow에 설치하고 사용할 수 있습니다.

Nextflow는 sub workflow와도 작동할 수 있는데, 이는 모듈과 비슷하지만 여러 process를 가지고 있습니다. 이것은 이 교육의 범위를 벗어나지만 기본적으로 같은 방식으로 작동합니다.

좋습니다. 살펴보겠습니다.

평소와 같이 training.nextflow.io로 이동하여 시작하십시오.

사이드바에서 "Hello Nextflow"로 가서 part 4: "Hello Modules"를 진행합니다.

이제 GitHub Code Spaces 환경으로 이동하여 "hello-modules" 파일을 살펴보겠습니다.

이전과 마찬가지로 이전 장의 끝 지점에서 시작하므로 이 스크립트가 익숙할 것입니다. 세 개의 process가 있습니다: say hello, convert to upper, collect greetings가 있고, 간단한 workflow에서 이 세 가지 명령을 실행하고 마지막에 메시지를 emit합니다. greeting과 batch라는 두 개의 매개변수가 있으며, 마지막에 수집된 출력 파일에 사용되는 이름을 지정합니다.

## 0. 워밍업: hello-modules.nf 실행

nextflow run hello, modules를 실행하여 이 workflow가 예상대로 작동하는지 확인할 수 있습니다.

좋습니다. 각 process에서 세 개의 작업, 하나의 수집 작업이 실행되었으며, 이 batch에 세 개의 greetings가 있다고 알려줍니다. results로 이동하면 수집된 test, batch 출력을 포함하여 다양한 출력 파일이 있습니다.

## 1. 모듈을 저장할 디렉토리 생성

좋습니다. 모듈화를 해봅시다.

일반적으로 깔끔하게 정리하기 위해 pipeline 저장소의 하위 폴더에 모듈을 넣는 것이 좋습니다. 원하는 이름으로 부를 수 있지만 관례상 보통 modules라고 부릅니다.

그럼 터미널로 가서 make the modules를 실행하겠습니다. VS Code의 사이드바에 나타나는 것을 볼 수 있습니다.

## 2. sayHello() 모듈 생성

그런 다음 첫 번째 모듈을 위한 새 파일을 만들겠습니다. "touch"나 "code"를 사용할 수 있고 사이드바에서 할 수도 있습니다. 정말 중요하지 않습니다. 그래서 code modules를 실행하고 process 이름을 따서 이름을 지정하겠습니다. 그래서 sayHello.nf라고 하겠습니다. NF는 Nextflow 파일의 전통적인 파일 확장자입니다.

여기서 저장하면 새 모듈 파일이 나타나는 것을 볼 수 있습니다.

## 2.2. sayHello process 코드를 모듈 파일로 이동

좋습니다. 다음으로 workflow에서 모듈 코드를 가져오겠습니다. 또한 여기 hash bang도 가져와서 먼저 복사하여 명확하게 Nextflow 파일임을 표시하겠습니다. 그런 다음 이 process를 가져와서 잘라내겠습니다. 메인 workflow 스크립트에서 제거하고 이 새 모듈에 붙여넣겠습니다.

이것이 이 모듈 파일에 포함될 모든 내용입니다. 단일 process만 있고, workflow도 없고, 로직도 없고, process만 있습니다.

이제 이 파일을 닫을 수 있습니다.

## 2.3. workflow 블록 앞에 import 선언 추가

이제 workflow에 첫 번째 process가 없으므로 import하여 다시 가져와야 합니다. 이에 대한 구문은 다른 프로그래밍 언어와 매우 유사하므로 익숙하게 느껴질 수 있습니다. include 중괄호를 사용하고, process 이름인 say hello를 넣은 다음 파일 경로 modules, say hello, nf에서 가져옵니다. 훌륭합니다.

여기 몇 가지 트릭이 있습니다. VS Code 확장 프로그램은 이것에 대해 영리합니다. 이 파일 경로를 인식하고 마우스를 올려놓고 링크 따라가기를 할 수 있습니다. 또는 Mac을 사용 중이라면 option을 누르고 클릭하면 이 파일이 열립니다. 그래서 빠르게 이동할 수 있습니다.

이 process 이름은 이제 아래 workflow에서 사용되고 있으며, 여기서도 같은 작업을 할 수 있습니다. 해당 process에 대한 약간의 정보를 보여주고, 다시 option을 누르고 클릭하면 편집기에서 열립니다.

따라서 다양한 process를 위한 많은 파일이 있을 때 VS Code에서 코드베이스를 빠르게 탐색할 수 있는 정말 빠른 방법입니다.

좋습니다. 기본적으로 이것이 이 장의 전부입니다. 이제 다른 process에 대해서도 같은 작업을 반복하면 됩니다.

## 3. convertToUpper() process 모듈화

그럼 여기에 새 파일을 만들어 봅시다. Convert to upper nf라고 부르겠습니다. 다시 hash bang을 복사합니다. 그런 다음 process를 잘라냅니다.

거기에 process 이름을 복사하고, 새 process 이름으로 새 include 문을 포함합니다.

## 4. collectGreetings() process 모듈화

그리고 세 번째 process에 대해서도 같은 작업을 합니다. 새 파일, connect. Greetings,

hash bang을 작성합니다. process를 잘라내고, process를 붙여넣고, 새 include 문을 작성합니다.

이제 여기에 잘못된 include 소스라고 하는 오류 밑줄이 있는 것을 볼 수 있습니다. 그리고 이것은 제가 너무 빨리 이동하여 실제로 발생한 진짜 오류입니다. 자세히 보면 T를 빠뜨리고 convert to upper라고 한 것을 볼 수 있습니다.

따라서 VS Code는 매우 유용하게도 제가 거기서 실수했다고 알려주었습니다. 파일 이름을 수정하면 오류가 사라집니다. 이것은 VS Code 내의 오류 검사가 Nextflow 코드를 작성하는 데 왜 그렇게 유용한지를 보여주는 좋은 예입니다. 그렇지 않았다면 이것을 발견하지 못했을 것이고 workflow를 실행하려고 할 때 훨씬 나중에야 알게 되었을 것입니다.

이제 메인 pipeline 스크립트가 훨씬 더 간단해 보입니다. process가 하나도 없고, 세 개의 include 문과 workflow만 있습니다. workflow의 로직을 변경하지 않았습니다. process 코드를 변경하지 않았으므로 정확히 같은 방식으로 작동해야 합니다.

## 4.4. 이전과 같은 작업을 수행하는지 확인하기 위해 workflow 실행

확인해 봅시다. 터미널을 열고 이전과 정확히 같은 명령을 실행하겠습니다.

확실히 process인 say hello, convert to upper collect greetings가 실행되었고, 다시 세 개의 greetings를 주었습니다.

그래서 코드를 이동했지만, workflow가 실행되는 방식에 대해서는 아무것도 변경하지 않았으며 완전히 변경되지 않았습니다. 유일한 차이점은 이제 더 깨끗한 코드를 가지고 있고, 유지 관리가 더 쉽고, 다른 사람과 공유하기가 더 쉽다는 것입니다.

그리고 그게 전부입니다. 짧은 장이었습니다. 간단한 개념이지만 매우 강력하고 더 복잡한 Nextflow workflow를 작성하는 방법의 핵심입니다. 따라서 이를 이해하고 사용하는 습관을 들이는 것이 중요합니다.

다음 장에서는 약간 속도를 바꿔서 Nextflow 코드를 작성하는 구문에 대해 생각하는 것을 멈추고, process 자체에서 소프트웨어를 사용하는 방법에 대해 조금 생각해 보겠습니다. part 5의 Hello Containers에 참여하십시오.

[다음 동영상 대본 :octicons-arrow-right-24:](05_hello_containers.md)
