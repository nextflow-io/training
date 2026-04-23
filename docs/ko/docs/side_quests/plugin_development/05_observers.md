# 파트 5: Trace Observer

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Trace observer를 사용하면 작업 완료, 파일 게시, 파이프라인 종료 등의 워크플로우 이벤트에 플러그인이 응답할 수 있습니다.
이를 통해 사용자 정의 보고서, Slack 알림, 메트릭 수집, 외부 모니터링 시스템과의 연동 등 다양한 활용이 가능합니다.
이 섹션에서는 완료된 작업 수를 세고 요약을 출력하는 observer를 구현합니다.

!!! tip "여기서부터 시작하시나요?"

    파트 5부터 시작하는 경우, 파트 4의 솔루션을 시작점으로 복사하세요:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. 기존 trace observer 이해하기

파이프라인을 실행했을 때 출력된 "Pipeline is starting!" 메시지는 플러그인의 `GreetingObserver` 클래스에서 나온 것입니다.

observer 코드를 확인합니다:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingObserver.groovy
```

```groovy title="Output" hl_lines="30 32-34 37-39"
/*
 * Copyright 2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package training.plugin

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.trace.TraceObserver

/**
 * Implements an observer that allows implementing custom
 * logic on nextflow execution events.
 */
@Slf4j
@CompileStatic
class GreetingObserver implements TraceObserver {    // (1)!

    @Override
    void onFlowCreate(Session session) {            // (2)!
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {                         // (3)!
        println "Pipeline complete! 👋"
    }
}
```

1. 워크플로우 생명주기 이벤트에 연결하기 위한 인터페이스
2. 워크플로우가 시작될 때 호출되며, 설정에 접근하기 위한 세션을 받습니다
3. 워크플로우가 성공적으로 완료될 때 호출됩니다

여기서 두 가지를 주목해야 합니다:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver`는 Nextflow가 정의한 인터페이스입니다. 클래스가 이 인터페이스를 구현하면, Nextflow가 이벤트 발생 시 해당 메서드를 호출할 수 있습니다.
2. **`@Override`**: `TraceObserver` 인터페이스는 `onFlowCreate`, `onFlowComplete` 등의 메서드를 정의합니다. 이 이름으로 메서드를 작성하고 `@Override` 어노테이션을 추가하면, Nextflow가 적절한 시점에 해당 메서드를 호출합니다. 오버라이드하지 않은 메서드는 무시됩니다.

현재 시점에서 연결할 수 있는 생명주기 이벤트의 전체 목록은 다음과 같습니다:

| 메서드              | 호출 시점          |
| ------------------- | ------------------ |
| `onFlowCreate`      | 워크플로우 시작    |
| `onFlowComplete`    | 워크플로우 완료    |
| `onProcessStart`    | 작업 실행 시작     |
| `onProcessComplete` | 작업 완료          |
| `onProcessCached`   | 캐시된 작업 재사용 |
| `onFilePublish`     | 파일 게시          |

전체 목록은 Nextflow 소스의 [TraceObserver 인터페이스](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy)를 참조하세요.

---

## 2. 작업 카운터 observer 추가하기

목표는 완료된 작업 수를 세고 마지막에 요약을 출력하는 observer를 구현하는 것입니다.
플러그인에 새 observer를 추가하려면 두 가지가 필요합니다: observer 클래스 작성과 Nextflow가 로드할 수 있도록 팩토리에 등록하는 것입니다.

### 2.1. 최소한의 observer 생성하기

새 파일을 생성합니다:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

작업이 완료될 때 메시지를 출력하는 가장 단순한 observer부터 시작합니다:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler       // (1)!
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer that responds to task completion
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {  // (2)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {  // (3)!
        println "✓ Task completed!"
    }
}
```

1. 필요한 클래스 임포트: `TraceObserver`, `TaskHandler`, `TraceRecord`
2. `TraceObserver`를 구현하는 클래스 생성
3. 작업이 완료될 때 실행할 코드를 위해 `onProcessComplete` 오버라이드

최소 요구사항은 다음과 같습니다:

- 필요한 클래스 임포트 (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- `TraceObserver`를 구현하는 클래스 생성
- 작업 완료 시 동작을 위해 `onProcessComplete` 오버라이드

### 2.2. Observer 등록하기

`GreetingFactory`가 observer를 생성합니다.
내용을 확인합니다:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingFactory.groovy
```

```groovy title="Output" hl_lines="25 27-29"
/*
 * Copyright 2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

@CompileStatic
class GreetingFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }

}
```

새 observer를 추가하기 위해 `GreetingFactory.groovy`를 수정합니다:

=== "후"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "전"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Groovy 리스트 문법"

    Java 스타일의 `List.<TraceObserver>of(...)`를 Groovy의 더 간단한 리스트 리터럴 `[...]`로 교체했습니다.
    둘 다 `Collection`을 반환하지만, 여러 항목을 추가할 때는 Groovy 문법이 더 읽기 쉽습니다.

### 2.3. 빌드, 설치, 테스트

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

!!! tip "왜 `-ansi-log false`를 사용하나요?"

    기본적으로 Nextflow의 ANSI 진행 표시는 이전 줄을 덮어써서 진행 상황을 깔끔하게 업데이트된 형태로 보여줍니다.
    이 경우 중간 메시지는 보이지 않고 *최종* 작업 수만 표시됩니다.

    `-ansi-log false`를 사용하면 이 동작이 비활성화되어 모든 출력이 순차적으로 표시됩니다. 실행 중에 메시지를 출력하는 observer를 테스트할 때 필수적입니다.

"✓ Task completed!"가 다섯 번(작업당 한 번) 출력되며, 기존 파이프라인 출력과 섞여서 표시됩니다:

```console title="Output (partial)"
...
[9b/df7630] Submitted process > SAY_HELLO (4)
Decorated: *** Hello ***
✓ Task completed!
✓ Task completed!
Decorated: *** Holà ***
✓ Task completed!
...
Pipeline complete! 👋
```

Observer가 정상적으로 동작하고 있습니다.
작업이 완료될 때마다 Nextflow가 `onProcessComplete`를 호출하고, 구현된 코드가 메시지를 출력합니다.

??? exercise "메시지 변경해보기"

    `onProcessComplete`의 메시지를 원하는 내용으로 변경하고, 다시 빌드하여 실행해 보세요.
    이를 통해 observer에 대한 전체 편집-빌드-실행 사이클이 정상적으로 동작하는지 확인할 수 있습니다.

### 2.4. 카운팅 로직 추가하기

최소한의 observer는 훅이 동작함을 확인해 주지만, 아무것도 추적하지 않습니다.

클래스는 객체의 생명주기 동안 유지되는 변수(필드 또는 인스턴스 변수라고 함)를 가질 수 있습니다.
즉, observer는 파이프라인 실행 중 여러 이벤트에 걸쳐 상태를 누적할 수 있습니다.

다음 버전에서는 0에서 시작하는 카운터 변수(`taskCount`)를 추가합니다.
작업이 완료될 때마다 카운터가 1씩 증가합니다.
전체 워크플로우가 완료되면 observer가 최종 합계를 출력합니다.

강조 표시된 변경 사항으로 `TaskCounterObserver.groovy`를 업데이트합니다:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="14 18-19 22-24"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer that counts completed tasks
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0                // (1)!

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++                          // (2)!
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFlowComplete() {                  // (3)!
        println "📈 Final task count: ${taskCount}"
    }
}
```

1. `taskCount`는 observer 객체에 속하는 변수입니다. 메서드 호출 간에 값을 유지하므로 전체 워크플로우 실행에 걸쳐 카운트를 누적할 수 있습니다. `private`은 이 클래스만 접근할 수 있음을 의미합니다.
2. `taskCount++`는 카운터에 1을 더합니다. 이 줄은 작업이 완료될 때마다 실행되므로 워크플로우가 진행됨에 따라 카운트가 증가합니다.
3. `onFlowComplete`는 두 번째 생명주기 훅입니다. 워크플로우가 완료될 때 한 번 실행되므로 요약을 출력하기에 적합합니다.

정리하면:

- `taskCount`는 메서드 호출 간에 유지되며 전체 실행에 걸쳐 카운트를 누적합니다
- `onProcessComplete`는 작업이 완료될 때마다 카운터를 증가시키고 누적 합계를 출력합니다
- `onFlowComplete`는 마지막에 한 번 실행되어 최종 카운트를 출력합니다

다시 빌드하고 테스트합니다:

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "출력"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `greet.nf` [pensive_engelbart] DSL2 - revision: 85fefd90d0
    Pipeline is starting! 🚀
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    [be/bd8e72] Submitted process > SAY_HELLO (2)
    [5b/d24c2b] Submitted process > SAY_HELLO (1)
    [14/1f9dbe] Submitted process > SAY_HELLO (3)
    Decorated: *** Bonjour ***
    Decorated: *** Hello ***
    [85/a6b3ad] Submitted process > SAY_HELLO (4)
    📊 Tasks completed so far: 1
    📊 Tasks completed so far: 2
    Decorated: *** Holà ***
    📊 Tasks completed so far: 3
    Decorated: *** Ciao ***
    [3c/be6686] Submitted process > SAY_HELLO (5)
    📊 Tasks completed so far: 4
    Decorated: *** Hallo ***
    📊 Tasks completed so far: 5
    Pipeline complete! 👋
    📈 Final task count: 5
    ```

    카운터 메시지는 observer가 작업 완료 시 실행되기 때문에 작업 제출 메시지와 섞여서 출력됩니다.

---

## 3. 게시된 파일 추적하기

Observer는 파일이 게시될 때도 응답할 수 있습니다.
`onFilePublish` 메서드는 대상 경로와 소스 경로를 받으며, 이를 사용하여 게시된 출력을 로깅, 검증, 또는 처리할 수 있습니다.

### 3.1. 게시 디렉토리 추가하기

먼저, `SAY_HELLO` 프로세스가 출력 파일을 게시하도록 `greet.nf`를 업데이트합니다:

=== "후"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // 커스텀 플러그인 함수를 사용하여 인사말을 꾸밉니다
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "전"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // 커스텀 플러그인 함수를 사용하여 인사말을 꾸밉니다
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. onFilePublish 메서드 추가하기

`TaskCounterObserver.groovy`에 `onFilePublish` 메서드와 필요한 임포트를 추가합니다:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1" hl_lines="5 23-26"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer that counts completed tasks
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private int taskCount = 0

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "📁 Published: ${destination.fileName}"
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount}"
    }
}
```

### 3.3. 빌드 및 테스트

```bash
cd nf-greeting && make install && cd ..
nextflow run greet.nf -ansi-log false
```

각 출력 파일에 대해 "Published:" 메시지가 작업 카운터 출력과 함께 표시됩니다:

```console title="Output (partial)"
...
📊 Tasks completed so far: 1
📁 Published: greeting.txt
📊 Tasks completed so far: 2
📁 Published: greeting.txt
...
📈 Final task count: 5
Pipeline complete! 👋
```

`onFilePublish` 메서드는 Nextflow가 `results` 디렉토리에 파일을 게시할 때마다 실행됩니다.
이 패턴은 감사 로그 구축, 다운스트림 작업 트리거, 또는 출력이 생성될 때 검증하는 데 유용합니다.

---

## 핵심 정리

이 섹션에서 학습한 내용:

- Trace observer는 `onFlowCreate`, `onProcessComplete`, `onFilePublish`, `onFlowComplete` 등의 워크플로우 생명주기 이벤트에 연결됩니다
- `TraceObserver`를 구현하고 팩토리에 등록하여 observer를 생성합니다
- Observer는 인스턴스 변수를 통해 이벤트 간 상태를 누적할 수 있습니다
- Observer는 사용자 정의 로깅, 메트릭 수집, 알림, 보고에 유용합니다

---

## 다음 단계

작업 카운터는 동작하지만 항상 활성화되어 있습니다.
실제 플러그인에서는 플러그인 소스 코드를 수정하지 않고 `nextflow.config`에서 기능을 활성화하거나 비활성화하거나 동작을 조정할 수 있어야 합니다.
다음 섹션에서는 observer를 설정 가능하게 만드는 방법과 완성된 플러그인을 다른 사람과 공유하는 방법을 다룹니다.

[파트 6으로 계속 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
