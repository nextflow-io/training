# 파트 3: 사용자 정의 함수

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

이 섹션을 마치면 플러그인에 사용자 정의 함수를 추가하고, 로컬에서 빌드 및 설치하여 실제 워크플로우에서 실행할 수 있게 됩니다.

!!! tip "여기서부터 시작하시나요?"

    파트 3부터 시작하는 경우, 파트 2의 해결책을 시작점으로 복사하세요:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. 템플릿이 생성한 내용 확인

직접 함수를 작성하기 전에, 템플릿이 생성한 예제 함수를 살펴보고 패턴을 이해합니다.

플러그인 디렉토리로 이동합니다:

```bash
cd nf-greeting
```

템플릿은 플러그인 함수가 정의되는 `GreetingExtension.groovy` 파일을 생성했습니다.
시작점을 확인하기 위해 파일을 열어봅니다:

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

```groovy title="Output" hl_lines="29 40-43"
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
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

/**
 * Nextflow 스크립트에서 가져올 수 있는
 * 사용자 정의 함수를 구현합니다.
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * 지정된 대상에게 인사합니다.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. 확장이 기반으로 하는 클래스입니다. Nextflow가 함수를 인식하려면 이 클래스가 필요합니다.
2. 플러그인이 로드될 때 호출됩니다. 초기화에 사용합니다.
3. 이 메서드를 `include`를 통해 워크플로우에서 호출 가능하게 만듭니다.

템플릿에는 `sayHello` 함수 예제가 포함되어 있습니다.
`@Function` 어노테이션은 메서드를 Nextflow 워크플로우에서 호출 가능하게 만드는 역할을 합니다.
이 어노테이션이 없으면 해당 메서드는 플러그인 코드 내부에서만 존재합니다.

Groovy(및 Java)에서는 메서드가 반환하는 타입과 매개변수의 타입을 선언합니다.
예를 들어, `String reverseGreeting(String greeting)`은 `String` 매개변수를 받아 `String`을 반환하는 메서드를 선언합니다.
`void` 키워드는 위의 `sayHello`처럼 메서드가 아무것도 반환하지 않음을 의미합니다.
이는 타입을 명시적으로 선언할 필요가 없는 Python이나 R과는 다릅니다.

---

## 2. sayHello를 reverseGreeting으로 교체

템플릿의 `sayHello` 함수는 플레이스홀더입니다.
함수를 직접 작성하고, 빌드하고, 플러그인 함수를 사용하는 전체 과정을 확인하기 위해 이를 교체합니다.

`src/main/groovy/training/plugin/GreetingExtension.groovy`를 편집하여 `sayHello` 메서드를 교체합니다:

=== "후"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * 인사말 문자열을 역순으로 변환합니다
         */
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. 메서드를 Nextflow 워크플로우에서 호출 가능하게 만듭니다.
    2. String을 받아 String을 반환합니다.
    3. Groovy의 내장 문자열 역순 변환 메서드입니다.

=== "전"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Nextflow 스크립트에서 가져올 수 있는
     * 사용자 정의 함수를 구현합니다.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * 지정된 대상에게 인사합니다.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

이 함수의 핵심 구성 요소:

- **`@Function`**: 메서드를 Nextflow 워크플로우에서 호출 가능하게 만듭니다.
- **`String reverseGreeting(String greeting)`**: String을 받아 String을 반환합니다.
- **`greeting.reverse()`**: Groovy의 내장 문자열 역순 변환 메서드입니다.

!!! tip "public 메서드와 private 메서드"

    `@Function`이 없는 메서드는 Nextflow 워크플로우에 노출되지 않습니다.
    워크플로우 네임스페이스에 노출될 걱정 없이 클래스에 보조 메서드를 추가할 수 있습니다.

---

## 3. 플러그인 빌드 및 설치

플러그인을 빌드하고 설치합니다:

```bash
make install
```

!!! tip "빌드가 실패하는 경우"

    오류 메시지를 주의 깊게 읽으세요. 보통 줄 번호와 문제 설명이 포함되어 있습니다.
    일반적인 원인으로는 문법 오류(괄호나 따옴표 누락), 클래스 이름 오타, 타입 불일치 등이 있습니다.
    해결이 어렵다면 예제와 코드를 한 글자씩 비교해 보세요.

---

## 4. 워크플로우에서 함수 사용

플러그인이 빌드되고 설치되었습니다.
다음 단계는 `reverseGreeting`을 워크플로우에서 사용하여 전체 동작을 검증하는 것입니다.

파이프라인 디렉토리로 돌아갑니다:

```bash
cd ..
```

`greet.nf`를 편집하여 `reverseGreeting`을 가져오고 사용합니다:

=== "후"

    ```groovy title="greet.nf" hl_lines="4 23-25" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "전"

    ```groovy title="greet.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

파이프라인을 실행합니다:

```bash
nextflow run greet.nf
```

??? example "출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Output: Hello
    Output: Bonjour
    Output: Holà
    Output: Ciao
    Output: Hallo
    Pipeline complete! 👋
    ```

첫 번째 사용자 정의 플러그인 함수가 실제 워크플로우에서 동작하고 있습니다.
파트 1에서 nf-hello 및 nf-schema와 함께 사용했던 `include { ... } from 'plugin/...'` 패턴이 직접 만든 플러그인에서도 동일하게 작동합니다.

---

## 5. decorateGreeting 추가

플러그인은 여러 함수를 제공할 수 있습니다.
인사말을 장식 마커로 감싸는 두 번째 함수를 추가합니다. 파트 6에서 이를 설정 가능하게 만들 예정입니다.

`GreetingExtension.groovy`를 편집하여 `reverseGreeting` 뒤, 클래스의 닫는 중괄호 앞에 `decorateGreeting`을 추가합니다:

=== "후"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * 인사말 문자열을 역순으로 변환합니다
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

        /**
         * 인사말에 축하 마커를 추가합니다
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Groovy 문자열 보간: `#!groovy ${...}`는 변수의 값을 문자열에 삽입합니다.

=== "전"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * 인사말 문자열을 역순으로 변환합니다
         */
        @Function
        String reverseGreeting(String greeting) {
            return greeting.reverse()
        }

    }
    ```

이 함수는 Groovy 문자열 보간(`"*** ${greeting} ***"`)을 사용하여 인사말 변수를 문자열 안에 삽입합니다.

빌드, 설치 후 워크플로우를 업데이트합니다:

```bash
cd nf-greeting && make install && cd ..
```

`greet.nf`를 업데이트하여 `decorateGreeting`도 가져오고 사용합니다:

=== "후"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // 플러그인에서 사용자 정의 함수를 가져옵니다
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // 사용자 정의 플러그인 함수를 사용하여 인사말을 꾸밉니다
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        // reverseGreeting 함수 사용을 시연합니다
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. 동일한 플러그인에서 여러 함수를 가져오려면 각각 별도의 `include` 구문이 필요합니다.
    2. 플러그인 함수는 process `script:` 블록 내부에서도 동작합니다.

=== "전"

    ```groovy title="greet.nf" linenums="1" hl_lines="4 12 15 28"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

```bash
nextflow run greet.nf
```

??? example "출력"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Decorated: *** Hello ***
    Decorated: *** Bonjour ***
    Decorated: *** Holà ***
    Decorated: *** Ciao ***
    Decorated: *** Hallo ***
    Pipeline complete! 👋
    ```

플러그인 함수는 process 스크립트(`SAY_HELLO` 내부의 `decorateGreeting`)와 워크플로우 연산(`map` 내부의 `reverseGreeting`) 모두에서 동작합니다.

---

## 핵심 정리

이 섹션에서 학습한 내용:

- 함수는 `PluginExtensionPoint` 서브클래스에서 `@Function` 어노테이션으로 정의합니다.
- `include`로 가져온 플러그인 함수는 직접 만든 플러그인이든 기존 플러그인이든 동일하게 동작합니다.
- 플러그인 함수는 process 스크립트와 워크플로우 연산 모두에서 사용할 수 있습니다.

---

## 다음 단계

함수가 동작하지만, 지금까지는 전체 파이프라인을 실행하고 출력을 눈으로 확인하는 방식으로만 검증했습니다.
이 방식은 확장성이 없습니다. 함수가 늘어날수록, 특히 변경 후에도 각 함수가 올바르게 동작하는지 더 빠르게 확인할 방법이 필요합니다.
다음 섹션에서는 파이프라인을 실행하지 않고도 개별 함수를 자동으로 검증할 수 있는 단위 테스트를 소개합니다.

[파트 4로 계속 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
