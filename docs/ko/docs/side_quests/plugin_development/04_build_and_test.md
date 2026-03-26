# 파트 4: 테스트

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

플러그인은 파이프라인 개발자가 신뢰해야 하는 단독 실행형 소프트웨어입니다.
각 기능을 파이프라인 외부에서 독립적으로 테스트하면, 누군가 워크플로우에 통합하기 전에 플러그인이 올바르게 작동하는지 확인할 수 있습니다.
이 섹션에서는 Spock 테스트 프레임워크를 사용하여 테스트를 작성하고 실행합니다.

!!! tip "여기서부터 시작하시나요?"

    파트 3의 해결책을 시작점으로 사용하려면 다음과 같이 복사하세요:

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    그런 다음 플러그인 디렉토리로 이동합니다:

    ```bash
    cd nf-greeting
    ```

플러그인 디렉토리에 있는지 확인합니다:

```bash
cd nf-greeting
```

---

## 1. 테스트가 필요한 이유

빌드가 성공했다는 것은 코드가 컴파일되었다는 의미이지, 예상대로 작동한다는 의미는 아닙니다.
단위 테스트는 함수가 주어진 입력에 대해 올바른 출력을 생성하는지 자동으로 확인하는 작은 코드 조각입니다.
예를 들어, `#!groovy reverseGreeting("Hello")`가 `"olleH"`를 반환하는지 확인하는 테스트를 작성할 수 있습니다.

테스트가 유용한 이유는 다음과 같습니다:

- 사용자보다 먼저 버그를 발견할 수 있습니다
- 기존 기능을 손상시키지 않고 변경할 수 있다는 자신감을 줍니다
- 함수를 어떻게 사용해야 하는지 보여주는 문서 역할을 합니다

---

## 2. Spock 테스트 이해하기

플러그인 템플릿은 Groovy용 테스트 프레임워크인 [Spock](https://spockframework.org/)을 사용합니다.
Spock은 이미 프로젝트에 설정되어 있으므로(`build.gradle`을 통해) 별도로 추가할 필요가 없습니다.

Python의 `pytest`나 R의 `testthat`과 같은 테스트 도구를 사용해 본 적이 있다면, Spock도 동일한 역할을 합니다. 알려진 입력으로 코드를 호출하고 출력을 확인하는 작은 함수를 작성합니다.
차이점은 Spock이 Nextflow process나 workflow와 유사한 레이블이 붙은 블록(`given:`, `expect:`, `when:`, `then:`)을 사용한다는 것입니다.

기본 구조는 다음과 같습니다:

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **따옴표로 감싼 테스트 이름**: 테스트가 무엇을 확인하는지 설명합니다. 평이한 영어로 작성합니다.
2. **`given:` 블록**: 테스트에 필요한 것을 설정합니다(객체 생성, 데이터 준비).
3. **`expect:` 블록**: 실제 검증 내용입니다. 테스트가 통과하려면 각 줄이 `true`여야 합니다.

이 구조는 테스트를 읽기 쉽게 만들어 줍니다: "extension 객체가 주어졌을 때, `reverseGreeting('Hello')`가 `'olleH'`와 같을 것으로 예상합니다."

---

## 3. 테스트 작성하기

파트 3에서 만든 두 함수 `reverseGreeting`과 `decorateGreeting`에 대한 테스트를 작성합니다.

### 3.1. 테스트 클래스 생성

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

편집기에서 파일을 열고 빈 테스트 클래스 골격을 추가합니다:

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * greeting extension 함수에 대한 테스트
 */
class GreetingExtensionTest extends Specification {  // (1)!

}
```

1. 모든 Spock 테스트 클래스는 `Specification`을 상속합니다. 이것이 모든 Spock 테스트 파일의 시작점입니다.

### 3.2. reverseGreeting 테스트

클래스 본문 안에 테스트 메서드를 추가합니다.
`given:` 블록은 `GreetingExtension` 인스턴스를 생성하고, `expect:` 블록은 `reverseGreeting`이 두 가지 다른 입력을 올바르게 뒤집는지 확인합니다.
이 테스트는 파이프라인을 실행하지 않고 함수를 직접 테스트합니다.

=== "후"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="10-17"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension 함수에 대한 테스트
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()            // (1)!

            expect:
            ext.reverseGreeting('Hello') == 'olleH'     // (2)!
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

    1. 파이프라인을 실행하지 않고 직접 테스트하기 위해 extension 인스턴스를 생성합니다.
    2. `expect:`의 각 줄은 단언(assertion)입니다. 모두 `true`일 때만 테스트가 통과합니다.

=== "전"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension 함수에 대한 테스트
     */
    class GreetingExtensionTest extends Specification {

    }
    ```

### 3.3. decorateGreeting 테스트

첫 번째 테스트 메서드 다음에 두 번째 테스트 메서드를 추가합니다.
이 테스트는 `decorateGreeting`이 입력 문자열의 양쪽에 `***`를 추가하는지 확인합니다.

=== "후"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18-25"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension 함수에 대한 테스트
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }

        def 'should decorate a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.decorateGreeting('Hello') == '*** Hello ***'
        }
    }
    ```

=== "전"

    ```groovy title="GreetingExtensionTest.groovy" linenums="1" hl_lines="18"
    package training.plugin

    import spock.lang.Specification

    /**
     * greeting extension 함수에 대한 테스트
     */
    class GreetingExtensionTest extends Specification {

        def 'should reverse a greeting'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('Hello') == 'olleH'
            ext.reverseGreeting('Bonjour') == 'ruojnoB'
        }
    }
    ```

---

## 4. 테스트 실행

```bash
make test
```

??? example "테스트 출력"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 6 executed
    ```

    **테스트 결과는 어디에 있나요?** 모든 테스트가 통과하면 Gradle은 상세 출력을 숨깁니다.
    "BUILD SUCCESSFUL"은 모든 것이 정상적으로 작동했다는 의미입니다.
    테스트가 실패하면 상세한 오류 메시지가 표시됩니다.

??? exercise "엣지 케이스 테스트 추가"

    `reverseGreeting`이 빈 문자열을 처리하는지 확인하는 테스트를 추가합니다.
    `reverseGreeting('')`은 무엇을 반환해야 할까요?
    테스트를 추가하고 `make test`를 실행하여 통과하는지 확인합니다.

    ??? solution "해결책"

        `GreetingExtensionTest.groovy`에 다음 테스트 메서드를 추가합니다:

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        빈 문자열을 뒤집어도 빈 문자열입니다.

---

## 5. 테스트 리포트 확인

Gradle은 각 테스트의 상세 결과를 담은 HTML 테스트 리포트를 생성합니다.
리포트 디렉토리에서 웹 서버를 시작합니다:

```bash
pushd build/reports/tests/test
python -m http.server
```

VS Code가 브라우저에서 애플리케이션을 열도록 안내합니다.
테스트 클래스를 클릭하면 개별 테스트 결과를 확인할 수 있습니다:

![모든 테스트가 통과된 테스트 리포트](./img/test_report.png)

리포트에는 각 테스트 메서드와 통과 또는 실패 여부가 표시됩니다.

++ctrl+c++를 눌러 서버를 중지한 후 이전 디렉토리로 돌아갑니다:

```bash
popd
```

메인 프로젝트 디렉토리로 돌아갑니다:

```bash
cd ..
```

---

## 핵심 정리

이 섹션에서 학습한 내용은 다음과 같습니다:

- Spock 테스트는 읽기 쉬운 `given:`/`expect:` 구조를 사용합니다
- `make test`로 테스트를 실행하고, `build/reports/tests/test/`에서 HTML 리포트를 확인할 수 있습니다
- 테스트는 동작을 검증하고 함수 사용 방법을 보여주는 문서 역할을 합니다

---

## 다음 단계

지금까지 플러그인에 파이프라인이 호출할 수 있는 사용자 정의 함수를 추가했습니다.
플러그인은 trace observer를 사용하여 워크플로우 이벤트(작업 완료, 파일 게시, 파이프라인 종료 등)에 반응할 수도 있습니다.
다음 섹션에서는 완료된 작업 수를 세고 파이프라인이 종료될 때 요약을 출력하는 observer를 만들어 봅니다.

[파트 5로 계속하기 :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
