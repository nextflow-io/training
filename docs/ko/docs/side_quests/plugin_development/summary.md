# 요약

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } AI 지원 번역 - [자세히 알아보기 및 개선 제안](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

플러그인 개발 교육 과정을 완료하셨습니다.
이 페이지에서는 각 파트에서 구현한 내용을 정리하고, 배포 방법과 다음 단계에 대한 안내를 제공합니다.

---

## 학습 내용 정리

### 파트 1: 플러그인 사용

사용자 관점에서 Nextflow 플러그인의 동작 방식을 학습했습니다.
nf-schema와 nf-co2footprint를 설치하고 `nextflow.config`를 통해 설정하였으며, 플러그인이 입력값을 검증하고, 함수를 추가하고, 파이프라인 생명주기 이벤트에 연결되는 방식을 확인했습니다.

### 파트 2: 환경 설정

Java 21+ 환경에서 플러그인 개발 환경을 구성하고, `nextflow plugin create` 명령어로 새 플러그인 프로젝트를 생성했습니다. 또한 Nextflow가 요구하는 프로젝트 구조인 소스 파일, 빌드 설정, Makefile 워크플로우를 학습했습니다.

### 파트 3: 커스텀 함수

`PluginExtensionPoint` 클래스에 `@Function` 어노테이션이 적용된 메서드를 작성하여 첫 번째 확장 포인트를 구현했습니다.
`reverseGreeting`과 `decorateGreeting`을 구현한 후, 파이프라인 스크립트에서 가져와 호출했습니다.

### 파트 4: 테스트

Groovy 테스트 프레임워크를 사용하여 커스텀 함수에 대한 단위 테스트를 작성했습니다.
`make test`로 테스트를 실행하고, 플러그인을 설치하기 전에 올바르게 동작하는지 검증하는 방법을 학습했습니다.

### 파트 5: 옵저버

`TraceObserver` 인터페이스를 구현하여 파이프라인 생명주기 이벤트에 연결했습니다.
파이프라인 시작 및 완료 시 반응하는 `GreetingObserver`와 완료된 작업을 카운트하는 `TaskCounterObserver`를 구현한 후, `TraceObserverFactory`를 통해 등록했습니다.

### 파트 6: 설정

`session.config.navigate()`를 사용하여 런타임에 값을 읽어 `nextflow.config`를 통해 플러그인을 설정 가능하게 만들었습니다.
`@ConfigScope` 클래스를 추가하여 플러그인 옵션을 공식적으로 선언함으로써 "Unrecognized config option" 경고를 제거하고 IDE 지원을 활성화했습니다.

---

## 배포

플러그인이 로컬에서 정상적으로 동작하면 Nextflow 플러그인 레지스트리를 통해 다른 사람들과 공유할 수 있습니다.

### 버전 관리

릴리스 시 [시맨틱 버전 관리](https://semver.org/)를 따르세요:

| 버전 변경                 | 사용 시점                     | 예시                      |
| ------------------------- | ----------------------------- | ------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | 호환성이 깨지는 변경          | 함수 제거, 반환 타입 변경 |
| **MINOR** (1.0.0 → 1.1.0) | 하위 호환 가능한 새 기능 추가 | 새 함수 추가              |
| **PATCH** (1.0.0 → 1.0.1) | 하위 호환 가능한 버그 수정    | 기존 함수의 버그 수정     |

각 릴리스 전에 `build.gradle`의 버전을 업데이트하세요:

```groovy title="build.gradle"
version = '1.0.0'  // 시맨틱 버전 관리 사용: MAJOR.MINOR.PATCH
```

### 레지스트리에 게시

[Nextflow 플러그인 레지스트리](https://registry.nextflow.io/)는 커뮤니티와 플러그인을 공유하는 공식 방법입니다.

게시 절차:

1. [레지스트리](https://registry.nextflow.io/)에서 플러그인 이름을 등록합니다 (GitHub 계정으로 로그인)
2. `~/.gradle/gradle.properties`에 API 자격 증명을 설정합니다
3. 테스트를 실행하여 모든 것이 정상 동작하는지 확인합니다: `make test`
4. `make release`로 게시합니다

단계별 안내는 [공식 게시 문서](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin)를 참조하세요.

게시 후에는 로컬 설정 없이 플러그인을 설치할 수 있습니다:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow는 첫 사용 시 레지스트리에서 플러그인을 자동으로 다운로드합니다.

---

## 플러그인 개발 체크리스트

- [ ] Java 21+ 설치
- [ ] `nextflow plugin create <name> <org>`로 프로젝트 생성
- [ ] `@Function` 메서드가 포함된 확장 클래스 구현
- [ ] 단위 테스트 작성 및 `make test`로 실행
- [ ] `make install`로 빌드 및 설치
- [ ] 필요 시 워크플로우 이벤트를 위한 `TraceObserver` 구현 추가
- [ ] 필요 시 플러그인 설정을 위한 `ConfigScope` 추가
- [ ] `nextflow.config`에서 `plugins { id 'plugin-id' }`로 활성화
- [ ] `include { fn } from 'plugin/plugin-id'`로 함수 가져오기
- [ ] 버전 관리 후 레지스트리에 게시

---

## 주요 코드 패턴

**함수 정의:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**플러그인 설정:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**워크플로우에서 사용:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## 확장 포인트 요약

| 타입                | 클래스/어노테이션 | 목적                                   |
| ------------------- | ----------------- | -------------------------------------- |
| Function            | `@Function`       | 워크플로우에서 호출 가능               |
| Trace Observer      | `TraceObserver`   | 워크플로우 생명주기 이벤트에 연결      |
| Configuration Scope | `@ScopeName`      | nextflow.config에서 플러그인 설정 정의 |

---

## 다음 단계

플러그인 개발을 계속하기 위한 실용적인 다음 단계를 소개합니다.

**실제 문제를 해결해 보세요.**
직접 업무에서 활용할 수 있는 사례를 선택하세요. 팀에서 반복적으로 사용하는 커스텀 함수, 파이프라인 완료 시 Slack 알림을 보내는 옵저버, 또는 조직의 파이프라인 전반에 걸쳐 옵션을 표준화하는 config scope 등이 좋은 예입니다.
실제 문제에서 출발하는 것이 이해를 심화하는 가장 빠른 방법입니다.

**nf-hello를 참조 자료로 활용하세요.**
[nf-hello](https://github.com/nextflow-io/nf-hello) 저장소는 공식 최소 플러그인 예제입니다.
새 프로젝트의 시작점으로 적합하며, 구조를 확인할 때 유용한 참조 자료입니다.

**공식 문서를 읽어보세요.**
Nextflow 문서에는 이 교육 과정 외에도 채널 팩토리, 연산자 오버로딩, 고급 옵저버 패턴 등 다양한 주제가 포함되어 있습니다.
[플러그인 개발](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) 가이드가 가장 포괄적인 참조 자료입니다.

**기존 플러그인을 분석해 보세요.**
[Nextflow 플러그인 저장소](https://github.com/nextflow-io/plugins)에는 nf-schema, nf-wave, nf-tower 등 공식 플러그인의 소스 코드가 포함되어 있습니다.
실제 운영 중인 플러그인 코드를 읽는 것은 입문 예제를 넘어서는 패턴과 관례를 학습하는 가장 좋은 방법 중 하나입니다.

---

## 추가 자료

**공식 문서:**

- [플러그인 사용](https://www.nextflow.io/docs/latest/plugins/plugins.html): 플러그인 설치 및 설정에 대한 종합 가이드
- [플러그인 개발](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): 상세한 플러그인 개발 참조 자료
- [Config 스코프](https://nextflow.io/docs/latest/developer/config-scopes.html): 플러그인을 위한 설정 스코프 생성

**플러그인 탐색:**

- [Nextflow 플러그인 레지스트리](https://registry.nextflow.io/): 사용 가능한 플러그인 탐색 및 검색
- [플러그인 레지스트리 문서](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): 레지스트리 문서

**예제 및 참조 자료:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): 간단한 예제 플러그인 (시작점으로 적합)
- [Nextflow 플러그인 저장소](https://github.com/nextflow-io/plugins): 참조용 공식 플러그인 모음
