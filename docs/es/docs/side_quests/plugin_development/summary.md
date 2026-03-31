# Resumen

<span class="ai-translation-notice">:material-information-outline:{ .ai-translation-notice-icon } Traducción asistida por IA - [más información y sugerencias](https://github.com/nextflow-io/training/blob/master/TRANSLATING.md)</span>

Has completado la capacitación de Desarrollo de Plugins.
Esta página resume lo que construiste en cada parte, cubre la distribución y proporciona orientación sobre los próximos pasos.

---

## Lo que aprendiste

### Parte 1: Uso de plugins

Descubriste cómo funcionan los plugins de Nextflow desde la perspectiva del usuario.
Instalaste nf-schema y nf-co2footprint, los configuraste mediante `nextflow.config`, y viste cómo los plugins pueden validar entradas, agregar funciones y conectarse a los eventos del ciclo de vida del pipeline.

### Parte 2: Configuración del entorno

Configuraste un entorno de desarrollo de plugins con Java 21+, creaste un nuevo proyecto de plugin usando el comando `nextflow plugin create`, y aprendiste la estructura de proyecto que Nextflow espera: archivos fuente, configuración de compilación y el workflow del Makefile.

### Parte 3: Funciones personalizadas

Implementaste tu primer punto de extensión creando métodos anotados con `@Function` en una clase `PluginExtensionPoint`.
Construiste `reverseGreeting` y `decorateGreeting`, luego los importaste y llamaste desde un script de pipeline.

### Parte 4: Pruebas

Escribiste pruebas unitarias para tus funciones personalizadas usando el framework de pruebas de Groovy.
Aprendiste a ejecutar las pruebas con `make test` y a verificar que tu plugin se comporta correctamente antes de instalarlo.

### Parte 5: Observadores

Implementaste la interfaz `TraceObserver` para conectarte a los eventos del ciclo de vida del pipeline.
Construiste `GreetingObserver` (que reacciona al inicio y la finalización del pipeline) y `TaskCounterObserver` (que cuenta las tareas completadas), luego los registraste a través de un `TraceObserverFactory`.

### Parte 6: Configuración

Hiciste tu plugin configurable mediante `nextflow.config` usando `session.config.navigate()` para leer valores en tiempo de ejecución.
Agregaste una clase `@ConfigScope` para declarar formalmente las opciones de tu plugin, eliminando las advertencias de "Unrecognized config option" y habilitando el soporte del IDE.

---

## Distribución

Una vez que tu plugin funciona localmente, puedes compartirlo con otros a través del registro de plugins de Nextflow.

### Versionado

Sigue el [versionado semántico](https://semver.org/) para tus versiones:

| Cambio de versión         | Cuándo usarlo                          | Ejemplo                                              |
| ------------------------- | -------------------------------------- | ---------------------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Cambios que rompen compatibilidad      | Eliminar una función, cambiar tipos de retorno       |
| **MINOR** (1.0.0 → 1.1.0) | Nuevas funcionalidades, retrocompatible | Agregar una nueva función                            |
| **PATCH** (1.0.0 → 1.0.1) | Corrección de errores, retrocompatible | Corregir un error en una función existente           |

Actualiza la versión en `build.gradle` antes de cada versión publicada:

```groovy title="build.gradle"
version = '1.0.0'  // Usa versionado semántico: MAJOR.MINOR.PATCH
```

### Publicación en el registro

El [registro de plugins de Nextflow](https://registry.nextflow.io/) es la forma oficial de compartir plugins con la comunidad.

El workflow de publicación:

1. Registra el nombre de tu plugin en el [registro](https://registry.nextflow.io/) (inicia sesión con tu cuenta de GitHub)
2. Configura tus credenciales de API en `~/.gradle/gradle.properties`
3. Ejecuta las pruebas para verificar que todo funciona: `make test`
4. Publica con `make release`

Para instrucciones paso a paso, consulta la [documentación oficial de publicación](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Una vez publicado, los usuarios instalan tu plugin sin ninguna configuración local:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow descarga automáticamente el plugin desde el registro en el primer uso.

---

## Lista de verificación para el desarrollo de plugins

- [ ] Java 21+ instalado
- [ ] Crear el proyecto con `nextflow plugin create <name> <org>`
- [ ] Implementar la clase de extensión con métodos `@Function`
- [ ] Escribir pruebas unitarias y ejecutarlas con `make test`
- [ ] Compilar e instalar con `make install`
- [ ] Opcionalmente agregar implementaciones de `TraceObserver` para eventos del workflow
- [ ] Opcionalmente agregar `ConfigScope` para la configuración del plugin
- [ ] Habilitar en `nextflow.config` con `plugins { id 'plugin-id' }`
- [ ] Importar funciones con `include { fn } from 'plugin/plugin-id'`
- [ ] Versionar y publicar en el registro

---

## Patrones de código clave

**Definición de función:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Configuración del plugin:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Uso en workflows:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Resumen de puntos de extensión

| Tipo                    | Clase/Anotación  | Propósito                                                    |
| ----------------------- | ---------------- | ------------------------------------------------------------ |
| Función                 | `@Function`      | Invocable desde workflows                                    |
| Observador de traza     | `TraceObserver`  | Conectarse a los eventos del ciclo de vida del workflow      |
| Ámbito de configuración | `@ScopeName`     | Definir la configuración del plugin en nextflow.config       |

---

## ¿Qué sigue?

Aquí hay algunos próximos pasos prácticos para continuar tu camino en el desarrollo de plugins.

**Construye algo real.**
Elige un caso de uso de tu propio trabajo: una función personalizada que tu equipo usa repetidamente, un observador que envía notificaciones de Slack al completarse el pipeline, o un ámbito de configuración que estandariza opciones en los pipelines de tu organización.
Comenzar desde un problema real es la forma más rápida de profundizar tu comprensión.

**Usa nf-hello como referencia.**
El repositorio [nf-hello](https://github.com/nextflow-io/nf-hello) es el ejemplo oficial mínimo de plugin.
Es un buen punto de partida para nuevos proyectos y una referencia útil cuando necesitas verificar cómo está estructurado algo.

**Lee la documentación oficial.**
La documentación de Nextflow cubre temas más allá de esta capacitación, incluyendo fábricas de canales, sobrecarga de operadores y patrones avanzados de observadores.
La guía de [desarrollo de plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) es la referencia más completa.

**Estudia los plugins existentes.**
El [repositorio de plugins de Nextflow](https://github.com/nextflow-io/plugins) contiene el código fuente de los plugins oficiales como nf-schema, nf-wave y nf-tower.
Leer código de plugins en producción es una de las mejores formas de aprender patrones y convenciones que van más allá de los ejemplos introductorios.

---

## Recursos adicionales

**Documentación oficial:**

- [Uso de plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): guía completa para instalar y configurar plugins
- [Desarrollo de plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): referencia detallada para el desarrollo de plugins
- [Ámbitos de configuración](https://nextflow.io/docs/latest/developer/config-scopes.html): creación de ámbitos de configuración para plugins

**Descubrimiento de plugins:**

- [Registro de Plugins de Nextflow](https://registry.nextflow.io/): explora y descubre los plugins disponibles
- [Documentación del registro de plugins](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): documentación del registro

**Ejemplos y referencias:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): plugin de ejemplo simple (excelente punto de partida)
- [Repositorio de plugins de Nextflow](https://github.com/nextflow-io/plugins): colección de plugins oficiales como referencia
