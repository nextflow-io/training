# Part 5: Trace Observers

Trace observers let your plugin respond to workflow events, such as a task completing, a file being published, or the pipeline finishing.
This enables use cases like custom reports, Slack notifications, metrics collection, or integration with external monitoring systems.
In this section, you'll build an observer that counts completed tasks and prints a summary.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 4 to use as your starting point:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Understanding the existing trace observer

The "Pipeline is starting!" message when you ran the pipeline came from the `GreetingObserver` class in your plugin.

Look at the observer code:

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

1. Interface for hooking into workflow lifecycle events
2. Called when the workflow starts; receives the session for accessing config
3. Called when the workflow finishes successfully

There are two things to notice here:

1. **`class GreetingObserver implements TraceObserver`**: `TraceObserver` is an interface defined by Nextflow. If your class implements this interface, Nextflow can hook into it and call your methods when events happen.
2. **`@Override`**: The `TraceObserver` interface defines methods like `onFlowCreate` and `onFlowComplete`. When you write methods with these names and add the `@Override` annotation, Nextflow calls them at the appropriate time. Any methods you don't override are ignored.

The full set of lifecycle events you can hook into at the time of writing are:

| Method              | When it's called        |
| ------------------- | ----------------------- |
| `onFlowCreate`      | Workflow starts         |
| `onFlowComplete`    | Workflow finishes       |
| `onProcessStart`    | A task begins execution |
| `onProcessComplete` | A task finishes         |
| `onProcessCached`   | A cached task is reused |
| `onFilePublish`     | A file is published     |

For a complete list, see the [TraceObserver interface](https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/trace/TraceObserver.groovy) in the Nextflow source.

---

## 2. Add a task counter observer

The goal is to build an observer that counts completed tasks and prints a summary at the end.
Adding a new observer to a plugin requires two things: writing the observer class, and registering it in the factory so Nextflow loads it.

### 2.1. Create a minimal observer

Create a new file:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Start with the simplest possible observer that prints a message when any task completes:

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

1. Import the required classes: `TraceObserver`, `TaskHandler`, and `TraceRecord`
2. Create a class that `implements TraceObserver`
3. Override `onProcessComplete` to run code when a task finishes

This is the minimum needed:

- Import the required classes (`TraceObserver`, `TaskHandler`, `TraceRecord`)
- Create a class that `implements TraceObserver`
- Override `onProcessComplete` to do something when a task finishes

### 2.2. Register the observer

The `GreetingFactory` creates observers.
Take a look at it:

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

Edit `GreetingFactory.groovy` to add the new observer:

=== "After"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

=== "Before"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3"
    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
    ```

!!! note "Groovy list syntax"

    We've replaced the Java-style `List.<TraceObserver>of(...)` with Groovy's simpler list literal `[...]`.
    Both return a `Collection`, but the Groovy syntax is more readable when adding multiple items.

### 2.3. Build, install, and test

```bash
cd nf-greeting && make assemble && make install && cd ..
nextflow run greet.nf -ansi-log false
```

You should see "✓ Task completed!" printed five times (once per task), interleaved with the existing pipeline output:

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

The observer is working.
Each time a task finishes, Nextflow calls `onProcessComplete`, and our implementation prints a message.

??? exercise "Customize the message"

    Try changing the message in `onProcessComplete` to something of your own, rebuild, and rerun.
    This confirms the full edit-build-run cycle works for observers.

### 2.4. Add counting logic

The minimal observer proves the hook works, but it doesn't track anything.

A class can hold variables (called fields or instance variables) that persist for the lifetime of the object.
This means an observer can accumulate state across multiple events during a pipeline run.

The next version adds a counter variable (`taskCount`) that starts at zero.
Each time a task completes, the counter goes up by one.
When the entire workflow finishes, the observer prints the final total.

Update `TaskCounterObserver.groovy` with the highlighted changes:

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

1. `taskCount` is a variable that belongs to the observer object. It keeps its value between method calls, so it can accumulate a count across the whole workflow run. `private` means only this class can access it.
2. `taskCount++` adds one to the counter. This line runs every time a task completes, so the count grows as the workflow progresses.
3. `onFlowComplete` is a second lifecycle hook. It runs once when the workflow finishes, making it a good place to print a summary.

In summary:

- `taskCount` persists across method calls, accumulating a count over the whole run
- `onProcessComplete` increments the counter and prints the running total each time a task finishes
- `onFlowComplete` runs once at the end, printing the final count

Rebuild and test:

```bash
cd nf-greeting && make assemble && make install && cd ..
nextflow run greet.nf -ansi-log false
```

??? example "Output"

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

    The counter messages are interleaved with task submissions because observers run as tasks complete.

!!! tip "Why `-ansi-log false`?"

    By default, Nextflow's ANSI progress display overwrites previous lines to show a clean, updating view of progress.
    This means you'd only see the *final* task count, not the intermediate messages.

    Using `-ansi-log false` disables this behavior and shows all output sequentially, which is essential when testing observers that print messages during execution.

---

## 3. Track published files

The observer can also respond when files are published.
The `onFilePublish` method receives the destination and source paths, which you can use to log, validate, or process published outputs.

### 3.1. Add a publish directory

First, update `greet.nf` so the `SAY_HELLO` process publishes its output files:

=== "After"

    ```groovy title="greet.nf" linenums="10" hl_lines="2"
    process SAY_HELLO {
        publishDir 'results'
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Use our custom plugin function to decorate the greeting
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

=== "Before"

    ```groovy title="greet.nf" linenums="10"
    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Use our custom plugin function to decorate the greeting
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated' > greeting.txt
        """
    }
    ```

### 3.2. Add the onFilePublish method

Add an `onFilePublish` method and the required import to `TaskCounterObserver.groovy`:

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

### 3.3. Build and test

```bash
cd nf-greeting && make assemble && make install && cd ..
nextflow run greet.nf -ansi-log false
```

You should see "Published:" messages for each output file alongside the task counter output:

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

The `onFilePublish` method fires each time Nextflow publishes a file to the `results` directory.
This pattern is useful for building audit logs, triggering downstream actions, or validating outputs as they are produced.

---

## Takeaway

You learned that:

- Trace observers hook into workflow lifecycle events like `onFlowCreate`, `onProcessComplete`, `onFilePublish`, and `onFlowComplete`
- Create observers by implementing `TraceObserver` and registering them in a Factory
- Observers can hold instance variables to accumulate state across events
- Observers are useful for custom logging, metrics collection, notifications, and reporting

---

## What's next?

The task counter works, but it's always on.
In a real plugin, users should be able to enable or disable features, or adjust behavior, from `nextflow.config` without editing the plugin source code.
The next section shows how to make your observer configurable and how to share your finished plugin with others.

[Continue to Part 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
