# Part 5: Trace Observers

Trace observers let your plugin respond to workflow events.
In this section, you'll build an observer that counts completed tasks and prints a summary.

This enables use cases like custom reports, Slack notifications, metrics collection, or integration with external monitoring systems.

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
@Slf4j                                              // (1)!
@CompileStatic
class GreetingObserver implements TraceObserver {    // (2)!

    @Override
    void onFlowCreate(Session session) {            // (3)!
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {                         // (4)!
        println "Pipeline complete! 👋"
    }
}
```

1. Adds a logger for writing to Nextflow's log file
2. Interface for hooking into workflow lifecycle events
3. Called when the workflow starts; receives the session for accessing config
4. Called when the workflow finishes successfully

There are three things to notice here:

1. **The class implements `TraceObserver`**: `TraceObserver` is an interface defined by Nextflow. An interface is a contract: it defines a set of methods that your class promises to provide. By writing `implements TraceObserver`, you're telling Nextflow "this class can respond to workflow events."
2. **Override methods for specific events**: The `TraceObserver` interface defines methods like `onFlowCreate` and `onFlowComplete`. You override the ones you care about to provide your own behavior. Any methods you don't override are simply ignored.
3. **The factory registers the observer**: `GreetingObserver` doesn't register itself. The `GreetingFactory` class (which you'll see in section 2.2) creates observer instances and hands them to Nextflow. Without this registration step, Nextflow wouldn't know about the observer.

The full set of lifecycle events you can hook into:

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
Rather than modifying the existing `GreetingObserver`, you'll create a separate observer class for this.
Adding a new observer to a plugin requires two things: writing the observer class, and registering it in the factory so Nextflow loads it.

Here's the plan:

1. Write a minimal `TaskCounterObserver.groovy` that prints a message on task completion
2. Register it in `GreetingFactory.groovy` so Nextflow creates it
3. Build and verify it works
4. Enhance the observer with counting logic and a summary

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

Update `TaskCounterObserver.groovy` to track a count and report a summary:

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

1. A private instance variable that persists across method calls
2. Increment the counter and print the running total each time a task completes
3. `onFlowComplete` is called once when the workflow finishes, perfect for a summary

The key additions:

- A private instance variable `taskCount` persists across method calls
- `onProcessComplete` increments the counter and prints the running total
- `onFlowComplete` is called once when the workflow finishes, perfect for a summary

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

## Takeaway

You learned that:

- Trace observers hook into workflow lifecycle events like `onFlowCreate`, `onProcessComplete`, and `onFlowComplete`
- Create observers by implementing `TraceObserver` and registering them in a Factory
- Observers are useful for custom logging, metrics collection, notifications, and reporting

---

## What's next?

The next section shows how plugins can read configuration from `nextflow.config`, and how to share your plugin.

[Continue to Part 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
