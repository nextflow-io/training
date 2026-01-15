# Part 5: Trace Observers

In Part 1, we saw that plugins can provide many types of extensions.
So far we've implemented custom functions.
This part explores **trace observers**, which let you hook into workflow lifecycle events.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 4 to use as your starting point:

    ```bash
    cp -r solutions/4-build-and-test/* .
    ```

---

## 1. Understanding the existing trace observer

Remember the "Pipeline is starting!" message when you ran the pipeline?
That came from the `GreetingObserver` class in your plugin.

Look at the observer code:

```bash
cat nf-greeting/src/main/groovy/training/plugin/GreetingObserver.groovy
```

This observer hooks into workflow lifecycle events.
Trace observers can respond to many events:

| Method              | When it's called        |
| ------------------- | ----------------------- |
| `onFlowCreate`      | Workflow starts         |
| `onFlowComplete`    | Workflow finishes       |
| `onProcessStart`    | A task begins execution |
| `onProcessComplete` | A task finishes         |
| `onProcessCached`   | A cached task is reused |
| `onFilePublish`     | A file is published     |

This enables powerful use cases like custom reports, Slack notifications, or metrics collection.

---

## 2. Try it: Add a task counter observer

Rather than modifying the existing observer, create a new one that counts completed tasks.
We'll build it up progressively: first a minimal version, then add features.

### 2.1. Create a minimal observer

Create a new file:

```bash
touch nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy
```

Start with the simplest possible observer.
Just print a message when any task completes:

```groovy title="nf-greeting/src/main/groovy/training/plugin/TaskCounterObserver.groovy" linenums="1"
package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer that responds to task completion
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        println "✓ Task completed!"
    }
}
```

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

```groovy title="GreetingFactory.groovy (starting point)"
@CompileStatic
class GreetingFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        return List.<TraceObserver>of(new GreetingObserver())
    }
}
```

Edit `GreetingFactory.groovy` to add our new observer:

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
nextflow run main.nf -ansi-log false
```

You should see "✓ Task completed!" printed five times (once per task):

```console title="Expected output (partial)"
...
[be/bd8e72] Submitted process > SAY_HELLO (2)
✓ Task completed!
[5b/d24c2b] Submitted process > SAY_HELLO (1)
✓ Task completed!
...
```

Our observer is responding to task completion events.
The next step makes it more useful.

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

    private int taskCount = 0

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++
        println "📊 Tasks completed so far: ${taskCount}"
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount}"
    }
}
```

The key additions:

- **Line 14**: A private instance variable `taskCount` persists across method calls
- **Lines 18-19**: Increment the counter and print the running total
- **Lines 22-24**: `onFlowComplete` is called once when the workflow finishes, perfect for a summary

Rebuild and test:

```bash
cd nf-greeting && make assemble && make install && cd ..
nextflow run main.nf -ansi-log false
```

```console title="Expected output"
N E X T F L O W  ~  version 25.10.2
Launching `main.nf` [pensive_engelbart] DSL2 - revision: 85fefd90d0
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

!!! tip "Why `-ansi-log false`?"

    By default, Nextflow's ANSI progress display overwrites previous lines to show a clean, updating view of progress.
    This means you'd only see the *final* task count, not the intermediate "Tasks completed so far" messages.
    They'd be overwritten as new output arrives.

    Using `-ansi-log false` disables this behavior and shows all output sequentially, which is essential when testing observers that print messages during execution.
    Without this flag, you might think your observer isn't working when it actually is.
    The output is just being overwritten.

---

## Takeaway

You learned that:

- Trace observers hook into workflow lifecycle events like `onFlowCreate`, `onProcessComplete`, and `onFlowComplete`
- Create observers by implementing `TraceObserver` and registering them in a Factory
- Observers are useful for custom logging, metrics collection, notifications, and reporting

---

## What's next?

The next section shows how plugins can read configuration from `nextflow.config`.

[Continue to Part 6 :material-arrow-right:](06_configuration.md){ .md-button .md-button--primary }
