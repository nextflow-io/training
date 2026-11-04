package training.plugin

import groovy.transform.CompileStatic
import nextflow.processor.TaskHandler
import java.nio.file.Path
import java.util.concurrent.atomic.AtomicInteger
import nextflow.trace.TraceObserver
import nextflow.trace.TraceRecord

/**
 * Observer that counts completed tasks
 */
@CompileStatic
class TaskCounterObserver implements TraceObserver {

    private final AtomicInteger taskCount = new AtomicInteger()

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        final count = taskCount.incrementAndGet()
        println "📊 Tasks completed so far: ${count}"
    }

    @Override
    void onFilePublish(Path destination, Path source) {
        println "📁 Published: ${destination.fileName}"
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount.get()}"
    }
}
