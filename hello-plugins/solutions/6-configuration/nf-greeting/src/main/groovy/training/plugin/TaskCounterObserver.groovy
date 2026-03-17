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

    private final boolean verbose
    private int taskCount = 0

    TaskCounterObserver(boolean verbose) {
        this.verbose = verbose
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        taskCount++
        if (verbose) {
            println "📊 Tasks completed so far: ${taskCount}"
        }
    }

    @Override
    void onFlowComplete() {
        println "📈 Final task count: ${taskCount}"
    }
}
