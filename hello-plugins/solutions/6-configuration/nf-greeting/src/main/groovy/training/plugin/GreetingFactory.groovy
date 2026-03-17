package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.trace.TraceObserver
import nextflow.trace.TraceObserverFactory

@CompileStatic
class GreetingFactory implements TraceObserverFactory {

    @Override
    Collection<TraceObserver> create(Session session) {
        // Read from GreetingConfig scope (GreetingConfig class provides validation)
        final enabled = session.config.navigate('greeting.enabled', true)
        if (!enabled) return []

        final verbose = session.config.navigate('greeting.taskCounter.verbose', true) as boolean
        return [
            new GreetingObserver(),
            new TaskCounterObserver(verbose)
        ]
    }
}
