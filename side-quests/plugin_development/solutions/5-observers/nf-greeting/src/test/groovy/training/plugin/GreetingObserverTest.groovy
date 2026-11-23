package training.plugin

import nextflow.Session
import spock.lang.Specification

/**
 * Implements a basic factory test
 *
 */
class GreetingObserverTest extends Specification {

    def 'should create the observer instance' () {
        given:
        def factory = new GreetingFactory()
        when:
        def result = factory.create(Mock(Session))
        then:
        result.size() == 2
        result[0] instanceof GreetingObserver
        result[1] instanceof TaskCounterObserver
    }

}
