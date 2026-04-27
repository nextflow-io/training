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
        def session = Mock(Session) {
            getConfig() >> [greeting: [enabled: true]]
        }
        when:
        def result = factory.create(session)
        then:
        result.size() == 2
        result[0] instanceof GreetingObserver
        result[1] instanceof TaskCounterObserver
    }

    def 'should return empty list when disabled' () {
        given:
        def factory = new GreetingFactory()
        def session = Mock(Session) {
            getConfig() >> [greeting: [enabled: false]]
        }
        when:
        def result = factory.create(session)
        then:
        result.isEmpty()
    }

}
