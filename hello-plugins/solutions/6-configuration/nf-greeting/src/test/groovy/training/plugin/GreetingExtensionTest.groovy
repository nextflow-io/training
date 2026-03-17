package training.plugin

import spock.lang.Specification

/**
 * Tests for the greeting extension functions
 */
class GreetingExtensionTest extends Specification {

    def 'should reverse a greeting'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.reverseGreeting('Hello') == 'olleH'
        ext.reverseGreeting('Bonjour') == 'ruojnoB'
    }

    def 'should decorate a greeting'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.decorateGreeting('Hello') == '*** Hello ***'
    }

    def 'should create friendly greeting with default name'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.friendlyGreeting('Hello') == 'Hello, World!'
    }

    def 'should create friendly greeting with custom name'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.friendlyGreeting('Hello', 'Alice') == 'Hello, Alice!'
    }
}
