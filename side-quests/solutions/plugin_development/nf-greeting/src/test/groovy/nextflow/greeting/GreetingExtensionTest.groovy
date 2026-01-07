package nextflow.greeting

import org.junit.jupiter.api.Test
import static org.junit.jupiter.api.Assertions.*

/**
 * Unit tests for GreetingExtension functions
 */
class GreetingExtensionTest {

    @Test
    void testReverseGreeting() {
        def ext = new GreetingExtension()

        assertEquals('olleH', ext.reverseGreeting('Hello'))
        assertEquals('ruojnoB', ext.reverseGreeting('Bonjour'))
        assertEquals('', ext.reverseGreeting(''))
    }

    @Test
    void testDecorateGreeting() {
        def ext = new GreetingExtension()

        assertEquals('*** Hello ***', ext.decorateGreeting('Hello'))
        assertEquals('*** Bonjour ***', ext.decorateGreeting('Bonjour'))
    }

    @Test
    void testFriendlyGreetingDefault() {
        def ext = new GreetingExtension()

        assertEquals('Hello, World!', ext.friendlyGreeting('Hello'))
        assertEquals('Bonjour, World!', ext.friendlyGreeting('Bonjour'))
    }

    @Test
    void testFriendlyGreetingWithName() {
        def ext = new GreetingExtension()

        assertEquals('Hello, Alice!', ext.friendlyGreeting('Hello', 'Alice'))
        assertEquals('Bonjour, Bob!', ext.friendlyGreeting('Bonjour', 'Bob'))
    }
}
