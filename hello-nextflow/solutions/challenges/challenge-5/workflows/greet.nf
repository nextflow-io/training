// Include modules
include { sayHello } from '../modules/sayHello.nf'
include { translateGreeting } from '../modules/translateGreeting.nf'
include { collectGreetings } from '../modules/collectGreetings.nf'

workflow GREET {

    take:
    names_ch
    language
    batch_name

    main:
    // emit a greeting
    sayHello(names_ch)
    // translate the greeting
    translateGreeting(sayHello.out, language)
    // collect all the greetings into one file
    collectGreetings(translateGreeting.out.collect(), batch_name, language)

    emit:
    greetings = sayHello.out
    translated = translateGreeting.out
    collected = collectGreetings.out.outfile
    report = collectGreetings.out.report
}
