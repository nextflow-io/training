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
    greetings_ch = sayHello(names_ch)
    // translate the greeting
    translated_ch = translateGreeting(greetings_ch, language)
    // collect all the greetings into one file
    collected_out = collectGreetings(translated_ch.collect(), batch_name, language)

    emit:
    greetings = greetings_ch
    translated = translated_ch
    collected = collected_out.outfile
    report = collected_out.report
}
