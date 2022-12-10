/*
 * We can parse JSON formatted data using groovy JsonSlurper (!)
 */

import groovy.json.JsonSlurper

def f = file('data/meta/regions.json')
def records = new JsonSlurper().parse(f)


for( def entry : records ) {
  log.info "$entry.patient_id -- $entry.feature"
}


/*
 * This requires NF 20.04.0-edge
 *
 * When using older versions replace `parse(f)` with `parseText(f.text)`
 *
 */
