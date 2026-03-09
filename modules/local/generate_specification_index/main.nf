/**
 * Converts an item to its string representation.
 *
 * Handles Path objects by converting them to URI strings,
 * and all other types by calling toString().
 *
 * @param item the object to convert (Path or any other type)
 * @return String representation of the item
 */
 def convertString(item) {
     if (item instanceof Path) {
         return item.toUriString()
     } else if (item instanceof List && item.isEmpty()) {
         return ""
     } else {
         return item.toString()
     }
 }

/**
 * Recursively converts all values in a map to strings.
 *
 * Traverses nested maps and lists, converting all leaf values
 * (including Path objects) to their string representations.
 * Preserves the map structure and key names.
 *
 * @param m the Map to stringify
 * @return a new Map with all values converted to strings
 */
def stringifyMap(Map m) {
    m.collectEntries { k, v ->
        def stringifiedValue = null

        if (v instanceof Map) {
            stringifiedValue = stringifyMap(v)
        } else if (v instanceof List) {
            stringifiedValue = v.isEmpty() ? "" : v.collect { item ->
               convertString(item)
            }
        } else {
            stringifiedValue = convertString(v)
        }

        [k, stringifiedValue]
    }
}

process GENERATE_SPECIFICATION_INDEX {
    executor 'local'
    tag { "${spec.id}" }

    input:
    val(spec)
    val(versions)

    output:
    val(out_spec), emit: spec

    exec:
    def prefix = task.ext.prefix ?: "${spec.id}"

    // Define output JSON file
    def json_file = file("${task.workDir}/${prefix}.json")

    // Extract the data we want from the map and clean it up
    def json_spec = spec.subMap(["data", "params", "tools"])

    // Replace tools list with tools map
    json_spec.tools = json_spec.tools.collectEntries { stage, tools ->
        if (tools.isEmpty()) {
            return  // skip this stage
        }

        def stageVersions = versions.subMap(tools)
            .values()
            .inject([:]) { a, b -> a + b }
        [stage, stageVersions]
    }

    // Filter to remove unused datasets
    json_spec.data = json_spec.data.findAll { _datatype, datamap -> datamap.id }

    // Convert all path objects in the spec to a string
    json_spec = json_spec.collectEntries { k, v ->
        [k, (v instanceof Map) ? stringifyMap(v) : convertString(v)]
    }

    // Write JSON index file
    json_file.text = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson(json_spec)
    )

    // Return the spec with an index JSON file field
    out_spec = spec + [index: json_file]
}
