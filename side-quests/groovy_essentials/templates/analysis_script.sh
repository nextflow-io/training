#!/bin/bash

# Nextflow template file - accessed via template directive in process
# This template has access to all variables from the process input
# Groovy expressions are evaluated at runtime

echo "Generating report for sample: ${meta.id}"
echo "Organism: ${meta.organism}"
echo "Quality score: ${meta.quality}"

# Conditional logic in template
<% if (meta.organism == 'human') { %>
echo "Including human-specific quality metrics"
human_qc_script.py --input ${results} --output ${meta.id}_report.html
<% } else { %>
echo "Using standard quality metrics for ${meta.organism}"
generic_qc_script.py --input ${results} --output ${meta.id}_report.html
<% } %>

# Groovy variables can be used for calculations
<%
def priority_bonus = meta.priority == 'high' ? 0.1 : 0.0
def adjusted_score = (meta.quality + priority_bonus).round(2)
%>

echo "Adjusted quality score: ${adjusted_score}"
echo "Report generation complete"
