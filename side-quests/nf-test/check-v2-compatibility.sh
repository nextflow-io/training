#!/usr/bin/env bash
#
# Diagnostic script to verify nf-test works with NXF_SYNTAX_PARSER=v2
# Run this in a Codespace after container setup completes
#
# Usage:
#   ./check-v2-compatibility.sh          # Check if nf-test works with v2
#   ./check-v2-compatibility.sh --fix    # Install fixed nf-test and check
#

set -e

NFTEST_COMMIT="350bb147a23a7f0aa657c13342d9726c0e3edacc"

# Handle --fix flag
if [[ "$1" == "--fix" ]]; then
    echo "=== Installing nf-test with v2 parser support ==="
    echo ""

    # Install Maven if not present
    if ! command -v mvn &> /dev/null; then
        echo "Installing Maven..."
        sudo apt-get update -qq && sudo apt-get install -y -qq maven > /dev/null
    fi

    echo "Building nf-test from commit ${NFTEST_COMMIT}..."
    cd /tmp
    rm -rf nf-test
    git clone --quiet https://github.com/askimed/nf-test.git
    cd nf-test
    git checkout --quiet ${NFTEST_COMMIT}
    mvn install -DskipTests -q
    cp target/nf-test.jar ~/.nf-test/nf-test.jar
    rm -rf /tmp/nf-test
    echo "Done!"
    echo ""
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== nf-test v2 Parser Compatibility Check ==="
echo ""

# Check environment
echo "Environment:"
echo "  NXF_SYNTAX_PARSER=${NXF_SYNTAX_PARSER:-not set}"
echo "  Working directory: $(pwd)"
echo ""

# Check nf-test version
echo "nf-test version:"
nf-test version | head -5
echo ""

# Clean up any previous test artifacts
rm -rf .nf-test tests nf-test.config work results 2>/dev/null || true

# Initialize and generate a process test
echo "Generating process test..."
nf-test init >/dev/null
nf-test generate process main.nf >/dev/null

# Add input to the test
sed -i 's|// input\[0\] = file("test-file.txt")|input[0] = "hello"|' tests/main.sayhello.nf.test

# Run the test - this is the critical check
echo "Running process-level test with v2 parser..."
echo ""
if nf-test test tests/main.sayhello.nf.test; then
    echo ""
    echo "=========================================="
    echo "SUCCESS: nf-test works with v2 parser!"
    echo "=========================================="
    EXIT_CODE=0
else
    echo ""
    echo "=========================================="
    echo "FAILURE: nf-test process tests failed"
    echo ""
    echo "If you see 'Unexpected input: *' error, run:"
    echo ""
    echo "  ./check-v2-compatibility.sh --fix"
    echo ""
    echo "This will build nf-test from main with the fix."
    echo "=========================================="
    EXIT_CODE=1
fi

# Clean up
rm -rf .nf-test tests nf-test.config work results 2>/dev/null || true

exit $EXIT_CODE
